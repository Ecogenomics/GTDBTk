import logging
import multiprocessing as mp
import os
import subprocess
import tempfile
from collections import defaultdict

from gtdbtk.config.config import LOG_TASK
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.external.hmm_aligner import HmmAligner
from gtdbtk.files.marker.copy_number import CopyNumberFile
from gtdbtk.files.marker.tophit import TopHitPfamFile, TopHitTigrFile
from gtdbtk.files.marker_info import MarkerInfoFile
from gtdbtk.tools import tqdm_log


def get_single_copy_hits_worker(job):
    """For a given genome, obtain the PFAM and TIGRFAM tophit files. Use
    this information to determine what hits are single copy.

    Parameters
    ----------
    job : Tuple[str, str, CopyNumberFile]
        The genome id, path to called genes, and domain-specific copy number file object.

    Returns
    -------
    Dict[str, Dict[str, str]]
        dict[marker id][genome id] = sequence
    """
    gid, aa_path, copy_number_file = job

    # Load the marker top hit files.
    marker_genes_dir = os.path.dirname(os.path.dirname(aa_path))
    pfam_tophit_file = TopHitPfamFile(marker_genes_dir, gid)
    tigr_tophit_file = TopHitTigrFile(marker_genes_dir, gid)
    pfam_tophit_file.read()
    tigr_tophit_file.read()

    # Process each of the genes to determine if they are single copy.
    cnf = copy_number_file('/dev/null', None)
    cnf.add_genome(gid, aa_path, pfam_tophit_file, tigr_tophit_file)
    single_copy = cnf.get_single_copy_hits(gid)

    # Store the output
    out = defaultdict(dict)
    for marker_id, marker_d in single_copy.items():
        out[marker_id][gid] = marker_d['seq']
    return out


def get_single_copy_hits(gid_dict: dict, copy_number_file, cpus):
    """Collect all of the single copy hits (both domains) for each genome.

    Parameters
    ----------
    gid_dict : dict
        A dictionary containing information about the genome, indexed by the id.
    copy_number_file : CopyNumberFile
        A domain-specific subclass of the copy number file.
    cpus : int
        The number of CPUs to use in sub-processes.

    Returns
    -------
    Dict[str, Dict[str, str]]
        dict[marker id][genome id] = sequence
    """

    # Generate a queue job jobs.
    queue = list()
    for gid, gid_info in gid_dict.items():
        queue.append((gid,
                      gid_info['aa_gene_path'],
                      copy_number_file))

    # Process the queue.
    with mp.Pool(processes=cpus) as pool:
        results = list(tqdm_log(pool.imap_unordered(get_single_copy_hits_worker, queue),
                                total=len(queue), unit='genome'))

    # Re-format the results.
    out = defaultdict(dict)
    for result in results:
        for marker_id, marker_d in result.items():
            for gid, seq in marker_d.items():
                out[marker_id][gid] = seq
    return out


def read_hmmalign_output(output, expected_gids):
    """Reads the output of hmmalign and extracts the alignment.

    Parameters
    ----------
    output : str
        The string output from hmmalign.
    expected_gids : frozenset
        A set containing all of the gids expected in the output.

    Returns
    -------
    Dict[str, str]
        dict[genome id] = sequence
    """
    # In case gids have spaces.
    exp_mapping = {x.split(' ', 1)[0]: x for x in expected_gids}

    # Get the sequences and the mask.
    unmasked = dict()
    mask = None
    for line in output.splitlines():
        splitline = line.split(' ', 1)

        # Sequence
        if splitline[0] in exp_mapping:
            rsplitline = line.rsplit(" ", 1)
            hit_seq = rsplitline[-1]
            unmasked[exp_mapping[splitline[0]]] = hit_seq

        # Mask
        elif line[0:len("#=GC RF")] == "#=GC RF":
            mask = [x == 'x' for x in line.rsplit(' ', 1)[-1]]

    # Sanity check.
    if mask is None:
        raise GTDBTkExit(f'Unable to get mask from hmmalign result file: {output}')
    if len(unmasked) != len(expected_gids):
        raise GTDBTkExit(f'Not all genomes could be aligned: {output}')

    # Mask each of the sequences and return them.
    out = dict()
    for gid, seq in unmasked.items():
        out[gid] = ''.join([s for s, m in zip(seq, mask) if m is True])
    return out


def run_hmm_align_worker(job):
    """A worker process for running hmmalign.

    Parameters
    ----------
    job : Tuple[str, str, str, frozenset]
        The marker id, hmm marker path, called genes path, expected gids.

    Returns
    -------
    List[Tuple[str, str, str]]
        A list containing the (genome id, marker id, sequence).
    """
    marker_id, marker_path, marker_fa, expected_gids = job

    # Run the process and capture stdout.
    args = ["hmmalign", "--outformat", "Pfam", marker_path, marker_fa]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = proc.communicate()

    # Exit if an error was raised.
    if proc.returncode != 0:
        arg_str = ' '.join(args)
        raise GTDBTkExit(f'hmmalign returned a non-zero exit code: {arg_str}')

    # Process the output and return the sequences.
    seqs = read_hmmalign_output(stdout, expected_gids)
    return [(gid, marker_id, seq) for gid, seq in seqs.items()]


def create_concat_alignment(list_seqs, marker_info_file):
    """Create a concatenated alignment where missing markers are filled with gaps.

    Parameters
    ----------
    list_seqs : list[tuple[str, str, str]]
        A list containing the (genome id, marker id, sequence).
    marker_info_file : MarkerInfoFile
        A domain specific subclass of the marker info file.

    Returns
    -------
    Dict[str, str]
        dict[gid] = sequence
    """
    # Index the results.
    d_gid_marker = defaultdict(dict)
    for cur_result in list_seqs:
        for gid, marker_id, seq in cur_result:
            d_gid_marker[gid][marker_id] = seq

    # Create the alignment.
    out = defaultdict(list)
    for gid, cur_marker_d in d_gid_marker.items():
        for marker_id, marker_info in sorted(marker_info_file.markers.items()):
            out[gid].append(cur_marker_d.get(marker_id, '-' * marker_info['size']))
    return {k: ''.join(v) for k, v in out.items()}


def align_marker_set(gid_dict, marker_info_file: MarkerInfoFile, copy_number_file: CopyNumberFile, cpus):
    """Aligns the set of genomes for a specific domain.

    Parameters
    ----------
    gid_dict : dict
        A dictionary containing information about the genome, indexed by the id.
    marker_info_file : MarkerInfoFile
        A domain specific subclass of the marker info file.
    copy_number_file : CopyNumberFile
        A domain-specific subclass of the copy number file.
    cpus : int
        The maximum number of CPUs to use in subprocesses.

    Returns
    -------
    Dict[str, str]
        dict[gid] = sequence
    """
    logger = logging.getLogger('timestamp')

    logger.log(LOG_TASK, f'Generating concatenated alignment for each marker.')
    single_copy_hits = get_single_copy_hits(gid_dict, copy_number_file, cpus)

    with tempfile.TemporaryDirectory(prefix='gtdbtk_tmp_') as dir_tmp:
        # Write each of the markers to disk.
        marker_paths = dict()
        for marker_id, marker_d in single_copy_hits.items():
            cur_path = os.path.join(dir_tmp, f'{marker_id}.fa')
            marker_paths[marker_id] = cur_path
            with open(cur_path, 'w') as fh:
                for cur_gid, cur_seq in marker_d.items():
                    fh.write(f'>{cur_gid}\n{cur_seq}\n')

        # Run hmmalign on all of the markers (in order of largest)
        hmmer_v = HmmAligner.get_version()
        logger.log(LOG_TASK, f'Aligning {len(marker_paths)} identified markers using hmmalign {hmmer_v}.')
        queue = list()
        for marker_id, marker_path in sorted(marker_paths.items(),
                                             key=lambda z: -marker_info_file.markers[z[0]]['size']):
            queue.append((marker_id,
                          marker_info_file.markers[marker_id]['path'],
                          marker_path,
                          frozenset(single_copy_hits[marker_id])))
        with mp.Pool(processes=cpus) as pool:
            results = list(tqdm_log(pool.imap_unordered(run_hmm_align_worker, queue),
                                    total=len(queue), unit='marker'))

    # Create the concatenated alignment.
    return create_concat_alignment(results, marker_info_file)
