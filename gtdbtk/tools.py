import hashlib
import json
import logging
from typing import Optional

import math
import os
import random
import re
import time
import urllib.request
from itertools import islice

import dendropy
from tqdm import tqdm

import gtdbtk.config.config as Config
from gtdbtk.config.output import CHECKSUM_SUFFIX
from gtdbtk.exceptions import GTDBTkExit

order_rank = ["d__", "p__", "c__", "o__", 'f__', 'g__', 's__']
##################################################
############MISC UTILITIES########################
##################################################


def get_reference_ids():
    """Create a set of reference IDs using the config taxonomy file. This set
    contains the base id and the NCBI formatted ID (e.g. GB_GCA.. and  GCA_...)

    Returns
    -------
    frozenset
        An immutable set with short and long accessions (e.g. GB_GCA_ and GCA_).
    """
    results = set()
    with open(Config.TAXONOMY_FILE) as tf:
        for line in tf:
            raw_id = line.split('\t')[0]
            results.add(raw_id)
            if raw_id[0:4] in ['GCF_', 'GCA_']:
                results.add(add_ncbi_prefix(raw_id))
            elif raw_id[0:3] in ['RS_', 'GB_']:
                results.add(raw_id[3:])
    return frozenset(results)

def get_ref_genomes():
    """Returns a dictionary of genome accession to genome path.

    Returns
    -------
    dict[str, str]
        Dict[genome_id] = fasta_path
    """
    ref_genomes = dict()
    with open(Config.FASTANI_GENOME_LIST) as g_path_file:
        for line in g_path_file:
            (full_name, path) = line.strip().split()
            if full_name.endswith(Config.FASTANI_GENOMES_EXT):
                accession = full_name.split(Config.FASTANI_GENOMES_EXT)[0]
            ref_genomes[accession] = os.path.join(Config.FASTANI_DIR, path, full_name)
    return ref_genomes

def aa_percent_msa(aa_string):
        aa_len = sum([1 for c in aa_string if c.isalpha()])
        aa_perc = float(aa_len) / len(aa_string)
        return round(aa_perc * 100, 2)

def truncate_taxonomy(taxonomy, rank):
    taxonomy_list = taxonomy.split(';')
    taxonomy_list = taxonomy_list[0:order_rank.index(rank)]
    taxonomy = standardise_taxonomy(';'.join(taxonomy_list), 'bac120')
    return taxonomy

def limit_rank(taxa, rank_idx):
    return taxa[0:rank_idx+1] + order_rank[rank_idx+1:]

def standardise_taxonomy(taxstring, marker_set=None):
        """Create a 7 rank taxonomy string from an incomplete taxonomy string

        Parameters
        ----------
        taxstring : str
            incomplete taxonomy string
        marker_set : str
            The marker set to use.

        Returns
        -------
        string
            7 rank taxonomy string.
        """
        # return taxstring

        taxlist = taxstring.split(";")
        while '' in taxlist:
            taxlist.remove('')
        if marker_set == 'bac120':
            if not taxlist or taxlist[0] !='d__Bacteria' :
                taxlist.insert(0, 'd__Bacteria')
        if marker_set == 'ar53':
            if not taxlist or taxlist[0] !='d__Archaea' :
                taxlist.insert(0, 'd__Archaea')
        taxlist.extend(order_rank[len(taxlist):])
        new_taxstring = ";".join(taxlist)
        return new_taxstring


def add_ncbi_prefix(refname):
    if refname.startswith("GCF_"):
        return "RS_" + refname
    elif refname.startswith("GCA_"):
        return "GB_" + refname
    else:
        return refname


def splitchunks(d, n):
    chunksize = int(math.ceil(len(d) / float(n)))
    it = iter(d)
    for _ in range(0, len(d), chunksize):
        yield {k: d[k] for k in islice(it, chunksize)}


def splitchunks_list(l, n):
    """Yield successive n-sized chunks from l."""
    chunksize = int(math.ceil(len(l) / float(n)))
    for i in range(0, len(l), chunksize):
        yield l[i:i + chunksize]


def generateTempTableName():
    rng = random.SystemRandom()
    suffix = ''
    for _ in range(0, 10):
        suffix += rng.choice(
            'abcefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    return "TEMP" + suffix + str(int(time.time()))


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def confirm(msg):
    raw = input(msg + " (y/N): ")
    if raw.upper() in ["Y", "YES"]:
        return True
    return False


def sha1_dir(path, progress):
    """Recursively add files found within the path and output a SHA1 hash.

    Parameters
    ----------
    path : str
        The path to traverse.
    progress : bool
        True if progress should be displayed to stdout, False otherwise.

    Returns
    -------
    str
        The SHA1 hash representing this directory.
    """
    if progress:
        print('\r[{}]'.format(path), end='', flush=True)

    # Generate a queue of files to process
    queue = list()
    for root, dirs, files in os.walk(path):
        for file in files:
            path_file = os.path.join(root, file)
            queue.append(path_file)
    queue = sorted(queue)

    # Setup the hasher
    block_size = 65536
    hasher = hashlib.sha1()

    # Process the queue and obtain a single hash
    for idx, cur_path in enumerate(queue):

        if progress:
            print('\r[{}] - {}/{} files ({}%)'.format(path,
                                                      idx,
                                                      len(queue),
                                                      round(100 * (idx / len(queue)), 2)),
                  end='', flush=True)

        # Add the hash of the file
        with open(cur_path, 'rb') as fh:
            buf = fh.read(block_size)
            while len(buf) > 0:
                hasher.update(buf)
                buf = fh.read(block_size)

    if progress:
        print('\r' + ' ' * 100, end='', flush=True)
        print('\r', end='', flush=True)

    return hasher.hexdigest()


def sha256(input_file):
    """Determine SHA256 hash for file.

    Parameters
    ----------
    input_file : str
        Name of file.
    Returns
    -------
    str
        SHA256 hash.
    """
    block_size = 65536
    hasher = hashlib.sha1()
    with open(input_file, 'rb') as afile:
        buf = afile.read(block_size)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(block_size)
    return hasher.hexdigest()


def file_has_checksum(file_path, checksum_suffix=CHECKSUM_SUFFIX):
    """Check that the file contents match the checksum.

    Parameters
    ----------
    file_path : str
        Name of the file to check.
    checksum_suffix : str
        Suffix used to denote the file checksum.

    Returns
    -------
    bool
        True if the file has a checksum and it matches the original contents,
        False otherwise.

    """
    check_path = file_path + checksum_suffix
    if os.path.isfile(file_path) and os.path.isfile(check_path):
        with open(check_path, 'r') as check_f:
            return sha256(file_path) == check_f.read()
    return False


def symlink_f(src, dst, force=True):
    """Create a symbolic link pointing to src named dst.

    Parameters
    ----------
    src : str
        The source file.
    dst : str
        The destination file.
    force : bool
        Overwrite any file found with the same name as dst.

    """
    if force and os.path.isfile(dst):
        os.remove(dst)
    os.symlink(src, dst)


def get_memory_gb():
    try:
        with open('/proc/meminfo', 'r') as mem:
            hits = re.findall(r'([\w()]+):\s+(\d+) kB', mem.read())
        return {k: round(int(v) / 1e6, 2) for k, v in hits}
    except Exception:
        return None


def get_proc_memory_gb(pid):
    virt, res = None, None
    try:
        with open(f'/proc/{pid}/status', 'r') as fh:
            contents = fh.read()
        virt = int(re.search(r'VmSize:[^\d]+(\d+)', contents).group(1)) / 1e6
        res = int(re.search(r'VmRSS:[^\d]+(\d+)', contents).group(1)) / 1e6
    finally:
        return virt, res


def get_gtdbtk_latest_version():
    if not Config.GTDBTK_VER_CHECK:
        return None
    try:
        resp = json.loads(urllib.request.urlopen('https://pypi.org/pypi/gtdbtk/json',
                                                 timeout=Config.GTDBTK_VER_TIMEOUT).read().decode('utf-8'))
        return resp['info']['version']
    except Exception:
        return None


class TreeTraversal(object):
    """Efficiently calculates leaf nodes of a given node without re-computing
    any information.
    """

    def __init__(self):
        self.d_node_desc = dict()

    def get_leaf_nodes(self, node):
        """Efficiently return all leaf nodes under a given node.

        Parameters
        ----------
        node : dendropy.Node

        Returns
        -------
        frozenset
            The set of all leaf nodes under the given node.
        """

        # Leaf nodes will always return themselves.
        if node.is_leaf():
            self.d_node_desc[node] = frozenset({node})

        # Stop traversing down if the descendants are already known.
        if node in self.d_node_desc:
            return self.d_node_desc[node]

        # Descendants are not known, traverse down to find them.
        desc_nodes = set()
        for child_node in node.child_node_iter():

            # Already calculated, add it.
            if child_node in self.d_node_desc:
                desc_nodes = desc_nodes.union(self.d_node_desc[child_node])

            # Needs to be calculated, recurse.
            else:
                desc_nodes = desc_nodes.union(self.get_leaf_nodes(child_node))

        # Store the desc and exit.
        self.d_node_desc[node] = frozenset(desc_nodes)
        return self.d_node_desc[node]


def calculate_patristic_distance(qry_node, ref_nodes, tt=None):
    """Computes the patristic distance from the query node to all reference
    nodes. Note that all nodes must be a leaf nodes under max_node.

    Parameters
    ----------
    qry_node : dendropy.Node
        The query taxon node that the distance to all ref nodes will be found.
    ref_nodes : List[dendropy.Node]
        A list of reference nodes that the qry_node will be calculated to.
    tt : Optional[TreeTraversal]
        A TreeTraversal index, if absent a new one will be created.

    Returns
    -------
    Dict[dendropy.Node, float]
        A dictionary keyed by each reference taxon, valued by patristic dist.
    """
    tt = tt or TreeTraversal()

    # Iterate over each of the ref_nodes to find the MRCA to qry_node.
    d_ref_to_mrca = dict()
    for ref_node in ref_nodes:
        cur_dist_to_mrca = ref_node.edge_length

        # Go up the tree until the descendants include qry_node.
        parent_node = ref_node.parent_node
        while parent_node is not None:
            leaf_nodes = tt.get_leaf_nodes(parent_node)

            # Found the MRCA node.
            if qry_node in leaf_nodes:
                d_ref_to_mrca[ref_node] = (parent_node, cur_dist_to_mrca)
                break

            # Keep going up.
            cur_dist_to_mrca += parent_node.edge_length
            parent_node = parent_node.parent_node

        # If the loop did not break, raise an exception.
        else:
            raise GTDBTkExit(f'Unable to find MRCA: {qry_node.taxon.label} / '
                             f'{ref_node.taxon.label}')

    # Compute the distance from the qry_node to each of the MRCAs.
    out = dict()
    for ref_node, (mrca_node, ref_mrca_dist) in d_ref_to_mrca.items():

        # Go up the tree until the MRCA is found again.
        cur_dist_to_mrca = qry_node.edge_length
        cur_node = qry_node.parent_node
        while cur_node is not None:

            # Found the MRCA node.
            if cur_node == mrca_node:
                out[ref_node] = cur_dist_to_mrca + ref_mrca_dist
                break

            # Keep going up.
            cur_dist_to_mrca += cur_node.edge_length
            cur_node = cur_node.parent_node

        # Impossible case, but throw an exception anyway.
        else:
            raise GTDBTkExit(f'Tree is inconsistent: {qry_node.taxon.label} / '
                             f'{ref_node.taxon.label}')

    return out


class tqdm_log(object):
    """A basic wrapper for the tqdm progress bar. Automatically reports the
    runtime statistics after exit.
    """

    def __init__(self, iterable=None, **kwargs):
        # Setup reporting information.
        self.logger = logging.getLogger('timestamp')
        self.start_ts = None

        # Set default parameters.
        default = {'leave': False,
                   'smoothing': 0.1,
                   'bar_format': '==> Processed {n_fmt}/{total_fmt} {unit}s '
                                 '({percentage:.0f}%) |{bar:15}| [{rate_fmt}, ETA {remaining}]'}
        merged = {**default, **kwargs}
        self.args = merged

        # Instantiate tqdm
        self.tqdm = tqdm(iterable, **merged)

    def log(self):
        try:
            # Collect values.
            delta = time.time() - self.start_ts
            n = self.tqdm.n
            unit = self.args.get('unit', 'item')

            # Determine the scale.
            if delta > 60:
                time_unit = 'minute'
                scale = 1 / 60
            elif delta > 60 * 60:
                time_unit = 'hour'
                scale = 1 / (60 * 60)
            elif delta > 60 * 60 * 24:
                time_unit = 'day'
                scale = 1 / (60 * 60 * 24)
            else:
                time_unit = 'second'
                scale = 1

            # Scale units.
            value = delta * scale
            per = n / value

            # Determine what scale to use for the output.
            if per > 1:
                per_msg = f'{per:,.2f} {unit}s/{time_unit}'
            else:
                per_msg = f'{1 / per:,.2f} {time_unit}s/{unit}'

            # Output the message.
            s = 's' if n > 1 else ''
            msg = f'Completed {n:,} {unit}{s} in {value:,.2f} {time_unit}s ({per_msg}).'
            self.logger.info(msg)
        except Exception:
            pass

    def __enter__(self):
        self.start_ts = time.time()
        return self.tqdm

    def __iter__(self):
        try:
            self.start_ts = time.time()
            for item in self.tqdm:
                yield item
            self.log()
        finally:
            self.tqdm.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.log()
        self.tqdm.close()

    def __del__(self):
        self.tqdm.close()

def assert_outgroup_taxon_valid(outgroup_taxon: Optional[str]):
    """Check that the outgroup taxon is a valid string, no checks are made
    to ensure that it is a real taxon. If an invalid string is provided, an
    exception is raised."""
    if not outgroup_taxon:
        raise GTDBTkExit('The outgroup taxon cannot be empty, please specify a '
                         'valid taxon (e.g. p__Proteobacteria).')
    if len(outgroup_taxon) < 4:
        raise GTDBTkExit('The outgroup taxon seems suspiciously short and is '
                         'not a valid GTDB taxon, please specify a valid '
                         'taxon (e.g. p__Proteobacteria).')
    if outgroup_taxon.startswith('d__'):
        raise GTDBTkExit('You cannot specify a domain as an outgroup, please '
                         'select a lower rank (i.e. phylum, class, order, '
                         'family, genus, species).')
    if outgroup_taxon[:3] not in {'p__', 'c__', 'o__', 'f__', 'g__', 's__'}:
        raise GTDBTkExit('You must specify the outgroup taxon in '
                         'Greengenes-style format (e.g. p__Proteobacteria).')
    return
