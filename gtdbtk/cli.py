import argparse
import tempfile
from contextlib import contextmanager

from gtdbtk.biolib_lite.custom_help_formatter import ChangeTempAction
from gtdbtk.biolib_lite.custom_help_formatter import CustomHelpFormatter
from gtdbtk.config.common import CONFIG


@contextmanager
def subparser(parser, name, desc):
    yield parser.add_parser(name, conflict_handler='resolve', help=desc,
                            formatter_class=CustomHelpFormatter)


@contextmanager
def mutex_group(parser, required):
    group = parser.add_argument_group(f'mutually exclusive '
                                      f'{"required" if required else "optional"} '
                                      f'arguments')
    yield group.add_mutually_exclusive_group(required=required)


@contextmanager
def arg_group(parser, name):
    yield parser.add_argument_group(name)


def __temp_dir(group):
    group.add_argument('--tmpdir', action=ChangeTempAction, default=tempfile.gettempdir(),
                       help="specify alternative directory for temporary files")


def __genes(group):
    group.add_argument('--genes', action='store_true', default=False,
                       help='indicates input files contain predicted proteins as amino acids (skip gene calling).'
                            'Warning: This flag will skip the ANI comparison steps (ani_screen and classification).')


def __genome_dir(group):
    group.add_argument(
        '--genome_dir', help="directory containing genome files in FASTA format")


def __batchfile(group):
    group.add_argument('--batchfile', help="path to file describing genomes - tab "
                                           "separated in 2 or 3 columns (FASTA "
                                           "file, genome ID, translation "
                                           "table [optional])")


def __archaea(group):
    group.add_argument('--archaea', action='store_true', default=False,
                       help='process archaeal genomes')


def __bacteria(group):
    group.add_argument('--bacteria', action='store_true', default=False,
                       help='process bacterial genomes')


def __outgroup_taxon(group, required):
    group.add_argument('--outgroup_taxon', type=str, default=None, required=required,
                       help="taxon to use as outgroup (e.g., "
                            "``p__Patescibacteria`` or ``p__Altiarchaeota``)")


def __out_dir(group, required):
    group.add_argument('--out_dir', type=str, default=None, required=required,
                       help="directory to output files")


def __extension(group):
    group.add_argument('-x', '--extension', type=str, default='fna',
                       help='extension of files to process, ``gz`` = gzipped')

def __skip_ani_screen(group):
    group.add_argument('--skip_ani_screen', action="store_true", default=False,
                       help="Skip the ani_screening step to classify genomes using mash and skani.")

def __skip_gtdb_refs(group):
    group.add_argument('--skip_gtdb_refs', action="store_true", default=False,
                       help='do not include GTDB reference genomes in multiple sequence alignment.')


def __taxa_filter(group):
    group.add_argument('--taxa_filter', type=str, default=None,
                       help=('filter GTDB genomes to taxa (comma separated) within '
                             + 'specific taxonomic groups (e.g.: ``d__Bacteria`` '
                             + 'or ``p__Proteobacteria,p__Actinobacteria``)'))


def __min_perc_aa(group):
    group.add_argument('--min_perc_aa', type=float, default=10,
                       help='exclude genomes that do not have at least this '
                            'percentage of AA in the MSA (inclusive bound)')


def __custom_msa_filters(group):
    group.add_argument('--custom_msa_filters', action="store_true", default=False,
                       help=('perform custom filtering of MSA with ``cols_per_gene``, ``min_consensus`` '
                             + '``max_consensus``, and ``min_perc_taxa`` parameters instead of using canonical mask'))


def __cols_per_gene(group):
    group.add_argument('--cols_per_gene', type=int, default=42,
                       help='maximum number of columns to retain per gene '
                            'when generating the MSA')


def __min_consensus(group):
    group.add_argument('--min_consensus', type=float, default=25,
                       help='minimum percentage of the same amino acid '
                            'required to retain column (inclusive bound)')


def __max_consensus(group):
    group.add_argument('--max_consensus', type=float, default=95,
                       help='maximum percentage of the same amino acid '
                            'required to retain column (exclusive bound)')


def __min_perc_taxa(group):
    group.add_argument('--min_perc_taxa', type=float, default=50,
                       help='minimum percentage of taxa required to retain column (inclusive bound)')


def __prot_model(group):
    group.add_argument('--prot_model', choices=['JTT', 'WAG', 'LG'],
                       help='protein substitution model for tree inference', default='WAG')


def __rnd_seed(group):
    group.add_argument('--rnd_seed', type=int, default=None,
                       help='random seed to use for selecting columns, e.g. ``42``')


def __no_support(group):
    group.add_argument('--no_support', action="store_true", default=False,
                       help="do not compute local support values using the Shimodaira-Hasegawa test")


def __gamma(group):
    group.add_argument('--gamma', action="store_true", default=False,
                       help="rescale branch lengths to optimize the Gamma20 likelihood")


def __gtdbtk_classification_file(group):
    group.add_argument('--gtdbtk_classification_file', type=str, default=None,
                       help="file with GTDB-Tk classifications produced by the `classify` command")


def __custom_taxonomy_file(group):
    group.add_argument('--custom_taxonomy_file', type=str, default=None,
                       help="file indicating custom taxonomy strings for user "
                            "genomes, that should contain any genomes belonging to the outgroup. "
                            "Format: GENOME_ID<TAB>d__;p__;c__;o__;f__;g__;s__")


def __prefix(group):
    group.add_argument('--prefix', type=str, default='gtdbtk',
                       help='prefix for all output files')


def __cpus(group):
    group.add_argument('--cpus', default=1, type=int,
                       help='number of CPUs to use')


def __force(group):
    group.add_argument('--force', action='store_true', default=False,
                       help='continue processing if an error occurs on a single genome')


def __debug(group):
    group.add_argument('--debug', action="store_true", default=False,
                       help='create intermediate files for debugging purposes')

# This argument should be hidden from help
def __skip_pplacer(group):
    group.add_argument('--skip_pplacer', action="store_true", default=False,
                       help=argparse.SUPPRESS)


def __help(group):
    group.add_argument('-h', '--help', action="help", help="show help message")


def __pplacer_cpus(group):
    group.add_argument('--pplacer_cpus', type=int, default=None,
                       help='number of CPUs to use during pplacer placement')


def __scratch_dir(group):
    group.add_argument('--scratch_dir', type=str, default=None,
                       help='reduce pplacer memory usage by writing to disk (slower).')


# def __recalculate_red(group):
#     group.add_argument('-r', '--recalculate_red', default=False, action='store_true',
#                        help='recalculate RED values based on the reference tree and all added user genomes')


def __full_tree(group):
    group.add_argument('-f', '--full_tree', default=False, action='store_true',
                       help='use the unsplit bacterial tree for the classify step; this is the original GTDB-Tk '
                            f'approach (version < 2) and requires more than {CONFIG.PPLACER_MIN_RAM_BAC_FULL} GB of RAM to load the reference tree')


def __identify_dir(group, required):
    group.add_argument('--identify_dir', type=str, default=None, required=required,
                       help="output directory of 'identify' command")


def __skip_trimming(group):
    group.add_argument('--skip_trimming', default=False, action="store_true",
                       help='skip the trimming step and return the full MSAs')


def __msa_file(group, required):
    group.add_argument('--msa_file', type=str, default=None, required=required,
                       help="multiple sequence alignment in FASTA format")


def __align_dir(group, required):
    group.add_argument('--align_dir', type=str, default=None, required=required,
                       help="output directory of 'align' command")


def __input_tree(group, required):
    group.add_argument('--input_tree', type=str, default=None, required=required,
                       help="path to the unrooted tree in Newick format")


def __input_tree__rooted(group, required):
    group.add_argument('--input_tree', type=str, default=None, required=required,
                       help="rooted input tree with labelled ingroup taxon")


def __output_tree(group, required):
    group.add_argument('--output_tree', type=str, default=None, required=required,
                       help='path to output the tree')


def __ingroup_taxon(group, required):
    group.add_argument('--ingroup_taxon', type=str, default=None, required=required,
                       help="labelled ingroup taxon to use as root for "
                            "establishing RED values (e.g., c__Bacilli or f__Lactobacillaceae")


def __no_mash(group):
    group.add_argument('--no_mash', default=False, action='store_true',
                       help='skip pre-filtering of genomes using Mash')


def __mash_k(group):
    group.add_argument('--mash_k', default=CONFIG.MASH_K_VALUE, type=int,
                       help='k-mer size [1-32]')


def __mash_s(group):
    group.add_argument('--mash_s', default=CONFIG.MASH_S_VALUE, type=int,
                       help='maximum number of non-redundant hashes')


def __mash_d(group):
    group.add_argument('--mash_d', default=CONFIG.MASH_D_VALUE, type=float,
                       help='maximum distance to keep [0-1]')


def __mash_v(group):
    group.add_argument('--mash_v', default=CONFIG.MASH_V_VALUE, type=float,
                       help='maximum p-value to keep [0-1]')

def __mash_max_distance(group):
    group.add_argument('--mash_max_distance', default=CONFIG.MASH_MAX_DISTANCE, type=float,
                       help='Maximum Mash distance to select a potential GTDB genome as representative '
                            'of a user genome.')

def __mash_db(group):
    group.add_argument('--mash_db', default=None, type=str,
                       help='path to save/read (if exists) the Mash reference sketch database (.msh)')


def __min_af(group):
    group.add_argument('--min_af', type=float, default=CONFIG.AF_THRESHOLD,
                       help='minimum alignment fraction to assign genome to a species cluster')


def __untrimmed_msa(group, required):
    group.add_argument('--untrimmed_msa', type=str, default=None, required=required,
                       help="path to the untrimmed MSA file")


def __keep_intermediates(group):
    group.add_argument('--keep_intermediates', default=False, action='store_true',
                       help='keep intermediate files in the final directory')


def __output(group, required):
    group.add_argument('--output', type=str, default=None, required=required,
                       help='output file')


def __mask_file(group):
    group.add_argument('--mask_file', type=str, default=None,
                       help="path to a custom mask file for trimming the MSA")


def __reference_mask(group):
    group.add_argument('--reference_mask',
                       choices=['arc', 'bac'],
                       help="reference mask already present in GTDB-Tk")


def __domain(group, required):
    group.add_argument('--domain', required=required, choices=['arc', 'bac'],
                       help="domain to export")

def __all_ranks(group):
    group.add_argument('--all_ranks', default=False, action='store_true',
                       help='add all missing ranks to the leaf nodes if they are present in the reference tree.')

def __db_version(group):
    group.add_argument('--db_version', type = int, default = None,
                       help="GTDB-Tk version package to test for compatibility.")


def __write_single_copy_genes(group):
    group.add_argument('--write_single_copy_genes', default=False, action='store_true',
                       help='output unaligned single-copy marker genes')

def get_main_parser():
    # Setup the main, and sub parsers.
    main_parser = argparse.ArgumentParser(
        prog='gtdbtk', add_help=False, conflict_handler='resolve')
    sub_parsers = main_parser.add_subparsers(help="--", dest='subparser_name')

    # de novo workflow.
    with subparser(sub_parsers, 'de_novo_wf', 'Infer de novo tree and decorate with GTDB taxonomy.') as parser:
        with mutex_group(parser, required=True) as grp:
            __genome_dir(grp)
            __batchfile(grp)
        with mutex_group(parser, required=True) as grp:
            __bacteria(grp)
            __archaea(grp)
        with arg_group(parser, 'required named arguments') as grp:
            __outgroup_taxon(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __extension(grp)
            __skip_gtdb_refs(grp)
            __taxa_filter(grp)
            __min_perc_aa(grp)
            __custom_msa_filters(grp)
            __cols_per_gene(grp)
            __min_consensus(grp)
            __max_consensus(grp)
            __min_perc_taxa(grp)
            __rnd_seed(grp)
            __prot_model(grp)
            __no_support(grp)
            __gamma(grp)
            __gtdbtk_classification_file(grp)
            __custom_taxonomy_file(grp)
            __write_single_copy_genes(grp)
            __prefix(grp)
            __genes(grp)
            __cpus(grp)
            __force(grp)
            __temp_dir(grp)
            __keep_intermediates(grp)
            __debug(grp)
            __help(grp)

    # Classify workflow.
    with subparser(sub_parsers, 'classify_wf', 'Classify genomes by placement in GTDB reference tree.') as parser:
        with mutex_group(parser, required=True) as grp:
            __genome_dir(grp)
            __batchfile(grp)
        with arg_group(parser, 'required named arguments') as grp:
            __out_dir(grp, required=True)
        with mutex_group(parser, required=True) as grp:
            __skip_ani_screen(grp)
            __mash_db(grp)
        with arg_group(parser, 'optional Mash arguments') as grp:
            __no_mash(grp)
            __mash_k(grp)
            __mash_s(grp)
            __mash_v(grp)
            __mash_max_distance(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __full_tree(grp)
            __extension(grp)
            __min_perc_aa(grp)
            __prefix(grp)
            __genes(grp)
            __cpus(grp)
            __pplacer_cpus(grp)
            __force(grp)
            __scratch_dir(grp)
            __write_single_copy_genes(grp)
            __keep_intermediates(grp)
            __min_af(grp)
            __temp_dir(grp)
            __debug(grp)
            __skip_pplacer(grp)
            __help(grp)

    # Identify marker genes in genomes.
    with subparser(sub_parsers, 'identify', 'Identify marker genes in genomes.') as parser:
        with mutex_group(parser, required=True) as grp:
            __genome_dir(grp)
            __batchfile(grp)
        with arg_group(parser, 'required named arguments') as grp:
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __extension(grp)
            __prefix(grp)
            __genes(grp)
            __cpus(grp)
            __force(grp)
            __write_single_copy_genes(grp)
            __temp_dir(grp)
            __debug(grp)
            __help(grp)

    # Create a multiple sequence alignment.
    with subparser(sub_parsers, 'align', 'Create multiple sequence alignment.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __identify_dir(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __skip_gtdb_refs(grp)
            __taxa_filter(grp)
            __min_perc_aa(grp)
            __cols_per_gene(grp)
            __min_consensus(grp)
            __max_consensus(grp)
            __min_perc_taxa(grp)
            __rnd_seed(grp)
            __prefix(grp)
            __cpus(grp)
            __temp_dir(grp)
            __debug(grp)
            __help(grp)
        with mutex_group(parser, required=False) as grp:
            __custom_msa_filters(grp)
            __skip_trimming(grp)

    # Infer a de novo tree.
    with subparser(sub_parsers, 'infer', 'Infer tree from multiple sequence alignment.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __msa_file(grp, required=True)
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __prot_model(grp)
            __no_support(grp)
            __gamma(grp)
            __prefix(grp)
            __cpus(grp)
            __temp_dir(grp)
            __debug(grp)
            __help(grp)

    # Classify genomes via placement with pplacer.
    with subparser(sub_parsers, 'classify', 'Determine taxonomic classification of genomes.') as parser:
        with mutex_group(parser, required=True) as grp:
            __genome_dir(grp)
            __batchfile(grp)
        with arg_group(parser, 'required named arguments') as grp:
            __align_dir(grp, required=True)
            __out_dir(grp, required=True)
        with mutex_group(parser, required=True) as grp:
            __skip_ani_screen(grp)
            __mash_db(grp)
        with arg_group(parser, 'optional Mash arguments') as grp:
            __no_mash(grp)
            __mash_k(grp)
            __mash_s(grp)
            __mash_v(grp)
            __mash_max_distance(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __extension(grp)
            __prefix(grp)
            __cpus(grp)
            __pplacer_cpus(grp)
            __scratch_dir(grp)
            __genes(grp)
            __full_tree(grp)
            __min_af(grp)
            __temp_dir(grp)
            __debug(grp)
            __skip_pplacer(grp)
            __help(grp)

    # Root a tree using an outgroup.
    with subparser(sub_parsers, 'root', 'Root tree using an outgroup.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __input_tree(grp, required=True)
            __outgroup_taxon(grp, required=True)
            __output_tree(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __gtdbtk_classification_file(grp)
            __custom_taxonomy_file(grp)
            __temp_dir(grp)
            __debug(grp)
            __help(grp)

    # Decorate a tree.
    with subparser(sub_parsers, 'decorate', 'Decorate tree with GTDB taxonomy') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __input_tree(grp, required=True)
            __output_tree(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __gtdbtk_classification_file(grp)
            __custom_taxonomy_file(grp)
            __temp_dir(grp)
            __debug(grp)
            __help(grp)

    # Establish taxonomic ranks of internal nodes using RED
    with subparser(sub_parsers, 'infer_ranks', 'Establish taxonomic ranks of internal nodes using RED.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __input_tree__rooted(grp, required=True)
            __ingroup_taxon(grp, required=True)
            __output_tree(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __temp_dir(grp)
            __debug(grp)
            __help(grp)

    # Calculate ANI to representative genomes.
    with subparser(sub_parsers, 'ani_rep', 'Calculate ANI to GTDB representative genomes.') as parser:
        with mutex_group(parser, required=True) as grp:
            __genome_dir(grp)
            __batchfile(grp)
        with arg_group(parser, 'required named arguments') as grp:
            __out_dir(grp, required=True)
        with arg_group(parser, 'optional Mash arguments') as grp:
            __no_mash(grp)
            __mash_k(grp)
            __mash_s(grp)
            __mash_d(grp)
            __mash_v(grp)
            __mash_db(grp)
        with arg_group(parser, 'optional skani arguments') as grp:
            __min_af(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __extension(grp)
            __prefix(grp)
            __cpus(grp)
            __temp_dir(grp)
            __debug(grp)
            __help(grp)

    # Run a test.
    with subparser(sub_parsers, 'test', 'Test the classify_wf pipeline with 3 archaeal genomes.') as parser:
        with arg_group(parser, 'optional arguments') as grp:
            __out_dir(grp, required=False)
            __cpus(grp)
            __debug(grp)
            __help(grp)

    # Trim MSA.
    with subparser(sub_parsers, 'trim_msa', 'Trim an untrimmed MSA file based on a mask.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __untrimmed_msa(grp, required=True)
            __output(grp, required=True)
        with mutex_group(parser, required=True) as grp:
            __mask_file(grp)
            __reference_mask(grp)
        with arg_group(parser, 'optional arguments') as grp:
            __debug(grp)
            __help(grp)

    # Remove labels
    with subparser(sub_parsers, 'remove_labels', 'Remove labels (bootstrap values, node labels) from an Newick tree to '
                                                 'to improve compatibility with tree viewers.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __input_tree(grp, required=True)
            __output_tree(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __debug(grp)
            __help(grp)

    # Remove labels
    with subparser(sub_parsers, 'convert_to_itol', 'Reformat the GTDB-Tk tree to be iTOL compatible.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __input_tree(grp, required=True)
            __output_tree(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __debug(grp)
            __help(grp)

    # Convert genome ids to species names.
    with subparser(sub_parsers, 'convert_to_species', 'Replace GTDB genomes ids with GTDB Species name.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __input_tree(grp, required=True)
            __output_tree(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __custom_taxonomy_file(grp)
            __all_ranks(grp)
            __debug(grp)
            __help(grp)

    # Export MSA.
    with subparser(sub_parsers, 'export_msa', 'Export the untrimmed archaeal or bacterial MSA file.') as parser:
        with arg_group(parser, 'required named arguments') as grp:
            __domain(grp, required=True)
            __output(grp, required=True)
        with arg_group(parser, 'optional arguments') as grp:
            __debug(grp)
            __help(grp)

    # Verify install.
    with subparser(sub_parsers, 'check_install', 'Verify third party programs and '
                                                 'GTDB reference package.') as parser:
        with arg_group(parser, 'optional arguments') as grp:
            __db_version(grp)
            __debug(grp)
            __help(grp)

    return main_parser
