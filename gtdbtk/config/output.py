from os.path import join

# Command: identify
DIR_IDENTIFY = 'identify'
DIR_IDENTIFY_INTERMEDIATE = join(DIR_IDENTIFY, 'intermediate_results')
DIR_MARKER_GENE = join(DIR_IDENTIFY_INTERMEDIATE, 'marker_genes')
DIR_IDENTIFY_FASTA = join(DIR_IDENTIFY_INTERMEDIATE, 'single_copy_fasta')
PATH_BAC120_MARKER_SUMMARY = join(DIR_IDENTIFY, '{prefix}.bac120.markers_summary.tsv')
PATH_AR122_MARKER_SUMMARY = join(DIR_IDENTIFY, '{prefix}.ar122.markers_summary.tsv')
PATH_TLN_TABLE_SUMMARY = join(DIR_IDENTIFY, '{prefix}.translation_table_summary.tsv')

# Command: identify -> marker genes
GENOME_FILE_SUFFIX = "_genomic.fna"
PROTEIN_FILE_SUFFIX = "_protein.faa"
NT_GENE_FILE_SUFFIX = "_protein.fna"
GFF_FILE_SUFFIX = "_protein.gff"
TRANSLATION_TABLE_SUFFIX = "_translation_table.tsv"
CHECKSUM_SUFFIX = ".sha256"
TIGRFAM_SUFFIX = "_tigrfam.tsv"
TIGRFAM_TOP_HIT_SUFFIX = "_tigrfam_tophit.tsv"
PFAM_SUFFIX = "_pfam.tsv"
PFAM_TOP_HIT_SUFFIX = "_pfam_tophit.tsv"

# Command: align
DIR_ALIGN = 'align'
DIR_ALIGN_INTERMEDIATE = join(DIR_ALIGN, 'intermediate_results')
PATH_BAC120_FILTERED_GENOMES = join(DIR_ALIGN, '{prefix}.bac120.filtered.tsv')
PATH_AR122_FILTERED_GENOMES = join(DIR_ALIGN, '{prefix}.ar122.filtered.tsv')
PATH_BAC120_MSA = join(DIR_ALIGN, '{prefix}.bac120.msa.fasta')
PATH_AR122_MSA = join(DIR_ALIGN, '{prefix}.ar122.msa.fasta')
PATH_BAC120_USER_MSA = join(DIR_ALIGN, '{prefix}.bac120.user_msa.fasta')
PATH_AR122_USER_MSA = join(DIR_ALIGN, '{prefix}.ar122.user_msa.fasta')
PATH_BAC120_MARKER_INFO = join(DIR_ALIGN_INTERMEDIATE, '{prefix}.bac120.marker_info.tsv')
PATH_AR122_MARKER_INFO = join(DIR_ALIGN_INTERMEDIATE, '{prefix}.ar122.marker_info.tsv')
DIR_ALIGN_MARKERS = join(DIR_ALIGN_INTERMEDIATE, 'markers')

# Command: classify
DIR_CLASSIFY = 'classify'
PATH_BAC120_TREE_FILE = join(DIR_CLASSIFY, '{prefix}.bac120.classify.tree')
PATH_AR122_TREE_FILE = join(DIR_CLASSIFY, '{prefix}.ar122.classify.tree')
PATH_BAC120_SUMMARY_OUT = join(DIR_CLASSIFY, '{prefix}.bac120.summary.tsv')
PATH_AR122_SUMMARY_OUT = join(DIR_CLASSIFY, '{prefix}.ar122.summary.tsv')
PATH_HIGH_BAC120_TREE_FILE = join(DIR_CLASSIFY, '{prefix}.high.bac120.classify.tree')
PATH_LOW_BAC120_TREE_FILE = join(DIR_CLASSIFY, '{prefix}.bac120.classify.tree.{iter}.tree')
PATH_BAC120_CONFLICT = join(DIR_CLASSIFY, '{prefix}.bac120.conflict.tsv')

DIR_CLASSIFY_INTERMEDIATE = join(DIR_CLASSIFY, 'intermediate_results')
PATH_BAC120_RED_DICT = join(DIR_CLASSIFY_INTERMEDIATE, '{prefix}.bac120.red_dictionary.tsv')
PATH_AR122_RED_DICT = join(DIR_CLASSIFY_INTERMEDIATE, '{prefix}.ar122.red_dictionary.tsv')
PATH_BAC120_PPLACER_CLASS = join(DIR_CLASSIFY_INTERMEDIATE, '{prefix}.bac120.classification_pplacer.tsv')
PATH_BAC120_HIGH_PPLACER_CLASS = join(DIR_CLASSIFY_INTERMEDIATE, '{prefix}.bac120.high.classification_pplacer.tsv')
PATH_BAC120_LOW_PPLACER_CLASS = join(DIR_CLASSIFY_INTERMEDIATE, '{prefix}.bac120.low.classification_pplacer_tree_{iter}.tsv')


PATH_AR122_PPLACER_CLASS = join(DIR_CLASSIFY_INTERMEDIATE, '{prefix}.ar122.classification_pplacer.tsv')

DIR_PPLACER = join(DIR_CLASSIFY_INTERMEDIATE, 'pplacer')
PATH_BAC120_PPLACER_OUT = join(DIR_PPLACER, 'pplacer.bac120.out')
PATH_AR122_PPLACER_OUT = join(DIR_PPLACER, 'pplacer.ar122.out')
PATH_BAC120_PPLACER_JSON = join(DIR_PPLACER, 'pplacer.bac120.json')
PATH_AR122_PPLACER_JSON = join(DIR_PPLACER, 'pplacer.ar122.json')

# SPLIT TREE
PATH_HIGH_BAC120_PPLACER_OUT = join(DIR_PPLACER, 'pplacer.high.bac120.out')
PATH_HIGH_BAC120_PPLACER_JSON = join(DIR_PPLACER, 'pplacer.high.bac120.json')
DIR_LOW_PPLACER = join(DIR_PPLACER, 'tree_{iter}')
PATH_LOW_BAC120_SUBMSA = join(DIR_LOW_PPLACER, 'user_msa_file.fasta')
PATH_LOW_BAC120_PPLACER_OUT = join(DIR_LOW_PPLACER, 'pplacer.low.bac120.out')
PATH_LOW_BAC120_PPLACER_JSON = join(DIR_LOW_PPLACER, 'pplacer.low.bac120.json')

# Command: infer
DIR_INFER = 'infer'
DIR_INFER_INTERMEDIATE = join(DIR_INFER, 'intermediate_results')
PATH_MARKER_UNROOTED_TREE = join(DIR_INFER_INTERMEDIATE, '{prefix}.{marker}.unrooted.tree')
PATH_UNROOTED_TREE = join(DIR_INFER_INTERMEDIATE, '{prefix}.unrooted.tree')
PATH_MARKER_TREE_LOG = join(DIR_INFER_INTERMEDIATE, '{prefix}.{marker}.tree.log')
PATH_TREE_LOG = join(DIR_INFER_INTERMEDIATE, '{prefix}.tree.log')
PATH_MARKER_FASTTREE_LOG = join(DIR_INFER_INTERMEDIATE, '{prefix}.{marker}.fasttree.log')
PATH_FASTTREE_LOG = join(DIR_INFER_INTERMEDIATE, '{prefix}.fasttree.log')

PATH_BAC120_UNROOTED_TREE = join(DIR_INFER_INTERMEDIATE, '{prefix}.bac120.unrooted.tree')
PATH_AR122_UNROOTED_TREE = join(DIR_INFER_INTERMEDIATE, '{prefix}.ar122.unrooted.tree')
PATH_BAC120_ROOTED_TREE = join(DIR_INFER_INTERMEDIATE, '{prefix}.bac120.rooted.tree')
PATH_AR122_ROOTED_TREE = join(DIR_INFER_INTERMEDIATE, '{prefix}.ar122.rooted.tree')
PATH_BAC120_DECORATED_TREE = join(DIR_INFER, '{prefix}.bac120.decorated.tree')
PATH_AR122_DECORATED_TREE = join(DIR_INFER, '{prefix}.ar122.decorated.tree')

# Command: ani_rep
DIR_ANI_REP_INT = join('intermediate_results')
DIR_ANI_REP_INT_MASH = join(DIR_ANI_REP_INT, 'mash')

# General files
PATH_WARNINGS = '{prefix}.warnings.log'
