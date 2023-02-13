import json
import os
import sys

"""
Load the reference package. This will fail if the directory doesn't exist.
"""
try:
    GENERIC_PATH = os.environ['GTDBTK_DATA_PATH']
except KeyError:
    print('\n' + '=' * 80)
    print(' ERROR '.center(80))
    print('_' * 80 + '\n')
    print("The 'GTDBTK_DATA_PATH' environment variable is not defined.".center(80) + '\n')
    print('Please set this variable to your reference data package.'.center(80))
    print('https://github.com/Ecogenomics/GTDBTk#installation'.center(80))
    print('=' * 80)
    sys.exit(1)

"""
If the reference package sub-folders still exist in GTDBTK_DATA_PATH, then there
is no need to edit the variables below.
"""
MIN_REF_DATA_VERSION = 'r207'

MSA_FOLDER = os.path.join(GENERIC_PATH, "msa/")
MASK_DIR = os.path.join(GENERIC_PATH, "masks/")
PPLACER_DIR = os.path.join(GENERIC_PATH, "pplacer/")
FASTANI_DIR = os.path.join(GENERIC_PATH, "fastani/")
MASH_DIR = os.path.join(GENERIC_PATH, "mash/")

TAX_FOLDER = os.path.join(GENERIC_PATH, "taxonomy/")
RADII_DIR = os.path.join(GENERIC_PATH, "radii/")
METADATA_DIR = os.path.join(GENERIC_PATH, "metadata/")
RED_DIR = os.path.join(GENERIC_PATH, "mrca_red/")
MARKER_DIR = os.path.join(GENERIC_PATH, 'markers/')
TIGRFAM_HMMS = os.path.join(MARKER_DIR, 'tigrfam/tigrfam.hmm')
PFAM_HMM_DIR = os.path.join(MARKER_DIR, 'pfam/')

SPLIT_DIR = os.path.join(GENERIC_PATH, 'split')
BACKBONE_SPLIT_DIR = os.path.join(SPLIT_DIR, 'backbone')
CLASS_LEVEL_SPLIT_DIR = os.path.join(SPLIT_DIR, 'class_level')
BACKBONE_PPLACER_DIR = os.path.join(BACKBONE_SPLIT_DIR, 'pplacer')
CLASS_LEVEL_PPLACER_DIR = os.path.join(CLASS_LEVEL_SPLIT_DIR, 'pplacer')
BACKBONE_RED_DIR = os.path.join(BACKBONE_SPLIT_DIR, 'red')
CLASS_LEVEL_RED_DIR = os.path.join(CLASS_LEVEL_SPLIT_DIR, 'red')

CLASS_LEVEL_TREE_MAPPING_FILE = os.path.join(CLASS_LEVEL_SPLIT_DIR, 'tree_mapping.tsv')

BACKBONE_PPLACER_REF_PKG = 'gtdbtk_package_backbone.refpkg'
BACKBONE_RED_FILE = os.path.join(BACKBONE_RED_DIR, 'backbone_red_value.tsv')
CLASS_LEVEL_PPLACER_REF_PKG = 'gtdbtk.package.{iter}.refpkg'
CLASS_LEVEL_RED_FILE = os.path.join(CLASS_LEVEL_RED_DIR, 'red_value_{iter}.tsv')

RED_DIST_BAC_DICT = ''
RED_DIST_ARC_DICT = ''
VERSION_DATA = ''
try:
    with open(os.path.join(METADATA_DIR, "metadata.txt")) as metadataData:
        for line in metadataData:
            try:
                line_infos = line.strip().split('=')
                if line_infos[0] == 'RED_DIST_BAC_DICT':
                    RED_DIST_BAC_DICT = json.loads(line_infos[1])
                elif line_infos[0] == 'RED_DIST_ARC_DICT':
                    RED_DIST_ARC_DICT = json.loads(line_infos[1])
                elif line_infos[0] == 'VERSION_DATA':
                    VERSION_DATA = line_infos[1]
            except ValueError:
                print(f"Skipping invalid line {repr(line)}")
except IOError:
    print('\n' + '=' * 80)
    print(' ERROR '.center(80))
    print('_' * 80 + '\n')
    print('The GTDB-Tk reference data does not exist or is corrupted.'.center(80))
    print(('GTDBTK_DATA_PATH=%s' % GENERIC_PATH).center(80) + '\n')
    print('Please compare the checksum to those provided in the download repository.'.center(80))
    print('https://github.com/Ecogenomics/GTDBTk#gtdb-tk-reference-data'.center(80))
    print('=' * 80)
    sys.exit(1)

# Relative Evolution Distance
RED_INTERVAL = 0.1
RED_MIN_SUPPORT = 0.0
RED_MIN_CHILDREN = 2

# Marker information
BAC120_MARKERS = {"PFAM": ["PF00380.20.hmm", "PF00410.20.hmm", "PF00466.21.hmm",
                           "PF01025.20.hmm", "PF02576.18.hmm", "PF03726.15.hmm"],
                  "TIGRFAM": ["TIGR00006.HMM", "TIGR00019.HMM", "TIGR00020.HMM",
                              "TIGR00029.HMM", "TIGR00043.HMM", "TIGR00054.HMM",
                              "TIGR00059.HMM", "TIGR00061.HMM", "TIGR00064.HMM",
                              "TIGR00065.HMM", "TIGR00082.HMM", "TIGR00083.HMM",
                              "TIGR00084.HMM", "TIGR00086.HMM", "TIGR00088.HMM",
                              "TIGR00090.HMM", "TIGR00092.HMM", "TIGR00095.HMM",
                              "TIGR00115.HMM", "TIGR00116.HMM", "TIGR00138.HMM",
                              "TIGR00158.HMM", "TIGR00166.HMM", "TIGR00168.HMM",
                              "TIGR00186.HMM", "TIGR00194.HMM", "TIGR00250.HMM",
                              "TIGR00337.HMM", "TIGR00344.HMM", "TIGR00362.HMM",
                              "TIGR00382.HMM", "TIGR00392.HMM", "TIGR00396.HMM",
                              "TIGR00398.HMM", "TIGR00414.HMM", "TIGR00416.HMM",
                              "TIGR00420.HMM", "TIGR00431.HMM", "TIGR00435.HMM",
                              "TIGR00436.HMM", "TIGR00442.HMM", "TIGR00445.HMM",
                              "TIGR00456.HMM", "TIGR00459.HMM", "TIGR00460.HMM",
                              "TIGR00468.HMM", "TIGR00472.HMM", "TIGR00487.HMM",
                              "TIGR00496.HMM", "TIGR00539.HMM", "TIGR00580.HMM",
                              "TIGR00593.HMM", "TIGR00615.HMM", "TIGR00631.HMM",
                              "TIGR00634.HMM", "TIGR00635.HMM", "TIGR00643.HMM",
                              "TIGR00663.HMM", "TIGR00717.HMM", "TIGR00755.HMM",
                              "TIGR00810.HMM", "TIGR00922.HMM", "TIGR00928.HMM",
                              "TIGR00959.HMM", "TIGR00963.HMM", "TIGR00964.HMM",
                              "TIGR00967.HMM", "TIGR01009.HMM", "TIGR01011.HMM",
                              "TIGR01017.HMM", "TIGR01021.HMM", "TIGR01029.HMM",
                              "TIGR01032.HMM", "TIGR01039.HMM", "TIGR01044.HMM",
                              "TIGR01059.HMM", "TIGR01063.HMM", "TIGR01066.HMM",
                              "TIGR01071.HMM", "TIGR01079.HMM", "TIGR01082.HMM",
                              "TIGR01087.HMM", "TIGR01128.HMM", "TIGR01146.HMM",
                              "TIGR01164.HMM", "TIGR01169.HMM", "TIGR01171.HMM",
                              "TIGR01302.HMM", "TIGR01391.HMM", "TIGR01393.HMM",
                              "TIGR01394.HMM", "TIGR01510.HMM", "TIGR01632.HMM",
                              "TIGR01951.HMM", "TIGR01953.HMM", "TIGR02012.HMM",
                              "TIGR02013.HMM", "TIGR02027.HMM", "TIGR02075.HMM",
                              "TIGR02191.HMM", "TIGR02273.HMM", "TIGR02350.HMM",
                              "TIGR02386.HMM", "TIGR02397.HMM", "TIGR02432.HMM",
                              "TIGR02729.HMM", "TIGR03263.HMM", "TIGR03594.HMM",
                              "TIGR03625.HMM", "TIGR03632.HMM", "TIGR03654.HMM",
                              "TIGR03723.HMM", "TIGR03725.HMM", "TIGR03953.HMM"]}



#New Version of AR53_MARKERS
AR53_MARKERS = {"PFAM": ["PF04919.13.hmm","PF07541.13.hmm","PF01000.27.hmm",
                        "PF00687.22.hmm","PF00466.21.hmm","PF00827.18.hmm","PF01280.21.hmm","PF01090.20.hmm",
                        "PF01200.19.hmm","PF01015.19.hmm","PF00900.21.hmm","PF00410.20.hmm"],
                "TIGRFAM":["TIGR00037.HMM","TIGR00064.HMM","TIGR00111.HMM",
                            "TIGR00134.HMM","TIGR00279.HMM","TIGR00291.HMM","TIGR00323.HMM",
                            "TIGR00335.HMM","TIGR00373.HMM","TIGR00405.HMM","TIGR00448.HMM",
                            "TIGR00483.HMM","TIGR00491.HMM","TIGR00522.HMM","TIGR00967.HMM",
                            "TIGR00982.HMM","TIGR01008.HMM","TIGR01012.HMM","TIGR01018.HMM",
                            "TIGR01020.HMM","TIGR01028.HMM","TIGR01046.HMM","TIGR01052.HMM",
                            "TIGR01171.HMM","TIGR01213.HMM","TIGR01952.HMM","TIGR02236.HMM",
                            "TIGR02338.HMM","TIGR02389.HMM","TIGR02390.HMM","TIGR03626.HMM",
                            "TIGR03627.HMM","TIGR03628.HMM","TIGR03629.HMM","TIGR03670.HMM",
                            "TIGR03671.HMM","TIGR03672.HMM","TIGR03673.HMM","TIGR03674.HMM",
                            "TIGR03676.HMM","TIGR03680.HMM"]}


# Information for Multiple hits markers:
DEFAULT_MULTIHIT_THRESHOLD = 10.0

# Information for aligning genomes
DEFAULT_DOMAIN_THRESHOLD = 10.0
AR_MARKER_COUNT = 53
BAC_MARKER_COUNT = 120

# Information about alignment Fraction to resolve fastANI results
AF_THRESHOLD = 0.5

# MSA file names
CONCAT_BAC120 = os.path.join(MSA_FOLDER, f"gtdb_{VERSION_DATA}_bac120.faa")
CONCAT_AR53 = os.path.join(MSA_FOLDER, f"gtdb_{VERSION_DATA}_ar53.faa")

# Taxonomy file name
TAXONOMY_FILE = os.path.join(TAX_FOLDER, "gtdb_taxonomy.tsv")

# Type Strain radii file
RADII_FILE = os.path.join(RADII_DIR, "gtdb_radii.tsv")

# Mask file names
MASK_BAC120 = f"gtdb_{VERSION_DATA}_bac120.mask"
MASK_AR53 = f"gtdb_{VERSION_DATA}_ar53.mask"
MASK_RPS23 = f"gtdb_{VERSION_DATA}_rps23.mask"

# Pplacer configuration
PPLACER_BAC120_REF_PKG = f"gtdb_{VERSION_DATA}_bac120.refpkg"
PPLACER_AR53_REF_PKG = f"gtdb_{VERSION_DATA}_ar53.refpkg"
PPLACER_RPS23_REF_PKG = f"gtdb_{VERSION_DATA}_rps23.refpkg"
PPLACER_MIN_RAM_BAC_FULL = 320
PPLACER_MIN_RAM_BAC_SPLIT = 55
PPLACER_MIN_RAM_ARC = 40

# Fastani configuration
FASTANI_SPECIES_THRESHOLD = 95.0
FASTANI_GENOMES = os.path.join(FASTANI_DIR, "database/")
FASTANI_GENOME_LIST = os.path.join(FASTANI_DIR, "genome_paths.tsv")
FASTANI_GENOMES_EXT = "_genomic.fna.gz"

# Mash configuration
MASH_SKETCH_FILE = 'gtdb_ref_sketch.msh'
MASH_K_VALUE = 16
MASH_S_VALUE = 5000
MASH_MAX_DISTANCE = 0.1
MASH_D_VALUE = MASH_MAX_DISTANCE
MASH_V_VALUE = 1.0


# MRCA RED VALUE
MRCA_RED_BAC120 = os.path.join(RED_DIR, f"gtdbtk_{VERSION_DATA}_bac120.tsv")
MRCA_RED_AR53 = os.path.join(RED_DIR, f"gtdbtk_{VERSION_DATA}_ar53.tsv")

# Hashing information for validating the reference package.
REF_HASHES = {PPLACER_DIR: '20903925a856a58b102a7b0ce160c5cbd2cf675b',
              MASK_DIR: '50e414a9de18170e8cb97f990f89ff60a0fe29d5',
              MARKER_DIR: '163f542c3f0a40f59df45d453aa235b39aa96e27',
              RADII_DIR: '8fd13b1c5d7a7b073ba96fb628581613b293a374',
              MSA_FOLDER: '24f250d7cf0eb0bc65dccd2f3c9247e553ea322f',
              METADATA_DIR: '9772fbeac1311b31e10293fa610eb33aa1ec8e15',
              TAX_FOLDER: '6fb0233b05633242369b40c026fd1ee53e266afa',
              FASTANI_DIR: '973c456c02f55bb82908a6811c7076e207e9b206',
              RED_DIR: '7b8b67b3157204b470c9eb809d3c39c4effffabc'}

# Config values for checking GTDB-Tk on startup.
GTDBTK_VER_CHECK = True
GTDBTK_VER_TIMEOUT = 3  # seconds

# Internal settings used for logging.
LOG_TASK = 21




#
#New Version of AR53_MARKERS
# AR53_MARKERS = {"PFAM": ["PF01868.17.hmm", "PF01282.20.hmm", "PF01655.19.hmm",
#                           "PF01092.20.hmm", "PF01000.27.hmm", "PF00368.19.hmm",
#                           "PF00827.18.hmm", "PF01269.18.hmm", "PF00466.21.hmm",
#                           "PF01015.19.hmm", "PF13685.7.hmm", "PF02978.20.hmm",
#                           "PF04919.13.hmm", "PF01984.21.hmm", "PF04104.15.hmm",
#                           "PF00410.20.hmm", "PF01798.19.hmm", "PF01864.18.hmm",
#                           "PF01990.18.hmm", "PF07541.13.hmm", "PF04019.13.hmm",
#                           "PF00900.21.hmm", "PF01090.20.hmm", "PF02006.17.hmm",
#                           "PF01157.19.hmm", "PF01191.20.hmm", "PF01866.18.hmm",
#                           "PF01198.20.hmm", "PF01496.20.hmm", "PF00687.22.hmm",
#                           "PF03874.17.hmm", "PF01194.18.hmm", "PF01200.19.hmm",
#                           "PF13656.7.hmm", "PF01280.21.hmm"],
#                  "TIGRFAM": ["TIGR00468.HMM", "TIGR01060.HMM", "TIGR03627.HMM",
#                              "TIGR01020.HMM", "TIGR02258.HMM", "TIGR00293.HMM",
#                              "TIGR00389.HMM", "TIGR01012.HMM", "TIGR00490.HMM",
#                              "TIGR03677.HMM", "TIGR03636.HMM", "TIGR03722.HMM",
#                              "TIGR00458.HMM", "TIGR00291.HMM", "TIGR00670.HMM",
#                              "TIGR00064.HMM", "TIGR03629.HMM", "TIGR00021.HMM",
#                              "TIGR03672.HMM", "TIGR00111.HMM", "TIGR03684.HMM",
#                              "TIGR01077.HMM", "TIGR01213.HMM", "TIGR01080.HMM",
#                              "TIGR00501.HMM", "TIGR00729.HMM", "TIGR01038.HMM",
#                              "TIGR00270.HMM", "TIGR03628.HMM", "TIGR01028.HMM",
#                              "TIGR00521.HMM", "TIGR03671.HMM", "TIGR00240.HMM",
#                              "TIGR02390.HMM", "TIGR02338.HMM", "TIGR00037.HMM",
#                              "TIGR02076.HMM", "TIGR00335.HMM", "TIGR01025.HMM",
#                              "TIGR00471.HMM", "TIGR00336.HMM", "TIGR00522.HMM",
#                              "TIGR02153.HMM", "TIGR02651.HMM", "TIGR03674.HMM",
#                              "TIGR00323.HMM", "TIGR00134.HMM", "TIGR02236.HMM",
#                              "TIGR03683.HMM", "TIGR00491.HMM", "TIGR00658.HMM",
#                              "TIGR03680.HMM", "TIGR00392.HMM", "TIGR00422.HMM",
#                              "TIGR00279.HMM", "TIGR01052.HMM", "TIGR00442.HMM",
#                              "TIGR00308.HMM", "TIGR00398.HMM", "TIGR00456.HMM",
#                              "TIGR00549.HMM", "TIGR00408.HMM", "TIGR00432.HMM",
#                              "TIGR00264.HMM", "TIGR00982.HMM", "TIGR00324.HMM",
#                              "TIGR01952.HMM", "TIGR03626.HMM", "TIGR03670.HMM",
#                              "TIGR00337.HMM", "TIGR01046.HMM", "TIGR01018.HMM",
#                              "TIGR00936.HMM", "TIGR00463.HMM", "TIGR01309.HMM",
#                              "TIGR03653.HMM", "TIGR00042.HMM", "TIGR02389.HMM",
#                              "TIGR00307.HMM", "TIGR03673.HMM", "TIGR00373.HMM",
#                              "TIGR01008.HMM", "TIGR00283.HMM", "TIGR00425.HMM",
#                              "TIGR00405.HMM", "TIGR03665.HMM", "TIGR00448.HMM"]}
