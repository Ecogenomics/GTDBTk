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
MIN_REF_DATA_VERSION = 'r95'

MSA_FOLDER = os.path.join(GENERIC_PATH, "msa/")
MASK_DIR = os.path.join(GENERIC_PATH, "masks/")
PPLACER_DIR = os.path.join(GENERIC_PATH, "pplacer/")
FASTANI_DIR = os.path.join(GENERIC_PATH, "fastani/")
TAX_FOLDER = os.path.join(GENERIC_PATH, "taxonomy/")
RADII_DIR = os.path.join(GENERIC_PATH, "radii/")
METADATA_DIR = os.path.join(GENERIC_PATH, "metadata/")
RED_DIR = os.path.join(GENERIC_PATH, "mrca_red/")
MARKER_DIR = os.path.join(GENERIC_PATH, 'markers/')
TIGRFAM_HMMS = os.path.join(MARKER_DIR, 'tigrfam/tigrfam.hmm')
PFAM_HMM_DIR = os.path.join(MARKER_DIR, 'pfam/')

SPLIT_DIR = os.path.join(GENERIC_PATH, 'split')
HIGH_SPLIT_DIR = os.path.join(SPLIT_DIR, 'high')
LOW_SPLIT_DIR = os.path.join(SPLIT_DIR, 'low')
HIGH_PPLACER_DIR = os.path.join(HIGH_SPLIT_DIR, 'pplacer')
LOW_PPLACER_DIR = os.path.join(LOW_SPLIT_DIR, 'pplacer')
HIGH_RED_DIR = os.path.join(HIGH_SPLIT_DIR, 'red')
LOW_RED_DIR = os.path.join(LOW_SPLIT_DIR, 'red')

LOW_TREE_MAPPING_FILE = os.path.join(LOW_SPLIT_DIR, 'tree_mapping.tsv')

HIGH_PPLACER_REF_PKG = 'gtdbtk_package_high_level'
HIGH_RED_FILE = os.path.join(HIGH_RED_DIR, 'high_red_value.tsv')
LOW_PPLACER_REF_PKG = os.path.join(LOW_PPLACER_DIR, 'gtdbtk.package.{iter}.refpkg')
LOW_RED_FILE = os.path.join(LOW_RED_DIR, 'red_value_{iter}.tsv')

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
BAC120_MARKERS = {"PFAM": ["PF00380.14.hmm", "PF00410.14.hmm", "PF00466.15.hmm",
                           "PF01025.14.hmm", "PF02576.12.hmm", "PF03726.9.hmm"],
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

AR122_MARKERS = {"PFAM": ["PF01868.11.hmm", "PF01282.14.hmm", "PF01655.13.hmm",
                          "PF01092.14.hmm", "PF01000.21.hmm", "PF00368.13.hmm",
                          "PF00827.12.hmm", "PF01269.12.hmm", "PF00466.15.hmm",
                          "PF01015.13.hmm", "PF13685.1.hmm", "PF02978.14.hmm",
                          "PF04919.7.hmm", "PF01984.15.hmm", "PF04104.9.hmm",
                          "PF00410.14.hmm", "PF01798.13.hmm", "PF01864.12.hmm",
                          "PF01990.12.hmm", "PF07541.7.hmm", "PF04019.7.hmm",
                          "PF00900.15.hmm", "PF01090.14.hmm", "PF02006.11.hmm",
                          "PF01157.13.hmm", "PF01191.14.hmm", "PF01866.12.hmm",
                          "PF01198.14.hmm", "PF01496.14.hmm", "PF00687.16.hmm",
                          "PF03874.11.hmm", "PF01194.12.hmm", "PF01200.13.hmm",
                          "PF13656.1.hmm", "PF01280.15.hmm"],
                 "TIGRFAM": ["TIGR00468.HMM", "TIGR01060.HMM", "TIGR03627.HMM",
                             "TIGR01020.HMM", "TIGR02258.HMM", "TIGR00293.HMM",
                             "TIGR00389.HMM", "TIGR01012.HMM", "TIGR00490.HMM",
                             "TIGR03677.HMM", "TIGR03636.HMM", "TIGR03722.HMM",
                             "TIGR00458.HMM", "TIGR00291.HMM", "TIGR00670.HMM",
                             "TIGR00064.HMM", "TIGR03629.HMM", "TIGR00021.HMM",
                             "TIGR03672.HMM", "TIGR00111.HMM", "TIGR03684.HMM",
                             "TIGR01077.HMM", "TIGR01213.HMM", "TIGR01080.HMM",
                             "TIGR00501.HMM", "TIGR00729.HMM", "TIGR01038.HMM",
                             "TIGR00270.HMM", "TIGR03628.HMM", "TIGR01028.HMM",
                             "TIGR00521.HMM", "TIGR03671.HMM", "TIGR00240.HMM",
                             "TIGR02390.HMM", "TIGR02338.HMM", "TIGR00037.HMM",
                             "TIGR02076.HMM", "TIGR00335.HMM", "TIGR01025.HMM",
                             "TIGR00471.HMM", "TIGR00336.HMM", "TIGR00522.HMM",
                             "TIGR02153.HMM", "TIGR02651.HMM", "TIGR03674.HMM",
                             "TIGR00323.HMM", "TIGR00134.HMM", "TIGR02236.HMM",
                             "TIGR03683.HMM", "TIGR00491.HMM", "TIGR00658.HMM",
                             "TIGR03680.HMM", "TIGR00392.HMM", "TIGR00422.HMM",
                             "TIGR00279.HMM", "TIGR01052.HMM", "TIGR00442.HMM",
                             "TIGR00308.HMM", "TIGR00398.HMM", "TIGR00456.HMM",
                             "TIGR00549.HMM", "TIGR00408.HMM", "TIGR00432.HMM",
                             "TIGR00264.HMM", "TIGR00982.HMM", "TIGR00324.HMM",
                             "TIGR01952.HMM", "TIGR03626.HMM", "TIGR03670.HMM",
                             "TIGR00337.HMM", "TIGR01046.HMM", "TIGR01018.HMM",
                             "TIGR00936.HMM", "TIGR00463.HMM", "TIGR01309.HMM",
                             "TIGR03653.HMM", "TIGR00042.HMM", "TIGR02389.HMM",
                             "TIGR00307.HMM", "TIGR03673.HMM", "TIGR00373.HMM",
                             "TIGR01008.HMM", "TIGR00283.HMM", "TIGR00425.HMM",
                             "TIGR00405.HMM", "TIGR03665.HMM", "TIGR00448.HMM"]}

# Information for Multiple hits markers:
DEFAULT_MULTIHIT_THRESHOLD = 10.0

# Information for aligning genomes
DEFAULT_DOMAIN_THRESHOLD = 10.0
AR_MARKER_COUNT = 122
BAC_MARKER_COUNT = 120

# Information about alignment Fraction to resolve fastANI results
AF_THRESHOLD = 0.65

# MSA file names
CONCAT_BAC120 = os.path.join(MSA_FOLDER, f"gtdb_{VERSION_DATA}_bac120.faa")
CONCAT_AR122 = os.path.join(MSA_FOLDER, f"gtdb_{VERSION_DATA}_ar122.faa")

# Taxonomy file name
TAXONOMY_FILE = os.path.join(TAX_FOLDER, "gtdb_taxonomy.tsv")

# Type Strain radii file
RADII_FILE = os.path.join(RADII_DIR, "gtdb_radii.tsv")

# Mask file names
MASK_BAC120 = f"gtdb_{VERSION_DATA}_bac120.mask"
MASK_AR122 = f"gtdb_{VERSION_DATA}_ar122.mask"
MASK_RPS23 = f"gtdb_{VERSION_DATA}_rps23.mask"

# Pplacer configuration
PPLACER_BAC120_REF_PKG = f"gtdb_{VERSION_DATA}_bac120.refpkg"
PPLACER_AR122_REF_PKG = f"gtdb_{VERSION_DATA}_ar122.refpkg"
PPLACER_RPS23_REF_PKG = f"gtdb_{VERSION_DATA}_rps23.refpkg"

# Fastani configuration
FASTANI_SPECIES_THRESHOLD = 95.0
FASTANI_GENOMES = os.path.join(FASTANI_DIR, "database/")
FASTANI_GENOMES_EXT = "_genomic.fna.gz"

# MRCA RED VALUE
MRCA_RED_BAC120 = os.path.join(RED_DIR, f"gtdbtk_{VERSION_DATA}_bac120.tsv")
MRCA_RED_AR122 = os.path.join(RED_DIR, f"gtdbtk_{VERSION_DATA}_ar122.tsv")

# Hashing information for validating the reference package.
REF_HASHES = {PPLACER_DIR: 'f41cfe0284ebaca4485b42e054936190c6a88bd1',
              MASK_DIR: '63551a43333bc6cbc9abf139ce881847ca19240b',
              MARKER_DIR: 'a325720422d8348d7a934143cc86112b6c92ac98',
              RADII_DIR: '1092727925f38a8a2b3f4fb40e3316c0083671f5',
              MSA_FOLDER: 'cf91d712c733e7e2535a41e6153c12b3c37d1ede',
              METADATA_DIR: 'e003b4d5d48302e85c536751f663a70447de83d4',
              TAX_FOLDER: '30c5970b2eaf5df654b2e01bfa39265302c0be89',
              FASTANI_DIR: '6a3555bb61d9cc3163c26e65772b96b8f58a2d84',
              RED_DIR: '6f661eef8e172a8a7e78af2a74fe4d079a3f5b0f'}

# Config values for checking GTDB-Tk on startup.
GTDBTK_VER_CHECK = True
GTDBTK_VER_TIMEOUT = 3  # seconds

# Internal settings used for logging.
LOG_TASK = 21
