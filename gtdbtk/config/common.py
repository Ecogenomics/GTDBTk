import json
import os
import sys
from functools import lru_cache


class __GTDBTkCommonConfig:
    """
    This class encapsulates all configuration options. It will protect against
    importing code that depends on a specific value throwing an exception
    that requires the setting of the GTDB-Tk reference data path.
    """

    MIN_REF_DATA_VERSION = 'r226'
    COMPATIBLE_REF_DATA_VERSIONS = ['r220','r226']

    BACKBONE_PPLACER_REF_PKG = 'gtdbtk_package_backbone.refpkg'
    CLASS_LEVEL_PPLACER_REF_PKG = 'gtdbtk.package.{iter}.refpkg'

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

    # New Version of AR53_MARKERS
    AR53_MARKERS = {"PFAM": ["PF04919.13.hmm", "PF07541.13.hmm", "PF01000.27.hmm",
                             "PF00687.22.hmm", "PF00466.21.hmm", "PF00827.18.hmm", "PF01280.21.hmm", "PF01090.20.hmm",
                             "PF01200.19.hmm", "PF01015.19.hmm", "PF00900.21.hmm", "PF00410.20.hmm"],
                    "TIGRFAM": ["TIGR00037.HMM", "TIGR00064.HMM", "TIGR00111.HMM",
                                "TIGR00134.HMM", "TIGR00279.HMM", "TIGR00291.HMM", "TIGR00323.HMM",
                                "TIGR00335.HMM", "TIGR00373.HMM", "TIGR00405.HMM", "TIGR00448.HMM",
                                "TIGR00483.HMM", "TIGR00491.HMM", "TIGR00522.HMM", "TIGR00967.HMM",
                                "TIGR00982.HMM", "TIGR01008.HMM", "TIGR01012.HMM", "TIGR01018.HMM",
                                "TIGR01020.HMM", "TIGR01028.HMM", "TIGR01046.HMM", "TIGR01052.HMM",
                                "TIGR01171.HMM", "TIGR01213.HMM", "TIGR01952.HMM", "TIGR02236.HMM",
                                "TIGR02338.HMM", "TIGR02389.HMM", "TIGR02390.HMM", "TIGR03626.HMM",
                                "TIGR03627.HMM", "TIGR03628.HMM", "TIGR03629.HMM", "TIGR03670.HMM",
                                "TIGR03671.HMM", "TIGR03672.HMM", "TIGR03673.HMM", "TIGR03674.HMM",
                                "TIGR03676.HMM", "TIGR03680.HMM"]}

    # Information for Multiple hits markers:
    DEFAULT_MULTIHIT_THRESHOLD = 10.0

    # Information for aligning genomes
    DEFAULT_DOMAIN_THRESHOLD = 10.0
    AR_MARKER_COUNT = 53
    BAC_MARKER_COUNT = 120

    # Information about alignment Fraction to resolve skani results
    AF_THRESHOLD = 0.5

    PPLACER_MIN_RAM_BAC_FULL = 320
    PPLACER_MIN_RAM_BAC_SPLIT = 55
    PPLACER_MIN_RAM_ARC = 40

    SKANI_SPECIES_THRESHOLD = 95.0
    SKANI_IDENTITY_SKETCH_THRESHOLD = 85.0
    SKANI_GENOMES_EXT = "_genomic.fna.gz"
    SKANI_MIN_AF = 15.0

    # Mash configuration
    MASH_SKETCH_FILE = 'gtdb_ref_sketch.msh'
    MASH_K_VALUE = 16
    MASH_S_VALUE = 5000
    MASH_MAX_DISTANCE = 0.15
    MASH_D_VALUE = MASH_MAX_DISTANCE
    MASH_V_VALUE = 1.0

    # Config values for checking GTDB-Tk on startup.
    GTDBTK_VER_CHECK = True
    GTDBTK_VER_TIMEOUT = 3  # seconds

    # Internal settings used for logging.
    LOG_TASK = 21

    # To avoid multiple hits of parsing files
    _generic_path = None
    _red_dist_bac_dict = None
    _red_dist_arc_dict = None
    _version_data = None

    @property
    def GENERIC_PATH(self):
        if self._generic_path is None:
            try:
                # expandvars is required to transform things like $HOME
                out = os.path.expandvars(os.environ['GTDBTK_DATA_PATH'])
                self._generic_path = out
            except KeyError:
                print('\n' + '=' * 80)
                print(' ERROR '.center(80))
                print('_' * 80 + '\n')
                print("The 'GTDBTK_DATA_PATH' environment variable is not defined.".center(80) + '\n')
                print('Please set this variable to your reference data package.'.center(80))
                print('https://ecogenomics.github.io/GTDBTk/installing/index.html'.center(80))
                print('=' * 80)
                sys.exit(1)
        return self._generic_path

    @property
    def MSA_FOLDER(self):
        return os.path.join(self.GENERIC_PATH, 'msa/')

    @property
    def MASK_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'masks/')

    @property
    def PPLACER_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'pplacer/')

    @property
    def SKANI_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'skani/')

    @property
    def TAX_FOLDER(self):
        return os.path.join(self.GENERIC_PATH, 'taxonomy/')

    @property
    def RADII_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'radii/')

    @property
    def METADATA_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'metadata/')

    @property
    def RED_DIR(self):
        return os.path.join(self.GENERIC_PATH, "mrca_red/")

    @property
    def MARKER_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'markers/')

    @property
    def TIGRFAM_HMMS(self):
        return os.path.join(self.MARKER_DIR, 'tigrfam/tigrfam.hmm')

    @property
    def PFAM_HMM_DIR(self):
        return os.path.join(self.MARKER_DIR, 'pfam/')

    @property
    def SPLIT_DIR(self):
        return os.path.join(self.GENERIC_PATH, 'split')

    @property
    def BACKBONE_SPLIT_DIR(self):
        return os.path.join(self.SPLIT_DIR, 'backbone')

    @property
    def CLASS_LEVEL_SPLIT_DIR(self):
        return os.path.join(self.SPLIT_DIR, 'class_level')

    @property
    def BACKBONE_PPLACER_DIR(self):
        return os.path.join(self.BACKBONE_SPLIT_DIR, 'pplacer')

    @property
    def CLASS_LEVEL_PPLACER_DIR(self):
        return os.path.join(self.CLASS_LEVEL_SPLIT_DIR, 'pplacer')

    @property
    def BACKBONE_RED_DIR(self):
        return os.path.join(self.BACKBONE_SPLIT_DIR, 'red')

    @property
    def CLASS_LEVEL_RED_DIR(self):
        return os.path.join(self.CLASS_LEVEL_SPLIT_DIR, 'red')

    @property
    def CLASS_LEVEL_TREE_MAPPING_FILE(self):
        return os.path.join(self.CLASS_LEVEL_SPLIT_DIR, 'tree_mapping.tsv')

    @property
    def BACKBONE_RED_FILE(self):
        return os.path.join(self.BACKBONE_RED_DIR, 'backbone_red_value.tsv')

    @property
    def CLASS_LEVEL_RED_FILE(self):
        return os.path.join(self.CLASS_LEVEL_RED_DIR, 'red_value_{iter}.tsv')

    def _read_metadata_file(self):
        if self._red_dist_bac_dict is None or self._red_dist_arc_dict is None or self._version_data is None:
            try:
                with open(os.path.join(self.METADATA_DIR, "metadata.txt")) as metadataData:
                    for line in metadataData:
                        try:
                            line_infos = line.strip().split('=')
                            if line_infos[0] == 'RED_DIST_BAC_DICT':
                                self._red_dist_bac_dict = json.loads(line_infos[1])
                            elif line_infos[0] == 'RED_DIST_ARC_DICT':
                                self._red_dist_arc_dict = json.loads(line_infos[1])
                            elif line_infos[0] == 'VERSION_DATA':
                                self._version_data = line_infos[1]
                        except ValueError:
                            print(f"Skipping invalid line {repr(line)}")
            except IOError:
                print('\n' + '=' * 80)
                print(' ERROR '.center(80))
                print('_' * 80 + '\n')
                print('The GTDB-Tk reference data does not exist or is corrupted.'.center(80))
                print(('GTDBTK_DATA_PATH=%s' % self.GENERIC_PATH).center(80) + '\n')
                print('Please compare the checksum to those provided in the download repository.'.center(80))
                print('https://github.com/Ecogenomics/GTDBTk#gtdb-tk-reference-data'.center(80))
                print('=' * 80)
                sys.exit(1)
        return self._red_dist_bac_dict, self._red_dist_arc_dict, self._version_data

    @property
    def RED_DIST_BAC_DICT(self):
        return self._read_metadata_file()[0]

    @property
    def RED_DIST_ARC_DICT(self):
        return self._read_metadata_file()[1]

    @property
    def VERSION_DATA(self):
        return self._read_metadata_file()[2]

    """
    MSA file names
    """

    @property
    def CONCAT_BAC120(self):
        return os.path.join(self.MSA_FOLDER, f"gtdb_{self.VERSION_DATA}_bac120.faa")

    @property
    def CONCAT_AR53(self):
        return os.path.join(self.MSA_FOLDER, f"gtdb_{self.VERSION_DATA}_ar53.faa")

    @property
    def TAXONOMY_FILE(self):
        return os.path.join(self.TAX_FOLDER, "gtdb_taxonomy.tsv")

    @property
    def RADII_FILE(self):
        return os.path.join(self.RADII_DIR, "gtdb_radii.tsv")

    """
    Mask file names
    """

    @property
    def MASK_BAC120(self):
        return f"gtdb_{self.VERSION_DATA}_bac120.mask"

    @property
    def MASK_AR53(self):
        return f"gtdb_{self.VERSION_DATA}_ar53.mask"

    @property
    def MASK_RPS23(self):
        return f"gtdb_{self.VERSION_DATA}_rps23.mask"

    @property
    def PPLACER_BAC120_REF_PKG(self):
        return f"gtdb_{self.VERSION_DATA}_bac120.refpkg"

    @property
    def PPLACER_AR53_REF_PKG(self):
        return f"gtdb_{self.VERSION_DATA}_ar53.refpkg"

    @property
    def PPLACER_RPS23_REF_PKG(self):
        return f"gtdb_{self.VERSION_DATA}_rps23.refpkg"

    @property
    def SKANI_GENOMES(self):
        return os.path.join(self.SKANI_DIR, "database/")

    @property
    def SKANI_GENOME_LIST(self):
        return os.path.join(self.SKANI_DIR, "genome_paths.tsv")

    @property
    def MRCA_RED_BAC120(self):
        return os.path.join(self.RED_DIR, f"gtdbtk_{self.VERSION_DATA}_bac120.tsv")

    @property
    def MRCA_RED_AR53(self):
        return os.path.join(self.RED_DIR, f"gtdbtk_{self.VERSION_DATA}_ar53.tsv")

    def get_REF_HASHES(self,version=None):
        compatible_versions = [int(x.replace('r','')) for x in CONFIG.COMPATIBLE_REF_DATA_VERSIONS]
        if version is not None and version not in compatible_versions:
            raise ValueError(f"Version {version} is not compatible with this version of GTDB-Tk. Compatible versions are {compatible_versions}")
        if version is None or version==226:
            return {
                self.PPLACER_DIR: '40d242bd5b84d4b2218586390220e6f741191804',
                self.MASK_DIR: 'beb296ca29ffadf2d5b86f4160f2a525c5954185',
                self.MARKER_DIR: '163f542c3f0a40f59df45d453aa235b39aa96e27',
                self.RADII_DIR: '577663ffcd7892551a0f8cb2864fb80544c3bfcd',
                self.MSA_FOLDER: '1d47808337897c2965cfc8a101df5e4813ce1581',
                self.METADATA_DIR: '56b3a753b0cec1f0ac983fea4c7db32a160d00cb',
                self.TAX_FOLDER: 'ea99878828e0baffe8977c7371ce7b7d5109c1e1',
                self.SKANI_DIR: '8100f226743d74572fa726ac3ea358c34d5b8222',
                self.RED_DIR: 'bd8cfe5c4f5b8fc180cebf974713194395efbc00'
            }

        elif version==220:
            return {
                self.PPLACER_DIR: '75fdd0e093c9af6a73cb510c3d0cd2041265e093',
                self.MASK_DIR: 'f4b8ebfa59526a7a86f09752b47e8de1efc384c7',
                self.MARKER_DIR: '163f542c3f0a40f59df45d453aa235b39aa96e27',
                self.RADII_DIR: '63d06ecc8b4547addd22c5b06ada4a28c5332bcc',
                self.MSA_FOLDER: '3d5c1cf5346b244fcb0a9d48d2f1a9358a71cc7a',
                self.METADATA_DIR: '01b8c23253cef097b1bc233d609dae9eb84c98e2',
                self.TAX_FOLDER: '6758173fa61ae4a77f5588ec2874ea52ed345feb',
                self.SKANI_DIR: 'ff58a1d7e0584da324d140701ee12cead4f0df9d',
                self.RED_DIR: '206bd781997fffbac951b4437dd75e6543139fd6'
            }


    REF_HASHES = property(get_REF_HASHES)


# Export the class for import by other modules
@lru_cache(maxsize=1)
def __get_config():
    return __GTDBTkCommonConfig()


CONFIG = __get_config()
