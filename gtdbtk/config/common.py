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

    MIN_REF_DATA_VERSION = 'r232'
    COMPATIBLE_REF_DATA_VERSIONS = ['r232']

    BACKBONE_PPLACER_REF_PKG = 'gtdbtk_package_backbone.refpkg'
    CLASS_LEVEL_PPLACER_REF_PKG = 'gtdbtk.package.{iter}.refpkg'

    # Relative Evolution Distance
    RED_INTERVAL = 0.1
    RED_MIN_SUPPORT = 0.0
    RED_MIN_CHILDREN = 2

    # Marker information
    BAC120_MARKERS = {
                "PFAM": [],
                "TIGRFAM": [
                    "custom_1.HMM", "custom_2.HMM", "custom_3.HMM", "custom_4.HMM",
                    "custom_5.HMM", "custom_6.HMM", "custom_7.HMM", "custom_8.HMM",
                    "custom_9.HMM", "custom_10.HMM", "custom_11.HMM", "custom_12.HMM",
                    "custom_13.HMM", "custom_14.HMM", "custom_15.HMM", "custom_16.HMM",
                    "custom_17.HMM", "custom_18.HMM", "custom_19.HMM", "custom_20.HMM",
                    "custom_21.HMM", "custom_22.HMM", "custom_23.HMM", "custom_24.HMM",
                    "custom_25.HMM", "custom_26.HMM", "custom_27.HMM", "custom_28.HMM",
                    "custom_29.HMM", "custom_30.HMM", "custom_31.HMM", "custom_32.HMM",
                    "custom_33.HMM", "custom_34.HMM", "custom_35.HMM", "custom_36.HMM",
                    "custom_37.HMM", "custom_38.HMM", "custom_39.HMM", "custom_40.HMM",
                    "custom_41.HMM", "custom_42.HMM", "custom_43.HMM", "custom_44.HMM",
                    "custom_45.HMM", "custom_46.HMM", "custom_47.HMM", "custom_48.HMM",
                    "custom_49.HMM", "custom_50.HMM", "custom_51.HMM", "custom_52.HMM",
                    "custom_53.HMM", "custom_54.HMM", "custom_55.HMM", "custom_56.HMM",
                    "custom_57.HMM", "custom_58.HMM", "custom_59.HMM", "custom_60.HMM",
                    "custom_61.HMM", "custom_62.HMM", "custom_63.HMM", "custom_64.HMM",
                    "custom_65.HMM", "custom_66.HMM", "custom_67.HMM", "custom_68.HMM",
                    "custom_69.HMM", "custom_70.HMM", "custom_71.HMM", "custom_72.HMM",
                    "custom_73.HMM", "custom_74.HMM", "custom_75.HMM", "custom_76.HMM",
                    "custom_77.HMM", "custom_78.HMM", "custom_79.HMM", "custom_80.HMM",
                    "custom_81.HMM", "custom_82.HMM", "custom_83.HMM", "custom_84.HMM",
                    "custom_85.HMM", "custom_86.HMM", "custom_87.HMM", "custom_88.HMM",
                    "custom_89.HMM", "custom_90.HMM", "custom_91.HMM", "custom_92.HMM",
                    "custom_93.HMM", "custom_94.HMM", "custom_95.HMM", "custom_96.HMM",
                    "custom_97.HMM", "custom_98.HMM", "custom_99.HMM", "custom_100.HMM",
                    "custom_101.HMM", "custom_102.HMM", "custom_103.HMM", "custom_104.HMM",
                    "custom_105.HMM", "custom_106.HMM", "custom_107.HMM", "custom_108.HMM",
                    "custom_109.HMM", "custom_110.HMM", "custom_111.HMM", "custom_112.HMM",
                    "custom_113.HMM", "custom_114.HMM", "custom_115.HMM", "custom_116.HMM",
                    "custom_117.HMM", "custom_118.HMM", "custom_119.HMM", "custom_120.HMM",
                    "custom_121.HMM", "custom_122.HMM", "custom_123.HMM", "custom_124.HMM",
                    "custom_125.HMM", "custom_126.HMM", "custom_127.HMM", "custom_128.HMM",
                    "custom_129.HMM", "custom_130.HMM", "custom_131.HMM", "custom_132.HMM",
                    "custom_133.HMM", "custom_134.HMM", "custom_135.HMM", "custom_136.HMM",
                    "custom_137.HMM", "custom_138.HMM", "custom_139.HMM", "custom_140.HMM",
                    "custom_141.HMM", "custom_142.HMM", "custom_143.HMM", "custom_144.HMM",
                    "custom_145.HMM", "custom_146.HMM", "custom_147.HMM", "custom_148.HMM",
                    "custom_149.HMM", "custom_150.HMM", "custom_151.HMM", "custom_152.HMM",
                    "custom_153.HMM", "custom_154.HMM", "custom_155.HMM", "custom_156.HMM",
                    "custom_157.HMM", "custom_158.HMM", "custom_159.HMM", "custom_160.HMM",
                    "custom_161.HMM", "custom_162.HMM", "custom_163.HMM", "custom_164.HMM",
                    "custom_165.HMM", "custom_166.HMM", "custom_167.HMM", "custom_168.HMM",
                    "custom_169.HMM", "custom_170.HMM", "custom_171.HMM", "custom_172.HMM",
                    "custom_173.HMM", "custom_174.HMM", "custom_175.HMM", "custom_176.HMM",
                    "custom_177.HMM", "custom_178.HMM", "custom_179.HMM", "custom_180.HMM",
                    "custom_181.HMM", "custom_182.HMM", "custom_183.HMM", "custom_184.HMM",
                    "custom_185.HMM", "custom_186.HMM", "custom_187.HMM", "custom_188.HMM",
                    "custom_189.HMM", "custom_190.HMM", "custom_191.HMM", "custom_192.HMM",
                    "custom_193.HMM", "custom_194.HMM", "custom_195.HMM", "custom_196.HMM",
                    "custom_197.HMM", "custom_198.HMM", "custom_199.HMM", "custom_200.HMM",
                    "custom_201.HMM", "custom_202.HMM"]}

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
    BAC_MARKER_COUNT = 202

    # Information about alignment Fraction to resolve skani results
    AF_THRESHOLD = 0.5

    PPLACER_MIN_RAM_BAC_FULL = 320
    PPLACER_MIN_RAM_BAC_SPLIT = 55
    PPLACER_MIN_RAM_ARC = 40

    SKANI_SPECIES_THRESHOLD = 95.0
    SKANI_IDENTITY_SKETCH_THRESHOLD = 85.0
    SKANI_GENOMES_EXT = "_genomic.fna.gz"
    SKANI_MIN_AF = 15.0


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
        return os.path.abspath("/data/pam/team162/sd28/scratch/prevotella_prj/full_isolate_tree/gene_alignments/hmms")

    @property
    def TIGRFAM_HMMS(self):
        return os.path.join(self.MARKER_DIR, 'custom_markers/custom_marker.hmm')

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
    def SKANI_SKETCHDB(self):
        return os.path.join(self.SKANI_DIR, "database/")

    @property
    def SKANI_REFERENCE_EXTENSION(self):
        return "_genomic.fna.gz"

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

        if version is None or version==232:
            return {
                self.PPLACER_DIR: '4e8ffa1133f10bde827bddb81413d06d62de341e',
                self.MASK_DIR: '84c7f7a17adc134db4161b15db10ae58389a90c1',
                #self.MARKER_DIR: '163f542c3f0a40f59df45d453aa235b39aa96e27',
                self.RADII_DIR: '546c1769ae19c946ba94b91482f32002a204f598',
                self.MSA_FOLDER: '8cf2ed4ea53f9201b127dbeed71bcc0fff27e204',
                self.METADATA_DIR: '1febdba7d2513a8f43c409423aee83ed99df7a78',
                self.TAX_FOLDER: 'c1e1766fa610229cb7cb45edc9a8fa0eafb44c96',
                self.SKANI_DIR: '00afba7b7e89a27e4d9aa9b4ce176550b56af09e',
                self.RED_DIR: '10a0fc1ca4199ff33e6c2fc4bc933318962ab212'
            }


    REF_HASHES = property(get_REF_HASHES)


# Export the class for import by other modules
@lru_cache(maxsize=1)
def __get_config():
    return __GTDBTkCommonConfig()


CONFIG = __get_config()
