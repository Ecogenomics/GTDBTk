import os
from shutil import copyfile

from gtdbtk.biolib_lite.common import make_sure_path_exists
from gtdbtk.config.config import CONCAT_AR53, CONCAT_BAC120
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.model.enum import Domain


def export_msa(domain: Domain, output_file: str):
    """Exports the GTDB MSA to the specified path.

    :param domain: The domain used to determine the marker set.
    :param output_file: The path to write the MSA.
    """
    if domain is Domain.ARCHAEA:
        file_to_export = CONCAT_AR53
    elif domain is Domain.BACTERIA:
        file_to_export = CONCAT_BAC120
    else:
        raise GTDBTkExit(f'Unknown domain: "{domain}"')

    make_sure_path_exists(os.path.dirname(output_file))
    copyfile(file_to_export, output_file)
