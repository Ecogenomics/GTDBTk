#!/usr/bin/env python3

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


__prog_name__ = 'generate_skani_sketchdb.py'
__prog_desc__ = 'Regenerate the SKANI sketch database from fasta files.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2026'
__credits__ = ['Pierre Chaumeil', 'Donovan Parks','Aaron Mussig']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

# !/usr/bin/env python3

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
###############################################################################

__prog_name__ = 'generate_skani_sketchdb.py'
__prog_desc__ = 'Regenerate the SKANI sketch database recursively from a directory.'
__version__ = '0.2.0'

import hashlib
import os
import sys
import argparse
import subprocess
import shutil
import re
import logging
import tempfile
import datetime
from tqdm import tqdm

# Configure Logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Constants
METADATA_FILE = 'metadata.tsv'
CHECKSUM_FILES = {'index.db', 'markers.bin', 'sketches.db'}
REQUIRED_FILES = CHECKSUM_FILES.union({METADATA_FILE})



def get_skani_version():
    """Returns the version of skani on the system path.

    Returns
    -------
    str
        The string containing the skani version.
    """
    try:
        proc = subprocess.Popen(['skani', '-V'], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = proc.communicate()
        version = re.search(r'skani (.+)', stdout)
        return version.group(1)
    except Exception as e:
        return 'unknown'

def md5sum(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def write_metadata(output_dir, skani_version,genomes, args):
    """
    Writes a metadata file to the output directory to track versioning and parameters.
    """
    meta_path = os.path.join(output_dir, METADATA_FILE)
    with open(meta_path, 'w') as f:
        f.write(f"skani_version:{skani_version}\n")
        f.write(f"date_created:{datetime.datetime.now().isoformat()}\n")
        f.write(f"input_directory:{os.path.abspath(args.input)}\n")
        f.write(f"extension_filter:{args.extension}\n")
        f.write(f"num_genomes_sketches:{len(genomes)}\n")
        # Add Checksum values for integrity verification
        for file_name in CHECKSUM_FILES:
            file_path = os.path.join(output_dir, file_name)
            if os.path.exists(file_path):
                checksum = md5sum(file_path)
                f.write(f"{file_name}_md5:{checksum}\n")


def check_output_dir(sketching_dir):
    """
    Checks if the output directory is valid and complete.
    """
    if os.path.exists(sketching_dir):
        existing_files = set(os.listdir(sketching_dir))

        # Directory exists — inspect contents
        missing = REQUIRED_FILES - existing_files
        extra = existing_files - REQUIRED_FILES

        if not missing and not extra:
            logger.info(f"Valid sketch database found at {sketching_dir} (includes metadata).")
            logger.warning("Please remove this directory manually or choose a new output location.")
            sys.exit(0)

        if extra:
            logger.error(f"Directory {sketching_dir} contains unexpected files: {extra}")
            logger.error("Please remove this directory manually or choose a new output location.")
            sys.exit(1)

        if missing:
            logger.warning(f"Directory {sketching_dir} is incomplete. Missing required files: {missing}")
            shutil.rmtree(sketching_dir)



def get_genome_files(input_dir, extension):
    """
    Recursively find genome files with a specific extension and non-zero size.
    """
    valid_genomes = []
    logger.info(f"Scanning {input_dir} for files ending in '{extension}'...")

    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(extension):
                file_path = os.path.join(root, file)
                if os.path.getsize(file_path) == 0:
                    logger.warning(f"Skipping empty file: {file_path}")
                    continue
                valid_genomes.append(os.path.abspath(file_path))

    return valid_genomes


def run_skani_sketch(genome_list, output_dir, cpus):
    """
    Runs skani sketch using a temporary list file.
    """
    total_genomes = len(genome_list)
    if total_genomes == 0:
        logger.error("No valid genomes found to sketch.")
        sys.exit(1)

    # Create a temporary file to list the genomes
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_list:
        for path in genome_list:
            tmp_list.write(path + '\n')
        list_file_path = tmp_list.name

    try:
        args = ['skani', 'sketch']

        args += ['-o', output_dir]
        args += ['-t', str(cpus)]
        args += ['-l', list_file_path]

        logger.info(f"Starting skani sketch on {total_genomes} genomes using {cpus} threads.")

        proc = subprocess.Popen(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding='utf-8',
            bufsize=1,
            universal_newlines=True
        )

        progress_re = re.compile(r"(\d+)\s+sequences sketched")
        pbar = tqdm(total=total_genomes, desc="Sketching", unit="genomes", ncols=100)
        last_count = 0

        for line in proc.stdout:
            line = line.strip()
            # Pass through errors/warnings to console, but keep progress bar clean
            if "error" in line.lower() or "warning" in line.lower():
                if not progress_re.search(line):
                    tqdm.write(line)

            m = progress_re.search(line)
            if m:
                count = int(m.group(1))
                delta = count - last_count
                if delta > 0:
                    pbar.update(delta)
                    last_count = count

        if last_count < total_genomes:
            pbar.update(total_genomes - last_count)

        pbar.close()
        return_code = proc.wait()

        if return_code != 0:
            logger.error(f"Skani failed with exit code {return_code}")
            sys.exit(return_code)

    finally:
        if os.path.exists(list_file_path):
            os.remove(list_file_path)


def main():

    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')


    parser = argparse.ArgumentParser(description=__prog_desc__)

    parser.add_argument('-i', '--input', required=True, help="Input directory containing genome files")
    parser.add_argument('-o', '--output', required=True, help="Output directory for the sketch database")
    parser.add_argument('-x', '--extension', default='.gz', help="File extension to filter by (default: .gz)")
    parser.add_argument('-t', '--cpus', type=int, default=os.cpu_count(), help="Number of threads to use")

    args = parser.parse_args()

    skani_version = get_skani_version()
    logger.info(f"Running skani v{skani_version}")

    check_output_dir(args.output)

    #Find Genomes
    genomes = get_genome_files(args.input, args.extension)

    #Run Sketch
    run_skani_sketch(genomes, args.output, args.cpus)

    #Write Metadata (Only on success)
    write_metadata(args.output, skani_version,genomes, args)
    logger.info(f"Sketching complete. Database and metadata saved to: {args.output}")


if __name__ == '__main__':
    main()

