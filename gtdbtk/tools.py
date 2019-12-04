import hashlib
import os
import random
import re
import time
from itertools import islice

import math

from gtdbtk.config.output import CHECKSUM_SUFFIX


##################################################
############MISC UTILITIES########################
##################################################


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
