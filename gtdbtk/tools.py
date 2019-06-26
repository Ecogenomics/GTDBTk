import hashlib
import math
import os
import random
import time
from itertools import islice

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
