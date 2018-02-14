# import prodigal
import math
import time
import random
import os

from itertools import islice
from biolib.seq_io import read_fasta
from biolib.common import remove_extension


##################################################
############MISC UTILITIES########################
##################################################

def add_ncbi_prefix(refname):
    if refname.startswith("GCF_"):
        return "RS_"+refname
    elif refname.startswith("GCA_"):
        return "GB_"+refname
    else:
        return refname


def splitchunks(d, n):
    chunksize = int(math.ceil(len(d) / float(n)))
    it = iter(d)
    for _ in xrange(0, len(d), chunksize):
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


def list_genomes_dir(userdir):
    """List fasta files in a specified directory

    Parameters
    ----------
    userdir : str
        Directory path where all fasta files are

    Returns
    -------
    dict
        Dictionary indicating the genomic file for each genome.
    """
    if not os.path.exists(userdir):
        raise ValueError('{0} does not exist.'.format(userdir))
    else:
        onlygenomefiles = {f: os.path.join(userdir, f) for f in os.listdir(userdir) if os.path.isfile(os.path.join(userdir, f))}
        for potential_file in onlygenomefiles:
            try:
                read_fasta(os.path.join(userdir, potential_file))
            except:
                raise IOError("{0} is not a fasta file." .format(os.path.join(userdir, potential_file)))
        return onlygenomefiles


def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z
