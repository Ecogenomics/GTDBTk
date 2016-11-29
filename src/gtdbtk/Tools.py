# import prodigal
import math
import time
import random
import os

from itertools import islice
from biolib.seq_io import read_fasta

import Config


##################################################
############MISC UTILITIES########################
##################################################


def splitchunks(d, n):
    chunksize = int(math.ceil(len(d) / float(n)))
    it = iter(d)
    for _ in xrange(0, len(d), chunksize):
        yield {k: d[k] for k in islice(it, chunksize)}


def generateTempTableName():
    rng = random.SystemRandom()
    suffix = ''
    for _ in range(0, 10):
        suffix += rng.choice(
            'abcefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    return "TEMP" + suffix + str(int(time.time()))


def fastaPathGenerator(path=None, prefix=None):

    genomeUserDir = None
    if Config.GTDB_GENOME_USR_DIR:
        genomeUserDir = Config.GTDB_GENOME_USR_DIR

    genomeGBKDir = None
    if Config.GTDB_GENOME_GBK_DIR:
        genomeGBKDir = Config.GTDB_GENOME_GBK_DIR

    genomeRSQDir = None
    if Config.GTDB_GENOME_RSQ_DIR:
        genomeRSQDir = Config.GTDB_GENOME_RSQ_DIR

    if prefix == 'U':
        return os.path.join(genomeUserDir, path)
    elif prefix == "GB":
        return os.path.join(genomeGBKDir, path)
    elif prefix == "RS":
        return os.path.join(genomeRSQDir, path)
    else:
        print "prefix {0} is not existing".format(prefix)


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
