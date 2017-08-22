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


def genomes_to_process(genome_dir, batchfile):
        """Get genomes to process.

        Parameters
        ----------
        genome_dir : str
          Directory containing genomes.
        batchfile : str
          File describing genomes.
          
        Returns
        -------
        genomic_files : d[genome_id] -> FASTA file
            Map of genomes to their genomic FASTA files.
        """
        
        genomic_files = {}
        if genome_dir:
            
            for f in os.listdir(genome_dir):
                genome_id = remove_extension(f)
                genomic_files[genome_id] = os.path.join(genome_dir, f)
                
        elif batchfile:
            for line_no, line in enumerate(open(batchfile, "rb")):
                line_split = line.strip().split("\t")
                if line_split[0] == '':
                    continue # blank line
                    
                if len(line_split) != 2:
                    self.logger.error('Batch file must contain exactly 2 columns.')
                    sys.exit()

                genome_file, genome_id  = line_split
                
                if genome_file is None or genome_file == '':
                    self.logger.error('Missing genome file on line %d.' % line_no+1)
                    self.exit()
                elif genome_id is None or genome_id == '':
                    self.logger.error('Missing genome ID on line %d.' % line_no+1)
                    self.exit()
                elif genome_id in genomic_files:
                    self.logger.error('Genome ID %s appear multiple times.' % genome_id)
                    self.exit()

                genomic_files[genome_id] = genome_file
                
        for genome_key in genomic_files.iterkeys():
            if genome_key.startswith("RS_") or genome_key.startswith("GB_"):
                raise Exception("Submitted genomes start with similar prefix (RS_,GB_) as reference Genomes in GtdbTk. This may cause issues for downstream analysis.") 
            
            
        return genomic_files
    
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
