from gtdbtk.biolib_lite.common import canonical_gid
from gtdbtk.config.config import RADII_FILE


class GTDBRadiiFile(object):
    """A wrapper for the gtdb_radii.tsv file included in the reference data."""
    path = RADII_FILE

    def __init__(self):
        self._rep_idx = None
        self._species_idx = None
        self._read()

    def _read(self):
        """Read the file and create any data."""
        self._rep_idx, self._species_idx = dict(), dict()
        with open(self.path) as fh:
            for line in fh.readlines():
                species, genome, ani = line.strip().split('\t')
                genome = canonical_gid(genome)
                ani = float(ani)
                self._rep_idx[genome] = {'species': species, 'ani': ani}
                self._species_idx[species] = {'rep': genome, 'ani': ani}

    def get_species_ani(self, species):
        """Returns the ANI for a specific species.

        Parameters
        ----------
        species : str
            The name of the species.

        Returns
        -------
        float
            The ANI value for that species.
        """
        return self._species_idx[species]['ani']

    def get_species_rep(self, species):
        """Returns the representative genome for this species.

        Parameters
        ----------
        species : str
            The name of the species.

        Returns
        -------
        str
            The genome which represents this species.
        """
        return self._species_idx[species]['rep']

    def get_rep_ani(self, rep):
        """Returns the ANI for the species circumscribed by this rep.

        Parameters
        ----------
        rep : str
            The representative genome.

        Returns
        -------
        float
            The ANI value for that species.
        """
        return self._rep_idx[rep]['ani']

    def get_rep_species(self, rep):
        """Returns the species that this genome represents..

        Parameters
        ----------
        rep : str
            The representative genome.

        Returns
        -------
        str
            The species represented.
        """
        return self._rep_idx[rep]['species']
