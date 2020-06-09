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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import logging
import re
from collections import defaultdict

import dendropy

from gtdbtk.biolib_lite.common import canonical_gid, is_float

"""
To do:
 1. There is a serious hack in taxonomic_consistency which should be resolved, but
     requires the viral and plasmid phylogenies to be taxonomically consistent.
"""


class Taxonomy(object):
    """Manipulation of Greengenes-style taxonomy files and strings.

    This class currently assumes a Greengenes-style taxonomy
    string with the following 7 taxonomic ranks:
      d__; c__; o__; f__; g__; s__

    Spaces after the semi-colons are optional.
    """

    DOMAIN_IDX = 0
    PHYLUM_IDX = 1
    CLASS_IDX = 2
    ORDER_IDX = 3
    FAMILY_IDX = 4
    GENUS_IDX = 5
    SPECIES_IDX = 6

    rank_prefixes = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
    rank_labels = ('domain', 'phylum', 'class', 'order',
                   'family', 'genus', 'species')
    rank_index = {'d__': 0, 'p__': 1, 'c__': 2,
                  'o__': 3, 'f__': 4, 'g__': 5, 's__': 6}

    unclassified_rank = 'unclassified'

    unclassified_taxon = []
    for p in rank_prefixes:
        unclassified_taxon.append(p + unclassified_rank)
    unclassified_taxon = ';'.join(unclassified_taxon)

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

    def taxa(self, tax_str):
        """Taxa specified by taxonomy string.

        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.

        Returns
        -------
        list : [<domain>, <phylum>, ..., <species>]
            Rank order list of taxa.
        """

        taxa = [x.strip() for x in tax_str.split(';')]

        return taxa

    def taxa_at_ranks(self, tax_str):
        """Taxon at each taxonomic rank.

        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.

        Returns
        -------
        dict : d[rank_label] -> taxon
            Taxon at each taxonomic rank.
        """

        taxa = self.taxa(tax_str)

        d = {}
        for rank, taxon in enumerate(taxa):
            d[Taxonomy.rank_labels[rank]] = taxon

    def check_full(self, tax_str):
        """Check if taxonomy string specifies all expected ranks.

        Parameters
        ----------
        tax_str : str
            Greengenes-style taxonomy string.

        Returns
        -------
        boolean
            True if string contains all expected ranks, else False.
        """

        taxa = [x.strip() for x in tax_str.split(';')]
        if len(taxa) != len(Taxonomy.rank_prefixes):
            self.logger.error(
                '[Error] Taxonomy string contains too few ranks:')
            self.logger.error('[Error] %s' % str(taxa))
            return False

        for r, taxon in enumerate(taxa):
            if taxon[0:3] != Taxonomy.rank_prefixes[r]:
                self.logger.error(
                    '[Error] Taxon is not prefixed with the expected rank, %s.:' % Taxonomy.rank_prefixes[r])
                self.logger.error('[Error] %s' % str(taxa))
                return False

        return True

    def fill_trailing_ranks(self, taxa):
        """Fill in missing trailing ranks in a taxonomy string.

        Parameters
        ----------
        taxa : [d__<taxon>, ..., s__<taxon>]
            List of taxa.

        Returns
        -------
        list
            List of taxa with filled trailing ranks.
        """

        if not taxa:
            return ';'.join(Taxonomy.rank_prefixes)

        last_rank = Taxonomy.rank_prefixes.index(taxa[-1][0:3])

        for i in range(last_rank + 1, len(Taxonomy.rank_prefixes)):
            taxa.append(Taxonomy.rank_prefixes[i])

        return taxa

    def fill_missing_ranks(self, taxa, warning=True):
        """Fill in any missing ranks in a taxonomy string.

        This function assumes the taxonomic ranks are
        in the proper rank order, but that some ranks may
        be missing.

        e.g., [d__<taxon>, p__<taxon>, p__<taxon>, f__<taxon>, s__<taxon>]

        Parameters
        ----------
        taxa : [d__<taxon>, ..., s__<taxon>]
            List of taxa.

        Returns
        -------
        list
            List of taxa with all ranks.
        """

        new_taxa = []
        prev_rank_index = -1
        for t in taxa:
            rank_index = Taxonomy.rank_index[t[0:3]]
            if rank_index == prev_rank_index + 1:
                # rank is in proper order
                new_taxa.append(t)
            elif rank_index == prev_rank_index:
                # rank is repeated which is fine
                new_taxa.append(t)
            elif rank_index - prev_rank_index > 1:
                # fill in all missing ranks
                for r in range(prev_rank_index + 1, rank_index):
                    new_taxa.append(Taxonomy.rank_prefixes[r])
                new_taxa.append(t)
            elif warning:
                # current rank is more basal than previous rank
                print('Taxa have none canonical rank ordering: %s' % taxa)
                return taxa

            prev_rank_index = rank_index

        return new_taxa

    def taxonomic_consistency(self, taxonomy, report_errors=True):
        """Determine taxonomically consistent classification for taxa at each rank.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        report_errors : boolean
            Flag indicating if errors should be written to screen.

        Returns
        -------
        dict : d[taxa] -> expected parent
            Expected parent taxon for taxa at all taxonomic ranks, or
            None if the taxonomy is inconsistent.
        """

        expected_parent = {}
        for genome_id, taxa in taxonomy.items():
            if taxa[0] == 'd__Viruses' or '[P]' in taxa[0]:
                # *** This is a HACK. It would be far better to enforce
                # a taxonomically consistent taxonomy, but
                # the viral taxonomy at IMG is currently not consistent
                continue

            for r in range(1, len(taxa)):
                if taxa[r] == Taxonomy.rank_prefixes[r]:
                    break

                if taxa[r] in expected_parent:
                    if report_errors:
                        if taxa[r - 1] != expected_parent[taxa[r]]:
                            self.logger.warning(
                                'Provided taxonomy is not taxonomically consistent.')
                            self.logger.warning('Genome %s indicates the parent of %s is %s.' % (
                                genome_id, taxa[r], taxa[r - 1]))
                            self.logger.warning('The parent of this taxa was previously indicated as %s.' % (
                                expected_parent[taxa[r]]))
                            # return None

                expected_parent[taxa[r]] = taxa[r - 1]

        return expected_parent

    def extract_valid_species_name(self, taxon):
        """Try to extract a valid species name from a taxonomic label.

        A full species name should be  binomial and include a 'generic name' (genus) and
        a 'specific epithet' (species), i.e. Escherichia coli. This method
        assumes the two names should be separated by a space.

        Parameters
        ----------
        taxon : str
            Taxon label to process.

        Returns
        -------
        str
            Valid species name, or None.
        """

        if ' bacterium' in taxon.lower() or 'sp.' in taxon.lower():
            return None

        taxon = taxon.replace('s__', '')
        taxon = taxon.replace('Candidatus', '')
        taxon = taxon.replace('candidatus', '')

        if not taxon or taxon[0].islower():
            return None

        taxon_split = taxon.split(' ')
        if len(taxon_split) < 2:
            return None

        # sanity check
        taxon = 's__' + ' '.join(taxon_split[0:2])
        self.validate_species_name(taxon)

        return taxon

    def validate_species_name(self, species_name, require_full=True, require_prefix=True):
        """Validate species name.

        A full species name should be  binomial and include a 'generic name' (genus) and
        a 'specific epithet' (species), i.e. Escherichia coli. This method
        assumes the two names should be separated by a space.

        Parameters
        ----------
        species_name : str
            Species name to validate
        require_full : boolean
            Flag indicating if species name must include 'generic name and 'specific epithet'.
        require_prefix : boolean
            Flag indicating if name must start with the species prefix ('s__').

        Returns
        -------
        boolean
            True if species name is valid, otherwise False.
        str
            Reason for failing validation, otherwise None.
        """

        if species_name == 's__':
            return True, None

        # test for prefix
        if require_prefix:
            if not species_name.startswith('s__'):
                return False, 'name is missing the species prefix'

        # remove prefix before testing other properties
        test_name = species_name
        if test_name.startswith('s__'):
            test_name = test_name[3:]

        # test for full name
        if require_full:
            if 'candidatus' in test_name.lower():
                if len(test_name.split(' ')) <= 2:
                    return False, 'name appears to be missing the generic name'
            else:
                if len(test_name.split(' ')) <= 1:
                    return False, 'name appears to be missing the generic name'

        # check for tell-tale signs on invalid species names
        if " bacterium" in test_name.lower():
            return False, "name contains the word 'bacterium'"
        if " archaeon" in test_name.lower():
            return False, "name contains the word 'archaeon'"
        if " archeaon" in test_name.lower():
            return False, "name contains the word 'archeaon'"
        if "-like" in test_name.lower():
            return False, "name contains '-like'"
        if " group " in test_name.lower():
            return False, "name contains 'group'"
        if " symbiont" in test_name.lower():
            return False, "name contains 'symbiont'"
        if " endosymbiont" in test_name.lower():
            return False, "name contains 'endosymbiont'"
        if " taxon" in test_name.lower():
            return False, "name contains 'taxon'"
        if " cluster" in test_name.lower():
            return False, "name contains 'cluster'"
        if " of " in test_name.lower():
            return False, "name contains 'of'"
        if test_name[0].islower():
            return False, 'first letter of name is lowercase'
        if 'sp.' in test_name.lower():
            return False, "name contains 'sp.'"

        return True, None

    def duplicate_names(self, taxonomy):
        """Identify duplicate names in taxonomy.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxon] -> lineages
            List of lineages for duplicate taxa.
        """

        # get lineages for each taxon name
        taxon_lineages = defaultdict(set)
        for taxa in taxonomy.values():
            for i, taxon in enumerate(taxa):
                if taxon != Taxonomy.rank_prefixes[i]:
                    taxon_lineages[taxon].add(';'.join(taxa[0:i + 1]))

        # identify taxon belonging to multiple lineages
        duplicates = {}
        for taxon, lineages in taxon_lineages.items():
            if len(lineages) >= 2:
                duplicates[taxon] = lineages

        return duplicates

    def validate(self, taxonomy,
                 check_prefixes,
                 check_ranks,
                 check_hierarchy,
                 check_species,
                 check_group_names,
                 check_duplicate_names,
                 report_errors=True):
        """Check if taxonomy forms a strict hierarchy with all expected ranks.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        check_prefixes : boolean
            Flag indicating if prefix of taxon should be validated.
        check_ranks : boolean
            Flag indicating if the presence of all ranks should be validated.
        check_hierarchy : boolean
            Flag indicating if the taxonomic hierarchy should be validated.
        check_species : boolean
            Flag indicating if the taxonomic consistency of named species should be validated.
        check_group_names : boolean
            Flag indicating if group names should be checked for invalid characters.
        check_duplicate_names : boolean
            Flag indicating if group names should be checked for duplicates.
        report_errors : boolean
            Flag indicating if errors should be written to screen.

        Returns
        -------
        dict : d[taxon_id] -> taxonomy
            Taxa with invalid number of ranks.
        dict : d[taxon_id] -> [taxon, taxonomy]
            Taxa with invalid rank prefixes.
        dict: d[taxon_id] -> [species name, error message]
            Taxa with invalid species names.
        dict: d[child_taxon_id] -> two or more parent taxon ids
            Taxa with invalid hierarchies.
        """

        # check for incomplete taxonomy strings or unexpected rank prefixes
        invalid_ranks = {}
        invalid_prefixes = {}
        invalid_species_name = {}
        invalid_group_name = {}
        for taxon_id, taxa in taxonomy.items():
            if check_ranks:
                if len(taxa) != len(Taxonomy.rank_prefixes):
                    invalid_ranks[taxon_id] = ';'.join(taxa)
                    continue

            if check_prefixes:
                for r, taxon in enumerate(taxa):
                    if taxon[0:3] != Taxonomy.rank_prefixes[r]:
                        invalid_prefixes[taxon_id] = [taxon, ';'.join(taxa)]
                        break

            if check_group_names:
                for taxon in taxa:
                    canonical_taxon = ' '.join(
                        [t.strip() for t in re.split('_[A-Z]+(?= |$)', taxon[3:])]).strip()
                    if canonical_taxon and re.match('^[a-zA-Z0-9- ]+$', canonical_taxon) is None:
                        invalid_group_name[taxon_id] = [
                            taxon, 'Taxon contains invalid characters']

            if check_species:
                genus_index = Taxonomy.rank_index['g__']
                species_index = Taxonomy.rank_index['s__']
                if len(taxa) > species_index:
                    species_name = taxa[species_index]
                    valid, error_msg = self.validate_species_name(
                        species_name, require_full=True, require_prefix=True)
                    if not valid:
                        invalid_species_name[taxon_id] = [
                            species_name, error_msg]

                    if species_name != 's__':
                        genus_name = taxa[genus_index]
                        generic_name = species_name.split()[0]
                        if genus_name[3:] != generic_name[3:]:
                            invalid_species_name[taxon_id] = [
                                species_name, 'Genus and generic names do not match: %s' % genus_name]

        # check for duplicate names
        invalid_duplicate_name = []
        if check_duplicate_names:
            invalid_duplicate_name = self.duplicate_names(taxonomy)

        # check for inconsistencies in the taxonomic hierarchy
        invalid_hierarchies = defaultdict(set)
        missing_parent = set()
        if check_hierarchy:
            expected_parent = self.taxonomic_consistency(taxonomy, False)
            for taxon_id, taxa in taxonomy.items():
                for r in range(1, len(taxa)):
                    if taxa[r] == Taxonomy.rank_prefixes[r]:
                        continue

                    if r == self.rank_index['s__'] and not check_species:
                        continue

                    if taxa[r] not in expected_parent:
                        missing_parent.add(taxa[r])
                    elif taxa[r - 1] != expected_parent[taxa[r]]:
                        invalid_hierarchies[taxa[r]].add(taxa[r - 1])
                        invalid_hierarchies[taxa[r]].add(
                            expected_parent[taxa[r]])

        if report_errors:
            if len(invalid_ranks):
                print('')
                print('Taxonomy contains too few ranks:')
                for taxon_id, taxa_str in invalid_ranks.items():
                    print('%s\t%s' % (taxon_id, taxa_str))

            if len(invalid_prefixes):
                print('')
                print('Taxonomy contains an invalid rank prefix:')
                for taxon_id, info in invalid_prefixes.items():
                    print('%s\t%s\t%s' % (taxon_id, info[0], info[1]))

            if len(invalid_group_name):
                print('')
                print('Taxa containing invalid characters:')
                for taxon_id, err_msg in invalid_group_name.items():
                    print('%s\t%s\t%s' % (taxon_id, err_msg[0], err_msg[1]))

            if len(invalid_species_name):
                print('')
                print('Taxonomy contains invalid species names:')
                for taxon_id, info in invalid_species_name.items():
                    print('%s\t%s\t%s' % (taxon_id, info[0], info[1]))

            if len(invalid_duplicate_name):
                print('')
                print('Taxonomy contains identical taxon names in multiple lineages:')
                for duplicate_name in invalid_duplicate_name.keys():
                    print('%s' % duplicate_name)

            if len(missing_parent):
                print('')
                print('Taxonomy contains taxa with an undefined parent:')
                for taxon in missing_parent:
                    print('%s' % taxon)

            if len(invalid_hierarchies):
                print('')
                print('Taxonomy contains taxa with multiple parents:')
                for child_taxon, parent_taxa in invalid_hierarchies.items():
                    print('%s\t%s' % (child_taxon, ', '.join(parent_taxa)))

        return invalid_ranks, invalid_prefixes, invalid_species_name, invalid_hierarchies, invalid_group_name

    def taxon_children(self, taxonomy):
        """Get children taxa for each taxonomic group.

        For species, this is a list of extant taxa.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxon] -> list of children taxa
            All children taxa for of each named taxonomic group.
        """

        taxon_children = defaultdict(set)
        for taxon_id, taxa in taxonomy.items():
            for i, taxon in enumerate(taxa):
                if len(taxon) == 3:
                    continue  # just rank prefix

                if len(taxa) > i + 1 and len(taxa[i + 1]) != 3:
                    taxon_children[taxon].add(taxa[i + 1])

            if len(taxa) > self.rank_index['s__']:
                taxon = taxa[self.rank_index['s__']]
                if taxon != 's__':
                    taxon_children[taxon].add(taxon_id)

        return taxon_children

    def children(self, taxon, taxonomy):
        """Get children of taxon.

        For species, this is a list of extant taxa. For higher
        ranks, this is named groups and does not include the
        extant taxa.

        Parameters
        ----------
        taxon : str
            Named taxonomic group of interest.
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        set : {child1, child2, ..., childN}
            All children taxa for the named taxonomic group.
        """

        c = set()
        for taxon_id, taxa in taxonomy.items():
            if taxon in taxa:

                if taxon.startswith('s__'):
                    c.add(taxon_id)
                else:
                    taxon_index = taxa.index(taxon)
                    for child in taxa[taxon_index + 1:]:
                        if len(child) > 3:  # not just an empty prefix
                            c.add(child)

        return c

    def parents(self, taxonomy):
        """Get parents for all taxa.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        d[taxon] -> list of parent taxa in rank order
            Parent taxa for each taxon.
        """

        p = defaultdict(list)
        for taxon_id, taxa in taxonomy.items():
            p[taxon_id] = taxa
            for i, taxon in enumerate(taxa):
                if i != 0:
                    p[taxon] = taxa[0:i]

        return p

    def extant_taxa(self, taxonomy):
        """Get extant taxa for all taxa.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxon] -> set of extant taxa
            Extant taxa for named groups at the specified rank.
        """

        extant_taxa = {}
        for rank_label in Taxonomy.rank_labels:
            extant_taxa.update(self.extant_taxa_for_rank(rank_label, taxonomy))

        return extant_taxa

    def extant_taxa_for_rank(self, rank_label, taxonomy):
        """Get extant taxa for all named groups at the specified rank.

        Parameters
        ----------
        rank_label : str (e.g., class or order)
            Rank of interest
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxon] -> set of extant taxa
            Extant taxa for named groups at the specified rank.
        """

        assert (rank_label in Taxonomy.rank_labels)

        d = defaultdict(set)
        rank_index = Taxonomy.rank_labels.index(rank_label)
        for taxon_id, taxa in taxonomy.items():
            if taxa[rank_index] != Taxonomy.rank_prefixes[rank_index]:
                d[taxa[rank_index]].add(taxon_id)

        return d

    def named_lineages_at_rank(self, taxonomy):
        """Get named lineages at each taxonomic rank.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[rank] -> set of taxa
            Taxa at each taxonomic rank.
        """

        named_lineages = defaultdict(set)
        for taxa in taxonomy.values():
            for i, taxon in enumerate(taxa):
                if taxon != Taxonomy.rank_prefixes[i]:
                    named_lineages[i].add(taxon)

        return named_lineages

    def lineages(self, taxonomy):
        """Get lineages for all taxon.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.

        Returns
        -------
        dict : d[taxon] -> lineage
            Lineage information for each taxon.
        """

        lineages = defaultdict(set)
        for taxa in taxonomy.values():
            for i, taxon in enumerate(taxa):
                lineages[taxon] = taxa[0:i]

        return lineages

    def read_from_tree(self, tree, warnings=True):
        """Obtain the taxonomy for each extant taxa as specified by internal tree labels.

        Parameters
        ----------
        tree : str or dendropy.Tree
            Filename of newick tree or dendropy tree object.

        Returns
        -------
        dict : d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
            Taxa indexed by unique ids.
        """

        if isinstance(tree, str):
            tree = dendropy.Tree.get_from_path(tree,
                                               schema='newick',
                                               rooting="force-rooted",
                                               preserve_underscores=True)

        taxonomy = {}
        for leaf in tree.leaf_node_iter():
            taxa = []

            node = leaf.parent_node
            while node:
                if node.label:
                    taxa_str = node.label
                    if ':' in taxa_str:
                        taxa_str = taxa_str.split(':')[1]

                    if not is_float(taxa_str):
                        if taxa_str[-1] == ';':
                            taxa_str = taxa_str[:-1]

                        # check for concatenated ranks of the form:
                        # p__Crenarchaeota__c__Thermoprotei
                        for prefix in Taxonomy.rank_prefixes:
                            split_str = '__' + prefix
                            if split_str in taxa_str:
                                taxa_str = taxa_str.replace(
                                    split_str, ';' + prefix)

                        # appears to be an internal label and not simply a
                        # support value
                        taxa = [x.strip() for x in taxa_str.split(';')] + taxa
                node = node.parent_node

            if warnings and len(taxa) > 7:
                self.logger.warning('Invalid taxonomy string read from tree for taxon %s: %s' % (
                    leaf.taxon.label, taxa))
                # sys.exit(-1)

            # check if genus name should be appended to species label
            if len(taxa) == 7:
                genus = taxa[5][3:]
                species = taxa[6][3:]
                if genus not in species:
                    taxa[6] = 's__' + genus + ' ' + species

            taxa = self.fill_trailing_ranks(taxa)
            taxonomy[leaf.taxon.label] = taxa

        return taxonomy

    def read(self, taxonomy_file, canonical_ids=False):
        """Read Greengenes-style taxonomy file.

        Expected format is:
            <id>\t<taxonomy string>

        where the taxonomy string has the formats:
            d__; c__; o__; f__; g__; s__

        Parameters
        ----------
        taxonomy_file : str
            Path to a Greengenes-style taxonomy file.
        canonical_ids : bool
            True if to use the canonical ID format, False otherwise.

        Returns
        -------
        dict[str, tuple[str, str, str, str, str, str, str]]
            d[unique_id] -> [d__<taxon>, ..., s__<taxon>]
        """

        try:
            d = {}
            with open(taxonomy_file, 'r') as f:
                for row, line in enumerate(f.readlines()):
                    line_split = line.split('\t')
                    unique_id = line_split[0]
                    if canonical_ids:
                        unique_id = canonical_gid(unique_id)

                    tax_str = line_split[1].rstrip()
                    if tax_str[-1] == ';':
                        # remove trailing semicolons which sometimes
                        # appear in Greengenes-style taxonomy files
                        tax_str = tax_str[0:-1]

                    d[unique_id] = [x.strip() for x in tax_str.split(';')]
        except:
            self.logger.error(
                'Failed to parse taxonomy file on line %d' % (row + 1))
            raise

        return d

    def write(self, taxonomy, output_file):
        """Write Greengenes-style taxonomy file.

        Parameters
        ----------
        taxonomy : d[unique_id] -> [d__<taxon>; ...; s__<taxon>]
            Taxonomy strings indexed by unique ids.
        output_file : str
            Name of output file.
        """

        with open(output_file, 'w') as fout:
            for genome_id, taxa in taxonomy.items():
                fout.write(genome_id + '\t' + ';'.join(taxa) + '\n')

    def sort_taxa(self, taxa, reverse=False):
        """Sort taxa by rank and then alphabetically.

        Parameters
        ----------
        taxa : iterable
            Taxa with rank prefixes.

        Returns
        -------
        list
            Taxa sorted by rank and alphabetically within each rank.
        """

        ordered_taxa = []
        for rank_prefix in Taxonomy.rank_prefixes:
            rank_taxa = []
            for taxon in taxa:
                if taxon.startswith(rank_prefix):
                    rank_taxa.append(taxon)

            ordered_taxa.extend(sorted(rank_taxa))

        if reverse:
            ordered_taxa = ordered_taxa[::-1]

        return ordered_taxa
