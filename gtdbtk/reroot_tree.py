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

import logging
import random
from typing import Set

import dendropy

from gtdbtk.exceptions import GTDBTkExit


class RerootTree(object):
    """Reroot tree."""

    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')

    def root_with_outgroup(self, input_tree: str, output_tree: str, outgroup: Set[str]):
        """Reroot the tree using the given outgroup.

        Parameters
        ----------
        input_tree
          File containing Newick tree to rerooted.
        output_tree
          Name of file for rerooted tree.
        outgroup
          Labels of taxa in outgroup.
        """

        tree = dendropy.Tree.get_from_path(input_tree,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)

        outgroup = set(outgroup)
        outgroup_in_tree = set()
        ingroup_leaves = set()
        for n in tree.leaf_node_iter():
            if n.taxon.label in outgroup:
                outgroup_in_tree.add(n.taxon)
            else:
                ingroup_leaves.add(n)

        self.logger.info(f'Identified {len(outgroup_in_tree):,} outgroup taxa in the tree.')
        self.logger.info(f'Identified {len(ingroup_leaves):,} ingroup taxa in the tree.')

        if len(outgroup_in_tree) == 0:
            self.logger.error('No outgroup taxa identified in the tree.')
            raise GTDBTkExit('Tree was not rerooted.')

        # Since finding the MRCA is a rooted tree operation,
        # the tree is first rerooted on an ingroup taxa. This
        # ensures the MRCA of the outgroup can be identified
        # so long as the outgroup is monophyletic. If the
        # outgroup is polyphyletic trying to root on it
        # is ill defined. To try and pick a "good" root for
        # polyphyletic outgroups, random ingroup taxa are
        # selected until two of them give the same size
        # lineage. This will, likely, be the smallest
        # bipartition possible for the given outgroup though
        # this is not guaranteed.
        mrca = tree.mrca(taxa=outgroup_in_tree)
        mrca_leaves = len(mrca.leaf_nodes())
        while True:
            rnd_ingroup = random.sample(ingroup_leaves, 1)[0]
            tree.reroot_at_edge(rnd_ingroup.edge,
                                length1=0.5 * rnd_ingroup.edge_length,
                                length2=0.5 * rnd_ingroup.edge_length)

            mrca = tree.mrca(taxa=outgroup_in_tree)
            if len(mrca.leaf_nodes()) == mrca_leaves:
                break

            mrca_leaves = len(mrca.leaf_nodes())

        if len(mrca.leaf_nodes()) != len(outgroup_in_tree):
            self.logger.info('Outgroup is not monophyletic. Tree will be '
                             'rerooted at the MRCA of the outgroup.')
            self.logger.info(f'The outgroup consisted of '
                             f'{len(outgroup_in_tree):,} taxa, while the MRCA '
                             f'has {len(mrca.leaf_nodes()):,} leaf nodes.')
            if len(mrca.leaf_nodes()) == len(tree.leaf_nodes()):
                self.logger.warning('The MRCA spans all taxa in the tree.')
                self.logger.warning('This indicating the selected outgroup is '
                                    'likely polyphyletic in the current tree.')
                self.logger.warning('Polyphyletic outgroups are not suitable '
                                    'for rooting. Try another outgroup.')
        else:
            self.logger.info('Outgroup is monophyletic.')

        if mrca.edge_length is None:
            self.logger.info('Tree appears to already be rooted on this outgroup.')
        else:
            self.logger.info('Rerooting tree.')
            tree.reroot_at_edge(mrca.edge,
                                length1=0.5 * mrca.edge_length,
                                length2=0.5 * mrca.edge_length)
            tree.write_to_path(output_tree, schema='newick',
                               suppress_rooting=True, unquoted_underscores=True)
            self.logger.info(f'Rerooted tree written to: {output_tree}')

    def midpoint(self, input_tree, output_tree):
        """Reroot tree bat midpoint.

        Parameters
        ----------
        input_tree : str
          File containing Newick tree to rerooted.
        output_tree : str
          Name of file for rerooted tree.
        """

        tree = dendropy.Tree.get_from_path(
            input_tree, schema='newick', rooting='force-rooted', preserve_underscores=True)
        tree.reroot_at_midpoint()
        tree.write_to_path(output_tree, schema='newick',
                           suppress_rooting=True, unquoted_underscores=True)
