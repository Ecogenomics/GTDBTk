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

import gzip
import hashlib
import json
import re

import dendropy

from gtdbtk.biolib_lite.newick import parse_label
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.exceptions import GTDBTkExit

TABLE_HEADER = ['edge_num', 'red_below', 'red_above', 'path_below', 'path_above', 'edge_length',
                'placement_taxonomy', 'parent_taxon', 'parent_red', 'child_taxon', 'child_red',
                'n_refs', 'first_ref', 'rank_names']
FINGERPRINTS = ['tree_sha1', 'red_sha1', 'taxonomy_sha1']
ORDER_RANK = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
MIXED = '*'                                                     # rank_names: several names below


def tree_sha1(jplace_tree):
    """Fingerprint of a jplace reference tree (its edge numbering)."""
    return hashlib.sha1(jplace_tree.encode('utf-8')).hexdigest()


def file_sha1(path):
    """Fingerprint of a file's content (the RED file the table was built from)."""
    with open(path, 'rb') as fh:
        return hashlib.sha1(fh.read()).hexdigest()


def best_placements(jplace_file):
    """Read a jplace: (tree string, {genome: (edge_num, distal_length)}) for each genome's
    best placement, as gplace tog takes it: highest like_weight_ratio, the first one on ties."""
    with open(jplace_file) as fh:
        jplace = json.load(fh)
    fields = jplace['fields']
    best = {}
    for rec in jplace['placements']:
        places = [dict(zip(fields, p)) for p in rec['p']]
        top = max(places, key=lambda p: p['like_weight_ratio'])         # max keeps the first on ties
        names = rec['n'] if 'n' in rec else [nm[0] for nm in rec['nm']]
        for name in ([names] if isinstance(names, str) else names):
            best[name] = (top['edge_num'], top['distal_length'])
    return jplace['tree'], best


def refpkg_jplace_tree(refpkg):
    """The reference tree of a gplace reference package as gplace writes it in every jplace:
    same edge numbers {n}, same text (so the same fingerprint)."""
    from gplace.engine.newick import parse_newick, to_newick
    from gplace.engine.reference import Reference
    return to_newick(parse_newick(Reference.of(refpkg).tree_text()), with_node_numbers=True)


class JplaceTableBuilder(object):
    """Per-edge data of a reference tree, so placements can be classified from the jplace
    without the gplace tog tree. Built once per reference package, from the package itself
    (see scripts/build_jplace_table.py).

    For each edge, with C the node below it and P the node above:
    - RED: RED of the closest reference nodes below and above, and the path lengths to them
      (as Classify._assign_mrca_red);
    - placement_taxonomy: labels above the edge (as Classify._get_gplace_taxonomy);
    - what the RED / topology rules read (Classify._classify_red_topology): the first labelled
      node above the edge (parent_taxon, parent_red), the first node below labelled with the
      next rank (child_taxon, child_red), the number of reference genomes below, the first of
      them, and per rank their name if they all share it ('*' otherwise).
    """

    def __init__(self, jplace_tree, red_file, taxonomy_file):
        """jplace_tree: the reference tree with gplace's edge numbers (refpkg_jplace_tree);
        red_file, taxonomy_file: the RED and taxonomy files of the same release."""
        self.jplace_tree = jplace_tree
        self.red_file = red_file
        self.taxonomy_file = taxonomy_file
        # {n} -> a comment dendropy keeps as an annotation of the node below the edge
        self.tree = dendropy.Tree.get(data=re.sub(r'\{(\d+)\}', r'[&edge_num=\1]', jplace_tree),
                                      schema='newick', rooting='force-rooted',
                                      preserve_underscores=True, extract_comment_metadata=True)
        self.node_of_edge = {}
        for nd in self.tree.preorder_node_iter():
            self.node_of_edge[int(nd.annotations.get_value('edge_num'))] = nd
            nd.red = None

        leaf_of = {leaf.taxon.label: leaf for leaf in self.tree.leaf_node_iter()}
        with open(red_file) as fh:
            for line in fh:
                label_ids, red_value = line.strip().split('\t')
                labels = label_ids.split('|')
                if len(labels) == 2:
                    nd = self.tree.mrca(taxa=[leaf_of[x].taxon for x in labels])
                else:
                    nd = leaf_of[labels[0]]
                nd.red = float(red_value)

        self._reference_genomes(taxonomy_file)

    def _reference_genomes(self, taxonomy_file):
        """Bottom-up: reference genomes below each node (count, first in tree order, name per
        rank) and the first labelled node per rank in pre-order; top-down: the closest labelled
        node at or above each node."""
        taxonomy = Taxonomy().read(taxonomy_file)
        for nd in self.tree.postorder_node_iter():
            _support, nd.taxon_label, _aux_info = parse_label(nd.label)
            if nd.is_leaf():
                tax = taxonomy.get(nd.taxon.label.replace("'", ''))
                nd.n_refs = 1 if tax else 0
                nd.first_ref = nd.taxon.label.replace("'", '') if tax else None
                nd.rank_names = list(tax) if tax else None
                nd.first_labelled = {}
                continue
            kids = nd.child_nodes()
            with_refs = [k for k in kids if k.n_refs]
            nd.n_refs = sum(k.n_refs for k in kids)
            nd.first_ref = with_refs[0].first_ref if with_refs else None
            nd.rank_names = None
            if with_refs:
                nd.rank_names = [names[0] if len(set(names)) == 1 else MIXED
                                 for names in zip(*(k.rank_names for k in with_refs))]
            # first internal node in pre-order (this one, then the children in order)
            # whose label starts with the rank, as Classify's preorder_iter search
            nd.first_labelled = {}
            for rank in ORDER_RANK:
                if nd.taxon_label and nd.taxon_label.startswith(rank):
                    nd.first_labelled[rank] = nd
                else:
                    for k in kids:
                        if rank in k.first_labelled:
                            nd.first_labelled[rank] = k.first_labelled[rank]
                            break
        for nd in self.tree.preorder_node_iter():
            if nd.taxon_label:
                nd.labelled = nd
            else:
                nd.labelled = nd.parent_node.labelled if nd.parent_node else None

    def _red_row(self, node):
        """(RED of C, RED of P, path C -> edge, path edge -> P); the path below is None when C
        is the edge's own lower node."""
        # C: first node with a RED in level order from the edge's lower node (as
        # _assign_mrca_red's levelorder_iter; a query leaf never has one)
        low = next(nd for nd in node.levelorder_iter() if nd.red is not None)
        below = None
        if low is not node:
            below, cur = max(low.edge_length or 0.0, 1e-9), low.parent_node
            while cur is not node:
                below += cur.edge_length
                cur = cur.parent_node
        # P: first node with a RED above the edge
        above, cur = 0.0, node.parent_node
        while cur.red is None:
            above += cur.edge_length
            cur = cur.parent_node
        return low.red, cur.red, below, above

    @staticmethod
    def _taxonomy(node):
        """Labels of the nodes above the edge, root excluded (as _get_gplace_taxonomy)."""
        taxa, cur = [], node.parent_node
        while cur.parent_node:
            if cur.taxon_label:
                for t in cur.taxon_label.split(';')[::-1]:
                    taxa.append(t.strip())
            cur = cur.parent_node
        return ';'.join(taxa[::-1])

    @staticmethod
    def _topology(node):
        """(parent_taxon, parent_red, child_taxon, child_red) as _classify_red_topology finds
        them for a genome placed on the edge above node."""
        parent = node.parent_node.labelled
        if parent is None:
            return '', None, '', None
        parent_rank = parent.taxon_label.split(';')[-1][0:3]
        child = None
        if parent_rank in ORDER_RANK[:-2]:                      # a child rank exists (not g__ / s__)
            child = node.first_labelled.get(ORDER_RANK[ORDER_RANK.index(parent_rank) + 1])
        return (parent.taxon_label, parent.red,
                child.taxon_label if child else '', child.red if child else None)

    def write(self, path):
        """Write the table (gzipped tsv; first lines: fingerprints of the tree, RED and
        taxonomy files)."""
        def _str(v):
            return '' if v is None else str(v)

        with gzip.open(path, 'wt') as fh:
            fh.write(f'# tree_sha1\t{tree_sha1(self.jplace_tree)}\n')
            fh.write(f'# red_sha1\t{file_sha1(self.red_file)}\n')
            fh.write(f'# taxonomy_sha1\t{file_sha1(self.taxonomy_file)}\n')
            fh.write('\t'.join(TABLE_HEADER) + '\n')
            for edge_num, node in sorted(self.node_of_edge.items()):
                if node.parent_node is None:
                    continue                                      # the root has no edge to place on
                red_below, red_above, below, above = self._red_row(node)
                parent_taxon, parent_red, child_taxon, child_red = self._topology(node)
                fh.write('\t'.join([str(edge_num), str(red_below), str(red_above), _str(below), str(above),
                                    str(node.edge_length or 0.0), self._taxonomy(node),
                                    parent_taxon, _str(parent_red), child_taxon, _str(child_red),
                                    str(node.n_refs), _str(node.first_ref),
                                    ';'.join(node.rank_names) if node.rank_names else '']) + '\n')


class JplaceTable(object):
    """Read a table written by JplaceTableBuilder (rows are parsed when an edge is used)."""

    def __init__(self, path):
        self.path = path
        self.fingerprints = {}
        self.lines = {}
        with gzip.open(path, 'rt') as fh:
            line = fh.readline()
            while line.startswith('# '):                            # '# name<tab>sha1'
                name, value = line[2:].rstrip('\n').split('\t')
                self.fingerprints[name] = value
                line = fh.readline()
            if line.rstrip('\n').split('\t') != TABLE_HEADER or set(FINGERPRINTS) - set(self.fingerprints):
                raise GTDBTkExit(f'Unexpected header in {path}: rebuild it with scripts/build_jplace_table.py.')
            for line in fh:
                self.lines[int(line[:line.index('\t')])] = line
        self.rows = {}

    @property
    def edges(self):
        return self.lines

    def check(self, jplace_tree, red_file, taxonomy_file):
        """Stop if a jplace was not placed on the tree, or the RED / taxonomy file is not the
        one, this table was built from."""
        if tree_sha1(jplace_tree) != self.fingerprints['tree_sha1']:
            raise GTDBTkExit(f'The jplace reference tree does not match {self.path}: '
                             f'rebuild it with scripts/build_jplace_table.py.')
        for name, path in (('red_sha1', red_file), ('taxonomy_sha1', taxonomy_file)):
            if file_sha1(path) != self.fingerprints[name]:
                raise GTDBTkExit(f'{path} does not match {self.path}: '
                                 f'rebuild it with scripts/build_jplace_table.py.')

    def row(self, edge_num):
        """The parsed row of an edge, as a dict."""
        if edge_num not in self.rows:
            values = self.lines[edge_num].rstrip('\n').split('\t')
            row = dict(zip(TABLE_HEADER, values))
            for col in ('red_below', 'red_above', 'path_above', 'edge_length'):
                row[col] = float(row[col])
            for col in ('path_below', 'parent_red', 'child_red'):
                row[col] = None if row[col] == '' else float(row[col])
            row['n_refs'] = int(row['n_refs'])
            row['first_ref'] = row['first_ref'] or None
            row['rank_names'] = row['rank_names'].split(';') if row['rank_names'] else []
            self.rows[edge_num] = row
        return self.rows[edge_num]

    def red(self, edge_num, distal_length):
        """RED of the node a genome placed at (edge_num, distal_length) is attached to."""
        row = self.row(edge_num)
        length, above, below = row['edge_length'], row['path_above'], row['path_below']
        up = (length - distal_length) + above                       # attachment point -> P
        if below is None:
            # C -> P as _assign_mrca_red sums it on the tog tree: C's edge there is the part below
            # the genome, floored at 1e-9, so a genome placed on C (distal 0) gets a RED slightly
            # below C's instead of equal to it (ties in the RED rules go the same way as with the tree)
            total = max(distal_length, 1e-9) + (length - distal_length) + above
        else:
            total = below + length + above
        return row['red_above'] + (row['red_below'] - row['red_above']) * (up / total)

    def placement_taxonomy(self, edge_num):
        """Taxonomy from the labels above the edge (not standardised)."""
        return self.row(edge_num)['placement_taxonomy']


class TableEdge(object):
    """What the RED and topology rules read for one placed genome, from the jplace table
    (same interface as gtdbtk.classify.TreeEdge)."""

    def __init__(self, table, edge_num, distal_length, gtdb_taxonomy):
        self.row = table.row(edge_num)
        self.gtdb_taxonomy = gtdb_taxonomy
        self.red = table.red(edge_num, distal_length)
        self.placement_taxonomy = self.row['placement_taxonomy']
        self.parent_taxon = self.row['parent_taxon']
        self.parent_red = self.row['parent_red']

    def n_refs(self):
        return self.row['n_refs']

    def first_ref(self):
        return self.row['first_ref']

    def single_ref(self):
        """The reference genome below the edge if there is only one (Classify._get_gplace_reference)."""
        return self.row['first_ref'] if self.row['n_refs'] == 1 else None

    def unique_name(self, rank):
        """The name all reference genomes below have for a rank, or None."""
        name = self.row['rank_names'][ORDER_RANK.index(rank)]
        return None if name == MIXED else name

    def child_taxon(self, child_rk):
        """(taxon, RED) of the first internal node labelled with child_rk below the edge; the
        table holds it for the rank after the parent taxon's, the only one the rules ask for."""
        if self.row['child_taxon'] and self.row['child_taxon'].startswith(child_rk):
            return self.row['child_taxon'], self.row['child_red']
        return None, None

    def redtax(self, closest_rank):
        """As Classify._get_redtax on the reference genomes below: the first genome's ranks
        above closest_rank, then the name they all share at that rank, or '<rank>undefined'."""
        rank = closest_rank[0:3]
        subtax = list(self.gtdb_taxonomy.get(self.row['first_ref'])[:ORDER_RANK.index(rank)])
        name = self.unique_name(rank)
        if name is not None and name.startswith(closest_rank):
            subtax.append(name)
        else:
            subtax.append(closest_rank + "undefined")
        return subtax
