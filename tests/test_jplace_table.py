import gzip
import json
import os
import random
import shutil
import tempfile
import unittest
from unittest import mock

from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.jplace_table import JplaceTable, JplaceTableBuilder, best_placements, refpkg_jplace_tree

# ((R1,R2)C,R3)P: C (RED 0.6, label g__X) 0.4 below P (RED 0.0, label f__F); post-order edges
TREE = "(((R1:0.3{0},R2:0.2{1})'0.9:g__X':0.4{2},R3:0.5{3})'f__F':0.1{4},R4:0.2{5}){6};"
RED = ("R1\t1.0\nR2\t1.0\nR3\t1.0\nR4\t1.0\nR1|R2\t0.6\nR1|R3\t0.0\nR1|R4\t0.0\n")
TAXONOMY = ("R1\td__Bacteria;p__P1;c__C1;o__O1;f__F;g__X;s__X a\n"
            "R2\td__Bacteria;p__P1;c__C1;o__O1;f__F;g__X;s__X b\n"
            "R3\td__Bacteria;p__P1;c__C1;o__O1;f__F;g__Y;s__Y a\n"
            "R4\td__Bacteria;p__P2;c__C2;o__O2;f__G;g__Z;s__Z a\n")
FIELDS = ['distal_length', 'edge_num', 'like_weight_ratio', 'likelihood', 'pendant_length']


class TestJplaceTable(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        self.red_file = os.path.join(self.tmp, 'red.tsv')
        self.taxonomy_file = os.path.join(self.tmp, 'taxonomy.tsv')
        with open(self.taxonomy_file, 'w') as fh:
            fh.write(TAXONOMY)
        self.table_file = os.path.join(self.tmp, 'table.tsv.gz')

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _table(self, red=RED):
        with open(self.red_file, 'w') as fh:
            fh.write(red)
        JplaceTableBuilder(TREE, self.red_file, self.taxonomy_file).write(self.table_file)
        return JplaceTable(self.table_file)

    def test_reproducible(self):
        # same inputs, different build time: same bytes (the mrca_red hash depends on them)
        self._table()
        other = os.path.join(self.tmp, 'other.tsv.gz')
        os.rename(self.table_file, other)
        with mock.patch('time.time', return_value=1e9):
            JplaceTableBuilder(TREE, self.red_file, self.taxonomy_file).write(self.table_file)
        with open(self.table_file, 'rb') as a, open(other, 'rb') as b:
            self.assertEqual(a.read(), b.read())

    def test_red(self):
        # 0.1 above C on the 0.4 edge: 0.0 + (0.6 - 0.0) * 0.3 / 0.4
        table = self._table()
        self.assertAlmostEqual(table.red(2, 0.1), 0.45)
        self.assertAlmostEqual(table.red(2, 0.0), 0.6)                # at C
        self.assertAlmostEqual(table.red(2, 0.4), 0.0)                # at P

    def test_red_missing_node(self):
        # C has no RED: the first reference node below in level order (R1) is used
        table = self._table(RED.replace('R1|R2\t0.6\n', ''))
        self.assertAlmostEqual(table.red(2, 0.1), 0.0 + 1.0 * 0.3 / (0.3 + 0.4))

    def test_placement_taxonomy(self):
        table = self._table()
        self.assertEqual(table.placement_taxonomy(2), 'f__F')         # above C, root excluded
        self.assertEqual(table.placement_taxonomy(0), 'f__F;g__X')
        self.assertEqual(table.placement_taxonomy(5), '')

    def test_check(self):
        table = self._table()
        table.check(TREE, self.red_file, self.taxonomy_file)
        with self.assertRaises(GTDBTkExit):                         # another tree
            table.check(TREE.replace('0.3{0}', '0.31{0}'), self.red_file, self.taxonomy_file)
        with open(self.taxonomy_file, 'a') as fh:                  # the taxonomy changed since the build
            fh.write('R5\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G a\n')
        with self.assertRaises(GTDBTkExit):
            table.check(TREE, self.red_file, self.taxonomy_file)
        with open(self.red_file, 'a') as fh:                       # the RED file changed since the build
            fh.write('R2|R3\t0.1\n')
        with self.assertRaises(GTDBTkExit):
            table.check(TREE, self.red_file, self.taxonomy_file)

    def test_old_table(self):
        # a table without the RED fingerprint must be rebuilt
        with gzip.open(self.table_file, 'wt') as fh:
            fh.write('# tree_sha1\tabc\n# red_sha1\tabc\n' + '\t'.join(['edge_num', 'red_below', 'red_above',
                     'path_below', 'path_above', 'edge_length', 'placement_taxonomy']) + '\n')
        with self.assertRaises(GTDBTkExit):
            JplaceTable(self.table_file)

    def test_best_placements(self):
        # highest like_weight_ratio, the first one on ties (as gplace tog)
        path = os.path.join(self.tmp, 'g.jplace')
        with open(path, 'w') as fh:
            json.dump({'tree': TREE, 'version': 3, 'fields': FIELDS, 'metadata': {'invocation': 't'},
                       'placements': [{'p': [[0.1, 3, 0.4, -2, 0.1], [0.1, 2, 0.6, -1, 0.1],
                                             [0.2, 0, 0.6, -1, 0.1]], 'nm': [['QA', 1], ['QB', 1]]}]}, fh)
        tree, best = best_placements(path)
        self.assertEqual(tree, TREE)
        self.assertEqual(best, {'QA': (2, 0.1), 'QB': (2, 0.1)})


@unittest.skipUnless(shutil.which('gplace'), 'gplace not installed')
class TestRefpkgTree(unittest.TestCase):
    """The tree built from the reference package is the tree gplace writes in its jplace."""

    def test_same_as_jplace(self):
        import subprocess
        import gplace
        data = os.path.join(os.path.dirname(gplace.__file__), 'tests', 'data')
        refpkg = os.path.join(data, 'mini.refpkg')
        tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        try:
            for args in ([], ['--tile_edges', '100']):
                jplace = os.path.join(tmp, 'q.jplace')
                subprocess.run(['gplace', 'place', '-c', refpkg, '-o', jplace, *args,
                                os.path.join(data, 'queries.fasta')], check=True, capture_output=True)
                tree, _best = best_placements(jplace)
                self.assertEqual(refpkg_jplace_tree(refpkg), tree)
        finally:
            shutil.rmtree(tmp)


@unittest.skipUnless(shutil.which('gplace'), 'gplace not installed')
class TestJplaceTableMatchesTree(unittest.TestCase):
    """Same RED and placement taxonomy as gplace tog + Classify._assign_mrca_red and
    Classify._get_gplace_taxonomy, on random labelled trees with a complete RED file
    (RED up to tog's 6 significant digits)."""

    def test_random_trees(self):
        import logging
        import subprocess
        import types
        import gtdbtk.classify as classify_mod
        from gtdbtk.classify import Classify
        from gtdbtk.config.common import CONFIG
        from gtdbtk.tools import standardise_taxonomy

        class _ClassifyFile(object):
            def __init__(self):
                self.data = {}

            def add_genome(self, gid, tax):
                self.data[gid] = tax

            def write(self):
                pass

        tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        stub = types.SimpleNamespace(logger=logging.getLogger('test'))
        try:
            for seed in range(20):
                rng = random.Random(seed)
                tree, lengths, red_lines = _random_jplace_tree(rng, rng.randint(3, 40))
                red_file = os.path.join(tmp, 'red.tsv')
                with open(red_file, 'w') as fh:
                    fh.write('\n'.join(red_lines) + '\n')
                taxonomy_file = os.path.join(tmp, 'taxonomy.tsv')
                with open(taxonomy_file, 'w') as fh:
                    for x in red_lines:
                        if '|' not in x:
                            r = x.split('\t')[0]
                            fh.write(f'{r}\td__Bacteria;p__{r};c__{r};o__{r};f__{r};g__{r};s__{r} sp\n')
                edges = [e for e, length in lengths.items() if length is not None]
                placements = []
                for q in range(rng.randint(1, 15)):
                    e = rng.choice(edges)
                    d = rng.choice([rng.uniform(0, lengths[e]), 0.0, lengths[e]])
                    placements.append({'p': [[round(d, 9), e, 1.0, -1.0, 0.1]], 'nm': [[f'Q{q}', 1]]})
                    placements.append({'p': [[round(d, 9), e, 1.0, -1.0, 0.2]], 'nm': [[f'T{q}', 1]]})  # tie
                jplace = os.path.join(tmp, 't.jplace')
                with open(jplace, 'w') as fh:
                    json.dump({'tree': tree, 'version': 3, 'fields': FIELDS, 'placements': placements,
                               'metadata': {'invocation': 'test'}}, fh)
                tree_file = os.path.join(tmp, 't.tree')
                subprocess.run(['gplace', 'tog', '-o', tree_file, jplace], check=True, capture_output=True)

                # tree path
                classify_mod.CONFIG = types.SimpleNamespace(MRCA_RED_BAC120=red_file, MRCA_RED_AR53=red_file)
                t = Classify._assign_mrca_red(stub, tree_file, 'bac120')
                msa = os.path.join(tmp, 'user.fasta')
                with open(msa, 'w') as fh:
                    fh.write(''.join(f'>{p["nm"][0][0]}\nA\n' for p in placements))
                tree_tax = Classify._get_gplace_taxonomy(stub, _ClassifyFile(), 'bac120', msa, t)

                # jplace path
                table_file = os.path.join(tmp, 'table.tsv.gz')
                jplace_tree, best = best_placements(jplace)
                JplaceTableBuilder(jplace_tree, red_file, taxonomy_file).write(table_file)
                table = JplaceTable(table_file)
                table.check(jplace_tree, red_file, taxonomy_file)

                refs = {x.split('\t')[0] for x in red_lines if '|' not in x}
                for leaf in t.leaf_node_iter():
                    gid = leaf.taxon.label
                    if gid in refs:
                        continue
                    cur = leaf.parent_node
                    while not {x.taxon.label for x in cur.leaf_iter()} & refs:
                        cur = cur.parent_node
                    self.assertAlmostEqual(table.red(*best[gid]), cur.rel_dist, delta=1e-5)
                    self.assertEqual(standardise_taxonomy(table.placement_taxonomy(best[gid][0]), 'bac120'),
                                     tree_tax[gid])
        finally:
            shutil.rmtree(tmp)
            classify_mod.CONFIG = CONFIG


def _random_jplace_tree(rng, n_leaves):
    """Random rooted binary jplace tree (some zero-length edges, some labelled nodes, some with
    several ranks) and a complete RED file."""
    nodes = [(f'R{i}', []) for i in range(n_leaves)]
    while len(nodes) > 1:
        a, b = nodes.pop(rng.randrange(len(nodes))), nodes.pop(rng.randrange(len(nodes)))
        nodes.append((None, [a, b]))
    counter, lengths, red_lines = [0], {}, []
    ranks = ['p__', 'c__', 'o__', 'f__', 'g__']

    def rec(nd, is_root=False, depth=0):
        name, kids = nd
        if kids:
            parts = [rec(k, depth=depth + 1) for k in kids]
            s, leaves = '(' + ','.join(p[0] for p in parts) + ')', parts[0][1] + parts[1][1]
            red_lines.append(f'{parts[0][1][0]}|{parts[1][1][0]}\t{0.0 if is_root else round(rng.random(), 6)}')
            if depth < len(ranks) and rng.random() < 0.5:
                label = ranks[depth] + f'T{counter[0]}'
                if depth + 1 < len(ranks) and rng.random() < 0.3:
                    label += ';' + ranks[depth + 1] + f'U{counter[0]}'
                s += f"'{rng.choice(['', '0.95:'])}{label}'"
        else:
            s, leaves = name, [name]
            red_lines.append(f'{name}\t1.0')
        e, counter[0] = counter[0], counter[0] + 1
        length = None if is_root else (0.0 if rng.random() < 0.08 else round(rng.uniform(0.001, 0.5), 6))
        lengths[e] = length
        return s + ('' if is_root else f':{length}') + '{%d}' % e, leaves

    tree, _leaves = rec(nodes[0], True)
    return tree + ';', lengths, red_lines


if __name__ == '__main__':
    unittest.main()
