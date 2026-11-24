import json
import logging
import os
import random
import shutil
import subprocess
import tempfile
import types
import unittest

import gtdbtk.classify as classify_mod
from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.common import CONFIG
from gtdbtk.exceptions import GTDBTkExit
from gtdbtk.jplace_table import JplaceTable, JplaceTableBuilder, TableEdge, best_placements
from gtdbtk.tools import TreeTraversal

ORDER_RANK = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
RED_DICT = {'d__': 0.0, 'p__': 0.33, 'c__': 0.46, 'o__': 0.62, 'f__': 0.77, 'g__': 0.92, 's__': 1.0}
FIELDS = ['distal_length', 'edge_num', 'like_weight_ratio', 'likelihood', 'pendant_length']


class _SummaryFile(object):
    def __init__(self):
        self.rows = {}

    def add_row(self, row):
        self.rows[row.gid] = row


def random_gtdb_tree(rng, n_leaves):
    """Random rooted binary jplace tree labelled like a GTDB tree (single and multi-rank
    labels, ranks increasing along each path), a taxonomy consistent with the labels
    (unlabelled ranks are specific to each genome), REDs increasing with depth."""
    nodes = [[f'R{i}', [], None] for i in range(n_leaves)]
    while len(nodes) > 1:
        a, b = nodes.pop(rng.randrange(len(nodes))), nodes.pop(rng.randrange(len(nodes)))
        nodes.append([None, [a, b], None])
    root = nodes[0]

    # top-down: labels, inherited taxonomy, RED
    taxonomy, counter = {}, [0]
    stack = [(root, 0, {'d__': 'd__Bacteria'}, 0.0, True)]
    while stack:
        nd, depth, names, parent_red, is_root = stack.pop()
        name, kids, _ = nd
        if kids:
            label = None
            if is_root:
                label, red = 'd__Bacteria', 0.0
            else:
                red = round(parent_red + (1 - parent_red) * rng.uniform(0.05, 0.45), 6)
                deepest = max(ORDER_RANK.index(r) for r in names)
                if deepest < 5 and rng.random() < 0.5:
                    k = 2 if deepest < 4 and rng.random() < 0.25 else 1
                    new = [f'{ORDER_RANK[deepest + i + 1]}T{counter[0]}' for i in range(k)]
                    counter[0] += 1
                    names = dict(names, **{x[:3]: x for x in new})
                    label = ';'.join(new)
            nd[2] = (label, red)
            for k in kids:
                stack.append((k, depth + 1, names, red, False))
        else:
            taxonomy[name] = [names.get(r, f'{r}{name}') for r in ORDER_RANK[:-1]] + [f's__{name} sp']
            nd[2] = (None, 1.0)

    # newick with jplace edge numbers (post-order) + RED file
    counter, lengths, red_lines = [0], {}, []

    def rec(nd, is_root=False):
        name, kids, (label, red) = nd
        if kids:
            parts = [rec(k) for k in kids]
            s, leaves = '(' + ','.join(p[0] for p in parts) + ')', parts[0][1] + parts[1][1]
            red_lines.append(f'{parts[0][1][0]}|{parts[1][1][0]}\t{red}')
            if label:
                s += f"'{rng.choice(['', '0.95:'])}{label}'"
        else:
            s, leaves = name, [name]
            red_lines.append(f'{name}\t1.0')
        e, counter[0] = counter[0], counter[0] + 1
        length = None if is_root else (0.0 if rng.random() < 0.05 else round(rng.uniform(0.001, 0.4), 6))
        lengths[e] = length
        return s + ('' if is_root else f':{length}') + '{%d}' % e, leaves

    tree, _leaves = rec(root, True)
    return tree + ';', lengths, red_lines, taxonomy


@unittest.skipUnless(shutil.which('gplace'), 'gplace not installed')
class TestClassifyJplaceTable(unittest.TestCase):
    """_classify_red_topology and _get_gplace_reference give the same results from the
    jplace table (TableEdge) as from the gplace tog tree (TreeEdge)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix='gtdbtk_tmp_')

    def tearDown(self):
        shutil.rmtree(self.tmp)
        classify_mod.CONFIG = CONFIG

    def _classify(self, taxonomy_file):
        classify = Classify.__new__(Classify)                    # no reference data needed
        classify.order_rank = ORDER_RANK
        classify.gtdb_taxonomy = Taxonomy().read(taxonomy_file)
        classify.reference_ids = frozenset(classify.gtdb_taxonomy)
        classify.logger = logging.getLogger('timestamp')
        return classify

    def _run_rules(self, classify, edges, gids):
        summary = _SummaryFile()
        try:
            classify._classify_red_topology(edges, {g: 'A' for g in gids}, {}, {}, {}, RED_DICT, [], summary,
                                            {}, None, False, {})
        except GTDBTkExit as e:
            return 'GTDBTkExit'
        return {gid: (row.classification, row.classification_method, row.red_value)
                for gid, row in summary.rows.items()}

    def test_random_trees(self):
        tmp, stats = self.tmp, {'genomes': 0, 'methods': set(), 'red_diff': 0.0}
        for seed in range(60):
            rng = random.Random(seed)
            tree, lengths, red_lines, taxonomy = random_gtdb_tree(rng, rng.randint(4, 60))
            red_file, taxonomy_file = os.path.join(tmp, 'red.tsv'), os.path.join(tmp, 'taxonomy.tsv')
            with open(red_file, 'w') as fh:
                fh.write('\n'.join(red_lines) + '\n')
            with open(taxonomy_file, 'w') as fh:
                fh.write(''.join(f'{g}\t{";".join(t)}\n' for g, t in taxonomy.items()))

            edges = [e for e, length in lengths.items() if length is not None]
            placements = []
            for q in range(rng.randint(1, 20)):
                e = rng.choice(edges)
                d = rng.choice([rng.uniform(0, lengths[e]), 0.0, lengths[e]])
                placements.append({'p': [[round(d, 9), e, 1.0, -1.0, 0.1]], 'nm': [[f'Q{q}', 1]]})
                if rng.random() < 0.2:                               # another genome on the same spot
                    placements.append({'p': [[round(d, 9), e, 1.0, -1.0, 0.2]], 'nm': [[f'T{q}', 1]]})
            jplace = os.path.join(tmp, 't.jplace')
            with open(jplace, 'w') as fh:
                json.dump({'tree': tree, 'version': 3, 'fields': FIELDS, 'placements': placements,
                           'metadata': {'invocation': 'test'}}, fh)
            gids = [p['nm'][0][0] for p in placements]
            msa = os.path.join(tmp, 'user.fasta')
            with open(msa, 'w') as fh:
                fh.write(''.join(f'>{g}\nA\n' for g in gids))
            classify = self._classify(taxonomy_file)

            # tree path: tog tree + _assign_mrca_red + TreeEdge
            tree_file = os.path.join(tmp, 't.tree')
            subprocess.run(['gplace', 'tog', '-o', tree_file, jplace], check=True, capture_output=True)
            classify_mod.CONFIG = types.SimpleNamespace(MRCA_RED_BAC120=red_file, MRCA_RED_AR53=red_file)
            t = Classify._assign_mrca_red(classify, tree_file, 'bac120')
            tree_refs = Classify._get_gplace_reference(t, classify.reference_ids)
            tree_edges = classify._tree_edges(t, msa, {}, TreeTraversal())
            tree_result = self._run_rules(classify, tree_edges, gids)

            # table path: jplace + table + TableEdge
            table_file = os.path.join(tmp, 'table.tsv.gz')
            jplace_tree, best = best_placements(jplace)
            JplaceTableBuilder(jplace_tree, red_file, taxonomy_file).write(table_file)
            table = JplaceTable(table_file)
            table.check(jplace_tree, red_file, taxonomy_file)
            table_edges = {g: TableEdge(table, *best[g], classify.gtdb_taxonomy) for g in gids}
            table_refs = {g: {'gplace_g': e.single_ref()} for g, e in table_edges.items()}
            table_result = self._run_rules(classify, table_edges, gids)

            with self.subTest(seed=seed):
                self.assertEqual(table_refs, tree_refs)
                self.assertEqual(type(table_result), type(tree_result))
                if tree_result == 'GTDBTkExit':
                    self.assertEqual(table_result, tree_result)
                    continue
                self.assertEqual(set(table_result), set(tree_result))
                for gid in tree_result:
                    tax_t, method_t, red_t = tree_result[gid]
                    tax_j, method_j, red_j = table_result[gid]
                    self.assertEqual((tax_j, method_j), (tax_t, method_t), gid)
                    self.assertAlmostEqual(red_j, red_t, delta=1e-5)
                    stats['red_diff'] = max(stats['red_diff'], abs(red_j - red_t))
                    stats['methods'].add(method_t)
                stats['genomes'] += len(tree_result)
        print(f"\n  {stats['genomes']} genomes classified by RED/topology on 60 trees, "
              f"methods {sorted(stats['methods'])}, max RED diff {stats['red_diff']:.2g}")


if __name__ == '__main__':
    unittest.main()
