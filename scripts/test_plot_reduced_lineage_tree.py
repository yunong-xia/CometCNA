import tempfile
from pathlib import Path
import unittest

from plot_reduced_lineage_tree import read_lineages, reduce_tree, plot_tree


class ReducedTreeTests(unittest.TestCase):
    def test_suppresses_unary_nodes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lineages.tsv"
            path.write_text(
                "sample_id\tlineage_index\tid\tancestor\n"
                "6\t0\t6\t5\n6\t1\t5\t4\n6\t2\t4\t3\n6\t3\t3\t1\n6\t4\t1\t0\n"
                "7\t0\t7\t5\n7\t1\t5\t4\n7\t2\t4\t3\n7\t3\t3\t1\n7\t4\t1\t0\n"
                "8\t0\t8\t2\n8\t1\t2\t1\n8\t2\t1\t0\n"
            )
            parent, sampled = read_lineages(path)
            keep, children = reduce_tree(parent, sampled)
            self.assertEqual(keep, {1, 5, 6, 7, 8})
            self.assertEqual(children[1], {5, 8})
            self.assertEqual(children[5], {6, 7})
            self.assertNotIn(3, keep)
            self.assertNotIn(4, keep)
            output = Path(directory) / "tree.svg"
            nodes, branches, leaves = plot_tree(parent, sampled, output)
            self.assertTrue(output.exists())
            self.assertEqual((nodes, branches, leaves), (5, 2, 3))


if __name__ == "__main__":
    unittest.main()
