import io
import tempfile
from pathlib import Path
import unittest

from sample_extant_lineages import trace_lineages, read_sample_ids


class MRCATests(unittest.TestCase):
    def test_shared_ancestor_and_single_cell(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'cells.tsv'
            path.write_text('id\tancestor\tdeath\n7\t3\t0\n6\t3\t0\n3\t1\t2\n1\t0\t1\n')
            for samples, expected in [(['6', '7'], '3'), (['6'], '6')]:
                output, mrca = io.StringIO(), io.StringIO()
                count, unresolved = trace_lineages(str(path), samples, output, mrca)
                self.assertFalse(unresolved)
                self.assertEqual(count, 3 * len(samples))
                self.assertEqual(mrca.getvalue().splitlines()[1].split('\t')[1], expected)
            sample = Path(directory) / 'sample.tsv'
            sample.write_text('id\n6\n7\n')
            self.assertEqual(read_sample_ids(str(sample), str(path)), ['6', '7'])
            sample.write_text('id\n3\n')
            with self.assertRaises(SystemExit):
                read_sample_ids(str(sample), str(path))

    def test_missing_and_unsorted_ancestry(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'cells.tsv'
            for rows in ['7\t3\n', '3\t1\n7\t3\n1\t0\n', '7\t7\n']:
                path.write_text('id\tancestor\n' + rows)
                with self.assertRaises(SystemExit):
                    trace_lineages(str(path), ['7'], io.StringIO(), io.StringIO())


if __name__ == '__main__':
    unittest.main()
