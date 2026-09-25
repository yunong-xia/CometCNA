import argparse
import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from compare_bulk_seeding import run, lineage, load_parents


class ComparisonTests(unittest.TestCase):
    def test_pair_branches_exclude_shared_events(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            files = {
                'population': 'id\tancestor\tdeath\n8\t4\t0\n9\t4\t0\n4\t2\t3\n2\t1\t2\n1\t0\t1\n',
                'dead_cells': 'id\tancestor\tdeath\n7\t3\t4\n3\t2\t3\n2\t1\t2\n1\t0\t1\n',
                'bulk_ids': 'id\n8\n9\n',
                'seeding_ids': 'id\tNs\n7\t100\n2\t200\n',
                'cna': 'id\tcna_event\tchr\tarm\tstart\tend\n2\tshared\t1\tp\t0\t1\n4\tprimary\t1\tp\t0\t1\n3\tseed1\t1\tp\t0\t1\n7\tseed2\t1\tp\t0\t1\n7\tseed3\t2\tp\t0\t1\n8\tbulk_private\t1\tp\t0\t1\n',
            }
            args = {}
            for name, content in files.items():
                path = root / (name + '.tsv.gz')
                with gzip.open(path, 'wt') as handle:
                    handle.write(content)
                args[name] = str(path)
            args['outdir'] = str(root / 'out')
            run(argparse.Namespace(**args))
            def rows(name):
                with open(root / 'out' / name) as handle:
                    return list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(rows('bulk_mrca.tsv')[0]['id'], '4')
            pairs = rows('pair_mrcas.tsv')
            self.assertEqual([r['primary_metastasis_mrca_id'] for r in pairs], ['2', '2'])
            self.assertEqual(pairs[1]['seeding_branch_cells'], '0')
            events = rows('branch_cna_events.tsv')
            self.assertEqual([r['cna_event'] for r in events if r['seeding_id'] == '7'],
                             ['primary', 'seed1', 'seed2', 'seed3'])
            self.assertEqual([r['id'] for r in rows('bulk_mrca_lineage.tsv')], ['4', '2', '1'])

    def test_missing_ancestry(self):
        with self.assertRaisesRegex(ValueError, 'Missing ancestry'):
            lineage(7, {7: 3})

    def test_conflicting_parents(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'cells.tsv'
            path.write_text('id\tancestor\n4\t2\n4\t1\n')
            with self.assertRaisesRegex(ValueError, 'Conflicting'):
                load_parents([path])


if __name__ == '__main__':
    unittest.main()
