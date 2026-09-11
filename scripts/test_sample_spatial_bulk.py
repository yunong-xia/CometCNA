import csv
import gzip
from pathlib import Path
import tempfile
import unittest

import numpy as np

from sample_spatial_bulk import load_population, main, sample_bulk


class SpatialSamplingTests(unittest.TestCase):
    def test_matches_brute_force(self):
        coordinates = np.random.default_rng(12).normal(size=(200, 3))
        ids = np.arange(1, 201)
        center = [0.2, -0.4, 0.7]
        selected, distances, _, center_id = sample_bulk(ids, coordinates, 17, center)
        expected_distances = np.linalg.norm(coordinates - center, axis=1)
        expected = np.argsort(expected_distances)[:17]
        np.testing.assert_array_equal(selected, expected)
        np.testing.assert_allclose(distances, expected_distances[expected])
        self.assertIsNone(center_id)

    def test_seeded_center_and_full_population(self):
        ids = np.arange(1, 21)
        coordinates = np.column_stack((ids, ids * 0, ids * 0))
        first = sample_bulk(ids, coordinates, 20, seed=9)
        second = sample_bulk(ids, coordinates, 20, seed=9)
        np.testing.assert_array_equal(first[0], second[0])
        self.assertEqual(first[3], second[3])
        self.assertEqual(ids[first[0][0]], first[3])
        self.assertEqual(first[1][0], 0)
        self.assertEqual(set(ids[first[0]]), set(ids))

    def test_ties_have_correct_distances(self):
        coords = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, 0, 2]])
        indices, distances, _, _ = sample_bulk(np.arange(1, 5), coords, 2, [0, 0, 0])
        self.assertTrue(set(indices).issubset({0, 1, 2}))
        np.testing.assert_array_equal(distances, [1, 1])

    def test_plain_and_gzip_cli_excludes_ancestors(self):
        content = ('x\ty\tz\tid\tancestor\tbirth\tdeath\tomega\tNt\n'
                   '0\t0\t0\t1\t0\t0\t1\t10\t1\n'
                   '1\t0\t0\t2\t1\t1\t0\t10\t2\n'
                   '3\t0\t0\t3\t1\t1\t0.0\t10\t2\n'
                   '2\t0\t0\t4\t1\t1\t0\t10\t3\n')
        with tempfile.TemporaryDirectory() as directory:
            for suffix in ('.tsv', '.tsv.gz'):
                path = Path(directory) / ('population' + suffix)
                opener = gzip.open if suffix.endswith('.gz') else open
                with opener(path, 'wt') as handle:
                    handle.write(content)
                output = Path(directory) / 'bulk.tsv'
                main([str(path), '-k', '2', '--center', '0', '0', '0', '-o', str(output)])
                with output.open() as handle:
                    rows = list(csv.DictReader(handle, delimiter='\t'))
                self.assertEqual([row['id'] for row in rows], ['2', '4'])
                self.assertEqual([float(row['distance']) for row in rows], [1, 2])

    def test_invalid_inputs(self):
        for size in (0, -1, 3):
            with self.assertRaises(ValueError):
                sample_bulk(np.array([1, 2]), np.zeros((2, 3)), size)
        with self.assertRaises(ValueError):
            sample_bulk(np.array([1]), np.zeros((1, 3)), 1, [float('nan'), 0, 0])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'population.tsv'
            for content in ('id\n', 'id\tdeath\tx\ty\tz\n',
                            'id\tdeath\tx\ty\tz\n1\t0\tbad\t0\t0\n'):
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_population(path)


if __name__ == '__main__':
    unittest.main()
