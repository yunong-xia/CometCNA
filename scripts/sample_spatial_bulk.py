#!/usr/bin/env python3
"""Select a spatial bulk of extant cells using a Euclidean KD-tree."""

import argparse
import csv
import gzip
import math
from pathlib import Path
import sys

import numpy as np
from sklearn.neighbors import NearestNeighbors


def extant_cells(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'id', 'death', 'x', 'y', 'z'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Population header must contain id, death, x, y, z')
        for row in reader:
            try:
                death = float(row['death'])
                if not math.isfinite(death):
                    raise ValueError('nonfinite death time')
                if death != 0:
                    continue
                cell_id = int(row['id'])
                coords = tuple(float(row[axis]) for axis in ('x', 'y', 'z'))
                if cell_id <= 0 or not all(map(math.isfinite, coords)):
                    raise ValueError('invalid ID or coordinates')
            except (ValueError, TypeError, KeyError) as exc:
                raise ValueError(f'Invalid population row at line {reader.line_num}') from exc
            yield cell_id, coords


def load_population(path):
    # Two passes avoid retaining Python objects for millions of population rows.
    count = sum(1 for _ in extant_cells(path))
    if not count:
        raise ValueError('Population contains no extant cells')
    ids = np.empty(count, dtype=np.int64)
    coordinates = np.empty((count, 3), dtype=np.float64)
    loaded = 0
    for cell_id, coords in extant_cells(path):
        if loaded >= count:
            raise ValueError('Population changed between reads')
        ids[loaded] = cell_id
        coordinates[loaded] = coords
        loaded += 1
    if loaded != count:
        raise ValueError('Population changed between reads')
    return ids, coordinates


def sample_bulk(ids, coordinates, sample_size, center=None, seed=None):
    if not 1 <= sample_size <= len(ids):
        raise ValueError(f'Sample size must be between 1 and {len(ids)}')
    center_id = None
    if center is None:
        index = np.random.default_rng(seed).integers(len(ids))
        center = coordinates[index].copy()
        center_id = int(ids[index])
    center = np.asarray(center, dtype=np.float64)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError('Center must contain three finite coordinates')
    model = NearestNeighbors(algorithm='kd_tree', metric='euclidean', n_jobs=1)
    model.fit(coordinates)
    distances, indices = model.kneighbors(center.reshape(1, -1), n_neighbors=sample_size)
    # Order selected ties by ID; membership at a tied boundary is sklearn's choice.
    order = np.lexsort((ids[indices[0]], distances[0]))
    return indices[0][order], distances[0][order], center, center_id


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('population', help='Population TSV or TSV.GZ')
    parser.add_argument('-k', '--sample-size', type=int, default=1000)
    parser.add_argument('--center', nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                        help='Default: coordinates of a uniformly sampled extant cell')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('-o', '--output', required=True, help='Output TSV path')
    args = parser.parse_args(argv)
    if Path(args.population).resolve() == Path(args.output).resolve():
        parser.error('Output must differ from population input')
    if args.sample_size < 1:
        parser.error('Sample size must be positive')
    try:
        ids, coordinates = load_population(args.population)
        indices, distances, center, center_id = sample_bulk(
            ids, coordinates, args.sample_size, args.center, args.seed)
        with open(args.output, 'w', encoding='ascii', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['id', 'x', 'y', 'z', 'distance', 'center_id',
                             'center_x', 'center_y', 'center_z'])
            for index, distance in zip(indices, distances):
                writer.writerow([int(ids[index]), *coordinates[index], distance,
                                 '' if center_id is None else center_id, *center])
    except (ValueError, OSError, OverflowError) as exc:
        parser.error(str(exc))
    print(f'extant={len(ids)} sampled={len(indices)} center={center.tolist()} '
          f'radius={distances[-1]:.6g} output={args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
