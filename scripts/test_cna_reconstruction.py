#!/usr/bin/env python3
"""single cell lineage tracing."""
import argparse
import csv
import gzip
import pandas as pd
from collections import defaultdict
from pathlib import Path
from bisect import bisect_left
import numpy as np
import tempfile
import unittest


import single_cell_lineage_tracing as sclt


class ComparisonTests(unittest.TestCase):
    def test_cna_delta(self):
        cn_states = {"1p":[2]}
        bps = {"1p":[0,123400000]}
        trajectory_types = {"1p":["UDT"]}
        trajectory_cn_at_least_2 = {"1p":[True]}

        sclt.apply_delta(arm="1p", start=100, end=1000, delta=1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types, trajectory_cn_at_least_2=trajectory_cn_at_least_2)
        self.assertEqual(cn_states,
                         {"1p":[2,3,2]})
        self.assertEqual(bps,
                         {"1p":[0,100,1000,123400000]})
        self.assertEqual(trajectory_types,
                         {"1p":["UDT","UDT","UDT"]})
        self.assertEqual(trajectory_cn_at_least_2,
                         {"1p":[True,True,True]})

        sclt.apply_delta(arm="1p", start=50, end=500, delta=1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types, trajectory_cn_at_least_2=trajectory_cn_at_least_2)
        self.assertEqual(cn_states,
                            {"1p":[2,3,4,3,2]})
        self.assertEqual(bps,
                            {"1p":[0,50,100,500,1000,123400000]})
        self.assertEqual(trajectory_types,
                         {"1p":["UDT","UDT","UDT","UDT","UDT"]})
        self.assertEqual(trajectory_cn_at_least_2,
                         {"1p":[True,True,True,True,True]})

        sclt.apply_delta(arm="1p", start=100, end=500, delta=-1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types, trajectory_cn_at_least_2=trajectory_cn_at_least_2)
        self.assertEqual(cn_states,
                         {"1p":[2,3,3,3,2]})
        self.assertEqual(trajectory_types,
                         {"1p":["UDT","UDT","BDT","UDT","UDT"]})
        self.assertEqual(trajectory_cn_at_least_2,
                         {"1p":[True,True,True,True,True]})

    def test_cn_below_2_flag_persists_after_gain(self):
        cn_states = {"1p":[2]}
        bps = {"1p":[0,1000]}
        trajectory_types = {"1p":["UDT"]}
        trajectory_cn_at_least_2 = {"1p":[True]}

        sclt.apply_delta(arm="1p", start=100, end=500, delta=-1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types, trajectory_cn_at_least_2=trajectory_cn_at_least_2)

        self.assertEqual(cn_states["1p"], [2,1,2])
        self.assertEqual(trajectory_cn_at_least_2["1p"], [True,False,True])

        sclt.apply_delta(arm="1p", start=100, end=500, delta=1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types, trajectory_cn_at_least_2=trajectory_cn_at_least_2)

        self.assertEqual(cn_states["1p"], [2,2,2])
        self.assertEqual(trajectory_cn_at_least_2["1p"], [True,False,True])

    def test_cn_at_least_2_initializes_from_starting_cn(self):
        cn_states = {"1p":[1,2]}
        bps = {"1p":[0,100,1000]}

        _, _, _, trajectory_cn_at_least_2 = sclt.reconstruct_copy_numbers(
            [], starting_cn_states=cn_states, starting_bps=bps
        )

        self.assertEqual(trajectory_cn_at_least_2["1p"], [False,True])


if __name__ == '__main__':
    unittest.main()
