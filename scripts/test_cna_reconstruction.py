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

        sclt.apply_delta(arm="1p", start=100, end=1000, delta=1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types)
        self.assertEqual(cn_states,
                         {"1p":[2,3,2]})
        self.assertEqual(bps,
                         {"1p":[0,100,1000,123400000]})
        self.assertEqual(trajectory_types,
                         {"1p":["UDT","UDT","UDT"]})

        sclt.apply_delta(arm="1p", start=50, end=500, delta=1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types)
        self.assertEqual(cn_states,
                            {"1p":[2,3,4,3,2]})
        self.assertEqual(bps,
                            {"1p":[0,50,100,500,1000,123400000]})
        self.assertEqual(trajectory_types,
                         {"1p":["UDT","UDT","UDT","UDT","UDT"]})

        sclt.apply_delta(arm="1p", start=100, end=500, delta=-1, cn_states=cn_states, bps=bps, trajectory_types=trajectory_types)
        self.assertEqual(cn_states,
                         {"1p":[2,3,3,3,2]})
        self.assertEqual(trajectory_types,
                         {"1p":["UDT","UDT","BDT","UDT","UDT"]})


if __name__ == '__main__':
    unittest.main()
