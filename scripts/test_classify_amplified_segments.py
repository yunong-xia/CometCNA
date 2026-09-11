import tempfile
from pathlib import Path
import unittest

from classify_amplified_segments import classify, read_events, read_lineages


class AmplifiedSegmentTests(unittest.TestCase):
    def test_status_and_route_type(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            lineages = directory / "lineages.tsv"
            lineages.write_text("sample_id\tid\n10\t4\n10\t3\n10\t2\n10\t1\n")
            cna = directory / "cna.tsv"
            cna.write_text(
                "id\tcna_event\tchr\tarm\tstart\tend\n"
                "2\tarm_gain\t8\t8q\t45200000\t145138636\n"
                "3\tfocal_loss\t8\t8q\t50000000\t60000000\n"
                "4\tarm_gain\t8\t8q\t45200000\t145138636\n"
            )
            rows = classify(read_lineages(lineages), read_events(cna))
            amplified = [r for r in rows if r["amplification_status"] == "amplified"]
            non_amplified = [r for r in rows if r["amplification_status"] == "non-amplified"]
            self.assertTrue(any(r["route_type"] == "bidirectional" for r in amplified))
            self.assertTrue(all(r["route_type"] == "NA" for r in non_amplified))
            self.assertTrue(any(r["cn_trajectory"].startswith("2,3") for r in amplified))


if __name__ == "__main__":
    unittest.main()
