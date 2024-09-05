from io import StringIO
import json
from unittest import TestCase


from mex_gene_archive.kbpython import (
    make_quantification_name,
    parse_kallisto_run_info_stream,
)



class TestKBPython(TestCase):
    def test_parse_kallisto_run_info_stream(self):
        kallisto_version = "0.50.1"
        kallisto_args = "~/kallisto bus -i ref/kallisto/c57bl6j.idx -o kallisto/igvf_b01/Subpool_1/ -x SPLIT-SEQ -t 12 --fr-stranded fastq/igvf_b01_Subpool_1_1_L001_R1_001.fastq.gz fastq/igvf_b01_Subpool_1_1_L001_R2_001.fastq.gz fastq/igvf_b01_Subpool_1_2_L001_R1_001.fastq.gz fastq/igvf_b01_Subpool_1_2_L001_R2_001.fastq.gz"
        stream = StringIO(json.dumps({
            "n_targets": 206337,
            "n_bootstraps": 0,
            "n_processed": 318003931,
            "n_pseudoaligned": 234474089,
            "n_unique": 175798579,
            "p_pseudoaligned": 73.7,
            "p_unique": 55.3,
            "kallisto_version": kallisto_version,
            "index_version": 13,
            "start_time": "Fri Aug  9 13:31:15 2024",
            "call": kallisto_args,
        }))
        results = parse_kallisto_run_info_stream(stream)
        self.assertEqual(results["software"], "kallisto")
        self.assertEqual(results["software_version"], kallisto_version)
        self.assertEqual(results["arguments"], kallisto_args)

    def test_make_quantification_name(self):
        self.assertEqual(make_quantification_name(["cell"]), "cell")
        self.assertEqual(make_quantification_name(["nascent", "ambiguous", "mature"]), "core")
        self.assertEqual(make_quantification_name(["ambiguous", "nascent", "mature"]), "core")
        self.assertEqual(make_quantification_name(["nascent", "ambiguous", "mature", "total"]), "total+core")
        self.assertEqual(make_quantification_name(["ambiguous", "total", "nascent", "mature"]), "total+core")
        self.assertRaises(ValueError, make_quantification_name, ["too", "hot"])
