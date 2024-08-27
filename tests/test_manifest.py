from io import BytesIO, StringIO
import os
from pathlib import Path
from unittest import TestCase

from mex_gene_archive.manifest import (
    compute_md5sums,
    compute_stream_md5sum,
    ConfigError,
    logger as manifest_logger,
    read_manifest,
    validate_config_metadata,
    write_manifest,
)


class TestManifest(TestCase):
    def test_compute_stream_md5sum(self):
        zero_file = BytesIO(b"\x00" * 10000)
        digest = compute_stream_md5sum(zero_file).hexdigest()
        self.assertEqual("b85d6fb9ef4260dcf1ce0a1b0bff80d3", digest)

    def test_compute_md5sum(self):
        b = BytesIO()
        initfile = Path(__file__).parent / "__init__.py"
        files = [b, initfile]
        digests = compute_md5sums(files)
        for row in digests:
            # this is just the md5 hash of a 0 byte file
            self.assertEqual(
                "d41d8cd98f00b204e9800998ecf8427e",
                row[1],
                "unexpected hash for {}".format(row[0]),
            )

    def test_validate_config_metadata_valid(self):
        config = {
            "read1": ["ENCFF150FBF", "ENCFF385IAW"],
            "read2": ["ENCFF351VBS", "ENCFF503CCI"],
            "include_introns": True,
            "stranded": "Forward",
            "experiment_accession": "ENCSR724KET",
            "description": "snRNA on a human adrenal gland",
            "library_accession": "ENCLB002DZK",
            "umi_version": 12,
            "genome_dir": "genome_dir",
            "allow_list": "allow_list",
            "mem_mb": 65535,
            "disk_mb": 51200,
        }
        validate_config_metadata(config)

    def test_validate_config_metadata_pipeline_variables(self):
        config = {
            "read1": ["ENCFF150FBF", "ENCFF385IAW"],
            "read2": ["ENCFF351VBS", "ENCFF503CCI"],
            "software": "STAR",
            "software_version": "dev_EoI_2.7.9a_2021-09-30",
            "arguments": "--help",
            "include_introns": True,
            "stranded": "Forward",
            "experiment_accession": "ENCSR724KET",
            "description": "snRNA on a human adrenal gland",
            "library_accession": "ENCLB002DZK",
            "umi_version": 12,
            "genome_dir": "genome_dir",
            "allow_list": "allow_list",
            "mem_mb": 65535,
            "disk_mb": 51200,
        }
        with self.assertLogs(manifest_logger) as log:
            self.assertRaises(ConfigError, validate_config_metadata, config)

    def test_validate_config_metadata_missing_variables(self):
        config = {
            "read1": ["ENCFF150FBF", "ENCFF385IAW"],
            "read2": ["ENCFF351VBS", "ENCFF503CCI"],
            "include_introns": True,
            "stranded": "Forward",
            "experiment_accession": "ENCSR724KET",
            "library_accession": "ENCLB002DZK",
            "umi_version": 12,
            "genome_dir": "genome_dir",
            "allow_list": "allow_list",
            "mem_mb": 65535,
            "disk_mb": 51200,
        }
        with self.assertLogs(manifest_logger) as log:
            self.assertRaises(ConfigError, validate_config_metadata, config)

    def test_validate_config_metadata_missing_library_accession(self):
        config = {
            "read1": ["ENCFF150FBF", "ENCFF385IAW"],
            "read2": ["ENCFF351VBS", "ENCFF503CCI"],
            "include_introns": True,
            "stranded": "Forward",
            "experiment_accession": "ENCSR724KET",
            "library_accession": "ENCLB002DZK",
            "umi_version": 12,
            "genome_dir": "genome_dir",
            "allow_list": "allow_list",
            "mem_mb": 65535,
            "disk_mb": 51200,
        }
        with self.assertLogs(manifest_logger) as log:
            self.assertRaises(ConfigError, validate_config_metadata, config)

    def test_create_igvf_metadata(self):
        config = {
            "input_file_sets": ["TSTDS34582101", "TSTDS07432728"],
            "description": ["RNA-seq of fresh tomato with basil"],
        }

        validate_config_metadata(config)

        with StringIO(newline="") as buf:
            write_manifest(buf, config)
            manifest = buf.getvalue()

        expected = {
            0: "name\tvalue\r\n",
            1: "input_file_sets\tTSTDS34582101\r\n",
            2: "\tTSTDS07432728\r\n",
            3: "description\tRNA-seq of fresh tomato with basil\r\n",
        }
        for i, line in enumerate(StringIO(manifest, newline="")):
            self.assertEqual(line, expected[i])

        with StringIO(manifest, newline="") as buf:
            round_trip = read_manifest(buf)

            self.assertEqual(round_trip["input_file_sets"], config["input_file_sets"])
