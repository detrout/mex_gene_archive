from unittest import TestCase

from mex_gene_archive.utils import cast_to_list

class TestUtils(TestCase):
    def test_cast_to_list(self):
        self.assertEqual(cast_to_list("hello"), ["hello"])
        self.assertEqual(cast_to_list(["hello"]), ["hello"])

    
