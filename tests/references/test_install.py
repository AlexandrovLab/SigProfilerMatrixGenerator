import pathlib
import unittest
from importlib import resources

import SigProfilerMatrixGenerator
from SigProfilerMatrixGenerator.references import ref_install


class TestReferenceDir(unittest.TestCase):
    def test_refdir_default(self):
        refdir = ref_install.reference_dir(None)
        expected = resources.files(SigProfilerMatrixGenerator)
        observed = refdir.path
        self.assertEqual(expected, observed)
        self.assertEqual(observed, observed.resolve())

    def test_refdir_custom_absolute(self):
        paths = ["/somewhere", pathlib.Path("/somewhere")]
        for custom_dir in paths:
            refdir = ref_install.reference_dir(custom_dir)
            expected = pathlib.Path("/somewhere")
            observed = refdir.path
            self.assertEqual(expected, observed)
            self.assertEqual(observed, observed.resolve())

    def test_refdir_custom_relative(self):
        refdir = ref_install.reference_dir("relative/to/here")
        here = pathlib.Path.cwd()
        expected = here / "relative/to/here"
        observed = refdir.path
        self.assertEqual(expected, observed)
