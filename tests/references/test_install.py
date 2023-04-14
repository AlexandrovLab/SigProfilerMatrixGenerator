import pathlib
import unittest

import pkg_resources

import SigProfilerMatrixGenerator
from SigProfilerMatrixGenerator.references import ref_install


class TestReferenceDir(unittest.TestCase):
    def test_reference_dir_default(self):
        refdir = ref_install.reference_dir(None)
        # using deprecated pkg_resources for compatibility with python 3.8
        expected = pathlib.Path(
            pkg_resources.resource_filename(SigProfilerMatrixGenerator.__name__, "")
        )
        observed = refdir.path
        self.assertEqual(expected, observed)
        # check path is absolute
        self.assertEqual(observed, observed.resolve())

    def test_reference_dir_custom(self):
        paths = ["/somewhere", pathlib.Path("/somewhere")]
        for provided in paths:
            expected = pathlib.Path("/somewhere")
            refdir = ref_install.reference_dir(provided)
            observed = refdir.path
            self.assertEqual(expected, observed)

    def test_path_relative(self):
        provided = pathlib.Path("relative/to/here")
        expected = pathlib.Path.cwd() / "relative/to/here"

        refdir = ref_install.ReferenceDir(provided)
        observed = refdir.path
        self.assertEqual(expected, observed)
