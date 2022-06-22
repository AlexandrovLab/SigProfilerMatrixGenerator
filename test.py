from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from pandas import read_csv
from pandas.testing import assert_frame_equal

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import (
    SigProfilerMatrixGeneratorFunc as matGen,
    CNVMatrixGenerator as scna,
)


REFERENCE_CNV = Path().resolve() / "SigProfilerMatrixGenerator" / "references" / "CNV"


class TestPurpleCnvMatrix(TestCase):
    def test_purple_cnv_file(self):
        """Test parse a CNV file in PURPLE format.

        We have transformed the first records in
        "all.breast.ascat.summary.sample.tsv" from ASCAT format to PURPLE
        format. This unit test verifies, irrespective of the data format, the
        same matrix is generated.
        """
        cnv_in_ascat = REFERENCE_CNV / "all.breast.ascat.summary.sample.tsv"
        cnv_in_purple = REFERENCE_CNV / "example.purple.sample.tsv"
        with TemporaryDirectory() as tmpdir:
            scna.generateCNVMatrix(
                file_type="ASCAT",
                input_file=cnv_in_ascat,
                project="ascat",
                output_path=tmpdir,
            )
            scna.generateCNVMatrix(
                file_type="PURPLE",
                input_file=cnv_in_purple,
                project="purple",
                output_path=tmpdir,
            )

            ascat_matrix_name = tmpdir + "ascat/ASCAT.CNV.matrix.tsv"
            df_ascat = read_csv(ascat_matrix_name, sep="\t", index_col=0)

            purple_matrix_name = tmpdir + "purple/PURPLE.CNV.matrix.tsv"
            df_purple = read_csv(purple_matrix_name, sep="\t", index_col=0)

            # Verify that the first two segments are identical.
            assert_frame_equal(df_ascat.iloc[:2], df_purple.iloc[:2])
