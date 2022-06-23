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
    def test_ascat_cnv_file(self):
        """Test parse a CNV file in ASCAT format."""
        cnv_in_ascat = REFERENCE_CNV / "all.breast.ascat.summary.sample.tsv"
        with TemporaryDirectory() as tmpdir:
            scna.generateCNVMatrix(
                file_type="ASCAT_NGS",
                input_file=cnv_in_ascat,
                project="ascat",
                output_path=tmpdir,
            )

            ascat_matrix_name = tmpdir + "ascat/ASCAT_NGS.CNV.matrix.tsv"
            df_ascat = read_csv(ascat_matrix_name, sep="\t", index_col=0).transpose()

            sample_name = "PD9067a"
            # Verify that there are 15 homozygous deletions (verified by counting).
            is_homdel = df_ascat.columns.map(lambda x: "homdel" in x)
            self.assertEqual(df_ascat.loc[sample_name, is_homdel].sum(), 15)
            # There are 158 records for sample "PD9067a".
            self.assertEqual(df_ascat.sum(axis=1)[sample_name], 158)

    def test_purple_cnv_file(self):
        """Test parse a CNV file in PURPLE format."""
        cnv_in_purple = REFERENCE_CNV / "example.purple.tsv"
        with TemporaryDirectory() as tmpdir:
            scna.generateCNVMatrix(
                file_type="PURPLE",
                input_file=str(cnv_in_purple),
                project="purple",
                output_path=tmpdir,
            )

            purple_matrix_name = tmpdir + "purple/PURPLE.CNV.matrix.tsv"
            df_purple = read_csv(purple_matrix_name, sep="\t", index_col=0).transpose()

            sample_name = "example.purple"
            # Verify that the right column is 1.
            self.assertEqual(df_purple.loc[sample_name, "3-4:het:>40Mb"], 1)
            # Check that the remaining columns are zero.
            self.assertEqual(df_purple.sum(axis=1)[sample_name], 1)
