import os

import pandas as pd
from pandas.testing import assert_frame_equal

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import (
    SigProfilerMatrixGeneratorFunc as matGen,
    ref_install,
)

reference_dir = ref_install.reference_dir()
TEST_INPUT_DIR = str(reference_dir.path / "references/tests/")
FILE_PREF = "test_example"
TEST_GENOMES = [
    "c_elegans",
    "dog",
    "GRCh37",
    "GRCh38",
    "mm9",
    "mm10",
    "mm39",
    "rn6",
    "yeast",
]


def load_and_compare(matrices, solution_dir):
    # loop through every key in matrices and load the corresponding solution file
    for key in matrices:
        solution_file = os.path.join(solution_dir, FILE_PREF + ".SBS" + key + ".all")
        solution_df = pd.read_csv(solution_file, sep="\t", index_col=0)
        # check that the dataframes are equal
        assert_frame_equal(matrices[key], solution_df)


def test_one_genome(genome, volume):
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        FILE_PREF,
        genome,
        os.path.join(TEST_INPUT_DIR, genome),
        plot=False,
        exome=False,
        bed_file=None,
        chrom_based=False,
        tsb_stat=False,
        seqInfo=False,
        cushion=100,
        volume=volume,
    )

    # load the results from the solution directory and compare them to the dataframes in matrices
    solution_dir = os.path.join(TEST_INPUT_DIR, "solutions/" + genome + "/")
    load_and_compare(matrices, solution_dir)


def install_genomes(genome_install_list):
    if genome_install_list is None:
        return
    # download all the genomes specified by user
    for tmp_genome in genome_install_list:
        if tmp_genome == "all":
            for genome in TEST_GENOMES:
                genInstall.install(genome)
        elif tmp_genome not in TEST_GENOMES:
            print(
                "Warning: Download for",
                tmp_genome,
                "is skipped because it is not a valid test genome.",
            )
        else:
            genInstall.install(tmp_genome)


def test_genomes(test_genome, volume=None):
    # Check for all test, otherwise test specific genomes
    if test_genome is None:
        print("No genomes specified for testing. Please specify a genome or all.")
    elif test_genome[0] == "all":
        for genome in TEST_GENOMES:
            test_one_genome(genome, volume=volume)
    else:
        for genome in test_genome:
            if genome not in TEST_GENOMES:
                print(
                    genome
                    + " is not a valid test genome. Please choose from: "
                    + str(TEST_GENOMES)
                )
                continue
            try:
                test_one_genome(genome, volume=volume)
                print("Completed test for " + genome)
            except Exception as e:
                assert False, "Test failed for " + genome + ":\n" + str(e)
