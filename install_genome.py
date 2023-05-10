#!/usr/bin/env python3
import sys

from SigProfilerMatrixGenerator import install as genInstall


def install_ref(ref_path, genome="GRCh37"):
    genInstall.install(genome, offline_files_path=ref_path)


if __name__ == "__main__":
    ref_path = sys.argv[1]
    genome = sys.argv[2]
    install_ref(ref_path, genome=genome)
