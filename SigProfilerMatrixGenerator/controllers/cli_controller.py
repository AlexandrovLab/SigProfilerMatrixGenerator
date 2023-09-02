import argparse
from typing import List

from SigProfilerMatrixGenerator import install, test_helpers
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as mg


def parse_arguments_test(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run tests for SigProfilerMatrixGenerator."
    )
    parser.add_argument(
        "-t",
        "--test_genome",
        help="Genomes to test (GRCh37, any test genome, or all)",
        nargs="+",
        default=None,
    )
    parser.add_argument(
        "-d",
        "--download_genomes",
        help="Download genomes for the test (GRCh37, any test genome, or all)",
        nargs="+",
        default=None,
    )
    parser.add_argument(
        "-v",
        "--volume",
        help="Specify a destination for the downloaded genomes (default: None, used for Docker)",
        default=None,
    )

    result = parser.parse_args(args)
    return result


def parse_arguments_install(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install reference genome files.")
    parser.add_argument(
        "genome",
        help="The reference genome to install. Supported genomes include {c_elegans, dog, ebv, GRCh37, GRCh38, mm9, mm10, mm39, rn6, yeast}.",
    )
    parser.add_argument(
        "-l",
        "--local_install_genome",
        help="""
            Install an offline reference genome downloaded from the Alexandrov Lab's FTP server.
            Provide the absolute path to the locally-stored genome file.
            For downloads, visit AlexandrovLab's server:
            ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator/
            """,
        default=None,
    )
    parser.add_argument(
        "-v",
        "--volume",
        help="Specify a destination for the downloaded genomes (default: None, used for Docker)",
        default=None,
    )
    result = parser.parse_args(args)
    return result


def parse_arguments_matrix_generator(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create mutational matrices for all types of somatic mutations"
    )

    # Mandatory arguments
    parser.add_argument("project", help="The name of the project.")
    parser.add_argument(
        "reference_genome",
        help="The name of the reference genome. Supported values {c_elegans, dog, ebv, GRCh37, GRCh38, mm9, mm10, mm39, rn6, yeast}.",
    )
    parser.add_argument("path_to_input_files", help="The path to the input files.")

    # Optional arguments
    parser.add_argument(
        "--exome",
        action="store_true",
        help="Downsamples mutational matrices to the exome regions of the genome. Default is False.",
    )
    parser.add_argument(
        "--bed_file",
        default=None,
        help="Downsamples mutational matrices to custom regions of the genome. Provide full path to BED file. Note: BED file header is required. Default is None.",
    )
    parser.add_argument(
        "--chrom_based",
        action="store_true",
        help="Outputs chromosome-based matrices. Default is False.",
    )
    parser.add_argument(
        "--plot",
        type=bool,
        default=False,
        help="Integrates with SigProfilerPlotting to output visualizations for each matrix. Default is False.",
    )
    parser.add_argument(
        "--tsb_stat",
        type=bool,
        default=False,
        help="Outputs the results of a transcriptional strand bias test for the respective matrices. Default is False.",
    )
    parser.add_argument(
        "--seqInfo",
        type=bool,
        default=True,
        help="Outputs original mutations into a text file with the SigProfilerMatrixGenerator classification for each mutation. Default is True.",
    )
    parser.add_argument(
        "--cushion",
        type=int,
        default=100,
        help="Adds an Xbp cushion to the exome/bed_file ranges for downsampling mutations. Default is 100.",
    )

    parser.add_argument(
        "-v",
        "--volume",
        help="Specify a destination for the downloaded genomes (default: None, used for Docker)",
        default=None,
    )

    result = parser.parse_args(args)
    return result


class CliController:
    def dispatch_install(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_install(user_args)
        install.install(
            parsed_args.genome,
            offline_files_path=parsed_args.local_install_genome,
            volume=parsed_args.volume,
        )

    def dispatch_test(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_test(user_args)
        test_helpers.install_genomes(parsed_args.download_genomes)
        test_helpers.test_genomes(parsed_args.test_genome, parsed_args.volume)

    def dispatch_matrix_generator(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_matrix_generator(user_args)

        mg.SigProfilerMatrixGeneratorFunc(
            project=parsed_args.project,
            reference_genome=parsed_args.reference_genome,
            path_to_input_files=parsed_args.path_to_input_files,
            exome=parsed_args.exome,
            bed_file=parsed_args.bed_file,
            chrom_based=parsed_args.chrom_based,
            plot=parsed_args.plot,
            tsb_stat=parsed_args.tsb_stat,
            seqInfo=parsed_args.seqInfo,
            cushion=parsed_args.cushion,
            volume=parsed_args.volume,
        )
