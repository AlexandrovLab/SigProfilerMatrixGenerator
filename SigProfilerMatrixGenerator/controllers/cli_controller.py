import argparse
from typing import List

from SigProfilerMatrixGenerator import install, test_helpers


def parse_arguments_test(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run tests for SigProfilerMatrixGenerator"
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
    result = parser.parse_args(args)
    return result


def parse_arguments_install(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install reference files")
    parser.add_argument("source_dir", help="Directory with local files")
    parser.add_argument("genome", help="Genome to install")
    result = parser.parse_args(args)
    return result


class CliController:
    def dispatch_install(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_install(user_args)
        install.install(parsed_args.genome, offline_files_path=parsed_args.source_dir)

    def dispatch_test(self, user_args: List[str]) -> None:
        parsed_args = parse_arguments_test(user_args)
        test_helpers.install_genomes(parsed_args.download_genomes)
        test_helpers.test_genomes(parsed_args.test_genome)
