from unittest import mock

import pytest

from SigProfilerMatrixGenerator import install, test_helpers
from SigProfilerMatrixGenerator.controllers import cli_controller


class TestController:
    @pytest.fixture
    def mock_install(self, monkeypatch):
        """Temporarily change install.install with a mock object for testing"""
        mock_install = mock.Mock()
        monkeypatch.setattr(install, "install", mock_install)
        return mock_install

    @pytest.fixture
    def mock_test(self, monkeypatch):
        """Temporarily change test_helpers.test_one_genome with a mock object for testing"""
        mock_test = mock.Mock()
        monkeypatch.setattr(test_helpers, "test_one_genome", mock_test)
        return mock_test

    def test_dispatch_install(self, mock_install):
        controller = cli_controller.CliController()
        controller.dispatch_install(["/somewhere", "GRCh37"])
        mock_install.assert_called_with("GRCh37", offline_files_path="/somewhere")

    genome_calls = [
        pytest.param(
            ["--test_genome", "GRCh37", "dog"],
            [mock.call("GRCh37"), mock.call("dog")],
            id="two genomes",
        ),
        pytest.param(
            ["--test_genome", "GRCh37"], [mock.call("GRCh37")], id="one genome"
        ),
        pytest.param(["--test_genome", "unknown_genome"], [], id="no known genomes"),
        pytest.param(
            ["--test_genome", "all"],
            [mock.call(genome) for genome in test_helpers.TEST_GENOMES],
            id="all genomes",
        ),
    ]

    @pytest.mark.parametrize("provided, expected", genome_calls)
    def test_dispatch_test_genome(self, mock_test, provided, expected):
        controller = cli_controller.CliController()
        controller.dispatch_test(provided)
        mock_test.assert_has_calls(expected)

    download_calls = [
        pytest.param(
            ["--download", "GRCh37", "dog"],
            [mock.call("GRCh37"), mock.call("dog")],
            id="two genomes",
        ),
        pytest.param(
            ["--download", "all"],
            [mock.call(genome) for genome in test_helpers.TEST_GENOMES],
            id="all genomes",
        ),
    ]

    @pytest.mark.parametrize("provided, expected", download_calls)
    def test_dispatch_test_download(self, mock_install, provided, expected):
        controller = cli_controller.CliController()
        controller.dispatch_test(provided)
        mock_install.assert_has_calls(expected)

    def test_dispatch_test_download_and_genome(self, mock_install, mock_test):
        provided = ["--download", "GRCh37", "--test_genome", "GRCh37"]
        controller = cli_controller.CliController()
        controller.dispatch_test(provided)
        mock_install.assert_has_calls([mock.call("GRCh37")])
        mock_test.assert_has_calls([mock.call("GRCh37")])
