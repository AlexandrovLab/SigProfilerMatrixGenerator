"""Supports installation of reference genomes"""
import pathlib
from importlib import resources
from typing import Optional, Union

import pkg_resources

import SigProfilerMatrixGenerator


def reference_dir(
    secondary_chromosome_install_dir: Optional[Union[str, pathlib.Path]] = None
) -> "ReferenceDir":
    """Small constructor for ReferenceDir"""
    if secondary_chromosome_install_dir is not None:
        install_dir = pathlib.Path(secondary_chromosome_install_dir)
    else:
        install_dir = None
    result = ReferenceDir(secondary_chromosome_install_dir=install_dir)
    return result


class ReferenceDir:
    """The directory where the reference genomes are installed"""

    _default_chromosome_dir = pathlib.Path("references", "chromosomes")

    def __init__(self, secondary_chromosome_install_dir: Optional[pathlib.Path] = None):
        self._secondary_chromosome_install_dir = secondary_chromosome_install_dir

    @property
    def path(self) -> pathlib.Path:
        """Where references that cannot be installed are"""
        root_path = self._get_package_installation_folder()
        return root_path.resolve()

    def get_fasta_dir(self) -> pathlib.Path:
        result = self.get_chromosomes_dir() / "fasta"
        return result

    def get_tsb_dir(self) -> pathlib.Path:
        result = self.get_chromosomes_dir() / "tsb"
        return result

    def get_chromosomes_dir(self) -> pathlib.Path:
        if self._secondary_chromosome_install_dir is None:
            root_dir = self.path / self._default_chromosome_dir
        else:
            root_dir = self._secondary_chromosome_install_dir.resolve()
        result = root_dir
        return result

    def _get_package_installation_folder(self):
        try:
            # python >= 3.9
            result = resources.files(SigProfilerMatrixGenerator)
        except AttributeError:
            # python 3.8, deprecated
            result = pathlib.Path(
                pkg_resources.resource_filename(SigProfilerMatrixGenerator.__name__, "")
            )
        return result
