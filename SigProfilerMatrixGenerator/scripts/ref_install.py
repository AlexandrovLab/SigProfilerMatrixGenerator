"""Supports installation of reference genomes"""
import pathlib
from importlib import resources
from typing import Optional, Union

import pkg_resources

import SigProfilerMatrixGenerator


def reference_dir(path: Optional[Union[str, pathlib.Path]] = None) -> "ReferenceDir":
    """Small constructor for ReferenceDir"""
    if path is None:
        try:
            # python >= 3.9
            path = resources.files(SigProfilerMatrixGenerator)
        except AttributeError:
            # python 3.8, deprecated
            path = pathlib.Path(
                pkg_resources.resource_filename(SigProfilerMatrixGenerator.__name__, "")
            )
    return ReferenceDir(pathlib.Path(path))


class ReferenceDir:
    """The directory where the reference genomes are installed"""

    def __init__(self, path: pathlib.Path):
        self._path = path

    @property
    def path(self) -> pathlib.Path:
        return self._path.resolve()
