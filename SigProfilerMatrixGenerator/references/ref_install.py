"""Supports installation of reference genomes"""
import pathlib
from importlib import resources
from typing import Optional, Union

import SigProfilerMatrixGenerator


def reference_dir(path: Optional[Union[str, pathlib.Path]] = None) -> "ReferenceDir":
    """Small constructor for ReferenceDir"""
    if path is None:
        path = resources.files(SigProfilerMatrixGenerator)
    return ReferenceDir(pathlib.Path(path))


class ReferenceDir:
    """The directory where the reference genome is installed"""

    def __init__(self, path: pathlib.Path):
        self._path = path

    @property
    def path(self) -> pathlib.Path:
        return self._path.resolve()
