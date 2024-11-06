# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.31] - 2024-11-6

### Fixed
- Resolved a bug where DBS and ID matrices were not being returned due to an indentation issue when the parameter 'chrom_based' was set as True.

### Notice
- Python 3.8 support will be phased out in future updates as we align with Python's ongoing release cycle. Users on Python 3.8 are encouraged to consider upgrading to Python 3.9 or newer.

## [1.2.30] - 2024-10-10

### Fixed
- Resolved a bug that caused the code to crash when running with `exome=True` during INDEL processing. The issue was due to incorrect indentation.

## [1.2.29] - 2024-09-24

### Added
- Added sort function to handle non-sorted input and BED files. This posed as an issue for when the BED file was not sorted by chromosome and start position in cases when certain reference genomes were used because of roman numerals and decimal (ie ChrX vs chr10).

### Changed
- Transcript files were modified to adopt the new sorting standard.


## [1.2.28] - 2024-08-06

### Added
- Added support for processing SV input for VCF versions 4.1, 4.2, and 4.3. The tool now supports both the previous input format (requiring the first six columns and either the "svclass" column or the "strand1" & "strand2" columns) and VCF files.

### Changed
- Updated the README command line examples to use CLI instead of calling the script directly.

## [1.2.27] - 2024-07-18

### Added
- Added environmental variable `SIGPROFILERMATRIXGENERATOR_VOLUME` to enhance configurability.
- Updated CLI to handle boolean arguments more effectively.

### Changed
- Fixed `binomtest` to use the parameter `k` instead of `x` for improved accuracy and clarity.
- Updated Dockerfile to require SigProfilerMatrixGenerator version 1.2.27 (previously 1.2.25).
- Updated `numpy` dependency to require `numpy>=1.18.5,<2.0.0` to prevent compatibility issues with `numpy` 2.0.0.
- Removed incorrect print statement for the SV plot file path.