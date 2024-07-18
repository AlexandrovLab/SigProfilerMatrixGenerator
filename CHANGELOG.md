# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.2.27] - 2024-07-18

### Added
- Added environmental variable `SIGPROFILERMATRIXGENERATOR_VOLUME` to enhance configurability.
- Updated CLI to handle boolean arguments more effectively.

### Changed
- Fixed `binom_test` to use the parameter `k` instead of `x` for improved accuracy and clarity.
- Updated Dockerfile to require SigProfilerMatrixGenerator version 1.2.27 (previously 1.2.25).
- Updated `numpy` dependency to require `numpy>=1.18.5,<2.0.0` to prevent compatibility issues with `numpy` 2.0.0.