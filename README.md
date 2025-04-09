<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/imgs/protlearn_logo.png" height="85" width="230">
</p>

<p align="center">
  <strong>A Python package for extracting protein sequence features (Fork)</strong>
  <br>
  <a href="https://protlearn.readthedocs.io/en/latest/">Original Documentation</a>
  ·
  <a href="https://github.com/ewijaya/protlearn/issues/new?assignees=&labels=&template=feature_request.md&title=%5BNEW+FEATURE%5D">Request a feature (Fork)</a>
  ·
  <a href="https://github.com/ewijaya/protlearn/issues/new?assignees=&labels=&template=bug_report.md&title=%5BBUG%5D">Report a bug (Fork)</a>
  <br><br>
  <!-- Travis CI and Codecov badges removed as they likely need setup for the fork -->
  <!-- <a href="https://travis-ci.org/tadorfer/protlearn"><img alt="Travis CI" src="https://img.shields.io/travis/tadorfer/protlearn"></a> -->
  <!-- <a href="https://codecov.io/gh/tadorfer/protlearn"><img alt="Codecov" src="https://codecov.io/gh/tadorfer/protlearn/branch/master/graph/badge.svg"></a> -->
  <a href="https://protlearn.readthedocs.io/en/latest/?badge=latest"><img alt="Docs" src="https://readthedocs.org/projects/protlearn/badge/?version=latest"></a>
  <a href="https://pypi.org/project/protlearn/"><img alt="Original PyPI" src="https://img.shields.io/pypi/v/protlearn"></a>
  <a href="https://anaconda.org/conda-forge/protlearn"><img alt="Original Conda version" src="https://img.shields.io/conda/vn/conda-forge/protlearn.svg"></a>
  <a href="https://img.shields.io/pypi/pyversions/protlearn"><img alt="Python versions" src="https://img.shields.io/pypi/pyversions/protlearn"></a>
  <a href="https://github.com/ewijaya/protlearn/blob/master/LICENSE"><img alt="License" src="https://img.shields.io/badge/License-MIT-blue.svg"></a>
</p>
<hr><br>

**Note:** This is a fork of the original `protlearn` package available at `https://github.com/tadorfer/protlearn`. This fork includes specific fixes or modifications.

*protlearn* is a Python package for the feature extraction of amino acid sequences.
It is comprised of three stages - preprocessing, feature computation, and
subsequent dimensionality reduction. This version has been tested with Python 3.12 and includes fixes for potential runtime warnings in autocorrelation calculations.

## Overview

<p align="center">
  <img src="https://raw.githubusercontent.com/tadorfer/protlearn/master/imgs/protlearn_summary.png" height="430" width="624">
</p>

For more information on how to use it, please refer to the [original documentation](https://protlearn.readthedocs.io/en/latest/).

## Installation

### Dependencies

- NumPy
- Pandas
- scikit-learn
- xgboost
- mlxtend
- biopython

### User Installation

#### Installing this Fork (Recommended for fixes)

To install this specific forked version including any fixes:

```bash
pip install git+https://github.com/ewijaya/protlearn.git
```
Or, after cloning the repository:
```bash
git clone https://github.com/ewijaya/protlearn.git
cd protlearn
pip install .
```

#### Installing the Original Package (Via PyPI/Conda)

**Warning:** The following commands will install the *original* package from PyPI or Conda Forge, which may *not* include the fixes present in this fork.

*   **PyPI:**
    ```bash
    pip install protlearn
    ```

*   **Conda:**
    ```bash
    conda install -c conda-forge protlearn
    ```

## Authors

This package was originally created by [Thomas Dorfer](https://github.com/tadorfer).
This fork is currently maintained by [ewijaya](https://github.com/ewijaya).

## License

This package is licensed under the [MIT License](https://github.com/ewijaya/protlearn/blob/master/LICENSE).
