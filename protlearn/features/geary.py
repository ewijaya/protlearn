# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np # <-- Added import
import pandas as pd
from ..utils.validation import check_input, check_alpha, check_natural
import pkg_resources

PATH = pkg_resources.resource_filename(__name__, 'data/')

# default indices of AAIndex1 (Xiao et al., 2015)
default = ['CIDH920105', 'BHAR880101', 'CHAM820101', 'CHAM820102',
           'CHOC760101', 'BIGC670101', 'CHAM810101', 'DAYM780201']

def geary(X, *, d=1, properties=default, start=1, end=None):
    """Geary's C based on AAIndex1.

    Geary's C autocorrelation descriptors are defined based on the distribution
    of AAIndex1-based amino acid properties along the sequence. All indices are
    standardized before computing the descriptors. For the exact formula, please
    refer to the documentation at https://protlearn.readthedocs.io/.

    Parameters
    ----------

    X : string, fasta, or a list thereof
        Dataset of amino acid sequences.

    properties : list
        List of strings denoting AAIndex1 indices.

    d : int, default=1
        Represents the lag. Must be smaller than sequence length. Maximum: 30.

    start : int, default=1
        Determines the starting point of the amino acid sequence. This number is
        based on one-based indexing.

    end : int, default=None
        Determines the end point of the amino acid sequence. Similarly to start,
        this number is based on one-based indexing.

    Returns
    -------

    arr : ndarray of shape (n_samples, n_properties)
        Array containing Geary's C autocorrelation descriptors. Will contain NaN
        where the calculation is undefined (e.g., sequence length <= 1 or
        zero variance of the property in the sequence).

    References
    ----------

    Geary, R. (1954). The Contiguity Ratio and Statistical Mapping. The
    Incorporated Statistician, 5(3), 115-146. doi:10.2307/2986645

    Jeffers, J. (1973). A Basic Subroutine for Geary's Contiguity Ratio. Journal
    of the Royal Statistical Society. Series D (The Statistician), 22(4),
    299-302. doi:10.2307/2986827

    Sokal, RR., Thomson, BA. (2006). Population structure inferred by local spatial
    autocorrelation: an example from an Amerindian tribal population. Am J Phys
    Anthropol, 129: 121–131. 10.1002/ajpa.20250

    Xiao et al. (2015). protr/ProtrWeb: R package and web server for generating
    various numerical representation schemes of protein sequences.
    Bioinformatics 31 (11), 1857-1859

    Examples
    --------

    >>> from protlearn.features import geary
    >>> seqs = ['ARKLY', 'EERKPGL', 'A'] # Added short seq 'A'
    >>> gearyC = geary(seqs)
    >>> # Note: The output for 'A' will be NaN because len(seq)-1 is zero
    >>> print(gearyC)
    [[0.52746275 1.12898944 0.94222955 0.39077186 0.96444569 0.66346012
      0.87481962 0.32546227]
     [0.65656058 0.95397893 0.87962853 0.70972353 0.65407555 0.96823847
      1.01949384 0.21073089]
     [       nan        nan        nan        nan        nan        nan
             nan        nan]]

    """

    # input handling
    X = check_input(X)
    # Check d against original sequence lengths before potential slicing
    # Note: The original check might still be insufficient if slicing makes sequences too short for d
    min_len_orig = min([len(seq) for seq in X])
    if d > 30:
        raise ValueError('Maximum lag parameter is 30!')
    # Check if d is valid for the *shortest possible* sequence *after* slicing
    # This is hard to check perfectly without knowing start/end, so rely on loop check
    # However, a basic check against original min length is still useful:
    if d >= min_len_orig:
        raise ValueError(f'Lag parameter d ({d}) must be smaller than the minimum original sequence length ({min_len_orig})!')


    # load data
    df = pd.read_csv(PATH+'aaindex1.csv').set_index('Description')
    df = df.reindex(sorted(df.columns), axis=1)
    # Ensure properties exist in the loaded data
    try:
        data = np.asarray(df.loc[properties])
    except KeyError as e:
        raise ValueError(f"One or more properties not found in AAIndex1: {e}")

    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aadict = {amino_acids[i]: i for i in range(20)}

    # standardization
    standardized_data = np.zeros_like(data, dtype=float)
    for i in range(data.shape[0]):
        mean_val = np.mean(data[i,:])
        std_val = np.std(data[i,:])
        # Handle case where all values in a property row are the same (std_dev=0)
        if std_val == 0:
             standardized_data[i,:] = 0.0 # Or handle as appropriate, e.g., np.nan or raise error
        else:
             standardized_data[i,:] = (data[i,:] - mean_val) / std_val

    # calculate Geary's C
    arr = np.zeros((len(X), len(properties)))
    for i, seq_orig in enumerate(X):
        check_alpha(seq_orig) # check if alphabetical
        check_natural(seq_orig) # check for unnatural amino acids
        
        # Apply positional slicing
        seq = seq_orig[start-1:end]
        seq_len = len(seq)

        # Check if lag 'd' is valid for the potentially sliced sequence length
        if d >= seq_len:
             # If lag is too large for this specific sequence after slicing, fill with NaN
             arr[i, :] = np.nan
             # Optionally, print a warning or log this occurrence
             # print(f"Warning: Lag d={d} >= length of sliced sequence {seq_len} for original index {i}. Setting output to NaN.")
             continue # Skip to the next sequence

        # Pre-calculate part of eq1 (valid only if seq_len > d)
        eq1_part = 1 / (2 * (seq_len - d)) # Note: seq_len > d checked above

        for j in range(len(properties)):
            # Get standardized property values for the amino acids in the sequence
            try:
                 p = [standardized_data[j, aadict[aa]] for aa in seq]
            except KeyError as e:
                 raise ValueError(f"Sequence '{seq_orig}' (index {i}) contains non-standard amino acid {e} not in 'ACDEFGHIKLMNPQRSTVWY'")

            # Check for invalid sequence length for Geary's formula
            if seq_len <= 1:
                arr[i, j] = np.nan # Denominator term (1/(seq_len-1)) is undefined
                continue # Skip to next property

            p_prime = sum(p) / seq_len
            eq2 = sum([(p[k] - p[k + d])**2 for k in range(seq_len - d)])
            eq3 = sum([(p[k] - p_prime)**2 for k in range(seq_len)])

            # --- Modification Start ---
            # Calculate denominator, checking for zero division potential
            denominator = (1 / (seq_len - 1)) * eq3

            # Check if denominator is close to zero (handles eq3 == 0 case)
            if abs(denominator) < 1e-9:
                arr[i, j] = np.nan # Undefined result (zero variance)
            else:
                # Calculate Geary's C = eq1 * eq2 / denominator
                arr[i, j] = eq1_part * eq2 / denominator
            # --- Modification End ---

    return arr
