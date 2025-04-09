# Author: Thomas Dorfer <thomas.a.dorfer@gmail.com>

import numpy as np
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
        Array containing Geary's C autocorrelation descriptors.

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
    >>> seqs = ['ARKLY', 'EERKPGL']
    >>> gearyC = geary(seqs)
    >>> gearyC
    array([[0.52746275, 1.12898944, 0.94222955, 0.39077186, 0.96444569,
            0.66346012, 0.87481962, 0.32546227],
           [0.65656058, 0.95397893, 0.87962853, 0.70972353, 0.65407555,
            0.96823847, 1.01949384, 0.21073089]])

    """

    # input handling
    X = check_input(X)
    min_len = min([len(seq) for seq in X])
    if d > 30:
        raise ValueError('Maximum lag parameter is 30!')
    if d >= min_len:
        raise ValueError('Lag parameter d must be smaller than sequence length!')

    # load data
    df = pd.read_csv(PATH+'aaindex1.csv').set_index('Description')
    df = df.reindex(sorted(df.columns), axis=1)
    data = np.asarray(df.loc[properties]) # Use specified or default properties

    # list of amino acids (IUPAC standard)
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    aadict = {amino_acids[i]: i for i in range(20)}

    # standardization
    for i in range(data.shape[0]):
        mean_val = np.mean(data[i, :])
        std_dev = np.std(data[i, :])
        if not np.isclose(std_dev, 0): # Avoid division by zero if std is zero
             data[i,:] = [(j - mean_val) / std_dev for j in data[i, :]]
        else:
             data[i,:] = [0.0 for j in data[i,:]] # Assign zero if std is zero

    # calculate Geary's C
    arr = np.zeros((len(X), len(properties)))
    for i, seq in enumerate(X):
        check_alpha(seq) # check if alphabetical
        check_natural(seq) # check for unnatural amino acids
        seq = seq[start-1:end] # positional information

        # Check for sequence length = 1 after slicing
        if len(seq) <= 1:
           arr[i, :] = 0.0 # Assign 0 or handle as appropriate for length 1
           continue # Skip calculation for this sequence

        eq1 = 1.0 / (2 * (len(seq) - d)) # Use 1.0 for float division
        for j in range(len(properties)):
            p = [data[j, aadict[aa]] for aa in seq]
            p_prime = sum(p) / len(seq)
            eq2 = sum([(p[k] - p[k+d])**2 for k in range(len(seq)-d)])
            eq3 = sum([(p[k] - p_prime)**2 for k in range(len(seq))])

            # --- Start Modification ---
            # Check for len(seq) - 1 == 0 explicitly, though d>=min_len check should cover it
            if len(seq) - 1 == 0:
                arr[i,j] = 0.0
                continue

            denominator_scale = (1.0 / (len(seq) - 1)) # Use 1.0 for float division
            denominator = denominator_scale * eq3

            # Check if the denominator is zero or very close to it
            if np.isclose(denominator, 0.0):
                arr[i, j] = 0.0  # Assign 0 if variance (eq3) is zero
            else:
                # Original calculation if denominator is non-zero
                arr[i, j] = (eq1 * eq2) / denominator
            # --- End Modification ---

    return arr
