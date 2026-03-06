import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
input is observed muts df (index =tirnuc context, col =alt base)
calculates cosine similarity between spectrums
plus runs a permutation test of the mutations
->permutates mutation with prob of being in branch1 =target1/(target1+target2)
the cosine diff between outputted spectrum
'''

sigs_dir='/home/maria/signatures/COSMIC_sigs_uncollapsed_sum1'

DEFAULT_SIGNATURE_LIST = [
    "SBS1","SBS2","SBS3","SBS4","SBS5","SBS6","SBS7a","SBS7b","SBS7c","SBS7d",
    "SBS8","SBS9","SBS10a","SBS10b","SBS10c","SBS10d","SBS11","SBS12","SBS13","SBS14",
    "SBS15","SBS16","SBS17a","SBS17b","SBS18","SBS19","SBS20","SBS21","SBS22a","SBS22b",
    "SBS23","SBS24","SBS25","SBS26","SBS28","SBS29","SBS30","SBS31","SBS32","SBS33",
    "SBS34","SBS35","SBS36","SBS37","SBS38","SBS39","SBS40a","SBS40b","SBS40c","SBS41",
    "SBS42","SBS44","SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90","SBS91",
    "SBS92","SBS93","SBS94","SBS96","SBS97","SBS98","SBS99"
]



def cosine_similarity_spectra(df1: pd.DataFrame,
                              df2: pd.DataFrame,
                              normalize: bool = True,
                              fill_value: float = 0.0) -> float:
    """
    Cosine similarity between two mutational spectra.

    Parameters
    ----------
    df1, df2 : pd.DataFrame
        Spectra with index = trinucleotide context and columns = alt base (e.g., A,C,G,T).
        Values can be counts, rates, or probabilities.
    normalize : bool
        If True, each spectrum is normalized to sum to 1 before comparison.
        Useful if inputs are counts with different totals.
    fill_value : float
        Value to fill missing contexts/columns after alignment (default 0).

    Returns
    -------
    float
        Cosine similarity in [0, 1] for non-negative spectra. (Can be [-1,1] in general.)
    """
    
    # Align to common shape (union) and fill missing entries
    a, b = df1.align(df2, join="outer", axis=None, fill_value=fill_value)
    
    # Flatten to vectors
    v1 = a.to_numpy(dtype=float).ravel()
    v2 = b.to_numpy(dtype=float).ravel()
    
    # Optional normalization
    if normalize:
        s1, s2 = v1.sum(), v2.sum()
        if s1 != 0:
            v1 = v1 / s1
        if s2 != 0:
            v2 = v2 / s2

    # Handle all-zero spectra
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return np.nan  # undefined if either spectrum is all zeros

    return float(np.dot(v1, v2) / (n1 * n2))

mismatch=[15,6,20,21,26,44]
epsiode11= [15,1,6,98,87,24,44,23,20,42,'10b',94,29,'7b',91,31,35,4,86,11,38,5,18,30,14,39,97,84,33,'40c',92]
sigs_dict={}
for sig in mismatch + [1,5]:
    sig=f'SBS{str(sig)}'
    df=pd.read_csv(f'{sigs_dir}/{sig}',index_col=0)
    sigs_dict[sig]=df



sig_names = list(sigs_dict.keys())
cosine_mat = pd.DataFrame(
    np.zeros((len(sig_names), len(sig_names))),
    index=sig_names,
    columns=sig_names
)

for sig1 in sig_names:
    for sig2 in sig_names:
        v1 = sigs_dict[sig1]
        v2 = sigs_dict[sig2]

        cosine_mat.loc[sig1, sig2] = cosine_similarity_spectra(v1, v2)

sns.clustermap(
    cosine_mat,
    cmap="viridis",
    vmin=0,
    vmax=1,
    figsize=(8,8),
    annot=False
)
plt.savefig('/home/maria/run_simulations_cactus_clean/plots_general/mismatch_corrs.png')
plt.show()
