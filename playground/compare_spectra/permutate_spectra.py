import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''
input is observed muts df (index =tirnuc context, col =alt base)
calculates cosine similarity between spectrums
plus runs a permutation test of the mutations
->permutates mutation with prob of being in branch1 =target1/(target1+target2)
the cosine diff between outputted spectrum


'''

branch2='hg38_Anc4'
#branch2='Anc4_Anc3'
#branch2='Anc1_Anc0'
#
branch1='fullTreeAnc212_fullTreeAnc213'
branch1_muts=f'/home/maria/run_simulations_cactus_clean/{branch1}/intron_matrix/muts'
branch1_trinucs=f'/home/maria/run_simulations_cactus_clean/{branch1}/intron_matrix/trinucs'
branch2_muts=f'/home/maria/run_simulations_cactus_clean/{branch2}/intron_matrix/muts'
branch2_trinucs=f'/home/maria/run_simulations_cactus_clean/{branch2}/intron_matrix/trinucs'

SBS5_file='/home/maria/signatures/TSB_signatures_scaled/SBS5'
SBS1_file='/home/maria/signatures/TSB_signatures_scaled/SBS1'
denovo_file = '/home/maria/find_intron_matrix/output_denovo/rates'

no_perms=100

BASES=['A','C','G','T']

muts1 =pd.read_csv(branch1_muts,index_col=0,sep='\t')
target1 =pd.read_csv(branch1_trinucs,index_col=0,sep='\t').rename(columns={'count':'target1'})
muts2 =pd.read_csv(branch2_muts,index_col=0,sep='\t')
target2 =pd.read_csv(branch2_trinucs,index_col=0,sep='\t').rename(columns={'count':'target2'})

rates1=muts1.div(target1['target1'],axis=0)
rates2=muts2.div(target2['target2'],axis=0)
target = target1.join(target2,how='inner')
target['prob'] = target['target1']/(target['target1']+target['target2'])

#sum muts to redistirbute
muts=muts1+muts2



def permutate(mut,target):
    perm1=mut.copy()
    for trinuc,row in mut.iterrows():
        for base in BASES:
            perm1.loc[trinuc,base]=np.random.binomial(row[base],target.loc[trinuc,'prob'])
    
    #remaining muts assigined to perm2
    perm2=mut-perm1
    rates1=perm1.div(target['target1'],axis=0)
    rates2=perm2.div(target['target2'],axis=0)
    return cosine_similarity_spectra(rates1,rates2)

#downsample total A to toal B
#WRONG-NOT PERSEVING TOTAL NUMBER OF COUNTS
#WANT TO MULTINOMIALLY SAMPLE ALL OF A with B_TOTAL INSTEAD (think same in expectation?)
def downsize(A,B_total):
    #assumes A has bigger total
    A_total=A.sum().sum()  
    A_downsampled=A.copy()
    A_prob=A/A_total
    for trinuc,row in A_prob.iterrows():
        for col in A.columns:
            A_downsampled.loc[trinuc,col]=np.random.binomial(B_total,row[col])
    return A_downsampled

def downsize_multinomial(A: pd.DataFrame, new_total: int, rng=None) -> pd.DataFrame:
    rng = np.random.default_rng(rng)
    probs = A.to_numpy(dtype=float).ravel()
    probs_sum = probs.sum()
    if probs_sum == 0:
        return A * 0
    probs = probs / probs_sum
    draws = rng.multinomial(new_total, probs)
    out = pd.DataFrame(draws.reshape(A.shape), index=A.index, columns=A.columns)
    return out

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


def down_sample_spectra(target1,muts1,target2,muts2):
    '''
    for two evolutionary brnaches have dfs of mut spectras and targets
    '''
    total_muts1=muts1.sum().sum()
    
    total_muts2=muts2.sum().sum()
    if total_muts1>total_muts2:
        
        target1_down=downsize_multinomial(target1,target2.sum()[0])
        muts1_down=downsize_multinomial(muts1,total_muts2)
        sim=cosine_similarity_spectra(muts1_down.div(target1_down['target1'],axis=0),muts2.div(target2['target2'],axis=0))
    return sim


real_diff=cosine_similarity_spectra(rates1,rates2,normalize=True)
print(f'cosine sim between {branch1} {branch2}: {real_diff}')

perms=[]
for i in range(no_perms):
    perms.append(down_sample_spectra(target1,muts1,target2,muts2))
plt.hist(x=perms)
plt.axvline(x=real_diff,color='red',label='observed cosine similarity')
plt.xlabel(f'expected cosine diff in {branch1} vs {branch2}')
plt.ylabel('frequency')
plt.show()


SBS1=pd.read_csv(SBS1_file,index_col=0)
SBS5=pd.read_csv(SBS5_file,index_col=0)
cosmic = 3*SBS5 + SBS1

cosmic_diff=cosine_similarity_spectra(rates1,cosmic,normalize=True)
print(f'cosine sim between {branch1}  3*SBS5 + SBS1: {cosmic_diff}')


denovo_df = pd.read_csv(denovo_file,sep='\t',index_col=0)
denovo_diff=cosine_similarity_spectra(rates1,denovo_df,normalize=True)
print(f'cosine sim between {branch1}  de novo spectrum: {denovo_diff}')

perms=[]
for i in range(no_perms):
    perms.append(permutate(muts,target))
plt.hist(x=perms)
plt.axvline(x=real_diff)
plt.xlabel(f'cosine diff in {branch1} vs {branch2}')
plt.ylabel('frequency')
plt.show()
