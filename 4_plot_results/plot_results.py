import pandas as pd
import glob
import matplotlib.pyplot as plt
pd.set_option("display.precision", 15)
import numpy as np
plt.rcParams.update({'font.size': 22})
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.patches import Rectangle
import os

TARGET='Anc4'
ANC='Anc3'
no_of_sims =100
lower_quant_cutoff = 0.025
upper_quant_cutoff = 0.975

cancer_genes= '/home/maria/run_simulations/data/Cosmic_CancerGeneCensus_names.txt'
gene_ids = '/home/maria/cactus_target_size/auxillary/gene_name_id.csv'
WEIGHTS_FILE = '/home/maria/data/signature_weights/calc_mut_weights.csv'
DIR = f'/home/maria/run_simulations_cactus_clean/{TARGET}_{ANC}'
TARGET_FILE= f'{DIR}/output/{TARGET}'
ANC_FILE= f'{DIR}/output/{ANC}'
PLOT_DIR= f'{DIR}/plots'

#make directory
if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)


#OPEN DFS TO BE SUED MULTIPLE TIMES
GENE_ID_DF = pd.read_csv(gene_ids)
FILTERED_SUM = pd.read_csv(WEIGHTS_FILE, index_col=0, names=['norm_weight'])
cancer_genes_df = pd.read_csv(cancer_genes, names=['gene_name'])
#convert cancer genes to ensembl gene ids
cancer_id_df = cancer_genes_df.merge(GENE_ID_DF, on=['gene_name'], how='inner')
CANCER_GENES_LIST = cancer_id_df['gene'].tolist()

#define functions
def extract_avg_l_per_sig(file_path, species, cancer_genes_list):
    #looking at current human 
    files_list = glob.glob(f'{file_path}/*')
    files_list.remove(f'{file_path}/M_n_non_syn_muts')
    l_c = []
    l_nc = []
    signature = []
    for l_file in files_list:
        
        l_df = pd.read_csv(l_file, sep=',', header=None, index_col = 0  )
        cg = set(cancer_genes_list)                      # faster membership
        mask = l_df.index.isin(cg)
        cancer_l_df = l_df[mask]
        mean_cancer_l = cancer_l_df.mean(numeric_only=True)
        l_c.append(mean_cancer_l.iat[0])

        non_cancer_l_df = l_df[~mask] 
        mean_non_cancer_l = non_cancer_l_df.mean(numeric_only=True)
        l_nc.append(mean_non_cancer_l.iat[0])

        signature_name = l_file[l_file.rindex('_')+1:]
        signature.append(signature_name) 
    print(f'number of genes is {len(l_df)}')
    print(f'number of cancer genes is {len(cancer_l_df)}')
    dict = {'signature': signature, f'l_c_{species}':l_c, f'l_nc_{species}': l_nc}
    df = pd.DataFrame(dict)
    df.set_index('signature',inplace=True)
    return df

def plot_c_diff_against_nc_diff(c_diff, nc_diff, signatures, sig_weights,title,c_errors=None, nc_errrors=None):
    
    # Scatter plot
    plt.scatter(nc_diff, c_diff, s=1000 * sig_weights, color='k', alpha=1)  # 'k' is black
    if c_errors:
        plt.errorbar(nc_diff, c_diff, xerr=nc_errrors,yerr=c_errors,fmt='none',alpha=0.4)  

    #add line y=x
    plt.plot(nc_diff, nc_diff, linewidth=1, alpha =0.3, label = 'line y=x')

    # Determine plot limits
    y_min = min(c_diff)*1.1
    y_max = max(c_diff)*1.1
    x_min = min(nc_diff)*1.1
    x_max = max(nc_diff)*1.1

    # Add shaded background rectangles
    ax = plt.gca()  # get current axes
    # Bottom-left 
    ax.add_patch(Rectangle((x_min, y_min), 0 - x_min, 0 - y_min, color='red', alpha=0.1,label='c increase in mutability'))
    # Bottom-right 
    ax.add_patch(Rectangle((0, y_min), x_max - 0, 0 - y_min, color='red', alpha=0.1))
    # Top-right (Red)
    ax.add_patch(Rectangle((0, 0), x_max - 0, y_max - 0, color='green', alpha=0.1, label='c decrease in mutability'))
    # Top-left (Yellow)
    ax.add_patch(Rectangle((x_min, 0), 0 - x_min, y_max - 0, color='green', alpha=0.1))

    # Axes lines
    plt.axvline(x=0, linestyle='-', color='black')
    plt.axhline(y=0, linestyle='-', color='black')

    # Labels
    plt.xlabel('$\Delta l^g, g \in \mathbf{non\ cancer\ genes}$')
    plt.ylabel('$\Delta l^g, g \in \mathbf{cancer\ genes}$')
    plt.title(title)

    # Signature labels for most impt sigs
    for i, sig in enumerate(signatures):
        if (sig_weights[i]>0.01) or (c_diff[i] > 0.003):
            plt.text(nc_diff[i], c_diff[i], sig[3:], fontsize=12, ha='right', va='bottom', alpha=1)

    #shaing region, x<y<0
    
    # Generate x values from the intersection of y=x and y=0 (i.e., where y = x < 0)
    x_fill = np.linspace(x_min, 0, 500)
    y1 = x_fill       # lower bound: y = x
    y2 = np.zeros_like(x_fill)  # upper bound: y = 0

    # Fill between y = x and y = 0, only where y < 0 and x < y
    plt.fill_between(x_fill, y1, y2, color='yellow', alpha=0.2, label='smaller increase in c mutability than nc')

    # Set axis limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.legend(fontsize=12)
    plt.tight_layout()

    return

#PLOT 1
title = 'Real Evolution'

human_df = extract_avg_l_per_sig(TARGET_FILE, TARGET, CANCER_GENES_LIST)
HCLCA_df = extract_avg_l_per_sig(ANC_FILE, ANC, CANCER_GENES_LIST)
merged_df = human_df.join(HCLCA_df)


weight_merged_df = merged_df.join(FILTERED_SUM, how='inner')
weight_merged_df['c_Delta'] = (weight_merged_df[f'l_c_{ANC}']-weight_merged_df[f'l_c_{TARGET}'])/( weight_merged_df[f'l_c_{ANC}']+weight_merged_df[f'l_c_{TARGET}'])
weight_merged_df['nc_Delta'] = (weight_merged_df[f'l_nc_{ANC}']-weight_merged_df[f'l_nc_{TARGET}'])/(weight_merged_df[f'l_nc_{ANC}']+weight_merged_df[f'l_nc_{TARGET}'])


plt.figure(figsize=(9,9))
plot_c_diff_against_nc_diff(weight_merged_df['c_Delta'], weight_merged_df['nc_Delta'], weight_merged_df.index, weight_merged_df['norm_weight'], title)
plt.savefig(f'{PLOT_DIR}/Real evolution_Delta.png')



#plot simulations
sim_file = f"{DIR}/output/sim_{{x}}"
title = f'{no_of_sims} sims, 95% CI'

HCLCA_df = extract_avg_l_per_sig(ANC_FILE, ANC, CANCER_GENES_LIST)

def safe_delta(A, B):
    denom = A + B
    out = (A - B) / denom
    out[denom == 0] = np.nan
    return out

c_delta_sims = []
nc_delta_sims = []
#for checking dist of l rather than delta
c_sims = []
nc_sims =[]

for i in range(no_of_sims):
    human_df = extract_avg_l_per_sig(sim_file.format(x=i), i, CANCER_GENES_LIST)
    # join the constants you need (HCLCA) just for this sim
    sim = human_df[[f'l_c_{i}',f'l_nc_{i}']].join(
        HCLCA_df[[f'l_c_{ANC}',f'l_nc_{ANC}']], how='inner'
    )

    c_delta_i  = safe_delta(sim[f'l_c_{ANC}'],  sim[f'l_c_{i}'])
    nc_delta_i = safe_delta(sim[f'l_nc_{ANC}'], sim[f'l_nc_{i}'])

    c_delta_sims.append(c_delta_i.rename(i))
    nc_delta_sims.append(nc_delta_i.rename(i))

    #to check l distn
    c_sims.append(sim[f'l_c_{i}'].rename(i))
    nc_sims.append(sim[f'l_nc_{i}'].rename(i))

#to see l dist, stack column wise
c_mat  = pd.concat(c_sims, axis=1)
nc_mat = pd.concat(nc_sims, axis=1)

# Stack sims column-wise
c_delta_mat  = pd.concat(c_delta_sims, axis=1)
nc_delta_mat = pd.concat(nc_delta_sims, axis=1)

# Mean across sims (rows = sigs)
c_mean  = c_delta_mat.mean(axis=1)
nc_mean = nc_delta_mat.mean(axis=1)

# Bring in weights (and ensure indices align)
plot_df = pd.DataFrame({'c_Delta': c_mean, 'nc_Delta': nc_mean}).join(
    FILTERED_SUM[['norm_weight']], how='inner'
)

# compute absolute mean deltas first (if not already in plot_df)
c_mean  = plot_df["c_Delta"]
nc_mean = plot_df["nc_Delta"]

#could also try running with percentiles:
c_quantile = c_delta_mat.quantile([lower_quant_cutoff,upper_quant_cutoff],axis=1)
nc_quantile = nc_delta_mat.quantile([lower_quant_cutoff,upper_quant_cutoff],axis=1)

#reindex to allign
cerr_df  = c_quantile.T.reindex(plot_df.index)
ncerr_df = nc_quantile.T.reindex(plot_df.index)

#convert absolute quantiles to distances from the mean
c_lower = c_mean - cerr_df[lower_quant_cutoff]
c_upper = cerr_df[upper_quant_cutoff] - c_mean
n_lower = nc_mean - ncerr_df[lower_quant_cutoff]
n_upper = ncerr_df[upper_quant_cutoff] - nc_mean

# build error arrays for plt.errorbar
cerr  = [c_lower, c_upper]
ncerr = [n_lower, n_upper]

plt.figure(figsize=(9,9))
plot_c_diff_against_nc_diff(plot_df['c_Delta'], plot_df['nc_Delta'], plot_df.index, plot_df['norm_weight'],
    title,cerr, ncerr)
plt.savefig(f'{PLOT_DIR}/Delta_sim.png')


#save lists of singatures reduced more thna expected by sims
c_quant = c_delta_mat.quantile([0.05,0.975],axis=1)
nc_quant = nc_delta_mat.quantile([0.05,0.975],axis=1)

#reindex to allign
cerr2_df  = c_quant.T.reindex(weight_merged_df.index)
ncerr2_df = nc_quant.T.reindex(weight_merged_df.index)
check =cerr2_df[0.975]<weight_merged_df['c_Delta']

# Save these genes to a file
with open(f'{PLOT_DIR}/signatures_reduced_more_than_expected.txt', 'w') as f:
    f.write("Signatures reduced more than expected (c_Delta):\n")
    for gene in check[check].sort_index().index.values:
        f.write(f"{gene}\n")
    f.write("\nSignatures reduced more than expected (nc_Delta):\n")
    for gene in ncerr2_df[0.975][ncerr2_df[0.975] < weight_merged_df['nc_Delta']].sort_index().index.values:
        f.write(f"{gene}\n")






