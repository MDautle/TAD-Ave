#! /usr/bin/env python3

'''
This script takes two Hi-C matrices (a control and a treatment) and calculates the average interactions within a TAD. 

Written by Madison Dautle 
Last updated 02.04.2026
Temple University, Department of Biology
'''

import pandas as pd
import numpy as np
import cooler
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from skimage.filters import threshold_li
from cooltools import insulation
from pathlib import Path
import argparse 

parser = argparse.ArgumentParser(description="Calculate average intra-TAD interactions for two Hi-C matrices") 

parser.add_argument("-hic1", '--matrix1', type=str, required=True, help="The path and file name of the first Hi-C matrix. This is the matrix that the TAD boundaries will be identified from") 
parser.add_argument("-hic2", '--matrix2', type=str, required=True, help="The path and file name of the second Hi-C matrix. This is the comparison condition") 
parser.add_argument("-o", '--output', type=str, required=True, help="Location where you want the results to be saved") 
parser.add_argument("-r", "--resolution", type=int, default=20000, help="The bin size for the TAD boundary identification") 

args = parser.parse_args()

def remove_outliers_iqr_multicol(df, columns):
    for column in columns:
        Q1 = df[column].quantile(0.25)
        Q3 = df[column].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        df = df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]
    return df

def calculateIntraTADS(region_TADs, clr): 
    stats_df = pd.DataFrame(index = ['Mean', 'Standard_Deviation','Median', 'Lower_Quartile','Upper_Quartile', 'Min_Val', 'Max_Val'])
    for k in range(0,len(region_TADs)): 
        region = region_TADs.iloc[k,0:3].tolist()
        region_array = clr.matrix(balance=False).fetch(region)
        zero_count = np.sum(region_array == 0) # how many zero elements there are to start
        #set number of boundary bins, This is where TADS from control are called (TADs_df.iloc[k,3] is the control intensity
        if region_TADs.iloc[k,3] <= 5: # Too small to remove diagonal and edges and have an average taken, skip this TAD
            continue
        elif 5 < region_TADs.iloc[k,3] <= 10: 
            boundary_bins = 1
        elif 10 < region_TADs.iloc[k,3] <= 20: 
            boundary_bins = 2
        elif 20 < region_TADs.iloc[k,3] <= 30:
            boundary_bins = 3
        elif region_TADs.iloc[k,3] > 30:
            boundary_bins = 4

        #remove the diagonal and "boundary_bins" from the diagonal
        n = region_array.shape[0]  
        for i in range(n):
            region_array[i, i] = 0
            for j in range(1,boundary_bins+1):
                if i + j < n:
                    region_array[i, i + j] = 0
        rows_to_remove = boundary_bins - 1
        columns_to_remove = -1 - boundary_bins 
        region_array[rows_to_remove, :] = 0 # remove boundary_bins as rows
        region_array[:, columns_to_remove] = 0 # remove boundary_bins as columns
        region_array[np.tril_indices(region_array.shape[0], -1)] = 0 # Make lower triangle all zero values
        region_array = region_array[region_array != 0] # Remove all zeroes from matrix
        region_vals = region_array.flatten().tolist() #change to list 
        for n in range(zero_count):
            region_vals.append(0) #add back appropriate number of zero values from original matrix read in
 
        #calculate statistics
        mean = np.mean(region_vals)
        std = np.std(region_vals)
        median = np.median(region_vals)
        lowerQ = np.percentile(region_vals, 25)
        upperQ = np.percentile(region_vals, 75)
        min_val = np.min(region_vals)
        max_val = np.max(region_vals)

        statistics = [mean, std, median, lowerQ, upperQ, min_val, max_val]
        colName = f'{region_TADs.iloc[k,0]}:{region_TADs.iloc[k,1]}-{region_TADs.iloc[k,2]}'
        stats_df[colName] = statistics
    
    return stats_df

directory_path = Path(args.output)
directory_path.mkdir(parents=True, exist_ok=True) 

conditions = {
    "matrix1": args.matrix1, 
    "matrix2": args.matrix2
}

clrs = {}
for label, path in conditions.items():
    clr = cooler.Cooler(f'{path}::/resolutions/{args.resolution}')
    clrs[label] = clr

windows = [3*args.resolution, 5*args.resolution, 10*args.resolution]
w = windows[2]

insulation_tables = {}
for label, path in conditions.items():
    clr = clrs[label]
    insulation_table = insulation(clr, windows, verbose=False)
    insulation_table = insulation_table[~np.isnan(insulation_table[f'boundary_strength_{w}'])]
    insulation_tables[label] = insulation_table

thresholds_li = {}
boundaries_li = {}
insulation_table = insulation_tables["matrix1"]

thresholds_li = threshold_li(insulation_table[f'boundary_strength_{w}'].dropna().values)
n_boundaries_li = (insulation_table[f'boundary_strength_{w}'].dropna()>=thresholds_li).sum()
boundaries_li = insulation_table[(insulation_table[f'boundary_strength_{w}']>=thresholds_li)]

pd.DataFrame(boundaries_li).to_csv(f"{args.output}/matrix1_TADboundaries.csv",index=False,sep='\t')

# Keep only "good" bins
boundaries_li = boundaries_li.query('is_bad_bin == False')

chromosomes = boundaries_li['chrom'].unique()
TADs = boundaries_li[['chrom']]
TADs_df = pd.DataFrame()

filein = f"{args.output}/matrix1_TADboundaries.csv"
boundaries_li = pd.read_csv(filein, sep='\t')
for chromosome in chromosomes: 
    chr_df = boundaries_li.query('chrom == @chromosome')
    result = pd.DataFrame({'chrom': chr_df['chrom'], 'pos1': chr_df['end'], 'pos2': chr_df['start'].shift(-1).fillna(0).astype('int')})
    result = result.iloc[:-1]
    TADs_df = pd.concat([TADs_df, result], ignore_index=True)
    del result
TADs_df.to_csv(f'{args.output}/matrix1_TADs_df.tsv', header=False, index=False, sep='\t')

TADs_df["matrix1"] = ((TADs_df['pos2']-TADs_df['pos1'])/args.resolution).astype(int) #Bins in each TAD
bins = len(TADs_df)

for label, path in conditions.items(): 
    clr = clrs[label] 
    overallStats_df=calculateIntraTADS(TADs_df, clr)
    fileout = f'{args.output}/{label}-statistics.tsv'
    overallStats_df.to_csv(fileout, header = True, index = True, sep='\t')


average_df = pd.DataFrame()
for label, path in conditions.items(): 
    filein = f'{args.output}/{label}-statistics.tsv'
    temp = pd.read_csv(filein, sep='\t', index_col=0)

    # Extract the 'Mean' row and make it a Series with TADs as index
    if 'Mean' in temp.index:
        mean_series = temp.loc['Mean']
        mean_series.name = label  # Rename series to match the condition
        average_df = pd.concat([average_df, mean_series], axis=1)
    else:
        print(f"'Mean' not found in {filein}, skipping...")

# Save the final merged DataFrame
filename = f'{args.output}/intraTAD-intensity_rawData.tsv'
average_df.to_csv(filename, sep='\t')


# Prepare stats
stat_full = []
p_value_full = []
comparison_labels = []

sub_df = remove_outliers_iqr_multicol(average_df, conditions.keys())

for label, path in conditions.items():
    if label == 'matrix1':
        continue
    stat, p_value = stats.ranksums(sub_df['matrix1'], sub_df[label])
    stat_full.append(stat)
    p_value_full.append(p_value)
    comparison_labels.append(label)

# Reshape for seaborn
long_df = sub_df.melt(var_name='Sample', value_name='Value')
long_df['Sample'] = long_df['Sample'].map(conditions)


# Calculate and print medians and means
summary_df = long_df.groupby('Sample')['Value'].agg(['median', 'mean']).reset_index()
summary_df.columns = ['Sample', 'Median', 'Mean']
print("Median and Mean per Sample:")
print(summary_df.to_string(index=False)) 

# Plot
linewidth = 1
plt.figure(figsize=(10, 6))
ax = sns.violinplot(
    x='Sample',
    y='Value',
    data=long_df,
    palette=['blue','red'],
    inner=None,
    linewidth=linewidth,
)
for violin in ax.collections:
    violin.set_edgecolor('black')
    violin.set_linewidth(linewidth)

sns.boxplot(
    x='Sample',
    y='Value',
    data=long_df,
    width=0.1,
    showmeans=False,
    meanprops={
        'marker': 'X',
        'markerfacecolor': 'red',
        'markeredgecolor': 'white',
        'markersize': 6,
        'zorder': 100
    },
    boxprops={'facecolor': 'none', 'edgecolor': 'black', 'linewidth': linewidth, 'zorder': 3},
    whiskerprops={'color': 'black', 'linewidth': linewidth, 'zorder': 3},
    capprops={'color': 'black', 'linewidth': linewidth, 'zorder': 3},
    medianprops={'color': 'black', 'zorder': linewidth, 'linewidth':linewidth},
    flierprops={'marker': 'o', 'markerfacecolor':'black', 'markeredgecolor':None, 'markersize': 2}
)

ax.set_ylabel('Short Range Interactions within TADs',fontweight='bold')
ax.set_xlabel('')

# Add p-values above each comparison 
y_max = long_df['Value'].max()
y_offset = (y_max - long_df['Value'].min()) * 0.02
for cond, p in zip(comparison_labels, p_value_full):
    xpos = list(conditions.keys()).index(cond)  # use the label order
    ax.text(
         xpos,
         y_max + y_offset,
         f"p = {p:.3e}",
         ha='center',
         va='bottom',
         fontsize=10,
         fontweight='bold'
     )
    
#plt.title("Intra-TAD Intensity Across Samples", fontsize=14)
plt.xlabel("")
plt.ylabel("Intra-TAD Intensity", fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{args.output}/intraTAD_violins.pdf', dpi=300, bbox_inches='tight')

stat_results = []
stat, p_value = stats.ranksums(sub_df['matrix1'], sub_df['matrix2'])
stat_results.append({
    'Condition 1': conditions['matrix1'],
    'Condition 2': conditions['matrix2'],
    'Statistic': stat,
    'P-value': p_value
})

# Convert to DataFrame
results_df = pd.DataFrame(stat_results)
results_df.to_csv(f'{args.output}/statistics.csv', header=True, index=False)