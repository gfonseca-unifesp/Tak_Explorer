# required libraries for TaK function
import pandas as pd
import numpy as np
# support libraries for visualization
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from adjustText import adjust_text

### The Taxonomic Knowledge - TaK function
def TaK(input_matrix, group_col="Phylum", dataset_col="Station", count_col="individualCount",
        taxonomic_ranks=None, simple_weight=True, 
        custom_weights=None):
    """
    Calculates Taxonomic Knowledge indices (TC and TR).
    
    Parameters:
    -----------

    input_matrix : pd.DataFrame
        Matrix containing taxonomic levels and abundance.
    group_col : str
        Higher-order group for analysis (e.g., "Phylum" or "Class").
        Defaults to "Phylum".
    dataset_col : str
        Column to group datasets (e.g., "Station", "Year").
        Defaults to "Station".
    count_col : str
        Abundance column. If NaN, assumes 1 individual.
        Defaults to "individualCount"
    taxonomic_ranks : list, optional
        Ranks to include in weight calculation.
    simple_weight : bool
        If True, all ranks have weight 1.
        If False, uses linear or custom weights.
    custom_weights : array-like, optional
        Specific weights for each rank in taxonomic_ranks.
        If empty and simple_weight is False a linear weight is used.
    """
    
    # --- 0. Setup ---
    df = input_matrix.copy()
    
    # --- 1. Set the taxonomic ranks ---
    # If not provided, it is identified as all columns that are not group_col, dataset_col or count_col
    if taxonomic_ranks is None:
        tax_cols = [c for c in df.columns if c not in [group_col,dataset_col,count_col]]
    else:
        tax_cols = taxonomic_ranks
    # Get the number of taxonomic levels
    n_levels = len(tax_cols)
    
    # --- 1. Set the weighting scheme ---
    # Default is simple, i.e. all taxonomic ranks have same weight
    # If not simple, set a provided custom weighting scheme
    # If custom_weights is not provided, set a linear weighting scheme
    if simple_weight == True:
        weight_vector = np.ones(n_levels)
    elif custom_weights is not None:
        weight_vector = np.array(custom_weights)
    else:
        weight_vector = np.arange(1, n_levels + 1)
    # calculate the value of integrated weight, i.e. the value to divide to obtain 1
    max_potential_weight = sum(weight_vector)

    # --- 3. Identify abundance column ---
    # In case the number has been read as str
    # Assumes that absence of count means only 1 individual
    ind_count = (df[count_col].values).astype(float)
    ind_count[np.isnan(ind_count)] = 1

    # --- 4. Row-level calculation ---
    # Non-empty cells are considered 1, multiplied by the respective weight, and then all columns of the row are added together
    row_sums_weighted = df[tax_cols].notna().astype(int).dot(weight_vector)
    df['Row_Sum_W'] = row_sums_weighted
    df['RSW*n'] = df['Row_Sum_W']*df[count_col]
    df['RSP*n'] = max_potential_weight*df[count_col]
    df['indv_TC'] = max_potential_weight*df[count_col]

    # --- 5. TC Calculation (Completeness) ---
    # Count number of unique entries for each group and station/dataset and then adds the weighted sums for each case
    tc_tr_data = df.groupby([dataset_col, group_col]).agg(
        n_lines=(group_col, 'size'),
        n_individuals=(count_col, 'sum'),
        Weighted_Sum=('Row_Sum_W', 'sum'),
        Weighted_Sum_Total=('RSW*n', 'sum'),
        Weighted_Sum_Potential=('RSP*n', 'sum') 
    ).reset_index()
    # TC considers the number of individuals, so it is the ratio (sum RSW*n)/(sum RSP*n)
    tc_tr_data['TC'] = tc_tr_data['Weighted_Sum_Total'] / (tc_tr_data['Weighted_Sum_Potential'])

    # --- 6. TR Calculation (Resolution) ---
    # Resolution considers just each lineage, regardless of relative abundance
    tc_tr_data['TR'] = tc_tr_data['Weighted_Sum']/(tc_tr_data['n_lines']*max_potential_weight) 
    
    # --- 7. Prepare the output
    # Return variables are the summary table, 
    summary_dataset = tc_tr_data.copy().drop(columns=['Weighted_Sum', 'Weighted_Sum_Total', 'Weighted_Sum_Potential'])
    summary_dataset = summary_dataset.rename(columns={dataset_col: 'Dataset', group_col: 'Group'})
    # Descriptive statistics variables
    stats_list = summary_dataset.groupby('Group').agg(
        Total_Lineages=('n_lines', 'sum'),
        Total_Individuals=('n_individuals', 'sum'),
        TR_Median=('TR', 'median'),
        TC_Median=('TC', 'median')
    ).reset_index()

    # Return the sumary     
    return {
        "summary_matrix": summary_dataset,
        "descriptive_stats": stats_list
    }

#### Example of the function usage

# --- Load Data ---
input_file = 'raw_data.csv'
taxa_data = pd.read_csv(input_file, sep=';')

# --- Run the TaK function ---
results = TaK(
    taxa_data,
    group_col="Phylum",
    dataset_col="Dataset",
    taxonomic_ranks=["Phylum","Class", "Order", "Family", "Genus", "Species"],
    simple_weight=False
)

# --- Verification Table ---
print("Performance Summary by Phylum Group:")
# Selecionando colunas específicas para o resumo
median_tak = results['descriptive_stats'][['Group', 'Total_Lineages', 'Total_Individuals','TR_Median', 'TC_Median']]
print(median_tak)

# --- Plotting ---
# Dumping relevant data for plotting in new variables
plot_df = results['summary_matrix']
tc_med = median_tak['TC_Median'].median()
tr_med = median_tak['TR_Median'].median()
# Create figure
fig,ax = plt.subplots(figsize=(12, 8))
sns.set_style("white")
# Sets the limits of background quadrants
x1 = [0, 0.5]
x2 = [0.5, 1]
y0 = [0,0]
ym = [0.5,0.5]
y1 = [1,1]
# Fills the background quadrants
ax.fill_between(x1,ym,y1, color='lightblue', alpha=0.1)  # Upper Left
ax.fill_between(x2,ym,y1, color='lightgreen', alpha=0.1) # Upper Right
ax.fill_between(x1,y0,ym, color='lightpink', alpha=0.1)  # Lower Left
ax.fill_between(x2,y0,ym, color='khaki', alpha=0.1)      # Lower Right
# The separation lines
ax.axvline(.5, linestyle=':', color='grey',lw=0.75)
ax.axhline(.5, linestyle=':', color='grey',lw=0.75)
# Scatterplot of the data
scatter = sns.scatterplot(
    data=plot_df, x='TC', y='TR', hue='Group', 
    size='n_individuals', sizes=(20, 300),
    alpha=0.7,
    palette='viridis',
    clip_on=False # the symbol will be drawn outside of axis region
)
# Text labels of the data
texts = []
for i in range(len(plot_df)):
    texts.append(plt.text(
        plot_df.TC[i], plot_df.TR[i], plot_df.Dataset[i], 
        fontsize=9, color='black'
    ))
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
# Title and labels
ax.set_title("Taxonomic Knowledge", fontsize=15,pad=15)
ax.set_xlabel("Taxonomic Completeness (TC)")
ax.set_ylabel("Taxonomic Resolution (TR)")
ax.set_xlim(0, 1.)
ax.set_ylim(0, 1.) 
# Diagonal line
ax.plot([0,1.],[0,1.],'--',color='gray',lw=0.75, zorder=0)
# Ticks and gridlines
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.25))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.125))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.125))
ax.grid(True, which='major', axis='both', linestyle='-', linewidth=0.5, color='gray', alpha=0.3)
ax.grid(True, which='minor', axis='both', linestyle='-', linewidth=0.5, color='gray', alpha=0.3)
# Hide the top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# Final adjustments
ax.set_aspect('equal', adjustable='box')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
fig.tight_layout()
fig.savefig('tak_taxa_test.pdf',bbox_inches='tight',pad_inches=.2)
plt.close('all')
