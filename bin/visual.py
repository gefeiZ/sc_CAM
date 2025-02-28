import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def QTL_visual(eqtl_result, expressionMatrix, snv_gene_pair_metadata, figure_output=None):
    SNVid = eqtl_result['SNVid']
    geneid = eqtl_result['Peakid']
    adjust_pvalue = eqtl_result['adjusted_pvalue']
    

    ref_cells = snv_gene_pair_metadata[snv_gene_pair_metadata['SNPid'] == SNVid]['Ref_cells'].values[0]
    ref_cells = ref_cells.split(',')
    alt_cells = snv_gene_pair_metadata[snv_gene_pair_metadata['SNPid'] == SNVid]['Alt_cells'].values[0]
    alt_cells = alt_cells.split(',')

    ref_cells = list(set(ref_cells) & set(expressionMatrix.columns))
    alt_cells = list(set(alt_cells) & set(expressionMatrix.columns))
    

    counts_Ref = []
    counts_Alt = []
    

    if len(ref_cells) > 0:
        counts_Ref = expressionMatrix.loc[geneid, ref_cells].values
    

    if len(alt_cells) > 0:
        counts_Alt = expressionMatrix.loc[geneid, alt_cells].values
    

    if len(counts_Ref) == 0 and len(counts_Alt) == 0:
        print(f"No valid cells found for SNPid: {SNVid} and geneid: {geneid}")
        return None

    mean_Ref = len(counts_Ref) > 0 if len(counts_Ref) > 0 else np.nan
    mean_Alt = np.mean(counts_Alt) if len(counts_Alt) > 0 else np.nan
    

    df = pd.DataFrame({
        'Genotype': ['REF', 'ALT'],
        'Mean Accessibility': [mean_Ref, mean_Alt]
    })
    

    plt.figure(figsize=(10, 6))
    sns.barplot(x='Genotype', y='Mean Accessibility', data=df, palette='Set2')
    plt.title(f'SNV: {SNVid}    Peak: {geneid}', fontsize=14)
    plt.xlabel('Genotype', fontsize=14)
    plt.ylabel('Mean Accessibility', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    

    if figure_output is not None:
        plt.savefig(figure_output)
    
    plt.show()
    return plt

