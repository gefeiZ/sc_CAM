import numpy as np
import pandas as pd
from scipy.stats import chi2
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import NegativeBinomial
import statsmodels.api as sm
from statsmodels.tools import add_constant
from statsmodels.stats.multitest import multipletests

def callQTL(peakMatrix, snv_peak_pair_metadata, p_adjust_method="fdr_bh", zero_threshold=0.5, overdispersion_threshold=0.1):
    
    """
    Perform caQTL analysis by testing association between SNP genotypes (Ref/Alt groups)
    and single-cell ATAC-seq peak accessibility using either Zero-Inflated Poisson (ZIP)
    or Negative Binomial (NB) regression, followed by multiple testing correction.

    Parameters
    ----------
    peakMatrix : pandas.DataFrame
        A peak-by-cell count matrix (accessibility), where rows are peak IDs and columns are cell barcodes.
        All values must be non-negative integers. Matrix must not contain NaNs.

    snv_peak_pair_metadata : pandas.DataFrame
        Metadata for SNP-peak pairs to test. Must contain the following columns:
            - 'SNVid': SNP identifier.
            - 'CellList': comma-separated string of all relevant cells for this SNP.
            - 'Ref_cells': comma-separated string of cell barcodes with reference allele.
            - 'Alt_cells': comma-separated string of cell barcodes with alternative allele.
            - 'PeakList': comma-separated string of peaks to be tested for this SNP.

    p_adjust_method : str, optional (default="fdr_bh")
        Method for multiple testing correction. Options include:
            - "fdr_bh" (Benjamini-Hochberg)
            - "bonferroni"
            - "fdr_by"
            - Any method supported by `statsmodels.stats.multitest.multipletests`.

    zero_threshold : float, optional (default=0.5)
        Threshold for proportion of zero counts. If either group has > `zero_threshold` zeros,
        a Zero-Inflated Poisson model will be used.

    overdispersion_threshold : float, optional (default=0.1)
        Threshold for overdispersion (variance > mean). If either group has overdispersion
        greater than this threshold, ZIP model will be used.

    Returns
    -------
    results_df : pandas.DataFrame
        A long-format DataFrame where each row corresponds to a SNP-peak pair test.
        Columns include:
            - 'SNVid', 'Peakid'
            - sample sizes: 'sample_size_Ref', 'sample_size_Alt'
            - means: 'total_mean_Ref', 'total_mean_Alt'
            - fitted parameters: ('theta', 'mu', or 'size') for both Ref and Alt groups
            - 'chi2LR1': likelihood ratio test statistic
            - 'pvalue': raw p-value from test
            - 'adjusted_pvalue': p-value corrected for multiple testing
            - 'Remark': error message if the test failed (e.g., due to fitting issues)

    Notes
    -----
    - The function chooses between ZIP and NB model for each test dynamically,
      based on zero inflation and overdispersion.
    - Model comparison is performed using likelihood ratio test (LRT) between group-specific models
      and a combined (null) model.
    - p-value adjustment is applied globally across all tested SNP-peak pairs.
    - Internally uses `statsmodels` for GLM/ZIP fitting and `scipy.stats.chi2.sf` for LRT p-values.

    """


    if not isinstance(peakMatrix, (pd.DataFrame, np.ndarray)):
        raise ValueError("Wrong data type of 'peakMatrix'")
    if pd.isna(peakMatrix).any().any():
        raise ValueError("NA detected in 'peakMatrix'")
    if (peakMatrix < 0).any().any():
        raise ValueError("Negative value detected in 'peakMatrix'")
    if (peakMatrix == 0).all().all():
        raise ValueError("All elements of 'peakMatrix' are zero")

    def CallcaQTL(i):
        snvid = snv_peak_pair_metadata.loc[i, "SNVid"]
        cells_character = snv_peak_pair_metadata.loc[i, "CellList"]
        cells = cells_character.split(",")
        ref_cells_character = snv_peak_pair_metadata.loc[i, "Ref_cells"]
        ref_cells = ref_cells_character.split(",")
        alt_cells_character = snv_peak_pair_metadata.loc[i, "Alt_cells"]
        alt_cells = alt_cells_character.split(",")
        peaks = snv_peak_pair_metadata.loc[i, "PeakList"].split(",")
        
        results_SNV = []

        for peak in peaks:
            counts_1 = peakMatrix.loc[peak, ref_cells].values.flatten()
            counts_2 = peakMatrix.loc[peak, alt_cells].values.flatten()
            
            totalMean_1 = np.mean(counts_1)
            totalMean_2 = np.mean(counts_2)
            
            results_peak = {
                "SNVid": snvid,
                "Peakid": peak,
                "sample_size_1": len(counts_1),
                "sample_size_2": len(counts_2),
                "theta_1": None,
                "theta_2": None,
                "mu_1": None,
                "mu_2": None,
                "size_1": None,
                "size_2": None,
                "prob_1": None,
                "prob_2": None,
                "total_mean_1": totalMean_1,
                "total_mean_2": totalMean_2,
                "chi2LR1": None,
                "pvalue": None,
                "adjusted_pvalue": None,
                "Remark": None
            }
            
            zero_proportion_1 = np.mean(counts_1 == 0)
            zero_proportion_2 = np.mean(counts_2 == 0)


            def calculate_overdispersion(counts):
                mean_counts = np.mean(counts)
                var_counts = np.var(counts)
                return (var_counts - mean_counts) / (mean_counts ** 2)

            overdispersion_1 = calculate_overdispersion(counts_1)
            overdispersion_2 = calculate_overdispersion(counts_2)

            if zero_proportion_1 > zero_threshold or zero_proportion_2 > zero_threshold or overdispersion_1 > overdispersion_threshold or overdispersion_2 > overdispersion_threshold:
        
                try:
                    fit_1 = sm.ZeroInflatedPoisson(counts_1, np.ones((len(counts_1), 1))).fit(disp=False)
                    fit_2 = sm.ZeroInflatedPoisson(counts_2, np.ones((len(counts_2), 1))).fit(disp=False)
                    theta_1 = fit_1.params[0]
                    lambda_1 = np.exp(fit_1.params[1])
                    theta_2 = fit_2.params[0]
                    lambda_2 = np.exp(fit_2.params[1])
                    
                    logL_full = fit_1.llf + fit_2.llf
                    combined_counts = np.concatenate([counts_1, counts_2])
                    combined_fit = sm.ZeroInflatedPoisson(combined_counts, np.ones_like(combined_counts)).fit(disp=False)
                    logL_null = combined_fit.llf
                    
                    chi2LR1 = 2 * (logL_full - logL_null)
                    pvalue = chi2.sf(chi2LR1, df=2)
                    
                    results_peak.update({
                        "theta_1": theta_1,
                        "theta_2": theta_2,
                        "mu_1": lambda_1,
                        "mu_2": lambda_2,
                        "chi2LR1": chi2LR1,
                        "pvalue": pvalue
                    })
                except Exception as e:
                    results_peak["Remark"] = str(e)
                    results_SNV.append(results_peak)
                    continue
            else:
               
                try:
                    fit_1 = NegativeBinomial(counts_1, np.ones((len(counts_1), 1))).fit(disp=False)
                    fit_2 = NegativeBinomial(counts_2, np.ones((len(counts_2), 1))).fit(disp=False)
                    size_1 = fit_1.params[0]
                    size_2 = fit_2.params[0]
                    mu_1 = np.exp(fit_1.params[1])
                    mu_2 = np.exp(fit_2.params[1])
                    
                    logL_full = fit_1.llf + fit_2.llf
                    combined_counts = np.concatenate([counts_1, counts_2])
                    combined_fit = NegativeBinomial(combined_counts, np.ones_like(combined_counts)).fit(disp=False)
                    logL_null = combined_fit.llf
                    
                    chi2LR1 = 2 * (logL_full - logL_null)
                    pvalue = chi2.sf(chi2LR1, df=2)
                    
                    results_peak.update({
                        "size_1": size_1,
                        "size_2": size_2,
                        "mu_1": mu_1,
                        "mu_2": mu_2,
                        "chi2LR1": chi2LR1,
                        "pvalue": pvalue
                    })
                except Exception as e:
                    results_peak["Remark"] = str(e)
                    results_SNV.append(results_peak)
                    continue
            
            results_SNV.append(results_peak)
        
        return results_SNV

    results = []
    for i in range(len(snv_peak_pair_metadata)):
        print(f"Processing {i + 1} / {len(snv_peak_pair_metadata)} term")
        results.extend(CallcaQTL(i))

    results_df = pd.DataFrame(results)
    

    results_df["adjusted_pvalue"] = multipletests(results_df["pvalue"], method=p_adjust_method)[1]

    results_df.rename(columns={
        'sample_size_1': 'sample_size_Ref',
        'sample_size_2': 'sample_size_Alt',
        'theta_1': 'theta_Ref',
        'theta_2': 'theta_Alt',
        'mu_1': 'mu_Ref',
        'mu_2': 'mu_Alt',
        'size_1': 'size_Ref',
        'size_2': 'size_Alt',
        'prob_1': 'prob_Ref',
        'prob_2': 'prob_Alt',
        'total_mean_1': 'total_mean_Ref',
        'total_mean_2': 'total_mean_Alt'
    }, inplace=True)
    
    return results_df
