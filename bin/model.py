import pandas as pd
import numpy as np
from scipy.stats import chi2, nbinom, poisson
from statsmodels.discrete.count_model import ZeroInflatedPoisson
from statsmodels.api import GLM, families
from statsmodels.discrete.discrete_model import NegativeBinomial

def caQTL(peakMatrix, snv_peak_pair_metadata, p_adjust_method="fdr_bh", zero_threshold=0.5):
 
    if not isinstance(peakMatrix, (pd.DataFrame, np.ndarray)):
        raise ValueError("Wrong data type of 'peakMatrix'")
    if pd.isna(peakMatrix).any().any():
        raise ValueError("NA detected in 'peakMatrix'")
    if (peakMatrix < 0).any().any():
        raise ValueError("Negative value detected in 'peakMatrix'")
    if (peakMatrix == 0).all().all():
        raise ValueError("All elements of 'peakMatrix' are zero")


    def CallcaQTL(i):
        snvid = snv_peak_pair_metadata.loc[i, "SNPid"]
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
            
   
            zero_proportion_1 = np.mean(counts_1 == 0)
            zero_proportion_2 = np.mean(counts_2 == 0)


            if zero_proportion_1 > zero_threshold or zero_proportion_2 > zero_threshold:
          
                fit_1 = ZeroInflatedPoisson(counts_1, np.ones_like(counts_1)).fit(disp=False)
                fit_2 = ZeroInflatedPoisson(counts_2, np.ones_like(counts_2)).fit(disp=False)

                theta_1 = fit_1.params['inflate']
                lambda_1 = np.exp(fit_1.params['const'])
                theta_2 = fit_2.params['inflate']
                lambda_2 = np.exp(fit_2.params['const'])

                logL_full = fit_1.llf + fit_2.llf
                combined_counts = np.concatenate([counts_1, counts_2])
                combined_fit = ZeroInflatedPoisson(combined_counts, np.ones_like(combined_counts)).fit(disp=False)
                logL_null = combined_fit.llf
            else:
          
                fit_1 = NegativeBinomial(counts_1, np.ones_like(counts_1)).fit(disp=False)
                fit_2 = NegativeBinomial(counts_2, np.ones_like(counts_2)).fit(disp=False)

                size_1 = fit_1.params['alpha']
                size_2 = fit_2.params['alpha']
                mu_1 = np.exp(fit_1.params['const'])
                mu_2 = np.exp(fit_2.params['const'])

                logL_full = fit_1.llf + fit_2.llf
                combined_counts = np.concatenate([counts_1, counts_2])
                combined_fit = NegativeBinomial(combined_counts, np.ones_like(combined_counts)).fit(disp=False)
                logL_null = combined_fit.llf
            

            chi2LR = 2 * (logL_full - logL_null)
            pvalue = chi2.sf(chi2LR, df=2)

        
            results_SNV.append({
                "SNVid": snvid,
                "Peakid": peak,
                "sample_size_1": len(counts_1),
                "sample_size_2": len(counts_2),
                "theta_1": theta_1 if zero_proportion_1 > zero_threshold else None,
                "theta_2": theta_2 if zero_proportion_2 > zero_threshold else None,
                "mu_1": lambda_1 if zero_proportion_1 > zero_threshold else mu_1,
                "mu_2": lambda_2 if zero_proportion_2 > zero_threshold else mu_2,
                "size_1": None if zero_proportion_1 > zero_threshold else size_1,
                "size_2": None if zero_proportion_2 > zero_threshold else size_2,
                "prob_1": 1 - theta_1 if zero_proportion_1 > zero_threshold else None,
                "prob_2": 1 - theta_2 if zero_proportion_2 > zero_threshold else None,
                "chi2LR1": chi2LR,
                "pvalue": pvalue,
                "adjusted_pvalue": None
            })

  
        pvalues = [result["pvalue"] for result in results_SNV]
        adjusted_pvalues = multipletests(pvalues, method=p_adjust_method)[1]
        for j, adj_p in enumerate(adjusted_pvalues):
            results_SNV[j]["adjusted_pvalue"] = adj_p

        return results_SNV


    results = []
    for i in range(len(snv_peak_pair_metadata)):
        results.extend(CallcaQTL(i))

    return pd.DataFrame(results)
