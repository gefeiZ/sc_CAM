import pandas as pd
import numpy as np

def scATAC_build_metadata(snv_matrix, peak_matrix, snp_position_df, peak_position_df,
                          snv_number_of_cells=5, peak_min=0, peak_number_of_cells=5,
                          distance_threshold=100000):

    common_cells = np.intersect1d(snv_matrix.columns, peak_matrix.columns)
    
    if len(common_cells) == 0:
        raise ValueError("No common cells between SNV Matrix and Peak Matrix!")
    

    snv_matrix = snv_matrix[common_cells]
    peak_matrix = peak_matrix[common_cells]

    snv_list = snv_matrix.index
    peak_list = peak_matrix.index
    cell_list = snv_matrix.columns

    snv_peak_metadata = pd.DataFrame({
        'SNPid': [''] * len(snv_list),
        'Num_cells_ref': [0] * len(snv_list),
        'Num_cells_alt': [0] * len(snv_list),
        'Ref_cells': [''] * len(snv_list),
        'Alt_cells': [''] * len(snv_list),
        'PeakList': [''] * len(snv_list),
        'Num_peak': [0] * len(snv_list),
        'CellList': [''] * len(snv_list)
    })
    
    useful_snv = 0

    for snvid in snv_list:
        cell_ref = snv_matrix.columns[snv_matrix.loc[snvid] == 0]  # REF group cells
        cell_alt = snv_matrix.columns[snv_matrix.loc[snvid] == 1]  # ALT group cells

        if len(cell_ref) >= snv_number_of_cells and len(cell_alt) >= snv_number_of_cells:
            
         
            valid_peaks = []
            for peak in peak_list:
                cell_ref_peak = peak_matrix.loc[peak, cell_ref]
                cell_alt_peak = peak_matrix.loc[peak, cell_alt]

                valid_ref_cells = np.sum(cell_ref_peak > peak_min)
                valid_alt_cells = np.sum(cell_alt_peak > peak_min)
                
             
                if valid_ref_cells >= peak_number_of_cells and valid_alt_cells >= peak_number_of_cells:
                    valid_peaks.append(peak)
                    print(f"Adding valid peak: {peak}")
                else:
                    print(f"Skipping peak: {peak}")
            
           
            if len(valid_peaks) > 0:
                snp_position = snp_position_df[snp_position_df['SNPid'] == snvid]
                chromosome = snp_position['chr'].values[0]
                snp_pos = snp_position['Position'].values[0]
                
                peak_position_df['Position'] = (peak_position_df['start'] + peak_position_df['end']) / 2
                
                nearby_peaks = peak_position_df[peak_position_df['chr'] == chromosome]
                nearby_peaks = nearby_peaks[np.abs(nearby_peaks['Position'] - snp_pos) <= distance_threshold]
                
                valid_nearby_peaks = np.intersect1d(valid_peaks, nearby_peaks['peak_ID'].values)
                
                if len(valid_nearby_peaks) > 0:
                    useful_snv += 1
                    snv_peak_metadata.loc[useful_snv-1, 'SNPid'] = snvid
                    snv_peak_metadata.loc[useful_snv-1, 'Num_cells_ref'] = len(cell_ref)
                    snv_peak_metadata.loc[useful_snv-1, 'Num_cells_alt'] = len(cell_alt)
                    snv_peak_metadata.loc[useful_snv-1, 'Ref_cells'] = ','.join(cell_ref)
                    snv_peak_metadata.loc[useful_snv-1, 'Alt_cells'] = ','.join(cell_alt)
                    snv_peak_metadata.loc[useful_snv-1, 'PeakList'] = ','.join(valid_nearby_peaks)
                    snv_peak_metadata.loc[useful_snv-1, 'Num_peak'] = len(valid_nearby_peaks)
                    snv_peak_metadata.loc[useful_snv-1, 'CellList'] = ','.join(np.concatenate([cell_ref, cell_alt]))
                else:
                    print("No valid peaks in valid_peaks!")

    snv_peak_metadata = snv_peak_metadata.iloc[:useful_snv]

    return snv_peak_metadata
