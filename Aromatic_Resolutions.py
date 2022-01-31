import numpy as np
import pandas as pd

root_path         = "/home/tnchevez/aromatic_res"
acbap_filepath    = f'{root_path}/Aromatic_AC-BAP.csv'
pattern_filepath  = f'{root_path}/Aromatic_Output.csv'
is_practice       = False
output_filepath   = f'{root_path}/resolution-Practice-output.csv'
columns           = ["Cl","Br","F","BrCl","ClF","BrF","BrClF","Other"]

if is_practice :
  acbap_filepath    = f'{root_path}/733_resolution-practice.csv'
  pattern_filepath  = f'{root_path}/pattern-practice.csv'

def readData(acbap_filepath, pattern_filepath):
  acbap_df   = pd.read_csv(acbap_filepath)
  pattern_df = pd.read_csv(pattern_filepath)
  return acbap_df, pattern_df

def comparePeaksArr(peaks1, peaks2):
  results = list()
  for p1 in peaks1:
    for p2 in peaks2:
      mass_diff = abs(p1-p2)
      if mass_diff < 1 and mass_diff > 0:
        mass2 = p1 if p1 > p2 else p2
        resolution = mass2 / mass_diff
        results.append(resolution)
  return results

def getLowestRes():
  acbap_df, pattern_df = readData(acbap_filepath, pattern_filepath)
  results_df = pd.DataFrame([], columns=(["Formula","Peaks","Group"] + columns)) 
  
  for i_acbap, acbap in acbap_df.iterrows():
    print(f"INITIATING ACBAP CYCLE #{i_acbap}")
    # Initialize arrays
    arrs = dict({
      "Other": list(),
      "Cl": list(),
      "Br":list(),
      "F":list(),
      "BrCl":list(),
      "ClF":list(),
      "BrF":list(),
      "BrClF":list()
    })

    peaks_str = acbap["Peaks"].split(",")
    peaks     = [float(x) for x in peaks_str ]
      
    for i_pat, pattern in pattern_df.iterrows():
      
      if acbap["Formula"] == pattern["Formula"]:
        continue
      
      peaks_pat_str = pattern["Peaks"].split(",")
      peaks_pat     = [float(x) for x in peaks_pat_str ]
      
      peak_comparison_arr = comparePeaksArr(peaks1=peaks, peaks2=peaks_pat)
      
      if len(peak_comparison_arr) == 0:
        continue
      
      arr_id = pattern["Group"]
      arrs[arr_id].append(min(peak_comparison_arr))

    output_data = dict({
      "Formula"  :acbap["Formula"],
      "Peaks"    :acbap["Peaks"],
      "Group"    :acbap["Group"]
    })
    
    for index in columns:
      if len(arrs[index]) > 0:
        output_data[index] = max(arrs[index])
    results_df = results_df.append(output_data, ignore_index=True)
    
  renamed_cols = dict()
  for index in columns:
    renamed_cols[index] = f"min_Res_{index}"
  results_df = results_df.rename(columns=renamed_cols)
  results_df.to_csv(output_filepath, index=False)

    
if __name__ == "__main__":
  getLowestRes()