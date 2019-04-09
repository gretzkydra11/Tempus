import vcf, pandas as pd, numpy as np, sys, subprocess, ast, os

merged_df4 = pd.read_csv("/Users/denisavey/Tempus_challenge/merged_df4.txt", sep = '\t')

exac_data =pd.DataFrame()
s= "-"
d={}
EXAC = "exac.hms.harvard.edu/rest/variant/variant/"
PATH = "/Users/denisavey/Tempus_challenge/"

for i in range(0,len(merged_df4)-1):
    #for i in (range(len(record.ALT)):
    variant = s.join([str(merged_df4['CHROM'][i]),str(merged_df4['POS'][i]),str(merged_df4['REF'][i]),str(merged_df4['ALT'][i])])
    site = EXAC+variant
    subprocess.call(('wget' + ' ' + site), shell=True)
    file_loc = PATH+variant
    with open(file_loc, "r") as file:
        d = ast.literal_eval(file.read())
    exac_data = exac_data.append(pd.DataFrame.from_dict(d,orient='index'))
    os.remove(file_loc)

exac_data.to_csv("/Users/denisavey/Tempus_challenge/exac_data.txt", sep = '\t')
