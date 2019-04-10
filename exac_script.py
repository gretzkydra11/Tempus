import pandas as pd, subprocess, ast, os

# Define local path to read/write files
current_wd = os.getcwd()
print("working directory:",current_wd)

#### 5-6. Extracting data from EXAC Database ###
print("Extracting EXAC data from exac.hms.harvard.edu/rest/variant/variant/")

# if running code in chunks, we can import the previously generated dataframe now
df = pd.read_csv(current_wd+"/VCF_DF.txt", sep='\t')

exac_df =pd.DataFrame()
exac_dict = {}
s= "-"
d={}
EXAC = "exac.hms.harvard.edu/rest/variant/variant/"

# For each variant, download EXAC data and write to a dataframe. Then make a dictionary of those dataframes
# where the key is the original dataframe's index.
for i in range(0, len(df)-1):
    variant: str = s.join([str(df['CHROM'][i]),str(df['POS'][i]),str(df['REF'][i]),str(df['ALT'][i])])
    site = EXAC+variant
    subprocess.call(('wget' + ' ' + site), shell=True)
    file_loc = current_wd+"/"+variant
    with open(file_loc, "r") as file:
        d = ast.literal_eval(file.read())
    exac_dict[i] = d
    os.remove(file_loc)

"""
Once, while running the code on my Mac (OS X), I encountered the following error:
...
BlockingIOError: [Errno 35] Resource temporarily unavailable

... presumably due to a slow/broken internet connection. In this case, the EXAC.py script can be run to avoid
needing to re-run all previous steps
"""
# Parse the dictionary of dataframes containing variant-specific exac data
# define temp objects/variables
temp_dict = {}
exac_df_parse = pd.DataFrame()
art_id = 0
print("Parsing EXAC data...")
# loop through exac_dict indices and store relevant information as a DataFrame, then finally combine with previous
# DataFrame to generate the final table and export to working directory
for i in exac_dict:
    if 'allele_freq' in exac_dict[i]:
        temp_dict = {art_id: [exac_dict[i]['allele_freq'],exac_dict[i]['rsid'],exac_dict[i]['vep_annotations']]}
    else:
        temp_dict = {art_id: ["NA","NA","NA"]}
    exac_df_parse = exac_df_parse.append(pd.DataFrame.from_dict(temp_dict, orient='index'))
    art_id += 1

exac_df_parse.columns = ['Allele frequency','rsid','VEP annotation']
print("Assembling final table...")
final_df = pd.concat([df,exac_df_parse],axis=1)
final_df.to_csv(current_wd+"/final_table.txt", sep='\t')