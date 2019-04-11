"""

Author: Denis Avey
Date:

This master script will generate a tab-delimited text file containing various info
parsed from the Challenge_data.vcf file, the EXAC database, and the CADD database.

Note 1: Before proceeding with the downstream steps, first complete steps outlined in CADD_instructions

Note 2: All below code was developed using Python 3.7.3

"""
import os, vcf, pandas as pd, numpy as np, subprocess, ast

# Define local path to read/write files
current_wd = os.getcwd()
print("working directory:",current_wd)

# Read in VCF file using PyVCF
vcf_rd = vcf.Reader(open(current_wd+"/Challenge_data.vcf", 'r'))
print("Parsing VCF file...")
### Loop through each row of the VCF file and store relevant information in pandas DataFrame:
temp_dict = {}
vcf_df = pd.DataFrame()
art_id = 0
for record in vcf_rd:
    temp_dict = {art_id:[record.CHROM,record.POS,record.REF,record.ALT,record.INFO['TYPE'],record.INFO['DPB'],
    record.INFO['DP'],record.INFO['RO'],record.INFO['AO']]}
    vcf_df = vcf_df.append(pd.DataFrame.from_dict(temp_dict,orient='index'))
    art_id+=1

columns = ['CHROM','POS','REF','ALT','Var_type','Seq_coverage','Read_depth','Ref_reads','Var_reads']
vcf_df.columns = columns
print("VCF file successfully parsed...")
# I tested many downstream steps using a subset of the vcf file, only records from Chromosome 1.
# see "Chr1_test_script.py"

# Check to make sure dataframe contains appropriate info
# >>>vcf_df.head()
#   CHROM      POS REF  ALT Var_type  Seq_coverage  Read_depth  Ref_reads Var_reads
# 0     1   931393   G  [T]    [snp]        4124.0        4124       4029      [95]
# 1     1   935222   C  [A]    [snp]        1134.0        1134        480     [652]
# 2     1  1277533   T  [C]    [snp]         786.0         786          0     [786]
# 3     1  1284490   G  [A]    [snp]         228.0         228          0     [228]
# 4     1  1571850   G  [A]    [snp]        4055.0        4055       3961      [94]

########1. (cont.) If there are multiple possibilities, annotate with the most deleterious possibility. ########
"""
For loci (rows) with multiple possible variants, I did not see a value in the VCF file that directly corresponds
to allele fitness, so I used CADD to determine the most deleterious mutations (see "CADD_instructions.txt")

"CADD can quantitatively prioritize functional, deleterious, and disease causal variants across a wide range
of functional categories, effect sizes and genetic architectures and can be used prioritize causal variation
in both research and clinical settings." (https://cadd.gs.washington.edu/)

CADD output is a tsv file with a new row and corresponding 'PHRED' score for every unique mutation.

Notably, by checking "Include Annotations", much of the resulting annotations are relevant.
I uploaded the annotated CADD output to GitHub but didn't use it for further analyses.
"""
# Expand series with multiple values per row to a new row, allowing direct comparison with CADD output
var_type_df = pd.DataFrame(vcf_df['Var_type'].tolist(), index=vcf_df.index).stack().reset_index(level=1, drop=False)
var_type_df = var_type_df.reset_index()
var_reads_df = pd.DataFrame(vcf_df['Var_reads'].tolist(), index=vcf_df.index).stack().reset_index(level=1, drop=False)
var_reads_df = var_reads_df.reset_index()
alt_allele_df = pd.DataFrame(vcf_df['ALT'].tolist(), index=vcf_df.index).stack().reset_index(level=1, drop=False)
alt_allele_df = alt_allele_df.reset_index()

# After completing the instructions in "CADD_instructions", read in the CADD file as a pandas DataFrame
# and sort by Chromosome and position (it should already be sorted in this way)

cadd = pd.read_csv(current_wd+"/cadd_output.tsv", header=1, sep='\t')
cadd = cadd.sort_values(['#CHROM','POS'], ascending=[True,True])
print("CADD file successfully parsed...")

# After expanding rows containing lists of variables, collapse back to a single DataFrame and combine with CADD PHRED scores
merged_df = pd.concat([var_type_df,var_reads_df[var_reads_df.columns[2]],alt_allele_df[alt_allele_df.columns[2]],cadd['PHRED']],axis=1)

# Next, we need to remove rows of loci for which there are multiple variants, keeping only the most deleterious
# I accomplished this using a for loop, and if/else to obtain a list of indices to remove from the dataframe...

temp_df = pd.DataFrame()
temp_df2 = pd.DataFrame()
for i in range(0, len(merged_df)-1):
    if merged_df['index'][i] == merged_df['index'][i+1]:                        # if multiple variants per loci
        temp_df = merged_df.loc[merged_df['index'] == merged_df['index'][i]]    # create subset df of these rows
        most_deleterious = temp_df['PHRED'].idxmax(axis=1)                      # determine which variant has highest 'PHRED' score
        temp_df = temp_df.loc[temp_df.index != most_deleterious]                # delete all other rows
        temp_df2 = temp_df2.append(temp_df)                                     # append to new df after each loop

# Drop indices (variants) from original dataframe that are not the most deleterious
filtered_df = merged_df.drop(temp_df2.index)
# Convert index back to original and rename columns
filtered_df = filtered_df.reset_index(drop=True)
column_names = ['index','level_1','Var_type','Var_reads','ALT','PHRED']
filtered_df.columns = column_names

print("Most deleterious possible variants detected...")
########### 2. Depth of sequence coverage at the site of variation. ###########

# Read depth info is contained within the original vcf_df dataframe ('Read_depth')
read_depth = vcf_df['Read_depth']

########### 3. Number of reads supporting the variant ###########

# Again, this data was parsed from the VCF file in the first step

# We can now combine all of this relevant information into a single DataFrame
filtered_df_new = pd.concat([vcf_df[['CHROM','POS','REF']],filtered_df[['ALT','Var_type','PHRED']],
                            vcf_df[['Seq_coverage','Read_depth','Ref_reads']],filtered_df['Var_reads']],axis=1)

#### 4.Percentage of reads supporting the variant versus those supporting reference reads. ####
# define objects
ref_r = filtered_df_new['Ref_reads']
var_r = filtered_df_new['Var_reads']
total_r = filtered_df_new['Read_depth']
ref_read_dict = {}
var_read_dict = {}
ref_read_df = pd.DataFrame()
var_read_df = pd.DataFrame()
art_id = 0
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
# iterate over the dataframe, returning a dictionary of percentages for either reference or variant reads,
# then converting these dictionaries to a dataframe
for i in range(0, len(filtered_df_new)):
    ref_read_percent = np.around((np.divide(ref_r[i],total_r[i])*100), decimals=3) # divide using numpy and return % with 3 decimals
    ref_read_dict = {art_id:[ref_read_percent]}
    var_read_percent = np.around((np.divide(var_r[i], total_r[i])*100), decimals=3)
    var_read_dict = {art_id: [var_read_percent]}
    ref_read_df = ref_read_df.append(pd.DataFrame.from_dict(ref_read_dict, orient='index'))
    var_read_df = var_read_df.append(pd.DataFrame.from_dict(var_read_dict, orient='index'))
    art_id += 1

ref_read_df.columns = ['% Reference Reads']
var_read_df.columns = ['% Variant Reads']

# combine % read data with previously generated DataFrame
df = pd.concat([filtered_df_new,ref_read_df,var_read_df],axis=1)

print("Values for questions #1-4 stored as pandas DataFrame")
print("df.head()")
df.head()
# Can save DataFrame to directory so downstream steps can be optimized w/o need for repeating upstream code #
df.to_csv(current_wd+"/VCF_DF.txt", sep='\t')

#### 5-6. Extracting data from EXAC Database ###
print("Extracting EXAC data from exac.hms.harvard.edu/rest/variant/variant/")

# if running code in chunks, we can import the previously generated dataframe now
# df = pd.read_csv(current_wd+"/VCF_DF.txt", sep='\t')

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
This is the most time-consuming step of the whole code, taking roughly 15 min to complete on my local machine
(w/o multi-threading)

Once, while running the code on my Mac (OS X), I encountered the following error:
...
BlockingIOError: [Errno 35] Resource temporarily unavailable

... Because of this, I created the "exac_script.py" script to circumvent the need to re-run all previous steps.
"""
# Parse the dictionary of dataframes containing variant-specific exac data
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
        temp_dict = {art_id: ["NA","NA","NA"]}  # EXAC data not available for every record
    exac_df_parse = exac_df_parse.append(pd.DataFrame.from_dict(temp_dict, orient='index'))
    art_id += 1

exac_df_parse.columns = ['Allele frequency','rsid','VEP annotation']
print("Assembling final table...")
final_df = pd.concat([df,exac_df_parse],axis=1)
final_df.to_csv(current_wd+"/final_table.txt", sep='\t')
print("Successfully completed!")
