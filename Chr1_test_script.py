"""

Author: Denis Avey
Date:

This master script will generate a tab-delimited text file containing various info
parsed from the Challenge_data.vcf file, the EXAC database, and the CADD database.

Note 1: Before proceeding with the downstream steps, first complete steps outlined in CADD_instructions

Note 2: All below code was developed using Python 3.7.3

"""
import os, vcf, pandas as pd, numpy as np, sys, subprocess, ast

# Define local path to read/write files
current_wd = os.getcwd()
print("working directory:",current_wd)

# Read in VCF file using PyVCF
#vcf_rd = vcf.Reader(open(current_wd+"/Challenge_data.vcf", 'r'))

# I tested many downstream steps using a subset of the vcf file, only records from Chromosome 1.
# I originally tried to subset the vcf_rd object, but to circumvent the need for the pysam module
# (as a required dependency for the "vcf.fetch" method), I just manually subsetted the vcf file in Excel
vcf_chr1 = vcf.Reader(open(current_wd+"/Chr1_data.vcf", 'r'))

### Loop through each row of the VCF file and store relevant information in pandas DataFrame:
temp_dict={}
vcf_chr1_df =pd.DataFrame()
art_id = 0
for record in vcf_chr1:
    temp_dict = {art_id:[record.CHROM,record.POS,record.REF,record.ALT,record.INFO['TYPE'],record.INFO['DPB'],
    record.INFO['DP'],record.INFO['RO'],record.INFO['AO']]}
    vcf_chr1_df = vcf_chr1_df.append(pd.DataFrame.from_dict(temp_dict,orient='index'))
    art_id+=1

columns = ['CHROM','POS','REF','ALT','Var_type','Seq_coverage','Read_depth','Ref_reads','Var_reads']
vcf_chr1_df.columns = columns

# Check to make sure dataframe contains appropriate info
# >>>vcf_chr1_df.head()
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

Notably, by checking "Include Annotations
"""
# Expand series with multiple values per row to a new row, allowing direct
# comparison with CADD output
var_type_df = pd.DataFrame(vcf_df['Var_type'].tolist(), index=vcf_df.index).stack().reset_index(level=1, drop=False)#.reset_index(name='Variant_type')[['Variant_type']]
var_type_df = var_type_df.reset_index()
var_reads_df = pd.DataFrame(vcf_df['Var_reads'].tolist(), index=vcf_df.index).stack().reset_index(level=1, drop=False)#.reset_index(name='Variant_type')[['Variant_type']]
var_reads_df = var_reads_df.reset_index()
alt_allele_df = pd.DataFrame(vcf_df['ALT'].tolist(), index=vcf_df.index).stack().reset_index(level=1, drop=False)#.reset_index(name='Variant_type')[['Variant_type']]
alt_allele_df = alt_allele_df.reset_index()

# After expanding rows containing lists of variables, collapse back to a single DataFrame
merged_df = pd.concat([var_type_df,var_reads_df[var_reads_df.columns[2]],
alt_allele_df[alt_allele_df.columns[2]]],axis=1)

# After submitting the VCF file, and downloading the CADD file,
# Change the filename to "CADD_output.tsv" and move to the working directory
# Read in the CADD file as a pandas DataFrame and sort by Chromosome and position

#convert CADD to cadd
CADD = pd.read_csv("/Users/denisavey/Tempus_challenge/CADD_output.tsv", header = 1, sep ='\t')
CADD = CADD.sort_values(['#CHROM','POS'], ascending =[True,True])

### Next, we need to remove rows of loci for which there are multiple variants, keeping only the most deleterious
### I accomplished this using a for loop, and if/else to obtain a list of indices to remove from the dataframe...

with open('/Users/denisavey/Tempus_challenge/CADD_output.tsv', 'w+') as cadd:


with open('/Users/denisavey/Tempus_challenge/drop.txt', 'w+') as drop:
    for i in range(3, len(merged_df)):
        if merged_df['index'][i] == merged_df['index'][i-1]:
            if CADD['PHRED'][i-1] <= CADD['PHRED'][i]:
                print(merged_df.index[i-1], file=drop)
            else:
                print(merged_df.index[i], file=drop)
        if merged_df['index'][i] == merged_df['index'][i-2]:
            if CADD['PHRED'][i-2] <= CADD['PHRED'][i]:
                print(merged_df.index[i-2], file=drop)
            else:
                print(merged_df.index[i], file=drop)
        if merged_df['index'][i] == merged_df['index'][i-3]:
            if CADD['PHRED'][i-3] <= CADD['PHRED'][i]:
                print(merged_df.index[i-3], file=drop)
            else:
                print(merged_df.index[i], file=drop)

drop_index = pd.read_csv("/Users/denisavey/Tempus_challenge/drop.txt", sep = '\t',header=None)
merged_df2 = merged_df.drop(merged_df.index[drop_index])
#merged_df2.to_csv("/Users/denisavey/Tempus_challenge/merged_df2.txt", sep = '\t')

### Can you think of a "tidier" way? Basically I just want to drop all the rows
### for which the 'level_1' column value is equivalent, and the corresponding
### PHRED column in the CADD file is not the max among them, to return only
### the most deleterious mutant for each locus.

###############################################################
### 2. Depth of sequence coverage at the site of variation. ###
###############################################################

# Read depth info is contained within the parsed VCF file (record.INFO['DPB']),
# and is a column of the vcf_df dataframe ('Read_depth')

read_depth = vcf_df['Read_depth']

###############################################################
########## 3. Number of reads supporting the variant ##########
###############################################################

### Again, this data was parsed from the VCF file in the first step
### We can now combine all of this relevant information into a single DataFrame
merged_df3 = merged_df2.reset_index(drop=True)
merged_df4 = pd.concat([vcf_df7,vcf_df8,vcf_df9,merged_df3,read_depth,vcf_df['Ref_reads']],axis=1)
new_col_names = ['CHROM','POS','REF','Index','Variant#','Var_type','Var_reads','ALT','Read_depth','Ref_reads']
merged_df4.columns = new_col_names

#### Step 4; just need to do some numpy calculation and iteration over dataframe values
#var_read_% =

# performed this step when I wanted to start from here and perform
# merged_df4.to_csv("/Users/denisavey/Tempus_challenge/merged_df4.txt", sep = '\t')



#    for i in range(0, len(merged_df)-1):
#### Steps 5-6; can't figure out subprocess...

#import vcf, pandas as pd, numpy as np, sys, subprocess, ast

merged_df4 = pd.read_csv("/Users/denisavey/Tempus_challenge/merged_df4.txt", sep = '\t')

exac_df =pd.DataFrame()
exac_dict = {}
s= "-"
d={}
EXAC = "exac.hms.harvard.edu/rest/variant/variant/"
PATH = "/Users/denisavey/Tempus_challenge/"

for i in range(0,len(merged_df4)-1):
    variant: str = s.join([str(merged_df4['CHROM'][i]),str(merged_df4['POS'][i]),str(merged_df4['REF'][i]),str(merged_df4['ALT'][i])])
    site = EXAC+variant
    subprocess.call(('wget' + ' ' + site), shell=True)
    file_loc = PATH+variant
    with open(file_loc, "r") as file:
        d = ast.literal_eval(file.read())
 #  vcf_df = vcf_df.append(pd.DataFrame.from_dict(d,orient='index'))
 #  append these info to the original vcf_df dataframe
 #  exac_dict[i] = exac_df
    os.remove(file)

# print head()
#merge with vcf_df
# exac_df.to_csv("/Users/denisavey/Tempus_challenge/exac_data.txt", sep = '\t')

#subprocess.run('wget exac.hms.harvard.edu/rest/variant/variant/14-21853913-T-C')
