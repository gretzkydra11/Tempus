import vcf, pandas as pd, numpy as np, sys, subprocess, ast, os

# Define local path to read/write files #
PATH = "/Users/denisavey/Tempus_challenge/"

### Read in VCF file using PyVCF ###
vcf_rd = vcf.Reader(open('/Users/denisavey/Tempus_challenge/Challenge_data.vcf', 'r'))

### Loop through each row of the VCF file and store relevant information:

vcf_df =pd.DataFrame()
vcf_df2 =pd.DataFrame()
vcf_df3 =pd.DataFrame()
vcf_df4 =pd.DataFrame()
vcf_df5 =pd.DataFrame()
vcf_df6 =pd.DataFrame()
vcf_df7 =pd.DataFrame()
vcf_df8 =pd.DataFrame()
vcf_df9 =pd.DataFrame()
art_id = 0

for record in vcf_rd:
    temp_type_dict = {art_id:[record.INFO['TYPE']]}
    temp_coverage_dict = {art_id:[record.INFO['DPB']]}
    temp_read_depth_dict = {art_id:[record.INFO['DP']]}
    temp_var_reads_dict = {art_id:[record.INFO['AO']]}
    temp_ref_reads_dict = {art_id:[record.INFO['RO']]}
    temp_alt_allele_dict = {art_id:[record.ALT]}
    chrom_dict = {art_id:[record.CHROM]}
    pos_dict = {art_id:[record.POS]}
    ref_dict = {art_id:[record.REF]}
    vcf_df = vcf_df.append(pd.DataFrame.from_dict(temp_type_dict,orient='index'))
    vcf_df2 = vcf_df2.append(pd.DataFrame.from_dict(temp_coverage_dict,orient='index'))
    vcf_df3 = vcf_df3.append(pd.DataFrame.from_dict(temp_read_depth_dict,orient='index'))
    vcf_df4 = vcf_df4.append(pd.DataFrame.from_dict(temp_var_reads_dict,orient='index'))
    vcf_df5 = vcf_df5.append(pd.DataFrame.from_dict(temp_ref_reads_dict,orient='index'))
    vcf_df6 = vcf_df6.append(pd.DataFrame.from_dict(temp_alt_allele_dict,orient='index'))
    vcf_df7 = vcf_df7.append(pd.DataFrame.from_dict(chrom_dict,orient='index'))
    vcf_df8 = vcf_df8.append(pd.DataFrame.from_dict(pos_dict,orient='index'))
    vcf_df9 = vcf_df9.append(pd.DataFrame.from_dict(ref_dict,orient='index'))

    art_id+=1

# Combine relevant data into pandas dataframe and name columns #
vcf_df = pd.concat([vcf_df,vcf_df2,vcf_df3,vcf_df4,vcf_df5,vcf_df6],axis=1)
columns = ['Var_type','Seq_coverage','Read_depth','Var_reads','Ref_reads','ALT']
vcf_df.columns = columns

#################################################################
######## 1. (cont.) If there are multiple possibilities, ########
######## annotate with the most deleterious possibility. ########
#################################################################

# For loci (rows) with multiple possible variants, I did not see a value in the
# VCF file that directly corresponds to allele fitness, so I used CADD to
# determine the most deleterious mutations (https://cadd.gs.washington.edu/)
# From their website: "CADD can quantitatively prioritize functional, deleterious, and disease causal variants across a wide range of functional categories, effect sizes and genetic architectures and can be used prioritize causal variation in both research and clinical settings."
# CADD output is a tsv file with a new row and corresponding 'PHRED' score
# for every unique mutation.

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
CADD = pd.read_csv("/Users/denisavey/Tempus_challenge/CADD_output.tsv", header = 1, sep ='\t')
CADD = CADD.sort_values(['#CHROM','POS'], ascending =[True,True])

### Next, we need to remove rows of loci for which there are multiple variants,
### keeping only the most deleterious one
### I accomplished this by a combination of a for loop, and if/and filters
### to obtain a list of indices to drop from the dataframe...

#sys.stdout = open('/Users/denisavey/Tempus_challenge/drop.txt', 'w+')
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
merged_df4.to_csv("/Users/denisavey/Tempus_challenge/merged_df4.txt", sep = '\t')

#### Step 4; just need to do some numpy calculation and iteration over dataframe values
#var_read_% =

#    for i in range(0, len(merged_df)-1):
#### Steps 5-6; can't figure out subprocess...

#import vcf, pandas as pd, numpy as np, sys, subprocess, ast

merged_df4 = pd.read_csv("/Users/denisavey/Tempus_challenge/merged_df4.txt", sep = '\t')

exac_data =pd.DataFrame()
s= "-"
d={}
temp =[]
EXAC = "exac.hms.harvard.edu/rest/variant/variant/"
art_id = 0
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
    os.remove(file)
    #art_id+=1

exac_data.to_csv("/Users/denisavey/Tempus_challenge/exac_data.txt", sep = '\t')

#subprocess.run('wget exac.hms.harvard.edu/rest/variant/variant/14-21853913-T-C')
