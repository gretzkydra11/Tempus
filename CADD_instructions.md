Description:

For loci (rows) in the VCF file containing multiple possible variants, I did not see a value that directly corresponds to allele fitness (Question #1, pt 2), so I used "CADD" to determine the most deleterious mutations.

"CADD can quantitatively prioritize functional, deleterious, and disease causal variants across a wide range of functional categories, effect sizes and genetic architectures and can be used prioritize causal variation in both research and clinical settings."

For more info, see the following manuscripts:

Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M.
CADD: predicting the deleteriousness of variants throughout the human genome.
Nucleic Acids Res. 2018 Oct 29. doi: 10.1093/nar/gky1016.
PubMed PMID: 30371827.

Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J.
A general framework for estimating the relative pathogenicity of human genetic variants.
Nat Genet. 2014 Feb 2. doi: 10.1038/ng.2892.
PubMed PMID: 24487276.

Instructions:

1) The CADD online tool requires input VCF files to be smaller than 2 MB, so first compress the VCF file. This can be accomplished in command line with the following command:

gzip &lt; ./Challenge_data.vcf &gt; ./Challenge_data.vcf.gz

2) Navigate to https://cadd.gs.washington.edu/score, and upload the compressed VCF file.

3) Ensure that 'GRCh37-v1.4' is selected in the drop-down menu. Optionally, click on 'Include Annotations', which will generate a larger, more detailed output file. For my downstream code, this step is not necessary.

4) Click on 'Upload Variants'. The online tool will take a few minutes to run (longer if annotations included).

5) Download the resultant CADD output file, unzip, rename to "cadd_output.tsv" and move to working directory.

Notes:

The CADD output is a tsv file with a corresponding 'PHRED' score for every unique mutation.

A higher PHRED score is indicative that a variant is expected to be more deleterious.

Because records containing multiple possible variants are expanded to separate rows in the CADD output file, it is helpful to re-structure some of the data from the VCF file in a similar format (see master_script.py).
