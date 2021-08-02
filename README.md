# General

A repository for some code related to genetics/ genomics under development.

## Introgression viewer

A program for quickly viewing a high level overview of introgressions or contrasting genotypes between two or more individuals. Designed for high density SNP data as an input. The program was designed for viewing introgressions in barley genotypes. It currently only works with fixed (ie non-heterozygous) individuals and will ignore region of heterozygozity, though there are plans to include het calls.

The genotype of a particular region is determoined from a sliding window and parameters set by the user. The larger the sliding window, the more accurate the prediction of the genotype for that particular region, but there is a trade off between accuracy and resolution. A histogram of identity displayed by the program should bbe used to fine tune the parameters. The user can also set a minimum introgression size to avoid false positive recombination events.

### Dependencies:

-  pandas
-  matplotlib

Introgression viewer requires python >=3.7


### Usage

python Introgression_viewer.py [-d] [-c] [-chrom] [-cutoff] [-step_size] [-window_size] [-colour] [-features] 

**Optional arguments:**

  -d               Input dataset
  -c               The comparisons you would like to make
  -chrom           Optional: The chromosomes you want to visualise
  -cutoff          Optional: Percentage similarity cutoff for determining alleles (default: 90)
  -step_size       Optional: Step size (default = 10)
  -window_size     Optional: Window size (default = 100)
  -colour          Optional: Colours for comparison (default = "purple,green")
  -features        Optional: Genes or other features to add
  
Default usage would look like this:

python Introgression_viewer.py -d my_genotypes.tab -c "Int_52-Int_17, Int_52-Barke"

Detailed explanation of arguments:

**-d Input dataset** This is a tab delimited text file with genotype information for each of your individuals of interest. The columns must be in the format Marker Name, chrom, position, line 1 allele, line 2 allele,.... line n allele. e.g:

    Marker	Chromosome	Barke position	Barke	Int_52	Int_33	Int_42	Int_17	Int_19	Int_56
    JHI-Hv50k-2016-7	chr1H	112201	T	T	T	T	T	T	T
    JHI-Hv50k-2016-24	chr1H	110264	T	T	T	T	T	T	T
    JHI-Hv50k-2016-64	chr1H	106615	A	A	A	A	A	A	A
    JHI-Hv50k-2016-66	chr1H	106465	T	T	T	T	T	T	T
    JHI-Hv50k-2016-72	chr1H	105691	G	G	G	G	G	G	G
    JHI-Hv50k-2016-73	chr1H	105623	NA	NA	NA	NA	NA	NA	NA
    JHI-Hv50k-2016-88	chr1H	104313	A	A	A	A	A	A	A
    JHI-Hv50k-2016-97	chr1H	103806	A	A	A	A	A	A	A
  
 Genotypes should be in the format A|C|G|T. Heterozygous genotypes (e.g A/T) and NAs will be removed automatically.
 
 **-c comparisons** The genotypes you want to compare. One or more comparisons in this format: genotype2-genotype1, genotype3-genotype1. Warning: Different comparisons may have different similarities, so you may want to run one at a time with paramters set for each one. Also make sure genotype names match genotype names in your -d input file, otherwise nothing will work!')
 
 **-chrom chromosomes** Optional: The chromosomes you want to visualise. e.g chr1H,chr2H'. Default will display all. 
 
 **-cutoff** Optional: Percentage similarity cutoff for determining alleles', default = 90. Use the histogram to set this value correctly
 
 **-step_size** Step size (int). Small value (1) can give more accurate introgression intervals, but will be fuzzy. Large value for sharp intervals, but loss of accuracy', default = 10
 
 **-window_size** Window size. Larger windows will mean less false positives/negatives, but could lead to more inaccurate boundry positions
 
 **-colour** Colours for comparison in format <colour1,colour2>. Default: purple and green.
 
 **-features** Genes or other features to add onto chromosomes. File must be in bed4 format
 
 
### Outputs

Outputs will appear in the same directory as the script

**Genotype2-Genotype1.bed** Introgressions (or regions of similarity/difference) in bed4 format. This is used as the input for the chromosome plot. Check here for any oddities.

**Genotype2-Genotype1_chromosomes.png**  The visualisation of chromosomes

Examples:

![Figure 1](https://github.com/maxecoulter/General/blob/main/figures/Int_19-Barke_chromosomes.png)

![Figure 2](https://github.com/maxecoulter/General/blob/main/figures/Int_52-Barke_chromosomes.png)

In these examples, those regions that have a different genotype between the two lines being compared are in yellow, whilst areas of the same genotype are in blue. If you compare a line from a population with its parent, the areas of difference will correspond to introgressions (as in these examples).


