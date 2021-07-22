"""
For generating ideograms showing introgressions, or similarities/differences between lines
Input: Marker Name, chrom, position, line 1 allele, line 2 allele,.... line n allele
"""
import numpy as np
from operator import neg
import pandas as pd
import statistics as st
import argparse
#from sklearn.neighbors import KernelDensity
#from sklearn.cluster import MeanShift, estimate_bandwidth
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection






def get_matches(df, comparison):
	"""Adds a column to dataframe with 1 for allele match and 2 for no match at each marker for the comparison specified"""
	line1, line2 = comparison.split("-")
	matches = pd.DataFrame({comparison:[1 if df.at[i,line1] == df.at[i,line2] else 0 for i in df.index]})
	df2 = pd.concat([df,matches],axis=1)
	return df2

def find_introgressions(df2,columns,comparison,window_size,step_size,cutoff,chromosomes):
	"""Predicts approximate position for introgressions for comparison"""
	chromosomes.sort()
	for c in range(0,len(chromosomes)):
		df3 = df2[df2[columns[1]] == chromosomes[c]]
		df3 = df3.reset_index(drop=True)
		matches = list(df3[comparison])
		counts = []
		counts_dict = {}
		for i in range(0,len(matches) - window_size,step_size):
			print(str(i))
			count = matches[i:i + window_size].count(1.0)
			counts.append((count/window_size) * 100)
			counts_dict[i + window_size/2] = (count/window_size) * 100
		#Cluster counts
		plt.hist(counts,bins=50)
		plt.savefig(f"{chromosomes[c]}_histogram.png")
		plt.clf()
		predictions = []
		allele = "NA"
		start_allele = "NA"
		for i in range(0,len(matches)):
			if i in counts_dict.keys():
				if counts_dict[i] >= cutoff:
					predictions.append("A")
					allele = "A"
					if start_allele == "NA":
						start_allele = "A"
				else:
					predictions.append("A")
					allele = "B"
					if start_allele == "NA":
						start_allele = "B"
			else:
				predictions.append(allele)
		predictions = [start_allele if a == "NA" else a for a in predictions]
		predictions_df = pd.DataFrame({f"predictions_{comparison}":predictions})
		df3 = pd.concat([df3,predictions_df],axis=1)
		if not c:
			df4 = df3
		else:
			df4 = pd.concat([df4,df3],axis=0)
	df4 = df4.reset_index(drop=True)
	return df4

def create_bed_file(df3,comparison):
	with open(f"{comparison}.bed","w") as out:
		columns = list(df3.columns)
		previous_chrom = "NA"
		previous_allele = "NA"
		previous_position = "NA"
		for i in df3.index:
			chromosome = df3.at[i,columns[1]]
			position = int(df3.at[i,columns[2]])
			allele = df3.at[i,f"predictions_{comparison}"]
			if chromosome != previous_chrom:
				if previous_chrom != "NA":
					out.write(f"{previous_position}\ttest\t{previous_allele}\n")
				out.write(f"{chromosome}\t0\t")
				previous_chrom = chromosome
				previous_allele = allele
				previous_position = position
			elif allele != previous_allele and chromosome == previous_chrom:
				out.write(f"{position - 1}\ttest\t{previous_allele}\n")
				out.write(f"{chromosome}\t{position}\t")
				previous_allele = allele
				previous_position = position
			else:
				previous_position = position
		out.write(f"{position}\ttest\t{previous_allele}\n")
	return



def drop_hets(df,columns):
	"""Remove any rows with heterozygous calls."""
	no_hets = []
	for i in df.index:
		hets = False
		for c in columns[3:]:
			if len(df.at[i,c].split("/")) == 2:
				hets = True
		no_hets.append(hets)
	filter = pd.Series(no_hets)
	df = df[~filter]
	df = df.reset_index(drop=True)
	return df

def chromosome_collections(df, y_positions, height,  **kwargs):
    """
	Taken from https://www.biostars.org/p/147364/
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']

def create_chromosome_plot(df3,comparison, chromosome_list, colour, features,**kwargs):
	"""Below adapted from https://www.biostars.org/p/147364/"""
	colour1, colour2 = [c.replace("%","#") if c.startswith("%") else c for c in colour.split(",")]

	# Height of each ideogram
	chrom_height = 1
	# Spacing between consecutive ideograms
	chrom_spacing = 1
	# Height of the gene track. Should be smaller than `chrom_spacing` in order to
	# fit correctly
	gene_height = 0.2
	# Padding between the top of a gene track and its corresponding ideogram
	gene_padding = 0.1
	# Width, height (in inches)
	length = len(chromosome_list) * 0.5 #More chromosomes, longer figure!
	if length < 2:
		length += 2
	figsize = (6, length)
	# Decide which chromosomes to use
	#chromosome_list = list(set(df3[columns[1]]))
	chromosome_list.sort()
	# Keep track of the y positions for ideograms and genes for each chromosome,
	# and the center of each ideogram (which is where we'll put the ytick labels)
	ybase = 0
	chrom_ybase = {}
	gene_ybase = {}
	chrom_centers = {}
	# Iterate in reverse so that items in the beginning of `chromosome_list` will
	# appear at the top of the plot
	for chrom in chromosome_list[::-1]:
		chrom_ybase[chrom] = ybase
		chrom_centers[chrom] = ybase + chrom_height / 2.
		gene_ybase[chrom] = ybase - gene_height - gene_padding
		ybase += chrom_height + chrom_spacing
	# Read in ideogram.txt, downloaded from UCSC Table Browser
	ideo = pd.read_table(f'{comparison}.bed',names=['chrom', 'start', 'end', 'name', 'gieStain'])
	# Filter out chromosomes not in our list
	ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]
	# Add a new column for width
	ideo['width'] = ideo.end - ideo.start
	# Colors for different chromosome stains
	color_lookup = {'A': colour2,'B': colour1}
	# Add a new column for colors
	ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])
	fig = plt.figure(figsize=figsize)
	ax = fig.add_subplot(111)
	# Now all we have to do is call our function for the ideogram data...
	print("adding ideograms...")
	for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors="black", linewidth=0.5):
		ax.add_collection(collection)
	# ...and the gene data
	if features:#If the user includes a feature file(bed)
		# Same thing for genes
		genes = pd.read_table('BaRT2v18.bed',names=['chrom', 'start', 'end', 'name'],usecols=range(4))
		genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
		genes['width'] = genes.end - genes.start
		genes['colors'] = '#2243a8'
		print("adding features...")
		for collection in chromosome_collections(genes, gene_ybase, gene_height, alpha=0.5, linewidths=0):
			ax.add_collection(collection)
	# Axes tweaking
	ax.set_yticks([chrom_centers[i] for i in chromosome_list])
	ax.set_yticklabels(chromosome_list)
	ax.axis('tight')
	#plt.show()
	plt.savefig(f"{comparison}_chromosomes.png")
	plt.clf()


	





def main():
	parser = argparse.ArgumentParser(description='Create plot to visualise introgressions on near isogenic lines')
	parser.add_argument('-d', dest = 'dataset', help = 'Dataset for analysis. Must be in format Marker Name, chrom, position, line 1 allele, line 2 allele,.... line n allele', type = str )
	parser.add_argument('-c', dest = 'comparisons', type = str, help = 'The genotypes you want to compare. One or more comparisons in this format: genotype2-genotype1, genotype3-genotype1. Warning: Different comparisons may have different similarities, so you may want to run one at a time. Also make sure genotype names match genotype names in columns, otherwise nothing will work!')
	parser.add_argument('-chrom', dest = 'chromosomes', type = str, help = 'Optional: The chromosomes you want to visualise. Default will display all. e.g chr1H,chr2H', default = '')
	parser.add_argument('-cutoff', dest = 'cutoff', type = int, help = 'Optional: Percentage similarity cutoff for determining alleles', default = 90)
	parser.add_argument('-step_size', dest = 'step_size', type = int, help = 'Step size. Small value (1) can give more accurate introgression intervals, but will be fuzzy. Large value for sharp intervals, but loss of accuracy', default = 10)
	parser.add_argument('-window_size', dest = 'window_size', type = int, help = 'WIndow size. Larger windows will lead to less false positives/negatives, but more inaccurate boundry positions', default = 100)
	parser.add_argument('-colour', dest = 'colour', type = str, help = 'Colours for comparison in format <colour1,colour2>. Default: purple and green',default = "purple,green")
	parser.add_argument('-features', dest = 'features', type = str, help = 'Genes or other features to add. File must be in bed4 format',default = "")
	parser.add_argument('--dash', dest = 'dash', type = str, help = 'Use plotly dash visualisation',action="store_true")
	args = parser.parse_args()
	#dataset = "introgression_viewer_test.txt"
	#comparisons = "Int_52-Int_17"
	
	
	df = pd.read_csv(args.dataset,delimiter="\t",engine="c")
	columns = list(df.columns)
	df = df.dropna(axis=0)#Remove rows with NANs
	df = df.reset_index(drop=True)
	df = df.drop_duplicates(subset=[columns[0]]) #Remove duplicate markers
	df_sorted = df.sort_values(by=[columns[1], columns[2]]) #Sort by chromosome,position
	df_sorted = df_sorted.reset_index(drop=True)
	df2 = drop_hets(df_sorted,columns)
	comparison_list = args.comparisons.split(",")
	if args.chromosomes:
		chromosomes = list(args.chromosomes.split(","))
	else:
		chromosomes = list(set(df2[columns[1]]))
	for comparison in comparison_list:
		df2 = get_matches(df2, comparison)
		df3 = find_introgressions(df2,columns,comparison,args.window_size,args.step_size,args.cutoff,chromosomes)
		#Output introgressions in bed4 format
		create_bed_file(df3,comparison)
		create_chromosome_plot(df3, comparison, chromosomes, args.colour,args.features)

if __name__ == "__main__":
	main()


