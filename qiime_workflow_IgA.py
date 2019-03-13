import sys, os, csv, math
import numpy as np
from scipy.optimize import minimize

def deconvolve_samples(sample_abundances_filename, FACS_percentages_filename, metadata_filename, taxonomy, output_abdundance_file, output_datafile, output_table, barplot_filename):
	# Read in two input files:
	#  1) Abundance file with OTU as columns, samples as rows
	#  2) IgA binding percentages for each sample
	# Output is a biom file with deconvolved OTUs

	FACS_percentages = open(FACS_percentages_filename, 'rU')
	FACS_reader = csv.reader(FACS_percentages, delimiter='\t')
	FACS_header = next(FACS_reader)

	#"B#_D#_G#" = [percent IgA+, percent IgA-]
	# Read in IgA binding percentages and store
	FACS_data = {}
	for line in FACS_reader:
		sample = line[0]
		sample_id = sample.split(" ")[0]
		sample_type = ""
		if "APC" in sample:
			sample_type = "US"
		elif "IgA-" in sample:
			sample_type = "NEG"
		elif "IgA+" in sample:
			sample_type = "POS"
		else:
			print("FACS file data error. Cannot determine type of sample for %s" % (sample))
			exit()
		sample_name = "B%i_D%i_G%s" % (int(sample_id.split("-")[0]), int(sample_id.split("-")[1]), sample_type)
		FACS_data[sample_name] = [float(line[1])/100., float(line[2])/100.]

	FACS_percentages.close()


	#Read in abundances from file
	df = open(sample_abundances_filename, 'rU')
	reader = csv.reader(df, delimiter='\t')
	sa_header = next(reader)
	bacteria = sa_header[1:]
	num_bacteria = len(bacteria)

	#Deconvolved output
	out_abdundance = "%s.txt" % (output_abdundance_file)
	out_a = open(out_abdundance, "w", newline="")
	writer_a = csv.writer(out_a, delimiter='\t')
	writer_a.writerow(sa_header)

	# Deconvolution
	for line in reader:
		parts = line[0].split("_")
		neg_name = line[0]
		neg_sample = "%s_%s_%s" % (parts[0], parts[1], parts[2])
		neg_counts = []
		neg_sum = np.sum(int(i) for i in line[1:])
		for i in range(num_bacteria):
			neg_counts.append(int(line[i+1]))

		line = next(reader)
		parts = line[0].split("_")
		pos_name = line[0]
		pos_sample = "%s_%s_%s" % (parts[0], parts[1], parts[2])
		pos_counts = []
		pos_sum = np.sum(int(i) for i in line[1:])
		for i in range(num_bacteria):
			pos_counts.append(int(line[i+1]))

		line = next(reader)
		parts = line[0].split("_")
		us_name = line[0]
		us_sample = "%s_%s_%s" % (parts[0], parts[1], parts[2])
		us_counts = []
		us_sum = np.sum(int(i) for i in line[1:])
		for i in range(num_bacteria):
			us_counts.append(int(line[i+1]))

		new_counts = {}
		for i in range(num_bacteria):
				# deconvolved_IgA_positive = (%IgA+ in IgA+ sample) * (num_reads in IgA+ sample) + (%IgA+ in IgA- sample) * (num_reads in IgA- sample)
			deconvolved_IgA_positive = (FACS_data[pos_sample][0] * pos_counts[i]) + (FACS_data[neg_sample][0] * neg_counts[i])
			deconvolved_IgA_negative = (FACS_data[pos_sample][1] * pos_counts[i]) + (FACS_data[neg_sample][1] * neg_counts[i])
			new_counts[bacteria[i]] = [int(deconvolved_IgA_positive), int(deconvolved_IgA_negative)]


		out_line = [neg_name] + [new_counts[b][1] for b in bacteria]
		writer_a.writerow(out_line)
		out_line = [pos_name] + [new_counts[b][0] for b in bacteria]
		writer_a.writerow(out_line)
		out_line = [us_name] + [us_counts[i] for i in range(num_bacteria)]
		writer_a.writerow(out_line)


	df.close()
	out.close()
	out_a.close()


	#Transpose abdundance file to convert to biom file for qiime
	output_file = "%s.txt" % (output_datafile)
	df = open(out_abdundance, 'rU')
	a = zip(*csv.reader(df, delimiter='\t'))
	header = next(a)
	output_header = ['OTU'] + [i for i in header[1:]] + ['taxonomy']
	output = open(output_file, 'w', newline="")
	writer = csv.writer(output, delimiter='\t')
	writer.writerow(output_header)
	i = 1
	for line in a:
		writer.writerow(['%i' %(i)] + [x for x in line[1:]] + [line[0]])
		i += 1
	
	df.close()
	output.close()

	#Make taxonomy file
	tax_file = '%s.txt' % (taxonomy)
	df = open(output_file, 'rU')
	reader = csv.reader(df, delimiter='\t')
	next(reader)
	output = open(tax_file, 'w', newline="")
	writer = csv.writer(output, delimiter='\t')
	for line in reader:
		writer.writerow([line[0], line[-1]])
	df.close()
	output.close()

	#Run qiime commands
	print("Convert %s to biom file:" % (output_file))
	bashCommand = 'biom convert -i %s.txt -o %s.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy' % (output_datafile, output_datafile)
	os.system(bashCommand)
	
	print("Create qiime artifact:")
	bashCommand = "qiime tools import --input-path %s.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path %s.qza" % (output_datafile, output_table)
	os.system(bashCommand)

	print("Create taxonomy file")
	bashCommand = "qiime tools import --input-path %s.txt --output-path %s.qza --type FeatureData[Taxonomy] --input-format HeaderlessTSVTaxonomyFormat" % (taxonomy, taxonomy)
	os.system(bashCommand)

	print("Create barplot")
	bashCommand = "qiime taxa barplot --i-table %s.qza --i-taxonomy %s.qza --m-metadata-file %s --o-visualization %s.qzv" % (output_table, taxonomy, metadata_filename, barplot_filename)
	os.system(bashCommand)


def main():
		#Input files
		#Ex: sample_abundances_filename = "Merged_abdundance_file.txt"

		#Abundances with OTU as columns, samples as rows
	sample_abundances_filename = ""
		#IgA binding file where each row has format: Sample name; %IgA Positive; %IgA Negative
	FACS_percentages_filename = ""
		#Sample metadata used in qiime commands
	metadata_filename = ""

		#Output files without extensions
		# Ex: taxonomy = "Deconvolution/taxonomy"

		#Taxonomy file
	taxonomy = ""
		#Abdundance file with OTUs as columns, samples as rows
	output_abdundance_file = ""
		#Biom file with same data as abundance file
	output_datafile = ""
		#Qiime table
	output_table = ""
		#Qiime barplot
	barplot_filename = ""

		#program call
	deconvolve_samples(sample_abundances_filename, FACS_percentages_filename, metadata_filename, taxonomy, output_abdundance_file, output_datafile, output_table, barplot_filename)

if __name__ == "__main__":
	main()