import logging
import time

logging.basicConfig(filename = 'kallisto_import.log', level=logging.INFO,
					format = '%(levelname)s:%(asctime)s:%(message)s')

def kallisto_import(path_to_file, output_file_ = None):
	'''Import Kallisto count from all run_info.json and abundance.tsv files into a matrix file for correlation analysis
	
	Args:
		path_to_file - Folder Directory with Kallisto Folders 
		output_file - output name
	#Code below is adopted from Qiao Wen
	'''

	logging.info('Importing Counts and Mapping Stats from Kallisto Counts')

	if output_file_:
		path_to_result_directory = os.path.join(path_to_file, output_file_, 'results')
	else:
		path_to_result_directory = os.path.join(path_to_file, 'results')
	os.makedirs(path_to_result_directory)

	map_stats = open(path_to_result_directory + '/mapping_stats.tsv', 'a+') #mapping statistics output
	if output_file_:
		output_file = open(path_to_result_directory + '/' + output_file_ +'.txt', 'a+') #full matrix file
	else: 
		output_file = open(path_to_result_directory + '/count_matrix.txt', 'a+') #full matrix file
	dict_log = {} #stores average read length from logfile.txt

	#Creates expression matrix and mapping statistics log

	dicto = {} #stores gene and corresponding tpm values from the experiment 
	output_header = 'gene\t'
	output_content = ''
	map_stats.write('Run ID\tProcessed reads\tMapped reads\t% reads mapped\t% genes mapped\tGene count\tTotal exon length (bp)\n') #formats mapping_stats.tsv
	for file in os.listdir(path_to_file): #go through /fastq/ directory for the kallisto outputs to create mapping_stats.tsv and matrix_raw.tsv
		gene_count = 0 #container for number of genes that  had non-zero tpm values
		gene_len = 0 #container for total exon length
		name = file.split("_")[1] #ID of run
		processed_reads = '' #total reads in experiment
		mapped_reads = 0 #total reads mapped
		percent_mapped = '' #percentage of reads mapped
		CDS_count = 0 #number of genes in genome
		for file2 in os.listdir(path_to_file+file): #looks through the files in the kallisto directory
			#extracting information from kallisto output file for mapping_stats.tsv
			logging.info('Importing from ' + file2)

			if 'run_info.json' in file2: 
				#print('In directory ' + file)
				kallisto_json = ast.literal_eval(open(path_to_file+file+'/'+file2, 'r').read())
				processed_reads = str(kallisto_json["n_processed"])
				mapped_reads = kallisto_json["n_pseudoaligned"]
				percent_mapped = str(round((kallisto_json["n_pseudoaligned"] / kallisto_json["n_processed"])*100, 2))
				CDS_count =  kallisto_json["n_targets"]
				#appends corresponding tpm values to genes into the dictionary <dicto>
			elif 'abundance.tsv' in file2:
				output_header += name + '\t'
				content = open(path_to_file+file+'/'+file2, 'r')
				content.readline()
				for item in content:
					values = item.rstrip().split('\t')
					item = values[0]
					tpm = values[-1]
					gene_len += int(values[1])
					if tpm != '0':
						gene_count += 1
					if item in dicto:
						dicto[item].append(tpm)
					else:
						dicto[item] = [tpm]
		per_gene_mapped = str(round((gene_count/int(CDS_count))*100, 2)) #percentage of genes mapped
		#est_depth = str(round((mapped_reads*float(dict_log[name]))/gene_len, 2)) #(total number of mapped reads * average read length) / total length of all the exons
		map_stats.write(name+'\t'+ processed_reads+'\t'+str(mapped_reads)+'\t'+percent_mapped + '\t'+ per_gene_mapped + '\t' + str(CDS_count) + '\t' + str(gene_len) + '\n')
	map_stats.close()

	if '' in dicto:
		dicto.pop('')

	#creates full matrix file	
	for key, value in dicto.items():
		line = ''
		line += key + '\t'
		for item in value[:-1]:
			line += item + '\t'
		line += value[-1] + '\n'
		output_content += line
	output_header += '\n'
	output_file.write(output_header)
	output_file.write(output_content)
	logging.info('Matrix has been created')
	
	output_file.close()
	logging.info('Run has been completed'
