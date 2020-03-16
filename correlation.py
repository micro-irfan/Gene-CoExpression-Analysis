import math

def correlatedGenes(targetGene, path_to_file, inputFile):
	'''Generate Correlation Values between Bait Gene and Other Genes
	Args:
		targetGene - Bait Gene (Gene of Comparison)
		path_to_file - path to folder
		inputFile - matrix files with values from kallisto		
	'''
	path_to_file_ = path_to_file + '_' + targetGene + '_correlation.txt'
	output_file = open(path_to_file_, 'a+')
	output_header = 'Bait Gene\tGene\tcorr\n'
	output_file.write(output_header)

	#input file -> matrix file with kallisto output 

	with open(inputFile, 'r') as f:
		geneLst = []
		valueLst = []
		targetValues  = []

		for line in f:
			line = line.strip('\n').split('\t')
			gene = line[0]
			try:
				value = [float(x) for x in line[:1]]
				valueLst.append(value)
			except: pass

			if gene not in geneLst:
				geneLst.append(gene)
			if gene == targetGene:
				targetValues = [float(x) for x in line[:1]]

		newlst = genelst.pop(0)

	k = len(genelst) #Total number of genes

	for i in range(k):
		gene = geneLst[i]
		value = valueLst[i]
		try:
			corr = str(pearson_def(targetValues, value))
			output_file.write(targetGene + '\t' + gene + '\t' + corr + '\n')
		except: 
			output_file.write(targetGene + '\t' + gene + '\t' + str(0) + '\n')

	output_file.close()

def average(x):
	assert len(x) > 0
	return float(sum(x)) / len(x)

def pearson_def(x, y):
	'''pearson Correlation'''

	assert len(x) == len(y)
	n = len(x)
	assert n > 0
	avg_x = average(x)
	avg_y = average(y)
	diffprod = 0
	xdiff2 = 0
	ydiff2 = 0
	for idx in range(n):
	    xdiff = x[idx] - avg_x
	    ydiff = y[idx] - avg_y
	    diffprod += xdiff * ydiff
	    xdiff2 += xdiff * xdiff
	    ydiff2 += ydiff * ydiff

	return diffprod / math.sqrt(xdiff2 * ydiff2) 

def sortGenes(targetGene, path_to_file, input_file, corr = 0, N = 100):
	'''Sort Genes correlation values in descending order
	Args:
		targetGene - Bait Gene (Gene of Comparison)
		path_to_file - path to folder
		input_file - File with Correlation Values
		corr - cut-off correlation value for analysis
		N - cut-off number of Genes for analysis 	
	'''
	path_to_file_ = path_to_file + '_' + targetGene + '_bestCorrelation.txt'
	output_file = open(path_to_file_, 'a+')

	correlationList = []

	#input_file -> correlated Genes File

	with open(input_file) as f:
		next(f)
		for line in f:
			line = line.strip('\n').split('\t')
			correlationList.append(line)
					
	#output_file.close()
	sorted_correlationLst = Sort(correlationList)
	topGenes = []
	if corr != 0 and corr > 0.5:
		for i in sorted_correlationLst:
			if corr < i[2]:
					topGenes.append(i[1])
	else:
		for i in range(0, N+1):
			topGenes.append(sorted_correlationLst[i][1])

	for i in topGenes:
		output_file.write(i + '\n')
		
def Sort(sub_li): 
	'''Takes a list and arrange according to the correlation values in descending order
	Args:
		sub_li - table of bait Gene, corresponding gene and correlation value
	
	# reverse = None (Sorts in Ascending order) 
	# key is set to sort using last element of  
	# sublist lambda has been used 
	'''

	sub_li = [[x[0], x[1], float(x[2])] for x in sub_li]
	sub_li.sort(reverse = True, key = lambda x: x[2]) 
	return sub_li 

def importMatrix(targetGene, path_to_file, inputMatrixFile):
	'''Import kallisto values from matrix file to a new matrix file for the bait gene
	Args:
		targetGene - Bait Gene (Gene of Comparison)
		path_to_file - path to folder
		inputMatrixFile - Original Matrix file
	'''

	path_to_bestCorr = path_to_file + '_' + targetGene + '_bestCorrelation.txt'
	path_to_file_ = path_to_file + '_' + targetGene + '_bestCorrelationMatrix.txt'
	#output_file = open(path_to_file_, 'a+')
	matrix = inputMatrixFile
	matrix = open(matrix, 'r')

	listOfGenes = []
	for line in path_to_bestCorr:
		line = line.strip('\n')
		if line not in listOfGenes:
			listOfGenes.append(line)

	k = len(listOfGenes)
	value = []
	with open(path_to_file_, 'a+') as f:
		for row in matrix:
			gene = row.strip('\n').split('\t')
			gene = gene[0]
			try:
				values = [float(x) for x in gene[1:]]
			except: pass
				
			for j in range(k):
				if listOfGenes[j] == gene:
					f.write(row)

def createMatrixFile(targetGene, path_to_file):
	'''Sort Genes correlation values in descending order
	Args:
		targetGene - Bait Gene (Gene of Comparison)
		path_to_file - path to folder
		input_file - File with Correlation Values
		corr - cut-off correlation value for analysis
		N - cut-off number of Genes for analysis 	
	'''

	path_to_bestCorrMatrix = path_to_file + '_' + targetGene + '_bestCorrelationMatrix.txt'
	path_to_file_ = path_to_file + '_' + targetGene + '_NetworkAnalysis.txt'
	output_file = open(path_to_file_, 'a+')

	with open(path_to_bestCorrMatrix, 'r') as f:
		lst = []
		value = []

	for line in f:
		line = line.strip('\n').split('\t')
		gene = (line[0])
		try:
			values = [float(x) for x in line[1:]]
			value.append(values)
		except: pass

		if gene not in lst:
			lst.append(gene)

		# loop in loop to generate
		# correlation between genes 
		# for network analysis

	for i in range(len(lst)):
		for j in range(i):
			gene1 = lst[i]
			gene2 = lst[j] 
			value1 = value[i]
			value2 = value[j]

			try:
				corr = str(pearson_def(value1, value2))
				output_file.write(gene1 + '\t' + gene2 + '\t' + corr + '\n')
			except: #some of the sum of values = 0 
				zero = 0 
				output_file.write(gene1 + '\t' + gene2 + '\t' + str(zero) + '\n')


def sortGenesAfterMatrixFileCreated(targetGene, path_to_file, input_file = None, corr = 0.7, N = 0):
	'''Sort Genes correlation values of bait Genes and top correlated genes in descending order
	Args:
		targetGene - Bait Gene (Gene of Comparison)
		path_to_file - path to folder
		input_file - File with Correlation Values
		corr - cut-off correlation value for analysis
		N - cut-off number of Genes for analysis 	
	'''
	path_to_file_ = path_to_file + '_' + targetGene + '_bestCorrelation_toCreateMatrix.txt'
	output_file = open(path_to_file_, 'a+')
	if not input_file:
		input_file =  path_to_file + '_' + targetGene + '_bestCorrelation.txt'
	else:
		input_file = path_to_file + '/' + input_file
	correlationList = []

	#input_file -> correlated Genes File

	with open(input_file) as f:
		next(f)
		for line in f:
			line = line.strip('\n').split('\t')
			correlationList.append(line)
					
	#output_file.close()
	sorted_correlationLst = Sort(correlationList)
	topGenes = []
	if N > 400 and N != 0:
		for i in range(0, N+1):
			topGenes.append(sorted_correlationLst[i][1])
	else:
		for i in sorted_correlationLst:
			if corr < i[2]:
				topGenes.append(i[1])

	for i in topGenes:
		output_file.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\n')
