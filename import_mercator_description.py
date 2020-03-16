#change text file name
input = open(' ', 'r')
matrix = open('mercator.txt', 'r')

lst = []
for line in input:
   	line = line.strip('\n')
   	if line not in lst:
   		lst.append(line)

k = len(lst)
#change text file name
with open(' ', 'w') as f:

    for row in matrix:
    	gene = row.strip('\n').split('\t')
    	gene = gene[2]
    	gene = gene.replace('\'',"")

    	for j in range(k):
    		if lst[j] == gene:
    			f.write(row)

print('Done')

