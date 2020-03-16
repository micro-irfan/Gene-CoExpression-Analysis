# -*- coding: utf-8 -*-
# Removes redundant information from Mercator

import openpyxl
import re

results = openpyxl.load_workbook('mercator.xlsx')
sheet = results.get_sheet_by_name('mercator')
pathtolocal = --PATH TO LOCAL--
output_file = open(pathtolocal + '/mercatorcleanupALL.txt', 'a+')
output_file.write('Identifier\tDescription\n')

for row in range(2, sheet.max_row + 1):
	identifier = sheet['C' + str(row)].value
	identifier1 = identifier.replace('\'',"")
	description = sheet['D' + str(row)].value
	description1 = description.replace('\'',"")
	if re.search(r'^\(.+?\)', description1):
		new = re.sub(r'^\(.+?\)', '', description1)
		
		output_file.write(identifier1+ '\t'+ new + '\n')
	else:
		output_file.write(identifier1+ '\t'+ description1 + '\n')

output_file.close()
print ("clean up executed")
