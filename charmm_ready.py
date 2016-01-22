#!/bin/usr/python
# $Id: charmm_ready.py 2099 2008-07-29 20:36:21Z denniej0 $

from string import strip
from string import strip
from math import sqrt
import re
import sys
import string

def fcn(mk_file):
            file = open(mk_file, 'r').readlines()
            new_file = []
            for i in range(len(file)):
		if file[i][:4] == 'ATOM':
	            	new_file.append((strip(file[i])))

            return new_file


''' Changes HIS to HSD and ILE CD1 to CD '''
def ca_file(pdb):
            file = fcn(pdb)
            new_file = []
            for i in range(len(file)):
               	if file[i][:4] == 'ATOM' and file[i][17:20] == 'HIS':
			new_file.append(file[i][:17]+'HSD  '+file[i][22:79])
			#print file[i][:17]+'HSD  '+file[i][22:79]
		elif file[i][:4] == 'ATOM' and file[i][17:20] == 'ILE' and file[i][13:16] == 'CD1':
			new_file.append(file[i][:13]+'CD  '+file[i][17:20]+'  '+file[i][22:79])
			#print file[i][:13]+'CD  '+file[i][17:20]+'  '+file[i][22:79]                    
		else:
			new_file.append(file[i][:20]+'  '+file[i][22:79])
        		#print file[i][:20]+'  '+file[i][22:79]
	    return new_file
	    

def renum(res, target):
	    file = ca_file(res)
	    #print target
	    new_file = open(target,'w')

	    x = 0
	    for i in range(len(file)):
		if file[i][23:26] != file[i-1][23:26]:
			x += 1
                        y=str(x)
			if len(str(x))==1:
				new_file.write(file[i][:22]+'   '+str(x)+file[i][26:79]+'\n')
				#print file[i][:22]+'   '+str(x)+file[i][26:79]
			elif len(str(x))==2:
				new_file.write(file[i][:22]+'  '+str(x)+file[i][26:79]+'\n')
				#print file[i][:22]+'  '+str(x)+file[i][26:79]
			else:
				new_file.write(file[i][:22]+' '+str(x)+file[i][26:79]+'\n')
        	                #print file[i][:22]+' '+str(x)+file[i][26:79]
		elif file[i][23:26] == file[i-1][23:26]:
			x = x
			y=str(x)
                        if len(str(x))==1:
                                new_file.write(file[i][:22]+'   '+str(x)+file[i][26:79]+'\n')
                                #print file[i][:22]+'   '+str(x)+file[i][26:79]
                        elif len(str(x))==2:
                                new_file.write(file[i][:22]+'  '+str(x)+file[i][26:79]+'\n')
                                #print file[i][:22]+'  '+str(x)+file[i][26:79]
                        else:
                                new_file.write(file[i][:22]+' '+str(x)+file[i][26:79]+'\n')
                                #print file[i][:22]+' '+str(x)+file[i][26:79]

	    new_file.close()

if __name__ == "__main__":
        protein = sys.argv[1]
	target = sys.argv[2]
	renum(protein, target)

