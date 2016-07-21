#!/bin/usr/python
# $Id: main.py 2182 2008-08-14 14:50:15Z denniej0 $

from Scientific.IO import ArrayIO
from string import split
from string import strip
from re import search
import re
import sys            
import string
import math 
import os
import arrayio
import numpy

# instaniating the file histo.py

import rota_alter
import charmm_ready

# Generating, analyzing data and 

m = open('all_rms_300-899.dat','w')

for man in ["co","oc"]:
	for post in range(300,900):
		pos = "%03d" % (post)
		try:
			command = '/share/apps/analysis/catdcd -o %s%s.dcd -first 0 -last 5000000 -stride 1 /nfs/greenwulf/xenon/denniej0/oli/project/%s/%s%s.dcd' % (man,pos,man,man,pos)
                	#print command
                	os.system(command)
			#rota_alter.fcn(man,pos)
			#command = '/share/apps/charmm/c32b1_large lig:%s <remove_lig%03d.inp > lig.out'
			#os.system(command)
	
			ids = [s[6:12] for s in open('./all_pdb_id')]
			for target in ids:
				x = 'AKeco_%s.pdb' % (target)
				y = '%s.pdb' % (target.lower())
				charmm_ready.renum(x,y)
	
				rota_alter.fcn(man, target.lower(), pos)
				command = '~/c35a1-dims-intel <rms_calc%s%s_%s.inp | grep "RMS DIFF IS" > rms.out ' % (man,target.lower(),pos)
				print command
				os.system(command)
	
				f = 'rms.out' ##% (man,target,pos)
	                        rms = arrayio.readFloatArrayWithSkip(f,[4])
				rms_min = numpy.min(rms)
				for i,j in enumerate(rms):
					if j == rms_min:
						if man == 'co':
							file_len = len(rms)
							tot_rms = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%s.out' % (man,man,pos),[4])
	                                        	tot_rms2 = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%s.out' % (man,man[1]+man[0],pos),[4])
							tot_rms = tot_rms[i]
							tot_rms2 = tot_rms2[i]
							deltaRMSD = tot_rms - tot_rms2
							m.write(man+'\t'+str(pos)+'\t'+target+'\t'+str(i)+'\t'+str(file_len)+'\t'+str(j)+'\t'+str(tot_rms)+'\t'+str(deltaRMSD)+'\n')
						if man == 'oc':
							file_len = len(rms)
                                                        tot_rms = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%s.out' % (man,man,pos),[4])
                                                        tot_rms2 = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%s.out' % (man,man[1]+man[0],pos),[4])
                                                        tot_rms = tot_rms[i]
                                                        tot_rms2 = tot_rms2[i]
                                                        deltaRMSD = tot_rms2 - tot_rms
                                                        m.write(man+'\t'+str(pos)+'\t'+target+'\t'+str(i)+'\t'+str(file_len)+'\t'+str(j)+'\t'+str(tot_rms2)+'\t'+str(deltaRMSD)+'\n')
	
		except:
			continue
m.close()
