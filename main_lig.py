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

import lig_alter
import remove_alter
import charmm_ready

# Generating, analyzing data and 


for man in ["adp","lig"]:
	m = open('all_rms_%s.dat' % (man),'w')
	for post in range(1,311):
		pos = "%03d" % (post)
		try:
			remove_alter.fcn(man,pos)
			command = '~/c35a1-dims-intel lig:%s <remove_lig%s.inp > lig.out' % (man,pos)
			os.system(command)
			command = 'mv co%s_%s_new.dcd co%s_%s.dcd' % (pos,man,pos,man)
			os.system(command) 
	
			ids = [s[6:12] for s in open('./all_pdb_id')]
			for target in ids:
				#x = 'AKeco_%s.pdb' % (target)
				#y = '%s.pdb' % (target.lower())
				#charmm_ready.renum(x,y)
				lig_alter.fcn(man, target.lower(), pos)
				command = '~/c35a1-dims-intel <rms_calc%s%s_%s.inp | grep "RMS DIFF IS" > rms.out ' % (man,target.lower(),pos)
				print command
				os.system(command)
	
				f = 'rms.out' ##% (man,target,pos)
	                        rms = arrayio.readFloatArrayWithSkip(f,[4])
				rms_min = numpy.min(rms)
				for i,j in enumerate(rms):
					if j == rms_min:
						file_len = len(rms)
						tot_rms = arrayio.readFloatArrayWithSkip('/home/denniej0/oli/project/rms/cooc_rms%s_%s.out' % (pos,man),[4])
	                                        tot_rms2 = arrayio.readFloatArrayWithSkip('/home/denniej0/oli/project/rms/coco_rms%s_%s.out' % (pos,man),[4])
						tot_rms = tot_rms[i]
						tot_rms2 = tot_rms2[i]
						deltaRMSD = tot_rms - tot_rms2
						m.write(man+'\t'+str(pos)+'\t'+target+'\t'+str(i)+'\t'+str(file_len)+'\t'+str(j)+'\t'+str(tot_rms2)+'\t'+str(deltaRMSD)+'\n')
		except:
			continue
	m.close()
