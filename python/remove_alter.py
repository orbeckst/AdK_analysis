#!/bin/usr/python
# $Id: rota_alter.py 2350 2008-10-07 15:56:14Z denniej0 $

from string import split
from string import strip
from math import sqrt
from re import search
import re
import sys            
import string
import math 

def fcn(man,pos):
	file = open('remove_lig.inp', 'r').readlines()
	m = open('remove_lig%s.inp' % (pos), 'w')
	for i in range(len(file)):
		#print file[i][:]
		if not file[i].find('FORCE') == -1: #in file[i]:
			p = re.compile('(FORCE)')
			m.write(p.sub("co"+str(pos)+"_"+str(man), file[i]))
		else:
			m.write(file[i])
	m.close()


if __name__ == "__main__":
	man = str(sys.argv[1])
	pos = int(sys.argv[2])
	fcn(man,pos)
