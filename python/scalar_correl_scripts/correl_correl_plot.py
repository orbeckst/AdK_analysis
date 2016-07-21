import numpy
from string import split
import math
from lizhist2D import *
from pylab import *


dcdfile = open('correl','r').readlines()
for i in range(0,9):
	for j in range(0,9):
		s = numpy.ones((156,75))
		s *= 10
		for man in range(len(dcdfile)):
			try:
				M = dcdfile[man]
				t = M.split()
				#print t[-1][5:7]
				if t[-1][5:6] == '9':
					scal_matrix = numpy.load("%s" % (t[-1]) )
					lid = int(t[-1][5:7])
					nmp = int(t[-1][8:10])
				else:
					scal_matrix = numpy.load("%s" % (t[-1]) )
					lid = int(t[-1][5:8])
					nmp = int(t[-1][9:11])
				s[lid][nmp] = scal_matrix[i][j]/math.sqrt(scal_matrix[i][i]*scal_matrix[j][j])
				#print s
			except:
				continue 
		liz = [-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0] 
		line = contourf(s,liz)
		colorbar()
		xlim(40,75)
		ylim(95,156)
		ylabel(" LID ANGLE ")
		xlabel(" NMP ANGLE ")
		savefig("correl_correl_%d_%d.eps" % (i,j))
		clf()
