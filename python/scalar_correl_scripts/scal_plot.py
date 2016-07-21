import numpy
from string import split
import math
from lizhist2D import *
from pylab import *

dcdfile = open('../../pca/id_dcd_reduced','r').readlines()
for man in range(len(dcdfile)):
	try:
		M = dcdfile[man]
		t = M.split()
		if t[-1][:1] == '9':
			scal_matrix = numpy.load("%s_scal.pickle" % (t[-1][:5]) )
		else:
			scal_matrix = numpy.load("%s_scal.pickle" % (t[-1][:6]) )
		s = numpy.zeros((len(scal_matrix[0]),len(scal_matrix[0])))	
		for i in range(len(scal_matrix[0])):
			for j in range(len(scal_matrix[0])):
				s[i][j] = scal_matrix[i][j]/math.sqrt(scal_matrix[i][i]*scal_matrix[j][j])
		liz = [-1,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1]
		line = contourf(s,liz)
		colorbar()
		ylabel(" Residue Number ")
		xlabel(" Residue Number ")
		if t[-1][:1] == '9':
			savefig("%s_scal_cont" % (t[-1][:5]))
			print "hi"
		else:
			savefig("%s_scal_cont" % (t[-1][:6]))
		clf()
	except:
		continue	
