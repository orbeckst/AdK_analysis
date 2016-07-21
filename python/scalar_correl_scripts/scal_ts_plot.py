import numpy
from string import split
import math
from lizhist2D import *
from pylab import *


s = numpy.ones((156,75))
s *= 10
dcdfile = open('ave','r').readlines()
for man in range(len(dcdfile)):
	try:
		M = dcdfile[man]
		t = M.split()
		if t[-1][:1] == '9':
			scal_matrix = numpy.load("%s" % (t[-1]) )
			lid = int(t[-1][:2])
			nmp = int(t[-1][3:5])
		else:
			scal_matrix = numpy.load("%s" % (t[-1]) )
			lid = int(t[-1][:3])
			nmp = int(t[-1][4:6])
		for i in [117]: ##range(len(scal_matrix[0])):
			for j in [135]:  ##range(len(scal_matrix[0])):
				s[lid][nmp] = scal_matrix[i][j]/math.sqrt(scal_matrix[i][i]*scal_matrix[j][j])
	except:
		continue  
liz = [-1,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1]
#lid_array = [95,105,115,125,135,145,155]
#nmp_array = [40,45,50,55,60,65,70,75]
line = contourf(s,liz)
colorbar()
xlim(40,75)
ylim(95,156)
ylabel(" LID ANGLE ")
xlabel(" NMP ANGLE ")
savefig("%d_%d_ts" % (i+1,j+1))
clf()
