from MDAnalysis import *
import numpy
from string import split

system =  AtomGroup.Universe("../pca/ca.psf" , "../pca/130_55.dcd")
prot = system.selectAtoms(" ( protein and name CA ) ")

dcdfile = open('../pca/id_dcd','r').readlines()
for man in range(len(dcdfile)):
	M = dcdfile[man]
	t = M.split()
	scal_matrix = numpy.zeros((len(prot),len(prot)))
	try:
		system.load_new_dcd("../pca/%s" % (t[-1]))
		ca = system.dcd.timeseries(prot)
		for i in range(len(ca[0])-1):
			y = i + 1
			x = ca[:,i]-ca[:,y]
			for j in range(len(x[:,0])):
				for k in range(len(x[:,0])):
					scal_matrix[j][k] += numpy.dot(x[:,0][j],x[:,0][k])
		scal_matrix /= system.dcd.numframes
		if t[-1][:1] == '9':
			scal_matrix.dump("%s_scal.pickle" % (t[-1][:5]) )
		else:
			scal_matrix.dump("%s_scal.pickle" % (t[-1][:6]) )
	except:
		continue	
