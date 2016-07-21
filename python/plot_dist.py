import numpy
from pylab import *
import arrayio

for x,man in enumerate(['co','oc']):
	subplot(2,1,x+1)
	for i in range(1,252):
		if x == 0:
			y = man[1]+man[0]
		else:
			y = man[1]+man[0]
		try:
			file = '../project/rms/%s%s_rms%03d.out' % (man,y,i)
			f = arrayio.readFloatArrayWithSkip(file,[4])
			line = numpy.load("%s_%d_dist.pickle" % (man,i))
			test = plot(f,line[0][1:])
		except:
			continue
	ylabel("LID/N distance (A) ")
	title("%s" % (man))
xlabel(" progress variable (RMS) ")
savefig("dist_lid_n.eps")
clf()
