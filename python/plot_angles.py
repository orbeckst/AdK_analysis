import numpy
from pylab import *
import arrayio

for x,man in enumerate(['co','oc']):
	subplot(2,1,x+1)
	for i in range(1,252):
		try:
			file = '../project/rms/%s%s_rms%03d.out' % (man,man,i)
			f = arrayio.readFloatArrayWithSkip(file,[4])
			line = numpy.load("%s_%d_angles.pickle" % (man,i))
			test = plot(f,line[:,0][:-1])
		except:
			continue
	ylabel("Angle (degrees) ")
	title("%s" % (man))
xlabel(" progress variable (RMS) ")
savefig("angle_lid.eps")
clf()

for x,man in enumerate(['co','oc']):
        subplot(2,1,x+1)
        for i in range(1,252):
                try:
                        file = '../project/rms/%s%s_rms%03d.out' % (man,man,i)
                        f = arrayio.readFloatArrayWithSkip(file,[4])
                        line = numpy.load("%s_%d_angles.pickle" % (man,i))
                        test = plot(f,line[:,1][:-1])
                except:
                        continue
        ylabel("Angle (degrees) ")
        title("%s" % (man))
xlabel(" progress variable (RMS) ")
savefig("angle_nmp.eps")
clf()

all = []
for x,man in enumerate(['co','oc']):
        for i in range(1,252):
                try:
			file = '../project/rms/%s%s_rms%03d.out' % (man,man,i)
                        f = arrayio.readFloatArrayWithSkip(file,[4])
                        line = numpy.load("%s_%d_angles.pickle" % (man,i))
			for j in range(len(f)):
				all.append([f[j],line[j][0],line[j][1]])
		except:
			continue
all = numpy.array(all)
all.dump('hinge_angles.pickle')
liz = numpy.histogram2d(all[:,1],all[:,2])
area = contourf(liz[1][:-1],liz[2][:-1],liz[0])
colorbar()
xlabel("LID domain (degrees)")
ylabel("NMP domain (degrees)")
savefig("hinge.eps")
clf()
