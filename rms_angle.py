import numpy
from pylab import *
import arrayio


for i in ['co','oc']:
	for j in range(1,999):
		try:
			tot_d = []
			tot_rms = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%03d.out' % (i,i,j),[4])
        		tot_rms2 = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%03d.out' % (i,i[1]+i[0],j),[4])
			salt = numpy.load('../../data/new_angles/DIMS/apo/standard/%s%d_angles.pickle' % (i,j))        	
			for man in range(len(salt[:,0])-1):
				rms = tot_rms[man]
        			rms2 = tot_rms2[man]
        			deltaRMSD = rms2 - rms
				print deltaRMSD
				tot_d.append([rms,rms2,deltaRMSD,salt[:,0][man],salt[:,1][man]])
			#print tot_d
			dist_profile = numpy.array([tot_d])
                	dist_profile.dump('%s_%d_angle_deltaRMSD.pickle' % (i,j) )
		except:
			continue

