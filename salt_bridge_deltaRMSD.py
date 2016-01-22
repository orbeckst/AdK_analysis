import numpy
from pylab import *
import arrayio

#file = open('../../data/salt_bridge.txt','r').readlines()

for i in ['oc']:
	for j in range(1,300):
		try:
			tot_d = []
			tot_rms = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%03d.out' % (i,i,j),[4])
        		tot_rms2 = arrayio.readFloatArrayWithSkip('../../data/rms/%s%s_rms%03d.out' % (i,i[1]+i[0],j),[4])
			salt = numpy.load('../../data/salt_bridge/%s_%d_salt.pickle' % (i,j))        	
			for man in range(len(salt[0][:,0])-1):
				rms = tot_rms[man]
        			rms2 = tot_rms2[man]
        			deltaRMSD = rms2 - rms
				#print deltaRMSD
				tot_d.append([deltaRMSD,salt[0][man][0],salt[0][man][1],salt[0][man][2],salt[0][man][3],salt[0][man][4],salt[0][man][5],salt[0][man][6],salt[0][man][7],salt[0][man][8]])
			#print tot_d
			dist_profile = numpy.array([tot_d])
                	dist_profile.dump('../../data/salt_bridge/%s_%d_salt_deltaRMSD.pickle' % (i,j) )
		except:
			continue

