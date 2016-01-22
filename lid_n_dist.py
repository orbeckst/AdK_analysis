from MDAnalysis_old import *
from numpy import *
import numpy

for man in ['co','oc']:
	path = '/home/denniej0/oli'
	universe = Universe('%s/open_final.psf' % (path), '%s/project/%s001.dcd' % (path,man) )
	
	for naso in range(230,252):
		tot_d = []
		try:
			universe.load_new_dcd('%s/project/%s%03d.dcd' % (path,man,naso))
			for ts in universe.dcd:
				ndom = universe.selectAtoms(" resid 30:60 and ( backbone or name CB ) ").centerOfGeometry().copy()
        			lidd = universe.selectAtoms(" resid 118:160 and ( backbone or name CB ) ").centerOfGeometry().copy()
	
				box = universe.dcd.ts.dimensions[:3]
				ndom = numpy.array([ndom])
				lidd = numpy.array([lidd])
				u1 = distances.distance_array(ndom,lidd, box)
				dist = numpy.array([min(u1[0])])	
				tot_d.append(dist)
			dist_profile = numpy.array([tot_d])
			dist_profile.dump('%s_%d_dist.pickle' % (man,naso) )
		except:
			continue
