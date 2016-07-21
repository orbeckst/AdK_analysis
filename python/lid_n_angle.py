from MDAnalysis_old import *
from numpy import *
import numpy
from Scientific.Geometry import Vector

for man in ['co','oc']:
	path = '/home/denniej0/oli'
	universe = Universe('%s/open_final.psf' % (path), '%s/project/%s001.dcd' % (path,man) )
	
	for naso in range(230,252):
		try:
			tot_d = []
			universe.load_new_dcd('%s/project/%s%03d.dcd' % (path,man,naso))
			for ts in universe.dcd:
				ndom = universe.selectAtoms(" resid 30:60 and ( backbone or name CB ) ").centerOfGeometry().copy()
        			lidd = universe.selectAtoms(" resid 118:160 and ( backbone or name CB ) ").centerOfGeometry().copy()
	
				box = universe.dcd.ts.dimensions[:3]
				ndom = numpy.array([ndom])
				lidd = numpy.array([lidd])
				u1 = distances.distance_array(ndom,lidd, box)
				dist = numpy.array([min(u1[0])])	
				new_coor = lidd - ndom
				coor = Vector(float((new_coor[0][0])),float((new_coor[0][1])), float((new_coor[0][2])))
				theta = 180.*math.atan(coor.x()/coor.y())/3.1416
				tot_d.append([dist[0],theta])
			dist_profile = numpy.array(tot_d)
			dist_profile.dump('%s_%d_angle-dist.pickle' % (man,naso) )
		except:
			continue
