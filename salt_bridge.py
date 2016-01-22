from MDAnalysis import *
from numpy import *
import numpy

for man in ['co','oc']:
	path = '/home/denniej0/oli'
	universe = Universe('%s/open_final.psf' % (path), '%s/project/%s001.dcd' % (path,man) )
	
	for naso in range(1,901):
		tot_d = []
		try:
			universe.load_new_dcd('%s/project/%s%03d.dcd' % (path,man,naso))
			for ts in universe.dcd:
				asp33 = universe.selectAtoms(" resid 33 and ( name OD1 or name OD2 ) ").centerOfGeometry().copy()
        			arg156 = universe.selectAtoms(" resid 156 and ( name NH1 or name NH2 ) ").centerOfGeometry().copy()
				asn190 = universe.selectAtoms(" resid 190 and ( name OD1 or name ND2 ) ").centerOfGeometry().copy()
                                lys97 = universe.selectAtoms(" resid 97 and ( name NZ ) ").centerOfGeometry().copy()
				asp54 = universe.selectAtoms(" resid 54 and ( name OD1 or name OD2 ) ").centerOfGeometry().copy()
				lys157 = universe.selectAtoms(" resid 157 and ( name NZ ) ").centerOfGeometry().copy()	
				arg167 = universe.selectAtoms(" resid 167 and ( name NH1 or name NH2 ) ").centerOfGeometry().copy()
				arg36 = universe.selectAtoms(" resid 36 and ( name NH1 or name NH2 ) ").centerOfGeometry().copy()
				asp158 = universe.selectAtoms(" resid 158 and ( name OD1 or name OD2 ) ").centerOfGeometry().copy()
				glu170 = universe.selectAtoms(" resid 170 and ( name OE1 or name OE2 ) ").centerOfGeometry().copy()
				lys57 = universe.selectAtoms(" resid 57 and ( name NZ ) ").centerOfGeometry().copy()
				box = universe.dcd.ts.dimensions[:3]
				#asp33 = numpy.array(asp33)
				#arg156 = numpy.array(arg156)
				#asn190 = numpy.array(asn190)
				#lys97 = numpy.array(lys97)
				#asp54 = numpy.array(asp54)
				#lys157 = numpy.array(lys157)
				#arg167 = numpy.array(arg167)
				#arg36 = numpy.array(arg36)
				#asp158 = numpy.array(asp158)
				#glu170 = numpy.array(glu170)
				#lys57 = numpy.array(lys57)
				asp33 = numpy.array([asp33])
                                arg156 = numpy.array([arg156])
                                asn190 = numpy.array([asn190])
                                lys97 = numpy.array([lys97])
                                asp54 = numpy.array([asp54])
                                lys157 = numpy.array([lys157])
                                arg167 = numpy.array([arg167])
                                arg36 = numpy.array([arg36])
                                asp158 = numpy.array([asp158])
                                glu170 = numpy.array([glu170])
                                lys57 = numpy.array([lys57])
				u1 = distances.distance_array(asp33,arg156, box)
				u2 = distances.distance_array(lys97,asn190,box)
				u3 = distances.distance_array(asp54,lys157,box)
				u4 = distances.distance_array(asp54,arg156,box)
				u5 = distances.distance_array(asp54,arg167,box)
				u6 = distances.distance_array(arg36,asp158,box)
				u7 = distances.distance_array(arg36,glu170,box)
				u8 = distances.distance_array(lys57,asp158,box)
				u9 = distances.distance_array(lys57,glu170,box)
				dist = numpy.array([min(u1[0]),min(u2[0]),min(u3[0]),min(u4[0]),min(u5[0]),min(u6[0]),min(u7[0]),min(u8[0]),min(u9[0])])	
				tot_d.append(dist)
			dist_profile = numpy.array([tot_d])
			#print tot_d
			dist_profile.dump('/home/denniej0/oli/AdK/data/salt_bridge/%s_%d_salt.pickle' % (man,naso) )
		except:
			continue
