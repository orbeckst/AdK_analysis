#!/usr/bin/env python
# $Id: angles_adk.py 2669 2008-12-22 16:23:50Z oliver $
"""Process DIMS AdK transitions: extract the angle between LID and CORE domain
and NMP and CORE domain and save the angle pairs as pickled numpy arrays.

Format of array
---------------

 shape (2,N)               N = number of saved steps in trajectory
 [LID_angle, NMP_angle]    angles in degree

Note on angle definitions
-------------------------

The definitions of the angles are described in
txt/observables.txt. They *must* match between here and the MMFP GEO
restraints in pmf/pmf_angles.inp and pmf/pmf_temp.inp.

The corresponding CHARMM code:

   MMFP
   GEO  MAXGEO 5000 -
        sphere RCM angle -
        harmonic symmetric force  @FORCECONSTANT tref @NMP dtoff 0.0 -
        select resid 115:125 .and. backbone end     select resid 90:100 .and. backbone end     select resid 35:55 .and. backbone end 
   GEO  sphere RCM angle -
        harmonic symmetric force  @FORCECONSTANT tref @LID dtoff 0.0 -
        select resid 179:185 .and. backbone end     select resid 115:125 .and. backbone end     select resid 125:153 .and. backbone end 
   END
"""

from MDAnalysis import *
from numpy import *
import numpy
from Scientific.Geometry import Vector

import sys

try:
	from numpy import degrees
except ImportError:
	_rad2deg = 180./numpy.pi
	def degrees(x):
		return _rad2deg * x

args = sys.argv[1:]
if len(args) < 2:
	raise ValueError("usage: angles_adk.py [oc|co] istart istop\n"
			 "  leaving out the direction sets both\n"
			 "  istart   first trajectory number\n"
			 "  istop    last (incl.) trajectory number\n")
if len(args) == 3:
	direction = args.pop(0)
	if direction not in ('oc','co'):
		raise ValueError("direction must be either 'oc' or 'co'.")
	directions = [direction]
else:
	directions = ['co','oc']

istart, istop = int(args[0]), int(args[1]) + 1


print "Analyzing trajectory numbers %d to %d inclusive for directions %r." % (istart,istop-1,directions) 

for man in directions:
	path = '/home/denniej0/oli'
	universe = Universe('%s/open_final.psf' % (path), '%s/project/%s001.dcd' % (path,man) )
	for naso in xrange(istart,istop):
		dcdfile = '%s/project/%s%03d.dcd' % (path,man,naso)
		savefile = '%s_%d_angles.pickle' % (man,naso)    
		print "dcdfile = %s" % dcdfile
		try:
			tot_d = []
			universe.load_new_dcd(dcdfile)
		except:
			print "Failed to load"
			continue
		try:
			for ts in universe.dcd:
				# angles: (top,core,bot)  
				# v0 = top - core, v1 = bot - core
				# cos alpha = v0*v1/v0v1

                                # These defintions match the ones in the MMFP GEO restraints
                                # of pmf_angles.inp or pmf_bat.inp                                
                                nbot = universe.selectAtoms("resid 115:125 and backbone").centerOfGeometry()
                                ntop = universe.selectAtoms("resid 35:55 and backbone").centerOfGeometry()
                                ncore = universe.selectAtoms("resid 90:100 and backbone").centerOfGeometry()

                                lidc = universe.selectAtoms("resid 115:125 and backbone").centerOfGeometry()
                                lidt = universe.selectAtoms("resid 125:153 and backbone").centerOfGeometry()
                                lidb = universe.selectAtoms("resid 179:185 and backbone").centerOfGeometry()

				l_v0 = Vector(lidb - lidc)
				l_v1 = Vector(lidt - lidc)
				n_v0 = Vector(ntop - ncore)
				n_v1 = Vector(nbot - ncore)

				LID_angle = degrees( l_v1.angle(l_v0) )
				NMP_angle = degrees( n_v1.angle(n_v0) )

				tot_d.append([LID_angle,NMP_angle])
			dist_profile = numpy.array(tot_d)
			dist_profile.dump(savefile)
			print "Saved %d angle-pairs to '%s'" % (len(tot_d),savefile)
			del dist_profile
			del tot_d
		except:
			print "Failed to process"
			continue
print "DONE"
