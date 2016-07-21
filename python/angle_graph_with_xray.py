# $Id: angle_graph_with_xray.py 3014 2009-02-23 18:17:50Z denniej0 $
"""Create LID-CORE angle vs NMP-CORE angle graph from pickle.

Reads the pickle file that contains the angles. By default it reads
Oli's 'directions_angles.pickle' but this can be changed in the
script. 

"""

import cPickle
import numpy
import pylab
#from pylab import *
from pylab import plot,contour,contourf,colorbar,clf,cm

datafile2 = 'xray_angles.pickle'
datafile = 'directions_angles.pickle'

angles = cPickle.load(open('../../data/new_angles/DIMS/apo/standard/'+datafile))
xangles =  cPickle.load(open('../../data/new_angles/DIMS/apo/standard/'+datafile2))

if datafile == 'directions_angles.pickle':
    # dict with oc and co arrays:
    # 0   1
    import sys
    run_ids = sys.argv[1:]
    print run_ids
    try:
        all = numpy.concatenate([angles[run_id] for run_id in run_ids])
    except (KeyError,ValueError):
        raise ValueError("Give run_ids on the commandline, one or more of %r." % angles.keys())
    ha_LID, ha_NMP = all[:,0], all[:,1]
elif datafile == 'hinge_angles.pickle':
    # file format of data file:
    # 0                1               2
    # progress(RMSD)   LID-CORE angle  CORE-NMP angle
    ha_LID, ha_NMP = angles[:,1], angles[:,2]
else:
    raise ValueError("No idea how to deal with data file named %s" % datafile)


d_LID = ha_LID.max() - ha_LID.min()  # ranges
d_NMP = ha_NMP.max() - ha_NMP.min()

ratio = d_NMP/d_LID
bins_LID = 50                   # fix bins in LID direction
bins_NMP = int(ratio*bins_LID)  # scale bins proportionally in NMP

# use different number of bins 
h,e1,e2 = numpy.histogram2d(ha_LID,ha_NMP,bins=(bins_LID,bins_NMP),range=([94,161],[36,76]),normed=True)

# midpoints of bins
mLID = 0.5*(e1[1:]+e1[:-1])  # LID
mNMP = 0.5*(e2[1:]+e2[:-1])  # NMP

# turn density into free energy
F = -numpy.log(h+1e-100)   # add 1e-100 to 'mask' the empty (0) entries
Fshifted = F - F.min()     # make the endpoints the zero points

clf()
max = 6.5
# note that plotting the array reverses the axes (we use C-order):
# array = (row, col) <--> (y, x) !!
# run twice to smooth the image (!)
cont = contourf(mNMP,mLID,Fshifted,numpy.arange(0,max,0.1),extend='max',cmap=cm.hot_r)
cont = contourf(mNMP,mLID,Fshifted,numpy.arange(0,max,0.1),extend='max',cmap=cm.hot_r)
cont.ax.set_aspect('equal')
#cont.ax.set_aspect(1/1.618)

colorbar(extend='max',ticks=numpy.arange(0,max))

# light lines at 0.5 kT intervals
contour(mNMP,mLID,Fshifted,numpy.arange(0,max,0.5),colors='k',linewidths=0.2,alpha=0.5)
# heavy lines at 1kT
contour(mNMP,mLID,Fshifted,numpy.arange(0,max,1),colors='k',linewidths=0.5)

pylab.xlabel(r'angle NMP-CORE') ## [$^\circ$]') # use degreeFormatter below
pylab.ylabel(r'angle LID-CORE') ## [$^\circ$]')
pylab.title(r'Conditional free energy of transition')
pylab.ylim(95,160)
pylab.xlim(37,75)
ax = pylab.gca()
degreeFormatter = pylab.matplotlib.ticker.FormatStrFormatter(r'%d$^\circ$')
ax.xaxis.set_major_formatter(degreeFormatter)
ax.yaxis.set_major_formatter(degreeFormatter)

#print xangles.keys()
for i in xangles.items():
	try:
		print i
		x = [i[1][1]]
		y = [i[1][0]]
		#print x,y
		if 60 < x[0] :
			if  y[0] < 110:
				print i
		if y[0] > 140:
			print i
		plot(x,y,'wo')
	except:
		print "liz"
	#points = plot(all_xray[:,1], all_xray[:,0],'wo')
#print xangles.items()
pylab.ylim(95,160)
pylab.xlim(37,75)
pylab.draw()
pylab.savefig('angles_with_xray.pdf')
