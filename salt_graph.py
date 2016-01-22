# $Id: salt_graph_with_xray.py 2526 2008-12-06 18:23:42Z www-data $
"""Create deltaRMSD-CORE salt vs salt1-CORE salt graph from pickle.

Reads the pickle file that contains the salt. By default it reads
Oli's 'directions_salt.pickle' but this can be changed in the
script. 

"""

import cPickle
import numpy
import pylab
#from pylab import *
from pylab import plot,contour,contourf,colorbar,clf,cm,subplot

datafile2 = 'xray_salt.pickle'
datafile = 'directions_salt_bridge.pickle'

salt = cPickle.load(open('../../data/salt_bridge/'+datafile))
#xsalt =  cPickle.load(open(datafile2))

if datafile == 'directions_salt_bridge.pickle':
    # dict with oc and co arrays:
    # 0   1
    import sys
    run_ids = sys.argv[1:]
    print run_ids
    try:
        all = numpy.concatenate([salt[run_id] for run_id in run_ids])
    except (KeyError,ValueError):
        raise ValueError("Give run_ids on the commandline, one or more of %r." % salt.keys())
    ha_deltaRMSD, ha_salt1, ha_salt2, ha_salt3, ha_salt4, ha_salt5, ha_salt6, ha_salt7, ha_salt8, ha_salt9 = all[:,0], all[:,1], all[:,2], all[:,3], all[:,4], all[:,5], all[:,6], all[:,7], all[:,8], all[:,9]
elif datafile == 'hinge_salt.pickle':
    # file format of data file:
    # 0                1               2
    # progress(RMSD)   deltaRMSD-CORE salt  CORE-salt1 salt
    ha_deltaRMSD, ha_salt1 = salt[:,1], salt[:,2]
else:
    raise ValueError("No idea how to deal with data file named %s" % datafile)

for x in range(1,2):
	subplot(3,3,x)
	d_deltaRMSD = ha_deltaRMSD.max() - ha_deltaRMSD.min()  # ranges
	d_salt1 = ha_salt9.max() - ha_salt9.min()

	ratio = d_salt1/d_deltaRMSD
	bins_deltaRMSD = 50                   # fix bins in deltaRMSD direction
	bins_salt1 = int(ratio*bins_deltaRMSD)  # scale bins proportionally in salt1

	# use different number of bins 
	h,e1,e2 = numpy.histogram2d(ha_deltaRMSD,ha_salt9,bins=(bins_deltaRMSD,bins_salt1),range=([-6.6,6.6],[-1,31]),normed=True)

	# midpoints of bins
	mdeltaRMSD = 0.5*(e1[1:]+e1[:-1])  # deltaRMSD
	msalt1 = 0.5*(e2[1:]+e2[:-1])  # salt1
	
	# turn density into free energy
	F = -numpy.log(h+1e-100)   # add 1e-100 to 'mask' the empty (0) entries
	Fshifted = F - F.min()     # make the endpoints the zero points
	Fshifted = numpy.transpose(Fshifted)	
	clf()
	max = 6.5
	# note that plotting the array reverses the axes (we use C-order):
	# array = (row, col) <--> (y, x) !!
	# run twice to smooth the image (!)
	pylab.ylim(0,30)
	cont = contourf(mdeltaRMSD,msalt1,Fshifted,numpy.arange(0,max,0.1),extend='max',cmap=cm.hot_r)
	cont = contourf(mdeltaRMSD,msalt1,Fshifted,numpy.arange(0,max,0.1),extend='max',cmap=cm.hot_r)
	cont.ax.set_aspect('equal')
	#cont.ax.set_aspect(1/1.618)

	colorbar(extend='max',ticks=numpy.arange(0,max))
	
	# light lines at 0.5 kT intervals
	contour(mdeltaRMSD,msalt1,Fshifted,numpy.arange(0,max,0.5),colors='k',linewidths=0.2,alpha=0.5)
	# heavy lines at 1kT
	contour(mdeltaRMSD,msalt1,Fshifted,numpy.arange(0,max,1),colors='k',linewidths=0.5)
	
	pylab.xlabel(r'deltaRMSD') ## [$^\circ$]') # use degreeFormatter below
	pylab.ylabel(r'salt bridge distance') ## [$^\circ$]')
	pylab.ylim(0,30)	
	pylab.xlim(-6.5,6.5)
	ax = pylab.gca()
	#degreeFormatter = pylab.matplotlib.ticker.FormatStrFormatter(r'%d$^\circ$')
	#ax.xaxis.set_major_formatter(degreeFormatter)
	#ax.yaxis.set_major_formatter(degreeFormatter)
'''	
#print xsalt.keys()
for i in xsalt.items():
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
#print xsalt.items()
'''
pylab.draw()
pylab.savefig('salt9_without_xray.png')
