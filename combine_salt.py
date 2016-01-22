#!/usr/bin/env python
# $Id: combine_salt_bridge.py 2555 2008-12-09 22:44:48Z denniej0 $
"""Combine the output of angle_adk.py into a single data structure."""

from glob import glob
import numpy
import cPickle

datafiles = {'oc': glob('../../data/salt_bridge/oc*deltaRMSD.pickle'),
             'co': glob('../../data/salt_bridge/co*deltaRMSD.pickle'),
             'adp':glob('../../data/salt_bridge/adp*deltaRMSD.pickle'),
             'lig':glob('../../data/salt_bridge/lig*deltaRMSD.pickle'),
             'T400K_co':glob('../../data/salt_bridge/temp_co*deltaRMSD.pickle'),
             'T400K_oc':glob('../../data/salt_bridge/temp_oc*deltaRMSD.pickle'),
             }
outfile = "../../data/salt_bridge/directions_salt_bridge.pickle"


datacounter = dict( [(run_id,0) for run_id in datafiles] )
trajcounter = dict( [(run_id,0) for run_id in datafiles] )
salt_bridge = dict( [(run_id,numpy.array([]).reshape(0,10)) for run_id in datafiles] )

def status(run_id):
    print "[%5s] total salt_bridge found:    %d" % (run_id,len(salt_bridge[run_id]))
    print "[%5s] salt_bridge in files:       %d" % (run_id,datacounter[run_id])
    print "[%5s] total data files:      %d" % (run_id,len(datafiles[run_id]))
    print "[%5s] data files with data:  %d" % (run_id,trajcounter[run_id])


for run_id,files in datafiles.items(): 
    for datafile in files:
        print "[%(run_id)5s] %(datafile)s" % locals(),
        new_salt_bridge = numpy.load(datafile)
        new_count = len(new_salt_bridge)
        datacounter[run_id] += new_count
        print "\tnumber of entries = %d" % new_count,
	print run_id
	#print new_salt_bridge.shape()
	#print salt_bridge.shape()
        if new_count > 0:
            salt_bridge[run_id] = numpy.concatenate( (salt_bridge[run_id], new_salt_bridge[0]) )
            trajcounter[run_id] += 1
        else:
            print "\t**** no data ****",
        print

cPickle.dump(salt_bridge,open(outfile,'wb'),cPickle.HIGHEST_PROTOCOL)

print "\n\nSummary\n"
for run_id in datafiles:
    status(run_id)
print "results:         salt_bridge"
print "saved file:      '%(outfile)s'" % locals()


