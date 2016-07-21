#!/usr/bin/env python
# $Id: combine_angles.py 3014 2009-02-23 18:17:50Z denniej0 $
"""Combine the output of angle_adk.py into a single data structure."""

from glob import glob
import numpy
import cPickle

datafiles = {'oc': glob('../../data/new_angles/DIMS/apo/standard/oc*.pickle'),
             'co': glob('../../data/new_angles/DIMS/apo/standard/co*.pickle'),
             'adp':glob('../../data/new_angles/adp*.pickle'),
             'lig':glob('../../data/new_angles/lig*.pickle'),
             'T400K_co':glob('../../data/new_angles/temp_co*.pickle'),
             'T400K_oc':glob('../../data/new_angles/temp_oc*.pickle'),
             }
outfile = "../../data/new_angles/DIMS/apo/standard/directions_angles.pickle"


datacounter = dict( [(run_id,0) for run_id in datafiles] )
trajcounter = dict( [(run_id,0) for run_id in datafiles] )
angles = dict( [(run_id,numpy.array([]).reshape(0,2)) for run_id in datafiles] )

def status(run_id):
    print "[%5s] total angles found:    %d" % (run_id,len(angles[run_id]))
    print "[%5s] angles in files:       %d" % (run_id,datacounter[run_id])
    print "[%5s] total data files:      %d" % (run_id,len(datafiles[run_id]))
    print "[%5s] data files with data:  %d" % (run_id,trajcounter[run_id])


for run_id,files in datafiles.items(): 
    for datafile in files:
        print "[%(run_id)5s] %(datafile)s" % locals(),
        new_angles = numpy.load(datafile)
        new_count = len(new_angles)
        datacounter[run_id] += new_count
        print "\tnumber of entries = %d" % new_count,
        if new_count > 0:
            angles[run_id] = numpy.concatenate( (angles[run_id], new_angles) )
            trajcounter[run_id] += 1
        else:
            print "\t**** no data ****",
        print

cPickle.dump(angles,open(outfile,'wb'),cPickle.HIGHEST_PROTOCOL)

print "\n\nSummary\n"
for run_id in datafiles:
    status(run_id)
print "results:         angles"
print "saved file:      '%(outfile)s'" % locals()


