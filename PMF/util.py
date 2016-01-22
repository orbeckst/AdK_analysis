# $Id: util.py 3001 2009-02-21 02:54:38Z oliver $
"""Unsorted helper functions and examples."""

import numpy
import os.path

def filename(name,ext=None,keep=False):
    """Return a new name that has suffix attached; replaces other extensions.

    name        filename; extension is replaced unless keep=True
    ext         extension
    keep        False: replace existing extension; True: keep if exists
    """ 
    name = str(name)
    if ext is None:
        return name
    if not ext.startswith(os.path.extsep):
        ext = os.path.extsep + ext
    #if name.find(ext) > 0:    # normally >= 0 but if we start with '.' then we keep it
    #    return name
    root, origext = os.path.splitext(name)
    if len(origext) == 0 or not keep:
        return root + ext
    return name


def endpoints(db):
    """Returns average endpoints of trajectories (actually, min/max)

    Plot on top of other graphs with

    plot(avgends[:,0], avgends[:,1], "wo", ms=10, alpha=0.8)
    xlim(30,90)
    ylim(80,150)    
    """

    xx = db.selection("""SELECT filename, MIN(NMP) AS Nclosed, MIN(LID) AS Lclosed,
                                          MAX(NMP) AS Nopen,   MAX(LID) AS Lopen
                           FROM __self__
                          WHERE GLOB("*/co_*",filename) GROUP BY filename""")

    endpoints = xx.recarray

    avgends = numpy.array(((endpoints.Nclosed.mean(), endpoints.Lclosed.mean()),
                           (endpoints.Nopen.mean(), endpoints.Lopen.mean())))

    # average closed and open state coordinates from the dims trajectories
    return avgends
    

def around(x,y,d):
    """Return coodinates of the 8 neighbours of (x,y) at lattice constant d."""
    return numpy.array([(x+dx, y+dy) for dx in (-d,0,d)
                        for dy in (-d,0,d)]).reshape((9,2))

def print_around(x,y,d,istart=1):
    for i,(_x,_y) in enumerate(around(x,y,d)):
        i += istart
        #123456789.12345
        #N   LID NMp   ; Liz uses pos 5-7 and 9-10 in pmf_bat.sge
        #  1              48       118  ; this script (hopefully makes pmf_bat fail)
        print "%(i)3d\t\t %(_x)g\t  %(_y)g" % vars()


class ColorRing(object):
    def __init__(self):
        self.count = 0
        self.colors = ['black','magenta','cyan','green','yellow','orange','red',]
    def get(self):
        i = self.count % len(self.colors)
        self.count += 1
        return self.colors[i]

