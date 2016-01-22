# $Id: pathfinder.py 2479 2008-11-19 01:41:08Z oliver $
# (c) 2008 Oliver Beckstein <oliver.beckstein@bioch.ox.ac.uk>
"""Read a PMF landscape. Store as a weighted graph with DeltaW along
the edges (i.e. we are really storing the gradient landscape). Find
Dijkstra least-cost path through the weighted graph (i.e. integrating
the gradient along the path).

ARGH... this is stupid: this will just get the free energy difference
between initial and final point when we take the gradient.

We really need to cheat and zero any downhill steps. So:

* connect in both directions (DiGraph)

* zero all gradients < 0.

Use this interactively from ipython.

  import PMF
  P = PMF.Landscape(datafile)
  P.show()
  P.minpath( (1,1), (12,34) )

Datafile format:

# comment lines (ignored)
# x,y: indices (1-based)
# W: PMF at these indices
x y W(x,y)
...

"""

import networkx as NX
import numpy

class Landscape(object):
    # start in one corner and progress in +1,+1 direction
    # (trivial in 2D, for 3D a list approach comes in handy)
    # This is a graph which is connected in both directions.
    _shifts = [numpy.array(x) for x in ([
                [0,+1],[+1,0],[-1,0],[0,-1],      # first shell
                [+1,+1],[+1,-1],[-1,+1],[-1,-1],  # second shell
                ])]

    def __init__(self,PMFdatafile):
        """Initialize the graph:
        
        P = Landscape('/Users/oliver/tmp/PMF.dat')

        """

        # setup graph
        self.graph = NX.XDiGraph()

        PMF = open(PMFdatafile,'r')

        # collect as a numpy array;
        # have to be careful because points maybe missing
        pointrecords = []
        for line in PMF:
            if line[0] == '#':
                continue
            x,y,W = line.strip().split()
            x = int(x)
            y = int(y)
            pointrecords.append((x,y,float(W)))
        PMF.close()
        # recarray for processing
        points = numpy.rec.fromrecords(pointrecords,names="x,y,W")
        del pointrecords  # not needed anymore, allow garbage collection
        

        # dimensions
        # (pad on upper and lower end; this is why we want 1-based indices)
        xmin,xmax = points.x.min() - 1, points.x.max() + 1
        ymin,ymax = points.y.min() - 1, points.y.max() + 1
        if xmin < 0 or ymin < 0:
            raise ValueError("Indices in data file must start at 1.")
        
        # final data structure (could used masked arrays but let's just use
        # nan for unknown points and screen later ourself)
        self.pmf_array = numpy.NaN * numpy.ones((xmax-xmin+1,ymax-ymin+1),dtype=float)

        # fill array
        self.pmf_array[points.x,points.y] = points.W
        self.pmf_masked = numpy.ma.MaskedArray(self.pmf_array,
                                               mask=numpy.isnan(self.pmf_array),
                                               fill_value=None)

        # neighbours on grid:
        # clunky, there's certainly a nicer way to do this...
        for x in xrange(xmax):
            for y in xrange(ymax):
                currentpos = numpy.array([x,y])        
                currentW = self.pmf_array[x,y]
                if numpy.isnan(currentW):
                    continue
                currentnode = tuple(currentpos)
                neighbours = [tuple(currentpos + shift) for shift in self._shifts]
                for n in neighbours:
                    W = self.pmf_array[n]
                    #print "W(%d, %d) = %r" % (n[0],n[1],W)
                    if numpy.isnan(W):
                        continue
                    DeltaW = W - currentW
                    if DeltaW < 0:
                        DeltaW = 0        # only count uphill steps
                    node = tuple(n)
                    edge = (currentnode, node, DeltaW)
                    self.graph.add_edge(edge)      # adds nodes and edge
                    #print "Added edge (%(currentnode)r, %(node)r, %(DeltaW)f)" % vars()
                    

    def minpath(self,start,end):
        """Return Dijkstra path on graph between start=(x1,y1) and end=(x2,y2)."""
        score,path = NX.path.bidirectional_dijkstra(self.graph,start,end)
        self.last_path = {'score':score,'path':path}
        return score,path
    
    def plot_path(self,start=None,end=None,**plotargs):
        """Plot last path on top of the current graph, or with the
        optional arguments, calculate a new graph."""

        import pylab
        if start is not None and end is not None:
            self.minpath(start,end)
        x,y = numpy.array(self.last_path['path']).T
        pylab.plot(x,y,**plotargs)

    def plot(self,**kwargs):
        """plot landscape (kwargs are passed on to imshow()

        Use interpolation='bilinear' or 'bicubic' for a smooth
        surface. Default is 'nearest', which shows exact bin
        boundaries.
        """
        import pylab

        kwargs.setdefault('interpolation','nearest')
        pylab.clf()
        pylab.xlabel('x')
        pylab.ylabel('y')
        pylab.imshow(self.pmf_masked.T,**kwargs)
        pylab.colorbar()
        pylab.show()

	
