#!/usr/bin/env python
# $Id: histogram2D.py 1112 2007-09-11 20:36:40Z denniej0 $
# (c) Oliver Beckstein 2007 <orbeckst@gmail.com>
# GPL

import numpy

class Histogram2D:
    def __init__(self,points,xlimits,ylimits):
        """Bin 2D array z

        h = Histogram2D(points, {'min':2,'max':12,'nbins':10},
                                {'min':-100,'max':100,'nbins':10})

        points   list of points (x,y)

        limits = {min: smallest value, max: largest value, nbins: # of bins}

        **** NO SANITY CHECKS ****
        
        Attributes:

          histogram      array (can be plotted)
          x              x axis
          y              y.axis
          points         input points as numpy array


        Methods:

          plot_points()    plot points and contour lines
          plot_contours()  plot filled contour lines


        Manually plot with contours:

        >>> import pylab
        >>> pylab.plot(h.points,".")
        >>> pylab.contour(x,y,h.histogram.transpose())

        """

        self.points = numpy.array(points)

        xlimits['delta'] = self._delta(xlimits)
        ylimits['delta'] = self._delta(ylimits)

        self.limits = [self._floatify(xlimits),self._floatify(ylimits)]

        # allocate initialised array
        self.histogram = numpy.zeros((xlimits['nbins'],ylimits['nbins']))

        # make the histogram (currently, weight = 1)
        # array: h[row,column]
        # so for plotting (y,x,h)
        
        for p in points:
            try:
                i,j = self._ij(p)         # i(x), j(y) 
                self.histogram[i,j] += 1 
            except TypeError:
                # out of limits (None returned)
                pass

        # make axes
        a = []
        for axis in [0,1]:
            a.append(
                [self.limits[axis]['min'] + n*self.limits[axis]['delta'] \
                 for n in range(0,int(self.limits[axis]['nbins'])) ]
                )
        self.x = numpy.array(a[0])
        self.y = numpy.array(a[1])

        del a

    def _floatify(self,d):
        """only works for dicts with numbers as values"""
        for k in d:
            d[k] = float(d[k])
        return d
        
    def _delta(self,limits):
        return (limits['max'] - limits['min'])/float(limits['nbins'])



        
    def _ij(self,point):
        """return index (i,j) in 2D histogram"""

        x,y = point[0],point[1]
        i,j = self._idx(x,0), self._idx(y,1)

        #print "_ij(): point=(%(x)f,%(y)f) i=%(i)d j=%(j)d" % locals()

        # discard points outside limits (or could use exceptions)
        if i<0 or i >= self.limits[0]['nbins'] or \
               j<0 or j >= self.limits[1]['nbins']:
            return None

        return i,j


    def _idx(self,val,axis):
        """find 1D index using the appropriate axis (0 = x, 1= y)"""
        return numpy.floor(
            (val - self.limits[axis]['min']) / self.limits[axis]['delta'])

    def plot_contours(self):
        """filled contour plot with lines"""
        try:
            from pylab import contour,contourf
        except ImportError:
            print "pylab is required to plot"
            return None
        contour(self.x,self.y,self.histogram.transpose())
        contourf(self.x,self.y,self.histogram.transpose())

    def plot_points(self):
        """plot input points with contour lines"""
        try:
            from pylab import plot,contour
        except ImportError:
            print "pylab is required to plot"
            return None
        
        x,y = self.points.transpose()
        plot(x,y,".")
        contour(self.x,self.y,self.histogram.transpose())        

    def __repr__(self):
        """Histogram representation"""
        return "<2D Histogram "+str(self.limits[0]["nbins"])+"x"+\
               str(self.limits[1]["nbins"])+" bins from "+\
               str(len(self.points))+" data points>"

    def __getitem__(self,i):
        """Act as the histogram array"""
        return self.histogram[i]

class WHistogram2D(Histogram2D):
    def __init__(self,points,xlimits,ylimits):
        """Input (x,y,z) and weight with z

        h = Histogram2D(points, {'min':2,'max':12,'nbins':10},
                                {'min':-100,'max':100,'nbins':10})

        points   list of points (x,y,z)

        limits = {min: smallest value, max: largest value, nbins: # of bins}

        **** NO SANITY CHECKS ****
        
        Attributes:

          histogram      array (can be plotted)
          x              x axis
          y              y.axis
          points         input points as numpy array

         whistogram      histogram wheighted with z value
         avg             averaged whistogram


        Methods:

          plot_points()    plot points and contour lines
          plot_contours()  plot filled contour lines


        Manually plot with contours:

        >>> import pylab
        >>> pylab.plot(h.points,".")
        >>> pylab.contour(x,y,h.histogram.transpose())

        """

        self.points = numpy.array(points)

        xlimits['delta'] = self._delta(xlimits)
        ylimits['delta'] = self._delta(ylimits)

        self.limits = [self._floatify(xlimits),self._floatify(ylimits)]

        # allocate initialised array
        self.histogram = numpy.zeros((xlimits['nbins'],ylimits['nbins']))
        self.whistogram = numpy.zeros((xlimits['nbins'],ylimits['nbins']))

        # make the histogram (currently, weight = 1)
        # array: h[row,column]
        # so for plotting (y,x,h)
        
        for p in points:
            try:
                i,j = self._ij(p)         # i(x), j(y) 
                self.histogram[i,j] += 1
		''' liz added code'''
		if self.whistogram[i,j] == 0:
                        self.whistogram[i,j] = p[2]
		elif self.whistogram[i,j] > p[2]:
			self.whistogram[i,j] = p[2]
		else: 
			continue
            except TypeError:
                # out of limits (None returned)
                pass

        # make axes
        a = []
        for axis in [0,1]:
            a.append(
                [self.limits[axis]['min'] + n*self.limits[axis]['delta'] \
                 for n in range(0,int(self.limits[axis]['nbins'])) ]
                )
        self.x = numpy.array(a[0])
        self.y = numpy.array(a[1])

        del a
        self.avg = self.whistogram/(self.histogram + numpy.where(self.histogram==0,1.0,0.0))



def _test(npoints=1000):
    """test Histogram2D by binning npoints random points which are
    Gaussian dsitributed"""

    xmin, xmax = -20, 20
    ymin, ymax = -3, 5

    from random import gauss
    from pylab import plot,contour,contourf

    xmu, xsigma = 0, 6.0
    ymu, ysigma = 0, 1.0

    points = [(gauss(xmu,xsigma), gauss(ymu,ysigma)) for i in range(0,npoints)]

    h = Histogram2D(points,xlimits={'min':xmin,'max':xmax,'nbins':20},
                   ylimits={'min':ymin,'max':ymax,'nbins':10})

    try:
        import pylab
        h.plot_points()
        h.plot_contours()        
        pylab.show
    except:
        pass
    
    return h
