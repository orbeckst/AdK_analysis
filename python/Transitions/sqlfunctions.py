# $Id: sqlfunctions.py 3379 2009-04-21 17:57:29Z oliver $
"""SQL functions to be added to a SQLite database

Example:

  from sqlfunctions import *
  self.connection.create_function("sqrt", 1, _sqrt)
  self.connection.create_function("fformat",2,_fformat)
  self.connection.create_aggregate("std",1,_Stdev)
  self.connection.create_aggregate("median",1,_Median)
  self.connection.create_aggregate("array",1,_NumpyArray)
  self.connection.create_aggregate("histogram",4,_NumpyHistogram)
  self.connection.create_aggregate("distribution",4,_NormedNumpyHistogram)
  self.connection.create_aggregate("meanhistogram",5,_MeanHistogram)
  self.connection.create_aggregate("stdhistogram",5,_StdHistogram)
  self.connection.create_aggregate("minhistogram",5,_MinHistogram)
  self.connection.create_aggregate("maxhistogram",5,_MaxHistogram)
  self.connection.create_aggregate("medianhistogram",5,_MedianHistogram)
  self.connection.create_aggregate("zscorehistogram",5,_ZscoreHistogram)

"""
import numpy
# compatibility check: we NEED consistent 1d histogram functions: we
# decided to use numpy 1.x style, which returns edges, NOT lower bin edges
_numpyversion = map(int, numpy.version.version.split('.'))
if _numpyversion[0] < 1:
    raise ImportError('Need at least numpy 1.x, only have %r' % numpy.version.version)
if _numpyversion[1] < 1:
    # we want a histogram that returns edges
    def histogram1d(*args,**kwargs):
        _range = kwargs.pop('range',None)
        if not _range is None:
            kwargs['range'] = (_range,)   # needs to be a sequence
        h,e = numpy.histogramdd(*args,**kwargs)
        return h,e[0]
    histogram1d.__doc__ = "1D histogram, based on numpy histogramdd; returns edges as in numpy 1.1.x\n"+\
                        numpy.histogram.__doc__    
else:
    # once deprecation for new=True sets in we can catch this here
    def histogram1d(*args,**kwargs):
        kwargs['new'] = True
        h,e = numpy.histogram(*args,**kwargs)
        return h,e
    histogram1d.__doc__ = numpy.histogram.__doc__
  

from sqlutil import adapt_numpyarray, convert_numpyarray,\
    adapt_object, convert_object


def _sqrt(x):
    try:
        x = float(x)
    except TypeError:
        return None
    return numpy.sqrt(x)        

def _fformat(format,x):
    return format % x

class _Stdev(object):
    """Implement standard deviation of the sample as SQL aggregate function.
    (Uses N-1 variance.)
    Do it in one pass (see eg
    http://smallcode.weblogs.us/2006/11/27/calculate-standard-deviation-in-one-pass/
    though we may run in an underflow by calculating N/N-1<X^2-<X>^2>.).

    Also, we don't check if our arguments are valid as numbers.
    """
    def __init__(self):
        self.x2 = 0
        self.x = 0
        self.n = 0
    def step(self,x):
        try:
            x = float(x)
            self.x2 += x*x
            self.x  += x
            self.n  += 1                    
        except TypeError:
            pass        # don't contribute to average
    def finalize(self):
        if self.n<2: return 0.0
        return numpy.sqrt((self.n*self.x2 - self.x*self.x)/(self.n*(self.n-1)))

class _Median(object):
    def __init__(self):
        self.data = []
    def step(self,x):
        try:
            x = float(x)
            self.data.append(x)
        except TypeError:
            pass        # don't contribute
    def finalize(self):
        return numpy.median(self.data)

class _NumpyArray(object):
    def __init__(self):
        self.data = []
    def step(self,x):
        self.data.append(x)
    def finalize(self):
        return adapt_numpyarray(numpy.array(self.data))

class _NumpyHistogram(object):
    def __init__(self):
        self.is_initialized = False
        self.data = []
    def step(self,x,bins,xmin,xmax):
        if not self.is_initialized:
            self.bins = bins
            self.range = (xmin,xmax)
            self.is_initialized = True
        self.data.append(x)
    def finalize(self):
        hist,edges = histogram1d(self.data,bins=self.bins,range=self.range,
                                 normed=False)
        return adapt_object((hist,edges))

class _NormedNumpyHistogram(_NumpyHistogram):
    def finalize(self):
        hist,edges = histogram1d(self.data,bins=self.bins,range=self.range,
                                 normed=True)
        return adapt_object((hist,edges))

class _FunctionHistogram(_NumpyHistogram):
    """Baseclass for histogrammed functions.

    A histogrammed function is created by applying a function
    to all values y that have been accumulated in a bin x.
    """
    def __init__(self):
        _NumpyHistogram.__init__(self)
        self.y = []
    def step(self,x,y,bins,xmin,xmax):
        _NumpyHistogram.step(self,x,bins,xmin,xmax)
        self.y.append(y)
    def finalize(self):
        raise NotImplementedError("_FunctionHistogram must be inherited from.")
        # return adapt_object( (...,...,...) )

class _MeanHistogram(_FunctionHistogram):
    """Mean of the weights in each bin. 
    Takes TWO column arguments: value and weight"""
    def finalize(self):
        return adapt_object(regularized_function(\
                self.data,self.y,numpy.mean,bins=self.bins,range=self.range))

class _StdHistogram(_FunctionHistogram):
    """Standard deviation of the weights in each bin. 
    Takes TWO column arguments: value and weight"""
    def finalize(self):
        return adapt_object(regularized_function(\
                self.data,self.y,numpy.std,bins=self.bins,range=self.range))

class _MinHistogram(_FunctionHistogram):
    """Min value of the weights in each bin. 
    Takes TWO column arguments: value and weight"""
    def _min(self,v):
        try:
            return numpy.min(v)
        except ValueError:  # empty array
            return numpy.nan

    def finalize(self):
        return adapt_object(regularized_function(\
                self.data,self.y,self._min,bins=self.bins,range=self.range))

class _MaxHistogram(_FunctionHistogram):
    """Max value of the weights in each bin. 
    Takes TWO column arguments: value and weight"""
    def _max(self,v):
        try:
            return numpy.max(v)
        except ValueError:  # empty array
            return numpy.nan

    def finalize(self):
        return adapt_object(regularized_function(\
                self.data,self.y,self._max,bins=self.bins,range=self.range))

class _MedianHistogram(_FunctionHistogram):
    """Median value of the weights in each bin. 
    Takes TWO column arguments: value and weight"""
    def finalize(self):
        return adapt_object(regularized_function(\
                self.data,self.y,numpy.median,bins=self.bins,range=self.range))

class _ZscoreHistogram(_FunctionHistogram):
    """Z-score of the weights in each bin abs(Y - <Y>)/std(Y). 
    Takes TWO column arguments: value and weight"""
    def Zscore(self,v):
        m = v.mean()
        s = v.std()
        return numpy.nan_to_num( numpy.mean(numpy.abs(v - m))/s )

    def finalize(self):
        return adapt_object(\
            regularized_function(self.data,self.y,self.Zscore,bins=self.bins,range=self.range))


# Helper functions

def regularized_function(x,y,func,bins=None,range=None):
    """Compute func() over data aggregated in bins.

    (x,y) --> (x', func(Y'))  with Y' = {y: y(x) where x in x' bin}

    First the data is collected in bins x' along x and then func is applied to
    all data points Y' that have been collected in the bin.

    :Arguments:
    x              abscissa values (for binning)
    y              ordinate values (func is applied)
    func           a numpy ufunc that takes one argument, func(Y')
    bins           number or array
    range          limits (used with number of bins)

    :Returns:
    F,edges        function and edges (midpoints = 0.5*(edges[:-1]+edges[1:]))
    """
    _x = numpy.asarray(x)
    _y = numpy.asarray(y)

    # setup of bins from numpy.histogram
    if (range is not None):
        mn, mx = range
        if (mn > mx):
            raise AttributeError('max must be larger than min in range parameter.')

    if not numpy.iterable(bins):
        if range is None:
            range = (_x.min(), _x.max())
        mn, mx = [float(mi) for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = numpy.linspace(mn, mx, bins+1, endpoint=True)
    else:
        bins = numpy.asarray(bins)
        if (numpy.diff(bins) < 0).any():
            raise AttributeError('bins must increase monotonically.')

    sorting_index = numpy.argsort(_x)
    sx = _x[sorting_index]
    sy = _y[sorting_index]

    # boundaries in SORTED data that demarcate bins; position in bin_index is the bin number
    bin_index = numpy.r_[sx.searchsorted(bins[:-1], 'left'),
                         sx.searchsorted(bins[-1], 'right')]

    # naive implementation: apply operator to each chunk = sy[start:stop] separately
    #
    # It's not clear to me how one could effectively block this procedure (cf
    # block = 65536 in numpy.histogram) because there does not seem to be a
    # general way to combine the chunks for different blocks, just think of
    # func=median
    F = numpy.zeros(len(bins)-1)  # final function
    F[:] = [func(sy[start:stop]) for start,stop in zip(bin_index[:-1],bin_index[1:])]
    return F,bins
