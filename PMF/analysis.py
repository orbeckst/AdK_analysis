# $Id: analysis.py 4107 2010-04-27 12:15:54Z root $
"""Analysis of observables in the database.

Some observables are directly entered into columns in PMF.angles.AngleDB.

This module provides for combining these observables, and creates new
tables that hold computed values.

The idea is that these derived tables can be used with the other
plotting and selection routines (maybe using SQL VIEWs and JOINs?)
"""


import os,errno
import warnings
import numpy
import config, angles


# WARNING: all observables are using the __data__ view instead of the
# __self__ table; if this is slow, use 'SELECT FROM __self__ WHERE TID
# NOTNULL' directly

class ObservableCreator(angles.Selection):
    observable_type = None   # descriptive string (eg 'FRET'), is used as a dirname
    qualifier = None         # additional descriptor (eg 'Kern'), is used in filenames
    init_SQL = None          # 'SELECT ... FROM __data__ ...' to create the observable
    base_columns = ('trajectory','frame','LID','NMP','oRMSD','cRMSD','DeltaRMSD')
    observables = None       # tuple containing the names of observable columns/names
    
    def __init__(self,db,**kwargs):        
        self.name = str(self.observable_type) + '_' + str(self.qualifier)
        if self.init_SQL is None:
            raise NotImplementedError
        super(ObservableCreator,self).__init__(db,sql=self.init_SQL)

    def __repr__(self):
        return "<%s observable_type=%r qualifier=%r>" % (self.__class__.__name__, 
                                                         self.observable_type,self.qualifier)

class SaltbridgeIndicatorCreator(ObservableCreator):
    """A saltbridge is said to exist if
    (1) the distance is smaller than a cutoff distance, and
    (2) the interaction energy is greather than a energy threshold.

        XXX: not sure how to define this: what is our reference, ie
        what is the energy of an intact salt bridge (this will depend
        on the nonbonded parameters!)? Or should we simply do a
        Schmitt-trigger to detect transitions between 'high' and 'low'
        (i.e. normalize the data to some scale and compare to a fixed
        switching distance).

    NOTE: (2) is not being used at the moment

    ATTENTION: The creator requires at the DeltaRMSD and the
               saltbridge data to be in the db.
    
        
    """
    observable_type = 'saltbridge'
    qualifier = 'indicator'
    # TODO: this should be configured in observables
    #       MUST BE UPDATED when new saltbridges are defined
    observables = ('ISB1','ISB2','ISB3','ISB4','ISB5',
                   'ISB6','ISB7','ISB8','ISB9','ISB10',
                   'ISB11','ISB12','ISB13','ISB14')  # column names
    
    def __init__(self,db,**kwargs):
        """Setup the SaltbridgeIndicator observables.

          SI = SaltbridgeIndicatorCreator(db, maxdistance=6.0)

        When this class is instantitated it computes the additional
        columns that are required for saltbridge indicator
        analysis.

        It returns a selection (a database-like object) that must be
        used for the saltbridge indicator analysis .... it's pretty
        crappy, I admit it.

        Arguments:

        maxdistance   maximum distance between amino/guanidino group and
                      carboxy group to be still considered a salt bridge [6 A]
        minenergy     minimum stabilization energy to be counted as a bridge
                      [NOT IMPLEMENTED]
        
        """
        self.maxdistance = kwargs.pop('maxdistance', 6.0)  # Angstroem
        self.minenergy = kwargs.pop('minenergy', 10.0)     # kcal/mol <-- ????
        maxD = self.maxdistance
        minE = self.minenergy
        columns = list(self.base_columns[:])  # copy so that the class variable is not changed
        # TODO: this should be less fragile ... it was a quick hack
        n_observables = len(self.observables)
        for n in xrange(1,n_observables+1):
            SB = 'SB%d' % n    # distance column
            ESB = 'ESB%d' % n  # energy column XXX
            ISB = 'ISB%d' % n
            #columns.append = '%(SB)s < %(maxD)f AND %(ESB)s > %(minE)s AS %(ISB)d' % vars()
            # ignore energy
            columns.append('%(SB)s < %(maxD)f AS %(ISB)s' % vars())
        selection = ", ".join(columns)
        # observe blacklisting, ie only use trajectories that are
        # listed in __meta__ by using __data__ instead of __self__        
        self.init_SQL = """SELECT %(selection)s FROM __data__ 
                           WHERE DeltaRMSD NOTNULL AND SB1 NOTNULL""" % vars()
        super(SaltbridgeIndicatorCreator,self).__init__(db,**kwargs)

        if len(self) == 0:
            warnings.warn("No ISB values computed --- did you load DeltaRMSD and saltbridges?")


class Observable(object):
    """A 1D observable y(x), a function of two observables X and Y.

    NOTE: In principle one needs to weigh each data point with the
    probability from the PMF. This is NOT implemented at the
    moment. See AngleprojectedObservable.
    """
    
    observable_type = None   # descriptive string (eg 'FRET'), is used as a dirname
    qualifier = None         # additional descriptor (eg 'Kern'), is used in filenames
    observable_column = None # column in __data__ (eg 'FRET_Kern')
                             # NOTE: if empty we use the class name!!!
    x_column = None          # column in __data__ that is used as the abscissa
    default_binwidth = 2.0

    # plots: default distance between major contour lines (minor are at 0.5*major)
    major_contour = 1.0
    
    # additional data to be written in legends/titles for the observable_type
    legend = '_nolegend_'
    xlabel = None

    # depending on plotting mode, construct different plot titles and filenames
    _title = {'observable': r'$d\ (\AA)$',        # is ok for 'FRET' and 'saltbridge'
              'stdev':      r'$\sigma_{d}\ (\AA)$',
              'counts':     r'counts $N$',
              }
              
    _suffix = {'observable': '',
               'stdev':      '_stdev',
               'counts':     '_counts',
               }

    def __init__(self,db,
                 nequilibration=0,binwidth=None,Nmin=10,
                 padding=0,blocksize=2**16,edges=None,**limits):
        """Analysis object for a 1D observable stored in the db or Selection.

        db          database or Selection with a __data__ table that contains the data
        edges       can be used instead of binwidth and limits to construct the histogram
        binwidth    bin size
        Nmin        minimum number of counts in a bin so that the bin
                    is included in the final histogram (must be at
                    least 1)
        **limits    limits xmin and xmax for the histogram [automatic]

        Notes:
        * nequilibration is ignored,
        * only padding == 0 works.
        """

        if self.observable_type is None and self.observable_column is None:
            raise NotImplementedError('This class should not be called itself but must be derived from.')
        if self.observable_column is None:
            # assume that the column is the class name
            self.observable_column = self.__class__.__name__

        self.name = str(self.observable_type) + '_' + str(self.qualifier)
            
        self.db = db        
        self.Nmin = Nmin
        assert Nmin > 1
        
        # NOTE: to make this class fairly universal simply we use the variable
        # column trick and alias the SQL output column with
        # '%(observable_column)s AS observable'; then use 'all.observable' in the
        # blockiterator further below
        self.selection = db.selection("""SELECT %s AS x,%s AS observable FROM __data__
                                         WHERE observable NOTNULL""" %
                                      (self.x_column,self.observable_column))

        self.nequilibration = nequilibration   # ignore --- should we honour this?
        self.padding = padding
        self.edges = edges
        self.binwidth = binwidth or self.default_binwidth

        # setting bins for the histogram
        if self.edges is not None:
            # overrides bins and limits arguments
            nbins = len(e)-1  # x
            self.limits = {'xmin': numpy.min(e), 'xmax': numpy.max(e)}
            self._histbinargs = {'bins': self.edges}
        else:            
            # limits
            [(xmin,xmax)] = self.selection.SQL('SELECT min(x), max(x) FROM __data__')
            limits.setdefault('xmin',xmin)  # use data limits unless something else
            limits.setdefault('xmax',xmax)  # was provided in **limits
            self.limits = limits
            xmin = self.limits['xmin']      # xmin and xmax to be used
            xmax = self.limits['xmax']

            nbins = int(numpy.ceil((xmax-xmin)/self.binwidth))
            x_range = [xmin, xmax]
            self._histbinargs = {'bins': nbins, 'range': x_range}                    
        # in order to reconstruct the exact histogram use **self._histbinargs

                                      
        # calculate observable
        # Do the histogramming in blocks for big selections because otherwise the memory blows up.
        h = numpy.zeros(nbins,dtype=numpy.float_) # final binned average
        h_i = h.copy()                            # tmp array for each block
        s = h.copy()                              # standard deviation
        s_i = s.copy()                            # tmp array for each block
        N = numpy.zeros(nbins,dtype=numpy.int_)   # total counts per bin
        N_i = N.copy()                            # tmp array for each block

        def histogram(*args,**kwargs):
            """New numpy.histogram with correct bins."""
            kwargs['normed'] = False   # never normalize for bin-averaging
            kwargs['new'] = True       # currently needed for proper edges and weights
            kwargs.update(self._histbinargs)
            return numpy.histogram(*args,**kwargs)
        
        # empty run to get edges e in case there are no data and the loop is skipped
        N_i[:],e = histogram([0])

        # TODO: The following the simplest aggregate statistic: the mean for one
        # bin. Its advantage is that it can be calculated in a single pass,
        # together with teh standard deviation. For more complicated ones we
        # would need to collect the distribution 'perpendicular' to the
        # abscissa. Technically that mean obtaining a histogram in the ordinate
        # for each data point. I could implement this as a 2D histogram.
        
        for all in self.selection.block_iterator(blocksize=blocksize,asrecarray=True):
            x, observable = all.x, all.observable
            del all                 # release memory (?)
            # number of values per bin
            N_i[:],e = histogram(x)
            # sum of values per bin
            h_i[:],e = histogram(x, weights=observable)
            # calculate variance in one pass: <(X - <X>)**2> = <X**2> - <X>**2
            s_i[:],e = histogram(x, weights=observable**2)
            
            del x  # release memory immediately

            h += h_i    # accumulate bin sums ...
            s += s_i    # ... sum of squares
            N += N_i    # ... and bin counts
            
        del h_i
        del s_i
        del N_i

        bins_with_data = (N > 0)                # only use bins that contain data
        h[bins_with_data] /= N[bins_with_data]  # average

        bins_with_data = (N > Nmin)             # N-1 standard deviation
        s[N <= Nmin] = 0                        # no values for N=0,1 data points (or mask array?)
        s[bins_with_data] /= N[bins_with_data]
        s[bins_with_data] -= h[bins_with_data]**2                         # VarX = <X**2> - <X>**2
        s[bins_with_data] *= N[bins_with_data] / (N[bins_with_data] - 1)  # N/(N-1) * VarX
        s[bins_with_data] = numpy.sqrt(s[bins_with_data])                 # final stddev

        self.observable       = numpy.ma.array(h, mask=numpy.logical_not(bins_with_data))
        self.observable_stdev = numpy.ma.array(s, mask=numpy.logical_not(bins_with_data))
        self.bincounts    = N

        # midpoints of bins
        self.xvalues = 0.5*(e[1:]+e[:-1])        

    def plot(self,**kwargs):
        """Plot observable with standard deviation.

          plot(**kwargs)

        Arguments are passed on to pylab.plot with the exception of
        'title'. The canvas is not cleared so one can plot multiple graphs on
        top of each other. The legend box can be shown with a final

          pylab.legend()

        command.
        
        :Arguments:
        title           title of the graph
        *               see pylab.plot
        """
        import pylab

        title = kwargs.pop('title',None)
        if title is not None:
            pylab.title(title)

        kwargs.setdefault('label',self.legend)
        kwargs.setdefault('linewidth',3)
        P, = pylab.plot(self.xvalues, self.observable, "-", **kwargs)

        # manual errorbars markers
        kwargs.setdefault('color',P.get_color())  # use same color as plot
        kwargs.pop('label')
        kwargs['linestyle'] = '--'
        kwargs['linewidth'] = 1.0
        kwargs.setdefault('markersize',3)
        kwargs['markeredgecolor'] = kwargs['color']
        kwargs['marker'] = 'v'
        kwargs.setdefault('alpha',0.5)
        pylab.plot(self.xvalues, self.observable + 0.5*self.observable_stdev,
                   **kwargs)
        kwargs['marker'] = '^'
        pylab.plot(self.xvalues, self.observable - 0.5*self.observable_stdev,
                   **kwargs)
        pylab.xlabel(self.xlabel)

    def plot_all(self,**kwargs):
        raise NotImplementedError

    def _auto_plots(self,mode,filebasename,figdir,plotargs):
        """Generate standard plots and write png and and pdf. Chooses filename and plot title."""
        import pylab

        try:
            os.makedirs(figdir)
        except OSError,err:
            if err.errno != errno.EEXIST:
                raise

        def figs(*args):
            return os.path.join(figdir,*args)
        
        modefilebasename = filebasename + self._suffix[mode]
        _plotargs = plotargs.copy()  # need a copy because of changing 'title'
        if plotargs.get('title') is None:  # None --> set automatic title
            _plotargs['title'] = self._title[mode]+' '+self.legend

        pylab.clf()
        self.plot(**_plotargs)
        pylab.savefig(figs(modefilebasename + '.png'))   # png
        pylab.savefig(figs(modefilebasename + '.pdf'))   # pdf
        
        print "--- Plotted %(modefilebasename)r (png,pdf)." % vars()

    def filebasename(self,runtype=None):
        """Create a filename that contains the RUNTYPE.

        filebasename(runtype=RUNTYPE)

        filebasename = observable_type + _ + RUNTYPE

        RUNTYPE is eg 'dims' or 'dims_co'.

        None      take the runtype from the observable
        ""        empty string leave the filename unchanged

        """
        if runtype is None:
            try:
                runtype = self.runtype
            except AttributeError:
                pass
        if runtype:
            runtype = '_' + runtype
        runtype = runtype or ""
        return self.observable_type + runtype
    


# NOTE: each observable requires a single column in the AngleDB
# TODO: change this to individual tables for each observable_type

class AngleprojectedObservable(object):
    """A 2D observable y(NMP,LID): Y is projected on the NMP-angle and LID-angle plane.

    The Observable is calculated by binning it on a NMP,LID grid and averaging
    over each bin. The correct bin average weights each data point by the
    probability from the PMF (FreeEnergy.P) and normalizes by the sum of
    probabilities. If one averages on the same bins as the PMF then one can
    ignore the weighting and it is probably still a decent approximation to
    ignore it in other cases.

    Derive a class from Observable and override definitions and
    legends.  Multiple inheritance should work with this one.
    """
    # AngleprojectedObservable class copied & pasted from wham.Energy (with modifications)
    
    # which definition of Observable to select from the db:
    # key <-- name; SQL column <--observable_column
    # Graphs are saved in figs/<name>/<qualifier> or figs/<name> if
    # qualifier == None
    observable_type = None   # descriptive string (eg 'FRET'), is used as a dirname
    qualifier = None         # additional descriptor (eg 'Kern'), is used in filenames
    observable_column = None # column in __data__ (eg 'FRET_Kern')
                             # NOTE: if empty we use the class name!!!

    # plots: default distance between major contour lines (minor are at 0.5*major)
    major_contour = 1.0
    
    # additional data to be written in legends/titles for the observable_type
    legend = None

    # depending on plotting mode, construct different plot titles and filenames
    _title = {'observable': r'$d\ (\AA)$',        # is ok for 'FRET' and 'saltbridge'
              'stdev':      r'$\sigma_{d}\ (\AA)$',
              'counts':     r'counts $N$',
              }
              
    _suffix = {'observable': '',
               'stdev':      '_stdev',
               'counts':     '_counts',
               }

    # default values for keyword args (not using a dict because this makes
    # inheritance more complicated; this way childrens can simply set a
    # new default_XXX variable and be done with it)
    default_shift = False
    default_binwidth = 2.0
    default_Nmin = 10

    def __init__(self,db, shift=None,P_PMF=None,
                 nequilibration=0,binwidth=None,Nmin=None,
                 reference=config.angles['reference_state'],
                 padding=0,blocksize=2**16,edges=None,**limits):
        """Analysis object for a observable projected on NMP and LID angles.

        db          database which includes metadata and energies
        P_PMF       probability interpolation function P_PMF(NMP,LID); if set to
                    None then all data points in a bin are weighted equally  [None]
        shift       False: display bin-averaged data
                    True:  subtract value of the observable in the reference bin
        edges       (NMP_edges, LID_edges) is used over binwidth and limits to construct
                    the histogram
        binwidth    bin size in degrees (same in NMP and LID)
        Nmin        minimum number of counts in a bin so that the bin
                    is included in the final histogram (must be at
                    least 1)
        reference   (NMP,LID) coordinates of state that is used as the reference state;
                    determines the value of Energy.reference.        
        **limits    limits for the histogram

        Notes:
        * nequilibration is ignored,
        * only padding == 0 works.
        * The default reference value is 1AKE. It is defined in 'lib/templates/config.py.in'.
        """
        
        if self.observable_type is None and self.observable_column is None:
            raise NotImplementedError('This class should not be called itself but must be derived from.')
        if padding != 0:
            raise NotImplementedError("Sorry, padding must be 0 at the moment.")

        if self.observable_column is None:
            # assume that the column is the class name
            self.observable_column = self.__class__.__name__
            
        self.db = db
        self.P_PMF = P_PMF   # spline interpolation function for probability, ie FreeEnergy.P
        if P_PMF is None:
            self.P_PMF = lambda *args: numpy.ones_like(args[0])    # equal weighting when no PMF
            warnings.warn("No PMF probability supplied; "
                          "bin averages will not be statistically correct.")

        # XXX crappy default handling
        if shift is None:
            shift = self.default_shift
        if Nmin is None:
            Nmin = self.default_Nmin
        if binwidth is None:
            binwidth = self.default_binwidth
            
        self.shift = shift
        self.Nmin = Nmin
        if Nmin < 1:
            raise ValueError('Nmin must be >= 1')
        
        # NOTE: to make this class fairly universal we simply use the variable
        # column trick and alias the SQL output column with
        # '%(observable_column)s AS observable'; then use 'all.observable' in the
        # blockiterator further below
        self.selection = db.selection("""SELECT NMP,LID,%s AS observable FROM __data__
                                         WHERE observable NOTNULL""" % self.observable_column)

        self.nequilibration = nequilibration   # ignore --- should we honour this?
        self.padding = padding                 # ignored
        self.edges = edges
        self.binwidth = binwidth

        if self.edges is not None:
            # overrides bins and limits arguments
            nbins = tuple([len(e)-1 for e in self.edges])  # NMP, LID
            anglerange = [[numpy.min(e),numpy.max(e)] for e in self.edges]
            self.limits = {'minNMP':anglerange[0][0], 'maxNMP':anglerange[0][1],
                           'minLID':anglerange[1][0], 'maxLID':anglerange[1][1]}
            self._histbinargs = {'bins': self.edges}
        else:            
            # limits
            # use angleranges wherever input contains None
            (minNMP,maxNMP),(minLID,maxLID) = self.db.angleranges(**limits)
            self.limits = {'minNMP':minNMP, 'maxNMP':maxNMP, 'minLID':minLID, 'maxLID':maxLID}

            # Bins in Project.wham (must be the same for consistency)
            nbins_NMP = int(numpy.ceil((maxNMP-minNMP)/self.binwidth))
            nbins_LID = int(numpy.ceil((maxLID-minLID)/self.binwidth))

            anglerange = [[minNMP,maxNMP],[minLID,maxLID]]
            nbins = (nbins_NMP,nbins_LID)
            self._histbinargs = {'bins': nbins, 'range': anglerange}                    
        # in order to reconstruct the exact histogram use **self._histbinargs
                
        # Do the histogramming in blocks for big selections because otherwise the memory blows up.
        h = numpy.zeros(nbins,dtype=numpy.float_) # final binned average
        h_i = h.copy()                           # tmp array for each block
        s = h.copy()                             # standard deviation
        s_i = s.copy()                           # tmp array for each block
        w = h.copy()                             # total prob. weight
        w_i = w.copy()                           
        N = numpy.zeros(nbins,dtype=numpy.int_)  # total counts per bin
        N_i = N.copy()

        def histogram2d(*args,**kwargs):
            """New numpy.histogram with correct bins."""
            kwargs['normed'] = False   # never normalize for bin-averaging
            kwargs.update(self._histbinargs)
            return numpy.histogram2d(*args,**kwargs)

        # empty run to get e1 and e2 in case there are no data and the loop is skipped
        N_i[:],e1,e2 = histogram2d([0],[0])
        
        for all in self.selection.block_iterator(blocksize=blocksize,asrecarray=True):
            NMP, LID, observable = all.NMP, all.LID, all.observable
            probabilities = self.P_PMF(NMP,LID)
            # number of values per bin or total probability weight
            N_i[:],e1,e2 = histogram2d(NMP,LID)
            w_i[:],e1,e2 = histogram2d(NMP,LID, weights=probabilities)
            # sum of values per bin
            h_i[:],e1,e2 = histogram2d(NMP,LID, weights=probabilities*observable)
            # calculate variance in one pass: <(E - <E>)**2> = <E**2> - <E>**2
            s_i[:],e1,e2 = histogram2d(NMP,LID, weights=probabilities*observable**2)

            del all
            del NMP  # release memory immediately
            del LID  # ... and note that they can be of len < blocksize, so no simple preallocation
            h += h_i    # accumulate bin sums ...
            s += s_i    # ... sum of squares
            N += N_i    # ... and bin counts
            w += w_i    # ... and sum of weights
            
        del h_i
        del s_i
        del N_i
        del w_i

        if P_PMF is None:
            normalization = N
        else:
            normalization = w
        bins_with_data = (N > 0)                # only use bins that contain data
        h[bins_with_data] /= normalization[bins_with_data]  # average

        bins_with_data = (N > Nmin)             # N-1 standard deviation
        s[N <= Nmin] = 0                        # no values for N=0,1 data points (or mask array?)
        s[bins_with_data] /= normalization[bins_with_data]
        s[bins_with_data] -= h[bins_with_data]**2
        if P_PMF is None:
            # when having equal weights then I know how to calculate the proper
            # standard deviation (with N-1):
            s[bins_with_data] *= N[bins_with_data] / (N[bins_with_data] - 1)  # N/(N-1) * VarX
        s[bins_with_data] = numpy.sqrt(s[bins_with_data])                 # final stddev

        self.observable       = numpy.ma.array(h, mask=numpy.logical_not(bins_with_data))
        self.observable_stdev = numpy.ma.array(s, mask=numpy.logical_not(bins_with_data))
        self.bincounts    = numpy.ma.array(N,fill_value=0)   # MA array for consistency
        self.binpartfunc  = numpy.ma.array(w, mask=numpy.logical_not(bins_with_data),
                                           fill_value=0.0)

        # midpoints of bins
        self.mNMP = 0.5*(e1[1:]+e1[:-1])  # NMP
        self.mLID = 0.5*(e2[1:]+e2[:-1])  # LID

        # shift observable to reference
        self.reference_state = reference
        try:
            self.reference = self.get_reference_observable(*reference)
        except IndexError:
            if self.shift:
                warnings.warn('Reference point not within data limits. Will not shift.')
            self.reference = None
            self.shift = False            
        if self.shift:
            self.observable -= self.reference

        # spline interpolation
        self.F = self.interpolationFunction()
        self.DeltaF = self.interpolationFunction(self.observable_stdev)
        self.DeltaF.__doc__ = "DeltaF(NMP,LID): Spline interpolation of the standard deviation"
        self.N = self.interpolationFunction(self.bincounts,cval=0)
        self.N.__doc__ = "N(NMP,LID): spline interpolation of the bin counts"
        self.Z = self.interpolationFunction(self.binpartfunc,cval=0.0)
        self.Z.__doc__ = "Z(NMP,LID): spline interpolation of the sum of bin probability weights"

    def get_reference_observable(self,refNMP,refLID):
        """Returns the value of the observable at the reference point."""
        # find bin of ref: the one with a single count
        ref,e1,e2 = numpy.histogram2d([refNMP],[refLID], normed=False, **self._histbinargs)        
        return self.observable[ref == 1][0]

    def interpolationFunction(self,data=None,spline_order=3,cval=None):
        """Returns a function F(x,y) that interpolates any values of the 2D observable.

        F = interpolationFunction(self,data=None,spline_order=3,cval=None)

        cval is set to data.mean() by default. cval cannot be chosen too large
        or too small or NaN because otherwise the spline interpolation breaks
        down near the region and produces wild oscillations.
        """
        # see http://www.scipy.org/Cookbook/Interpolation
        from scipy import ndimage

        if data is None:
            data = self.observable

        if cval is None:
            cval = data.mean()
        _data = data.filled(cval)   # fill with moderate value, not 1000 to keep spline happy

        coeffs = ndimage.spline_filter(_data,order=spline_order)
        x0,y0 = self.mNMP[0], self.mLID[0]
        dx,dy = self.mNMP[1] - x0, self.mLID[1] - y0
        def _transform(cnew, c0, dc):
            return (numpy.atleast_1d(cnew) - c0)/dc
        def interpolatedF(NMP,LID):
            """B-spline function over the 2D observable, F(NMP,LID).

            Example usage:
            >>> XX,YY = numpy.mgrid[40:75:0.5,96:150:0.5]
            >>> FF = Observable.F(XX,YY)
            >>> contourf(XX,YY,FF, numpy.linspace(-3,11,100),extend='both')
            """
            return ndimage.map_coordinates(coeffs,
                                           numpy.array([_transform(NMP,x0,dx),_transform(LID,y0,dy)]),
                                           prefilter=False,mode='constant',cval=cval)
        return interpolatedF            


    def contour_area(self,edges,resolution=0.1,**limits):
        """Return the area of each contour interval within the limits.

        edges           list containing the contour edges (1D array)
        resolution      resolution of the spline map used to calculate the area

        TODO: better handling of limits
        """

        lim = self.limits.copy()
        lim.update(limits)
        XX,YY = numpy.mgrid[lim['minNMP']:lim['maxNMP']:resolution,
                            lim['minLID']:lim['maxLID']:resolution]
        A = self.F(XX,YY)
        e = numpy.asarray(edges)
        area = numpy.zeros(len(e)-1,dtype=int)
        area[:] = [len(A[numpy.logical_and(lo <= A, A < hi)]) for lo,hi in zip(e[:-1],e[1:])]
        return area * resolution**2

        
    def plot(self,mode="observable",figname=None,cmap=None,title=None,alpha=1.0,
             min_contour=None,max_contour=None,n_contours=100,delta_contour=None,
             with_colorbar=True,with_contour=True,with_filled=True,clf=True,
             linewidth=1,colors='k',extend=None,fmt="%.1f",
             use_spline=True, resolution=0.5, maskfactor=0.1,
             overlay=None):
        """Plot data, defined by mode.

           plot(mode='observable',**kwargs)

        :Arguments:
        mode           observable   (average observable)
                       stdev        (standard deviation (N-1))
                       counts       (bin counts)
        figname        if not None, write figure to file [None]
        title          add title to figure [None]
                       
        min_contour    | boundaries for the contour plot
        max_contour    |
        delta_contour  distance of major contour lines [None]
        n_contours     number of filled contoiurs [100]
        cmap           colormap instance, eg cm.hot_r  [cm.jet]
        extend         None: automatic; 'neither','min','max','both'
        fmt            format string for the colorbar [%.1f]

        with_colorbar  | use these options to customize the plot
        with_contour   | when plotting over an existing one;
        with_filled    | all are True by default [True]
        clf            set to None when plotting on top of other graphs
        linewidth      base linewidth for contours [1]
        colors         color of contour lines ['black']

        use_spline    use spline-interpolated data instead of original unbiased histogram [True]
        resolution    when using splines, resolution in degrees along each axis [0.5]
        maskfactor    with splines, adjust to tune masking of areas without data 0...1 [0.1]

        overlay   function that is called with no arguments, eg
                  overlay = lambda : plot(....)        
        """
        import pylab

        # TODO: find a neater way to organize defaults
        _datasources = {}
        _datasources['spline'] = {'observable': self.F,
                                  'stdev':  self.DeltaF,
                                  'counts': self.N,
                                  'pweights': self.Z,
                                  }
        _datasources['raw'] = {'observable': self.observable,
                               'stdev':  self.observable_stdev,
                               'counts': self.bincounts,
                               'pweights': self.binpartfunc,
                               }
        if use_spline:
            _NMP,_LID = numpy.mgrid[self.mNMP.min():self.mNMP.max():resolution,
                                    self.mLID.min():self.mLID.max():resolution]
            NMP,LID = numpy.ogrid[self.mNMP.min():self.mNMP.max():resolution,
                                  self.mLID.min():self.mLID.max():resolution]
            NMP,LID = NMP.ravel(), LID.ravel()
            
            try:
                _data = _datasources['spline'][mode](_NMP,_LID)
            except KeyError:
                raise ValueError("mode must be one of %r." % _datasources.keys())            

            mask = _datasources['raw'][mode].mask
            if numpy.any(mask):
                # now mask regions that had data missing (i.e. filled with cval)
                # 1) must histogram mask on new grid
                X,Y = numpy.meshgrid(self.mNMP, self.mLID)
                X0 = X[mask.T]  # pick coordinates with True (masked)
                Y0 = Y[mask.T]  # ... need transposed here; check: plot(X1,Y1,'ks')
                # bin edges from midpoints (from PMF.wham.FreeEnergy)
                def _mid2edges(m):
                    # construct 2 additional points before first and and after last midpoint
                    # so that m[0] and m[-1] will end up at the centre of the first and last bin
                    # (simply subtract/add first/last binwidth: m' = m[0] - (m[1]-m[0])
                    x = numpy.concatenate(([2*m[0]-m[1]], m, [2*m[-1]-m[-2]]))
                    return 0.5*(x[1:] + x[:-1])  # edges are the midpoints of the extended midpoints
                edgesX = _mid2edges(NMP)
                edgesY = _mid2edges(LID)
                _mask,e1,e2 = numpy.histogram2d(X0,Y0,bins=[edgesX,edgesY])
                # not working yet: this only gives points with space in between; 
                # 2) we need to convolve with disks of radius ~ delta_old
                deltaX = (self.mNMP[1:] - self.mNMP[:-1]).mean()
                deltaY = (self.mLID[1:] - self.mLID[:-1]).mean()            

                import scipy.ndimage    # use a Gaussian filter
                _smoothed_mask = scipy.ndimage.gaussian_filter(_mask,[deltaX,deltaY])
                cutoff = maskfactor * _smoothed_mask.max()  # heuristic choice...
                new_mask = (_smoothed_mask > cutoff)
                # 3) mask out data
                data = numpy.ma.array(_data, mask=new_mask)
                del _mask
                del _smoothed_mask
            else:
                # no data points missing so ...
                data = _data
            del _data
        else:
            NMP,LID = self.mNMP,self.mLID
            try:
                data = _datasources['raw'][mode]
            except KeyError:
                raise ValueError("mode must be one of %r." % _datasources.keys())

        _formats = {'observable': fmt ,
                    'stdev': fmt ,
                    'counts': "%d",
                    'pweights': fmt,
                    }
        _extend = {'observable': extend or "both",  # 'neither' for no pointy ends
                   'stdev': "max",
                   'counts': "max",
                   'pweights': "max",
                   }
        
        _min_contour = {'observable': data.min(),
                        'stdev': 0.0,
                        'counts': 0.0,
                        'pweights': 0.0,
                        }                               
        
        was_interactive = pylab.matplotlib.is_interactive()
        pylab.matplotlib.interactive(False)

        cmap = cmap or pylab.cm.jet

        # make sure that min/max are aligned with contours
        major_contour = float(delta_contour or self.major_contour)
        minor_contour = 0.5*major_contour

        if min_contour is None:
            min_contour = numpy.floor(_min_contour[mode])
        if max_contour is None:
            max_contour = numpy.ceil(data.max())

        nbot = numpy.floor(min_contour/major_contour)
        ntop = numpy.ceil(max_contour/major_contour)
        min_contour = nbot * major_contour
        max_contour = ntop * major_contour
        nticks = nbot + ntop 

        if clf or with_colorbar:
            pylab.clf()

        # be selective for what we plot
        if not (with_filled or with_contour):
            raise ValueError('At least one of with_filled and with_contour must be True.')

        # Note that plotting the array reverses the axes (we use C-order):
        # array = (row, col) <--> (y, x) !!
        # and hence we use the TRANSPOSED of the array.
        XYZ = (NMP,LID, data.T)
        def contourf(*N,**kwargs):
            args = XYZ + N
            if with_filled:
                return pylab.contourf(*args,**kwargs)
        def contour(*N,**kwargs):
            kwargs['colors'] = colors
            args = XYZ + N
            if with_contour:
                return pylab.contour(*args,**kwargs)
        def colorbar(**kwargs):
            if with_colorbar:
                return pylab.colorbar(**kwargs)

        def ndiv(xmin,xmax,delta):
            return int(numpy.ceil((xmax-xmin)/delta)) + 1

        fudge = 1e-8  # must start a tad below eg 0; otherwise the exact value is left out
        contf = contourf(numpy.linspace(min_contour-fudge,max_contour,n_contours,endpoint=True),
                         extend=_extend[mode],cmap=cmap,alpha=alpha)        
        if mode in ('counts','pweights'):
            colorbar(extend=_extend[mode],format=_formats[mode])
            cont = contour(linewidths=0.5)        
        else:
            colorbar(extend=_extend[mode],
                     ticks=numpy.linspace(min_contour,max_contour,
                                          ndiv(min_contour,max_contour,major_contour)),
                     format=_formats[mode])
            cont = contour(numpy.linspace(min_contour,max_contour,
                                          ndiv(min_contour,max_contour,minor_contour)),
                           linewidths=0.2*linewidth,alpha=0.5)        
            cont = contour(numpy.linspace(min_contour,max_contour,
                                          ndiv(min_contour,max_contour,major_contour)),
                           linewidths=0.5*linewidth)        
        try:
            contf.ax.set_aspect('equal')
        except AttributeError:
            cont.ax.set_aspect('equal')

        pylab.xlabel(r'angle NMP-CORE')
        pylab.ylabel(r'angle LID-CORE')

        if title:
            pylab.title(title)

        ax = pylab.gca()
        degreeFormatter = pylab.matplotlib.ticker.FormatStrFormatter(r'%d$^\circ$')
        ax.xaxis.set_major_formatter(degreeFormatter)
        ax.yaxis.set_major_formatter(degreeFormatter)

        L = self.limits
        pylab.xlim(L['minNMP'],L['maxNMP'])
        pylab.ylim(L['minLID'],L['maxLID'])        

        if overlay:
            # this is very powerful/dangerous/stupid: user can execute anything here
            overlay.__call__()

        if was_interactive:
            pylab.draw()
        
        if figname:
            # http://banyan.usc.edu/log/python/plot : use str to work around
            # 'TypeError: cannot return std::string from Unicode object' with Agg
            pylab.savefig(str(figname))

        pylab.matplotlib.interactive(was_interactive)  # revert to previous state


    def filebasename(self,runtype=None):
        """Create a filename that contains the RUNTYPE.

        filebasename(runtype=RUNTYPE)

        filebasename = observable_type + _ + RUNTYPE

        RUNTYPE is eg 'dims' or 'dims_co'.

        None      take the runtype from the observable
        ""        empty string leaved the filename unchanged

        """
        if runtype is None:
            try:
                runtype = self.runtype
            except AttributeError:
                pass
        if runtype:
            runtype = '_' + runtype
        runtype = runtype or ""
        return self.observable_type + runtype


    def autoplot(self,mode='observable',runtype=None,figdir=None, **plotargs):
        """Make the standard plot of the data.

        autoplot(AngleprojectedObservable, mode=MODE,runtype='RUNTYPE',figdir=None,**plotargs)

        plotargs (see AngleprojectedObservable.plot()):

        title       None:   generate title for each individual plot
                            using the observables legend
                    string: use provided title for _all_ plots
                    "":     no title
        min_contour |
        max_contour +- boundaries for the contour plot                    
        cmap        pylab colormap instance or None for the default [cm.jet]
        alpha       alpha for the filled contours [1.0]

        overlay   function that is called with no arguments, eg
                  overlay = lambda : coRMSD.DeltaRMSD.plot(with_colorbar=False,clf=False,with_filled=False,linewidth=3,colors='white')

        """
        import pylab

        qualifier = self.qualifier or ""  # gives just trailing slash if no self.qualifier
        figdir = figdir or os.path.join(config.basedir,'figs',self.observable_type,qualifier)
        filebasename = self.filebasename(runtype)
        
        was_interactive = pylab.matplotlib.is_interactive()
        pylab.matplotlib.interactive(False)

        self._auto_plots(mode,filebasename,figdir,plotargs)

        pylab.matplotlib.interactive(was_interactive)  # revert to previous state
        return figdir

    
    def plot_all(self,runtype=None,figdir=None, **plotargs):
        """Make the 6 standard plots of the observable data.

        plot_all(runtype='RUNTYPE',figdir=None,**plotargs)

        plotargs (see AngleprojectedObservable.plot()):

        title       None:   generate title for each individual plot
                            using the observables legend
                    string: use provided title for _all_ plots
                    "":     no title
        min_contour |
        max_contour +- boundaries for the contour plot                    
        cmap        pylab colormap instance or None for the default [cm.jet]
        alpha       alpha for the filled contours [1.0]

        overlay   function that is called with no arguments, eg
                  overlay = lambda : coRMSD.DeltaRMSD.plot(with_colorbar=False,clf=False,with_filled=False,linewidth=3,colors='white')

        """

        kwargs = dict(runtype=runtype,figdir=figdir,**plotargs)
        self.autoplot('observable',**kwargs)
        self.autoplot('stdev',**kwargs)
        figdir = self.autoplot('counts',**kwargs)        

        print "-- Figures in %r." % figdir
        return figdir

    def _auto_plots(self,mode,filebasename,figdir,plotargs):
        """Generate standard plots and write png and and pdf. Chooses filename and plot title."""
        import pylab
        from angles import make_canonical_plot

        try:
            os.makedirs(figdir)
        except OSError,err:
            if err.errno != errno.EEXIST:
                raise

        def figs(*args):
            return os.path.join(figdir,*args)
        
        modefilebasename = filebasename + self._suffix[mode]
        _plotargs = plotargs.copy()  # need a copy because of changing 'title'
        if plotargs.get('title') is None:  # None --> set automatic title
            if self.legend is None:
                legend = ""
            else:
                legend = " "+self.legend
            _plotargs['title'] = self._title[mode]+legend
        _plotargs['figname'] = figs(modefilebasename + '.png')

        pylab.clf()
        self.plot(mode=mode,**_plotargs)
        pylab.savefig(figs(modefilebasename + '.pdf'))   # pdf
        
        make_canonical_plot(xray=True)
        pylab.savefig(figs(modefilebasename + '_canonical.png'))
        pylab.savefig(figs(modefilebasename + '_canonical.pdf'))
        print "--- Plotted plain and canonical %(modefilebasename)r (png,pdf)." % vars()


    def save(self,filename):
        """Write all data to a pickle file. (not used yet)"""
        import cPickle

        filename = os.path.splitext(filename)[0] + '.pickle'        

        # self.__dict__ simpler? no, need to filter db and selection
        data = {'observable_type': self.observable_type,
                'qualifier': self.qualifier,
                'observable_column': self.observable_column,
                'legend': self.legend,
                'observable': self.observable,
                'observable_stdev': self.observable_stdev,
                'bincounts': self.bincounts,
                'reference': self.reference, 'reference_state': self.reference_state,
                'mLID': self.mLID, 'mNMP': self.mNMP,
                'binwidth': self.binwidth, 'padding': self.padding,
                '_histbinargs': self._histbinargs,
                }
        f = open(filename,'wb')
        try:
            cPickle.dump(data, f, cPickle.HIGHEST_PROTOCOL)
        finally:
            f.close()

class Projector(Observable):
    """Project Observable A on observable X.

      p(X) = p(X(NMP),X(LID))
      A(X) dx = A(NMP(X),LID(X)) p(X) dx  (?)


    This requires the 2D PMF. By default A is 1 and hence one obtains the
    probability distribution of observable X.
    """
    default_binwidth = 0.05    # degrees, on the spline map
    default_resolution = 0.5   # unit of the observable X, bins of the 1D histogram
    
    def __init__(self,db=None,PMF=None,X=None,A=None,binwidth=None,
                 resolution=None,edges=None,debug=False):
        """
        PMF             free energy object
        X               Angle projected observable object
        A               Angle projected observable object (XXX) or 1 [1]
        binwidth        degrees, resolution
        resolution      binwidth of the final histogram along X;
                        the histogram stretched from min to max of the data.
        edges           explicitly set edges of the X-histogram [XXX]
        
        Data limits are set by looking at the domain of the PMF: use PMF.set_datalimits()!
        """
        
        if PMF is None or X is None:
            raise ValueError('PMF and X are required')
        if A is not None:
            raise NotImplemented

        # db is not needed and ignored; it's only here to fit in the ObservableRegistry framework
        # self.db = db
        self.PMF = PMF
        self.P = self.PMF.P
        self.observable_X = X
        self.X = self.observable_X.F     # spline function
        self.observable_A = A or 1
        try:
            self.A = self.observable_A.F # spline function
        except AttributeError:
            self.A = self.observable_A   # might work for numbers and proper user supplied funcs

        self.binwidth = binwidth or self.default_binwidth    # degrees, on the spline map
        self.resolution = resolution or self.default_resolution


        # limits are based on whatever the PMF has set
        lim = self.PMF.limits
        
        NMP,LID = numpy.mgrid[lim[0,0]:lim[1,0]:self.binwidth,
                              lim[0,1]:lim[1,1]:self.binwidth]  # hi-res grid
        X = self.X(NMP,LID)  # interpolated on grid with fine resolution
        P = self.P(NMP,LID)  # interpolated on grid with fine resolution
        #P[P<0] = 0           # fix spline errors: no negative probabilities
        assert numpy.any( P<0 ) == False  # new interpolated probabilites should be ok

        # X-histogram contours are defined by the edges
        if edges:
            e = numpy.asarray(edges)
        else:
            # get them from min/max of the data
            xmin,xmax = X.min(),X.max()
            e = numpy.linspace(xmin,xmax, numpy.ceil((xmax-xmin)/self.resolution)+1)

        # integration between contours defined by X
        #
        contours = zip(e[:-1],e[1:])
        
        h = numpy.zeros(len(e)-1,dtype=numpy.float_)
        s = h.copy()
        area = h.copy()
        # sum probability over the whole contour
        # NOTE: weighting by A not included yet;
        def _contour_integral(X,lo,hi,P,A):
            # NOT USED YET
            contour = numpy.logical_and(lo <= X, X < hi)
            return A[contour] * P[contour]

        h[:] = [ P[numpy.logical_and(lo <= X, X < hi)].sum()    for lo,hi in contours ]
        area[:] = [ len(P[numpy.logical_and(lo <= X, X < hi)])  for lo,hi in contours ]
        area *= self.binwidth**2

        self.normalization = h.sum()
        dx = e[1:] - e[:-1]
        self.observable = h / (self.normalization * dx)
        self.observable_stdev = self.observable / numpy.sqrt(area/self.binwidth**2) # ???
        self.bincounts = area      # leave raw area for diagnostics 

        # midpoints and edges of bins
        self.xvalues = 0.5*(e[1:]+e[:-1])
        self.xedges = e

        if debug:
            # for diagnostics
            self._debug = {'NMP':NMP,'LID':LID,'P':P,'X':X,
                           'deltaNorm': self.normalization - P.sum()}
        
# Not recommended/incomplete.
class PMFProjectedObservable(Observable):
    """Probability distribution of a 1D observable X.

    DO NOT USE. WRONG.

    Requires a PMF object and a angle-projected observable X.

    The probability from the PMF is summed over contours of the 1D observable X
    and normalized to the area of each contour. This should be the
    experimentally observed distribution of the observable X.

    THIS class calculates the sum of probabilities for all data points within
    the limits defined by the PMF.

    (Alternatively, use the Projector class instead, which does not require the
    db, only the spline functions.)


    NOTE: The data are restricted to strictly within the limits of the PMF and
     the contour area analysis; otherwise we would count data we would not
     normalize for. The only way to change the limits is to set different
     limits in the PMF, using

        PMF.set_datalimits(...)
    """
    
    default_binwidth = 0.2
    x_column = None         # HACK! will be set from observable_column
    
    def __init__(self,db,PMF=None,observable=None,
                 binwidth=None,Nmin=10,blocksize=2**16,edges=None,**xlimits):
        """P = PMF_X_Projector(F,X)

        db          AnglesDB
        PMF         PMF.wham.FreeEnergy object
        observable  AngleprojectedObservable object that represents X(NMP,LID)
        
        edges       can be used instead of binwidth and limits to construct the histogram
        binwidth    bin size
        Nmin        minimum number of counts in a bin so that the bin
                    is included in the final histogram (must be at
                    least 1)        
        **xlimits    limits ``xmin`` and ``xmax`` for the histogram [automatic]

        Note: choose correct data limits in the PMF using

           PMF.setdatalimits(...)

        before setting up the Projector; otherwise the spline interpolation
        might not be well defined near regions with no data.
        """
        warnings.warn("Do not use this class; the output is wrong as the normalization "
                      "is incomplete. See source for details.",category=DeprecationWarning)
        
        if self.observable_type is None and self.observable_column is None:
            raise NotImplementedError('This class should not be called itself but must be derived from.')
        if self.observable_column is None:
            # assume that the column is the class name
            self.observable_column = self.__class__.__name__

        # HACK: the ObservableRegistry machinery sets
        # observable_column (which we don't need here) but we do need variable x_column
        self.x_column = self.observable_column
        self.observable_column = None

        self.name = str(self.observable_type) + '_' + str(self.qualifier)
        
        if PMF is None or observable is None:
            raise ValueError('PMF and observable arguments must be provided')
        # we _could_ calculate the X-observable now, but that takes more time...
        
        self.db = db
        self.PMF = PMF
        self.X = observable
        self.Nmin = Nmin
        assert Nmin > 1
        self.P = self.PMF.P     # spline interpolation of the probability

        # restrict selection to PMF limits
        lim = self.PMF.limits
        self.selection = self.db.selection(
            """SELECT NMP,LID, %s AS x FROM __data__
               WHERE x NOTNULL AND NMP BETWEEN %f AND %f
                               AND LID BETWEEN %f AND %f
            """ % (self.x_column, lim[0,0], lim[1,0], lim[0,1], lim[1,1]) )

        self.edges = edges
        self.binwidth = binwidth or self.default_binwidth

        # setting bins for the histogram
        if self.edges is not None:
            # overrides bins and limits arguments
            nbins = len(e)-1  # x
            self.xlimits = {'xmin': numpy.min(e), 'xmax': numpy.max(e)}
            self._histbinargs = {'bins': self.edges}
        else:            
            # limits
            [(xmin,xmax)] = self.selection.SQL("SELECT min(x), max(x) FROM __data__")
            xlimits.setdefault('xmin',xmin)  # use data limits unless something else
            xlimits.setdefault('xmax',xmax)  # was provided in **xlimits
            self.xlimits = xlimits
            xmin = self.xlimits['xmin']      # xmin and xmax to be used
            xmax = self.xlimits['xmax']

            nbins = int(numpy.ceil((xmax-xmin)/self.binwidth))
            x_range = [xmin, xmax]
            self._histbinargs = {'bins': nbins, 'range': x_range}                    
        # in order to reconstruct the exact histogram use **self._histbinargs
        
        
        # calculate observable
        # Do the histogramming in blocks for big selections because otherwise the memory blows up.
        h = numpy.zeros(nbins,dtype=numpy.float_) # final binned average
        h_i = h.copy()                            # tmp array for each block
        N = numpy.zeros(nbins,dtype=numpy.int_)   # total counts per bin
        N_i = N.copy()                            # tmp array for each block

        def histogram(*args,**kwargs):
            """New numpy.histogram with correct bins."""
            kwargs['normed'] = False   # never normalize for bin-averaging
            kwargs['new'] = True       # currently needed for proper edges and weights
            kwargs.update(self._histbinargs)
            return numpy.histogram(*args,**kwargs)
        
        # empty run to get edges e in case there are no data and the loop is skipped
        N_i[:],e = histogram([0])

        # For the projection of the PMF we add the probability weights P'(X) =
        # P(NMP(X),LID(X)) in each bin X.
        for all in self.selection.block_iterator(blocksize=blocksize,asrecarray=True):
            NMP,LID,x, = all.NMP, all.LID, all.x
            del all                 # release memory (?)
            # number of values per bin (for error bars)
            N_i[:],e = histogram(x)
            # sum of values per bin (projection of the probability distribution)
            # note how simple the spline-interpolated probability fits in!!!
            h_i[:],e = histogram(x, weights=self.P(NMP,LID))
            
            del x  # release memory immediately

            h += h_i    # accumulate bin sums ...
            N += N_i    # ... and bin counts
            
        del h_i
        del N_i

        bins_with_data = (N > 0)                # only use bins that contain data

        self.observable       = numpy.ma.array(h, mask=numpy.logical_not(bins_with_data))
        # jacobian: divide each bin by the area of the contour
        self.normalization = self.X.contour_area(e,**self.PMF.limits_dict)
        self.observable /= self.normalization

        # simple error estimat --- not sure if this makes much sense
        self.observable_stdev = self.observable/numpy.sqrt(N)
        self.bincounts    = N

        # midpoints and edges of bins
        self.xvalues = 0.5*(e[1:]+e[:-1])
        self.xedges = e



# ------------------------------------------------------------
# testing/quick hacks

# outside so that i can use for existing objects, can be removed soon
# OB 2009-04-18
def autoplot(self,mode='observable',runtype=None,figdir=None, **plotargs):
    """Make the standard plot of the data.

    autoplot(AngleprojectedObservable, mode=MODE,runtype='RUNTYPE',figdir=None,**plotargs)

    plotargs (see AngleprojectedObservable.plot()):

    title       None:   generate title for each individual plot
                        using the observables legend
                string: use provided title for _all_ plots
                "":     no title
    min_contour |
    max_contour +- boundaries for the contour plot                    
    cmap        pylab colormap instance or None for the default [cm.jet]
    alpha       alpha for the filled contours [1.0]

    overlay   function that is called with no arguments, eg
              overlay = lambda : coRMSD.DeltaRMSD.plot(with_colorbar=False,clf=False,with_filled=False,linewidth=3,colors='white')

    """
    import pylab

    qualifier = self.qualifier or ""  # gives just trailing slash if no self.qualifier
    figdir = figdir or os.path.join(config.basedir,'figs',self.observable_type,qualifier)
    filebasename = self.filebasename(runtype)

    was_interactive = pylab.matplotlib.is_interactive()
    pylab.matplotlib.interactive(False)

    self._auto_plots(mode,filebasename,figdir,plotargs)

    pylab.matplotlib.interactive(was_interactive)  # revert to previous state
    return figdir


def plot(self,col):
    from pylab import plot,xlabel,title
    X = self.recarray   # quick hack... dont use for big dbs!
    Y = X.__getattribute__(col)
    hI,e = numpy.histogram(X.DeltaRMSD, bins=50, weights=Y, new=True)
    hN,e = numpy.histogram(X.DeltaRMSD, bins=50, new=True)

    havg = hI.astype(float)/hN

    progress = 0.5*(e[1:]+e[:-1])

    plot(progress, havg, label=col)
    xlabel(r'progress $\Delta$RMSD (\AA)$')
    # title('Saltbridge indicator SB1, c->o')

def plot_all(self):
    for N in xrange(1,11):
        col = 'ISB%d' % N
        print "-- Plotting %(col)r" % vars()
        plot(self,col)



def DIMS_trajectory_length(db, all=True, binwidth=1.0, direction='co', endpoint='ensemble',
                           **plotargs):
    """Histograms of DIMS run lengths.

    If all = True:
       Do statistics on whole db (use a selection if more control is needed.)     

       Use direction (co, oc) and endpoint (fixed, ensemble) for labels.
    
    If all = False:
       - Distinguishes between O->C and C->O runs.
       - Distinguishes between fixed endpoint runs (000 - 299) and
         ensemble endpoint runs (>=300).

         HARDCODED!!! -- see txt/trajectories.txt

         In order for this to work, it assumes that the trajectorynames are either

           ocNNN

         or

           coNNN

         Anything else is ignored.
    """

    # parts of the search pattern: in order to distinguish
    # trajectories we have to parse the trajectory name
    p = {'direction': {'co':'co', 'oc':'oc'},
         'endpoint':  {'fixed':    '[0-2][0-9][0-9]',
                       'ensemble': '[0-9]*[3-9][0-9][0-9]'}
         }

    # WARNING ... crappy code ahead (patched all == True in later...)

    counts = []
    mins = []
    maxs = []
    legends = []

    if all == True:
        n = numpy.array([x[0] for x in db.SQL("SELECT COUNT(*) FROM __self__ GROUP BY TID")])
        mins.append(n.min())
        maxs.append(n.max())
        counts.append(n)
        legends = [(direction, endpoint)]          
    else:        
        for endpoint in ('fixed', 'ensemble'):
            for direction in ('co', 'oc'):
                n = numpy.array([x[0] for x in
                                 db.SQL("SELECT COUNT(*) FROM __self__ WHERE REGEXP(?, trajectory) GROUP BY TID",
                                        p['direction'][direction]+p['endpoint'][endpoint])])
                mins.append(n.min())
                maxs.append(n.max())
                counts.append(n)
                legends.append((direction, endpoint))

    # global histogram boundaries
    hmin = numpy.min(mins) - 1
    hmax = numpy.max(maxs) + 1
    bins = numpy.arange(hmin, hmax, binwidth)

    histograms = [numpy.histogram(n, bins, new=False) for n in counts]

    try:
        from pylab import step, plot, xlabel, ylabel, legend, xlim
        for (direction, endpoint),(h,m) in zip(legends, histograms):
            plotargs['label'] = "%(direction)s (%(endpoint)s)" % vars()
            plot(m,h, **plotargs)
        xlabel("DIMS trajectory length (ps)")
        ylabel("count")
        legend(loc="best")
    except ImportError:
        print "Cannot import pylab, no picture :-("

    return histograms
    
