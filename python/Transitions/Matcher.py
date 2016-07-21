# $Id: Matcher.py 3425 2009-04-29 19:01:11Z oliver $
"""
Motivated by [Vonrhein1995] we try to match each pdb to a frame in
each trajectory by computing the RMSD. Associate the pdb with the time
step that has the smallest RMS. This will order the pdbs for a given
transition i.

More formally:

The i-th DIMS trajectory consists of frames X at time step t: X_i(t)

The static pdbs S are numbered sequentially; eg S_0 = 1AKE, S_1 =
4AKE, ... S_j, ... S_N

for each DIMS traj X_i:
	for each S_j:
		t_ij = {t: RMSD(X_i(t),S_j) == min_t' RMSD(X_i(t'),S_j)}
                       # t is the time step t' that minimizes the RMSD

Then for each trajectory i, sort the static pdbs j by increasing t_ij: 

   t_ij0 <= t_ij1 <= ... <= t_ijN   ==> (j0, j1, ..., jN) 

This gives us a 'movie' of static pdb snapshots

  (S_j0, S_j1, ..., S_jN) 


Now we want to know:
- By how much does the RMSD vary?
- Are the S_j always ordered in the same order (i.e. what's the
  correlation between the t_i ?).

IF the DIMS trajectories are different by OM or any other metric but
the ordering is very similar then I would say we're picking up some
sort of real transition.

Example:

  import Transitions.Matcher
  M = Transitions.Matcher.PDBMovie('data/all_rms1-251.dat')
  M.imshow('co',colorbar=True,vmax=1.2,
     title=r'matches against O$\rightarrow$C, with C$\rightarrow$O ordered accordingly')
  savefig('dist_combined_o-c.pdf')

  M.imshow('oc',colorbar=True,vmax=1.2,
     title=r'matches against C$\rightarrow$O, with O$\rightarrow$C ordered accordingly')
  savefig('dist_combined_c-o.pdf')

  # collect phis for each target
  M.data.select('target,array(phi) AS "phis [NumpyArray]"','WHERE direction = "oc" GROUP BY target')

  # distribution (histogram density)
  M.data.select('target,distribution(phi,50,0,1) as "phi_dist [Object]"','where direction = "oc" group by target')

import Transitions.Matcher
M = Transitions.Matcher.PDBMovie('data/all_rms1-899.dat')

M.imshow('oc',colorbar=True,vmax=1.5,title=r"AdK: crystal structure matches against DIMS open $\rightarrow$ closed")

M.imshow('co',colorbar=True,vmax=1.5,title=r"AdK: crystal structure matches against DIMS closed  $\rightarrow$ open")

M.imshow('oc','RMSD',colorbar=True,vmax=4,title=r"open $\rightarrow$ closed: Mean RMSD for matches ($\AA$)")

M.imshow('co','RMSD',colorbar=True,vmax=4,title=r"closed $\rightarrow$ open: Mean RMSD for matches ($\AA$)")

"""

import numpy
import Transitions.cm
import Transitions.sqlarray

class PDBMovie(object):
    """Assemble one (or an ensemble) of 'movies' (sequentially ordered
    images) from pdb files. The order is determined by MD
    trajectories. For each MD trajectory one movie is created. The pdb
    structures are sorted by RMSD to the MD trajectory.

        Datafile format (white space separated):

        +-transition direction (useful for book keeping)
        |       +-trajectory number (book keeping)
        |       |       +-target (filename of pdb file)
        |       |       |       +- frame # at which RMSD is minimal
        |       |       |       |       +-traj length (in sampled frames)
        |       |       |       |       |       +-RMSD at frame #
        |       |       |       |       |       |         +-progress at frame
        |       |       |       |       |       |         |        +-DeltaRMSD
        |       |       |       |       |       |         |        |
        cl      160     4AKE_B  0       83      0.827793  0.0      -6.2345
        cl      161     1AK2_A  31      87      3.667041  3.520789 -5.8888
        ...

        The combination (direction,number) must be unique as it is
        used as a key to uniquely identify the source trajectory.

        The position along the transition is calculated as

            phi = frame/trajectorylength
        i.e.
            phi(i) = i/N    with   1 <= i <= N

        This is different from RMSD-based progress variables as it
        contains time information; it is the time normalized to the
        length of the trajectory.

        DIMS uses by default the RMSD from the target state as progress variable:

           progress(i) = RMSD(X_i, X_final)

        The DeltaRMSD (Delta\rho) reaction coordinate is defined as

            DeltaRMSD(i) := RMSD(X_i,X_inital) - RMSD(X_i,X_final)

        DeltaRMSD < 0 means closer to initial state, DeltaRMSD > 0
        means closer to final state.
    """

    allowed_variables = ("progress","phi","DeltaRMSD")
    _variable_tex = {'progress':r'RMSD $\rho$ ($\AA$)',
                     'phi':r'$\phi$',
                     'DeltaRMSD':r'$\Delta\rho := \rho_{i} - \rho_{f}$ ($\AA$)'}


    def __init__(self,datafilename,bins=50,progressvariable="DeltaRMSD"):
        """Assemble the movie from the pre-computed datafile.

        M = PDBMovie(<datafilename>,bins=50,progressvariable="DeltaRMSD")

        datafilename      describes matches of x-ray structures to DIMS trajectories
                          (see class doc string for the definition of the format)
        bins              number of bins of the histogram/distribution along the
                          progress variable
        progressvariable  "DeltaRMSD", "progress", "phi"

                          DeltaRMSD(i) := RMSD(X_i,X_inital) - RMSD(X_i,X_final)
                          progress(i) := RMSD(X_i,X_final)     [used by standard DIMS]
                          phi(i) := i/N   with 1 <= i <= N     [relative time]

                          [X_i: configuration at frame i
                           X_initial: first frame, start of DIMS trajectory
                           X_final: last frame, target of DIMS trajectory]
        """
        self.datafilename = datafilename
        self.bins = bins
        self._observables = {}   # register observables with name and data structure
        matches = open(datafilename,'r')
        reclist = []
        for line in matches:
            f = line.strip().split()
            direction,trjno,target,frame,trjlen,rmsd,progress,DeltaRMSD = \
                str(f[0]),str(f[1]),str(f[2]),int(f[3]),int(f[4]),\
                float(f[5]),float(f[6]),float(f[7])
            phi = frame/(trjlen-1.0)  # approximation for real progressvariable
            record = (direction,trjno,target,frame,trjlen,phi,rmsd,progress,DeltaRMSD)
            reclist.append(record)
        matches.close()
        recarray = numpy.rec.fromrecords(\
            reclist,
            names='direction,number,target,frame,length,phi,rmsd,progress,DeltaRMSD')
        table_name = self.__class__.__name__ + "_" + str(id(self))  # unique name
        self.data = Transitions.sqlarray.SQLarray(table_name, recarray)
        self.data.sql_index('trajectory',('direction','number'),unique=False)
        self.data.sql_index('direction_target',('direction','target'),unique=False)
        del recarray
        self._directions = self.data.select('direction','GROUP BY direction').direction
        self._targets = self.data.select('target','GROUP BY target ORDER BY target ASC').target
        # calculate phi distribution and table of phi statistics
        self._make_distributions(variable=progressvariable,bins=self.bins)
        self._make_phis()
        self.change_progressvariable(progressvariable) # calculate all observables along progress
        
    def statistics(self):
        """Returns recarray of statistics aggregated across targets and directions."""
        return self.data.select(\
            """
            direction,target,
            median(progress) AS median_progress, avg(progress) AS mean_progress, 
               std(progress) AS std_progress,
            median(DeltaRMSD) AS median_DeltaRMSD, avg(DeltaRMSD) AS mean_DeltaRMSD, 
               std(DeltaRMSD) AS std_DeltaRMSD,
            median(phi) AS median_phi, avg(phi) AS mean_phi, std(phi) AS std_phi,
            median(rmsd) AS median_rmsd, avg(rmsd) AS mean_rmsd, std(rmsd) AS std_rmsd
            """,
            """
            GROUP BY direction,target ORDER BY direction,median_progress,median_phi ASC
            """)

    def change_progressvariable(self,progressvariable,bins=None,normed=True):
        """Recalculate distributions with a different progressvariable."""
        self._make_distributions(variable=progressvariable,bins=bins,normed=normed)
        self._make_RMSD_functions(variable=progressvariable,bins=bins)
        self._make_RMSD_std_functions(variable=progressvariable,bins=bins)
        self._make_RMSD_zscore_functions(variable=progressvariable,bins=bins)
        self._make_RMSD_max_functions(variable=progressvariable,bins=bins)
        self._make_RMSD_min_functions(variable=progressvariable,bins=bins)

    def _compute_observable(self,name,resultdictname,SELECT_template,
                            variable=None,bins=None):
        """Wrapper to compute observables from SQL histogram functions

        SELECT_template can contain %(variable)s,  %(bins)d, %(vmin)f, %(vmax)f
        """
        if variable is None:
            variable = self.progressvariable   # throws exception if called first time without arg
        if bins is None:
            bins = self.bins                   # ... but that's ok because it's not for the user
        # register
        self.__setattr__(resultdictname,{})
        resultdict = self.__getattribute__(resultdictname)  # aliased
        self._observables[name] = resultdict   # registry

        # choose (global) progress variable
        if variable not in self.allowed_variables:
            raise ValueError("Variable must be one of %r, not %r." % (self.allowed_variables,variable))
        self.progressvariable = variable

        # limits for the histogram across ALL data
        vmin,vmax = self.data.limits(variable)
        for direction in self._directions:
            # typically uses new aggregate functions (see SQLarray._init_sqlite_functions())
            result = self.data.select(SELECT_template % vars(),
                    """
                    WHERE direction = "%(direction)s" GROUP BY target
                    """ % vars(), 
                    asrecarray=False)
            resultdict[direction] = dict(result)    # sets self.__dict__[resultdictname] !

    def _make_distributions(self,variable="progress",bins=None,normed=True):
        if normed:
            HISTOGRAM = "distribution"
        else:
            HISTOGRAM = "histogram"
        self._compute_observable("distribution", "distribution",
                                 """target,"""+
                                 str(HISTOGRAM)+"""(%(variable)s,%(bins)d,%(vmin)f,%(vmax)f) 
                                                   AS "dist [Object]" """,
                                 variable=variable,bins=bins)

    def _make_RMSD_functions(self,variable="progress",bins=None):
        self._compute_observable("RMSD","RMSD_functions",
                                 """target, 
                                    meanhistogram(%(variable)s,rmsd,%(bins)d,%(vmin)f,%(vmax)f) 
                                    AS "func [Object]"
                                 """,
                                 variable=variable,bins=bins)

    def _make_RMSD_std_functions(self,variable="progress",bins=None):
        self._compute_observable("RMSD_std","RMSD_std_functions",
                                 """target, 
                                    stdhistogram(%(variable)s,rmsd,%(bins)d,%(vmin)f,%(vmax)f) 
                                    AS "func [Object]"
                                 """,
                                 variable=variable,bins=bins)

    def _make_RMSD_zscore_functions(self,variable="progress",bins=None):
        self._compute_observable("RMSD_zscore","RMSD_zscore_functions",
                                 """target, 
                                    zscorehistogram(%(variable)s,rmsd,%(bins)d,%(vmin)f,%(vmax)f) 
                                    AS "func [Object]"
                                 """,
                                 variable=variable,bins=bins)

    def _make_RMSD_max_functions(self,variable="progress",bins=None):
        self._compute_observable("RMSD_max","RMSD_max_functions",
                                 """target, 
                                    maxhistogram(%(variable)s,rmsd,%(bins)d,%(vmin)f,%(vmax)f) 
                                    AS "func [Object]"
                                 """,
                                 variable=variable,bins=bins)

    def _make_RMSD_min_functions(self,variable="progress",bins=None):
        self._compute_observable("RMSD_min","RMSD_min_functions",
                                 """target, 
                                    minhistogram(%(variable)s,rmsd,%(bins)d,%(vmin)f,%(vmax)f) 
                                    AS "func [Object]"
                                 """,
                                 variable=variable,bins=bins)

    def _make_phis(self):
        result = self.data.select(
            """direction,target,
               median(progress) AS median_progress, avg(progress) AS mean_progress, 
                 std(progress) AS std_progress,
               median(DeltaRMSD) AS median_DeltaRMSD, avg(DeltaRMSD) AS mean_DeltaRMSD, 
                 std(DeltaRMSD) AS std_DeltaRMSD,
               median(phi) AS median_phi, avg(phi) AS mean_phi, std(phi) AS std_phi,
               median(rmsd) AS median_rmsd, avg(rmsd) AS mean_rmsd, std(rmsd) AS std_rmsd
            """,
            """GROUP BY direction,target ORDER BY direction,median_progress,mean_progress""")
        self.phis = Transitions.sqlarray.SQLarray('phis',result,connection=self.data.connection)
        self.phis.sql_index('targetorder',('direction','target'),unique=True)
        self.phis.sql_index('sortorder',('median_progress','mean_progress'),unique=False)


    def get_ordered_targets(self,direction):
        """Return a list of target PDB IDs ordered according to the trajectory matches.

        pdbs = M.get_ordered_targets(direction)

        direction        'co' or 'oc' (is matched against direction in self.phis)
        """
        progress_variable = self.progressvariable
        return self.phis.select('target',
                         """WHERE direction = "%(direction)s" 
                         ORDER BY median_%(progress_variable)s,mean_%(progress_variable)s""" \
                                % vars())
        

    def plot(self,direction,legend=True,**kwargs):
        """Plot all distributions; colors are set automatically to kwargs[cmap]."""
        import pylab
        pnormalize = pylab.normalize(vmin=1,vmax=len(self.distribution[direction]))

        kwargs.setdefault('alpha',1.0)
        kwargs.setdefault('linewidth',4)
        fmt = kwargs.pop('fmt','-')
        cmap = kwargs.pop('cmap',pylab.cm.jet)

        count = 0
        for target,(hist,edges) in self.distribution[direction].items():
            count += 1
            midpoints = 0.5*(edges[1:] + edges[:-1])
            kwargs['color'] = cmap(pnormalize(count))
            pylab.plot(midpoints,hist,fmt,label="%s" % target,**kwargs)
        if legend:
            pylab.legend(loc='best',
                         prop=pylab.matplotlib.font_manager.FontProperties(size=6))

    def imshow(self,direction,observable="distribution",**kwargs):
        """Plot the distribution of structures being located closest to phi.

        imshow(direction,observable="distribution",
               colorbar=False,with_lines=True,
               fontsize=6,title=None,**kwargs)

        'observable' selects which data is to be plotted:

        "distribution"    Plot a density plot of p(s,phi) for each X-ray structure
                          s. The probability to find s located in the interval
                          [phi,phi+dphi[ in the transition trajectories is

                          p(s,phi) * dphi
        "RMSD"            average RMSD of each match (across the bin)
        "RMSD_std"        standard deviation of the average in the bin
        "RMSD_max"        max or min of the RMSD values in the bin
        "RMSD_min"       
        "RMSD_zscore"     absolute deviation of the RMSD x from the mean, divided by the
                          standard deviation:

                             z = |x - <x>|/std(x)

        On the y-axis the structures are listed, sorted by the median
        of their phi value with respect to all loaded transitions. The
        x-axis labels the histogram bins [phi,phi+dphi]. The colour
        shade measures the value of the distribution p.

        colorbar          True: add colorbar indicating density
        with_lines        True: add graph for median
        title             String that becomes the title of the graph. One can use 
                          a matplotlib raw string with TeX markup, eg 
                          r"transition open $\\rightarrow$ closed"
        cmap              colour map instance, eg cm.hot [Transitions.cm.jet2]
        auto_vmax         heuristic: cutoff colour scale so that auto_vmax of the total
                          data points are resolved [0.8]

        Other kwargs are directly passed on to pylab.imshow(). It is
        worthwhile adjusting vmax manually in order to increase
        resolution in the lower range.

        Example:        
        >>> M.imshow('oc',vmax=1.5,colorbar=True,
                title=r"distribution of pdb matches against closed $\\rightarrow$ open")
        >>> savefig('dist_oc.pdf')
        """
        import pylab

        # sanity check; if not caught it gives hard to diagnose errors
        if direction not in self._directions:
            raise ValueError("direction can only be one of %r." % self._directions)

        try:
            OBSERVABLES = self._observables[observable]
        except KeyError:
            raise ValueError("observable='%s' must be one of %r" % \
                                 (observable,self._observables.keys()))

        # PHI is either progress or phi
        progress_variable = self.progressvariable
        PHI = self.phis.select('*',"""WHERE direction = "%(direction)s" 
                                      ORDER BY median_%(progress_variable)s,mean_%(progress_variable)s""" \
                                   % vars())

        vmin,vmax = self.data.limits(progress_variable)
        kwargs.setdefault('aspect','auto')
        kwargs.setdefault('interpolation','nearest')
        kwargs.setdefault('origin','lower')
        kwargs.setdefault('cmap',Transitions.cm.jet2)
        kwargs['extent'] = (vmin,vmax,0,len(PHI))
        auto_vmax = kwargs.pop('auto_vmax',0.8)   # Switch off with 0 or None; then use vmax
        if auto_vmax and not 'vmax' in kwargs:
            kwargs['vmax'] = upperlimit(\
                numpy.array([a[0] for a in OBSERVABLES[direction].values()]),
                quantile=auto_vmax,nonzero=True,bins=50)
        with_colorbar = kwargs.pop('colorbar',False)
        with_lines = kwargs.pop('with_lines',True)        
        if with_lines:
            kwargs['hold'] = True
        lines_combined = kwargs.pop('lines_combined',True)
        fontsize = kwargs.pop('fontsize',8)
        title = kwargs.pop('title',None)

        sorted_dist_lists = {}
        distarrays = {}
        for xdirection,dist in OBSERVABLES.items():
            sorted_dist_lists[xdirection] = []
            for target in PHI.target:
                # sort all distributions along the selected direction
                sorted_dist_lists[xdirection].append(OBSERVABLES[xdirection][target][0])
            distarrays[xdirection] = numpy.array(sorted_dist_lists[xdirection])
        del sorted_dist_lists

        # now sort all phi (B) in the same way as phi[direction] (A)
        sorted_phis = {}
        for xdirection in self._directions:
            sorted_phis[xdirection] = self.phis.select(\
           """B.direction, B.target,
              B.median_DeltaRMSD, B.mean_DeltaRMSD, B.std_DeltaRMSD,
              B.median_progress, B.mean_progress, B.std_progress,
              B.median_phi, B.mean_phi, B.std_phi""",
           """AS A LEFT JOIN __self__ AS B 
                 ON (A.target = B.target 
                     AND A.direction = '%(direction)s'
                     AND B.direction = '%(xdirection)s')
                 WHERE B.direction NOTNULL AND B.target NOTNULL
              ORDER BY A.median_%(progress_variable)s,A.mean_%(progress_variable)s""" % vars())

        pylab.imshow(distarrays[direction],**kwargs)
        
        ax = pylab.gca()
        pylab.ylabel(r"crystal structure sorted by median(progress)")
        pylab.xlabel(r"progress %s" % self._variable_tex[self.progressvariable])
        def _idx2pdb(x,pos):
            try: return PHI.target[int(x)]
            except IndexError: return ''
        idx2pdbFormatter = pylab.matplotlib.ticker.FuncFormatter(_idx2pdb)
        integerLocator = pylab.matplotlib.ticker.IndexLocator(base=1,offset=0)
        ax.yaxis.set_major_locator(integerLocator)
        ax.yaxis.set_major_formatter(idx2pdbFormatter)
        major_labels = ax.get_yticklabels() 
        pylab.setp(major_labels,
                   horizontalalignment='right',
                   verticalalignment='bottom',
                   fontsize=fontsize)
        if title:
            pylab.title(title)
        if with_colorbar:
            cb = pylab.colorbar(aspect=30,extend='max')
            # manually clean up (bad interaction with hold/plot)
            cb.ax.set_xticks([])
            cb.ax.set_xticklabels('')
        if with_lines:
            ax.hold(True)
            median = "median_"+progress_variable     # plot median(progress) or median(phi)
            X = numpy.arange(len(PHI.target)) + 0.5  # shift to center of box
            if lines_combined:
                for d,phi in sorted_phis.items():
                    ax.plot(phi.field(median),X,'wo-',alpha=0.6,scalex=False,label=d)
            else:
                ax.plot(PHI.field(median),X,'wo-',alpha=0.6,scalex=False,label=direction)

            ax.hold(False)
        pylab.draw_if_interactive()

    def whiskerplot(self,direction,**kwargs):
        raise NotImplementedError
        kwargs.setdefault('mediancolor','k')     # for median in whiskerplot
        kwargs.setdefault('boxcolor','darkgray')

        # reformat
        for target, (edges,dist) in self.distribution[direction].items():
            pass # wait for SQL tables
        pylab.plot(x,ymean,'o-',color=kwargs['color'],ms=3)
        components = pylab.boxplot(all_y,positions=x,
                                   widths=logwidth(x,kwargs.setdefault('widths',0.15)))
        ax = pylab.gca()
        ax.set_xscale('log')
        ax.set_yscale('linear')
        # must modify components for customisation
        for l in components['boxes']:
            l.set_color(kwargs['boxcolor'])
        for l in components['whiskers']:
            l.set_color(kwargs['color'])
            l.set_linestyle('-')
        for l in components['caps']:
            l.set_color(kwargs['color'])            
        for l in components['medians']:
            l.set_color(kwargs['mediancolor'])
            l.set_linewidth(3.0)
        for l in components['fliers']:
            l.set_color(kwargs['color'])

        pylab.draw_if_interactive()


    def barplot(self,direction,alpha=0.4):
        import pylab
        pnormalize = pylab.normalize(vmin=1,vmax=len(self.distribution[direction]))
        count = 0
        # dx = 0.1/(2*len(self.distribution[direction]))
        for target,(edges,heights) in self.distribution[direction].items():
            count += 1
            lefts = edges[:-1]                # + count*dx
            widths = (edges[1:] - edges[:-1]) # - count*dx
            RGB = pylab.cm.jet(pnormalize(count))[:3]
            pylab.bar(lefts,heights,width=widths,
                      alpha=alpha,
                      color=RGB,
                      label="%s" % target)

    def __repr__(self):
        return "PDBMovie('%s')" % self.datafilename
        

def upperlimit(data,quantile=0.8,nonzero=True,bins=50):
    """Determine the upper value in the data that contains the lower quantile of the whole data."""
    data1d = numpy.ravel(data)
    if nonzero:
        data1d = data1d[data1d != 0]
    h,e = numpy.histogram(data1d,bins=bins,normed=True,new=True)
    m = 0.5*(e[1:]+e[:-1])
    cdf = numpy.cumsum(h * (e[1:]-e[:-1]))  # cumulative distribution
    return m[cdf >= quantile][0]

def vmd_load_targets_script(pdbmovie,direction,cmdfile='run_vmd'):
    import AdK.config as config
    import os
    modeldir = os.path.join(config.basedir,'Models','all')
    # psf = os.path.join(config.basedir, 'psf', 'open_final.psf')
    fmt = os.path.join(modeldir,'AKeco_%s.pdb')
    cmd = open(cmdfile,'w')
    try:
        # cmd.write("vmd %(psf)s \\\n" % vars())
        cmd.write("vmd  \\\n" % vars())
        cmd.write(" \\\n".join([fmt % p for p, in pdbmovie.get_ordered_targets(direction)]))
    finally:
        cmd.close()
    os.chmod(cmdfile,0755)
    print "Now run %(cmdfile)r." % vars()
    
