# $Id: wham.py 4107 2010-04-27 12:15:54Z root $
"""Functions to create input suitable for Alan Grossfield's wham-2d 2.0.2

Usage
=====

You need a database of the umbrella-sampled angle trajectories (pickle
files). The table 'trajectories' must have been filled (with
db.add_metadata()); this table contains the meta data.

Make a wham project

  P = Project(db)

By default this will place files under pmf/wham. Then write out input for Alan
Grossfield's wham-2d code:

  P.write_meta_file()
  P.write_data_files()

Then run wham-2d with appropriate options.

  P.run_wham(**kwargs)

See the doc string of the function for all kwargs parameters. If
successful, it write the free energy file and immediately loads it
into a FreeEnergy object

  P.F

One can plot the free energy with

  P.F.plot(**kwargs)

 
Any ``free.dat`` file can be analyzed and plotted independently:

 F = FreeEnergy('free.dat')
 F.plot()

If energies were loaded, build an energy map using

 E = Energy(db, edegs=F.edges)
 

Example session
===============

   db = PMF.angles.setup(pmfonly=True)
   P = PMF.wham.Project(db)
   P.write_data_files(nequilibration=2000)
   P.write_meta_file()
   P.run_wham(binwidth=2.5)
   P.F.plot()
   PMF.angles.make_canonical_plot()
   savefig('figs/pmf/all_bw2_5_equil2k.png')
   
TODO
====
* obtain good boundaries (using coverage as criterion)


WHAM
====

For a good explanation of the WHAM algorithm see B. Roux's 1995 paper.

wham-2d assumes that the bias has the harmonic form

  w(xi) = k/2 * (xi - xi_0)**2

However, Charmm uses

  w(xi) = k' * (xi - xi_0)**2

thus,

  k = 2 * k'

(The values in the database are the Charmm values k' but we write the
standard k's to the meta file.)


File formats
============

metadata file
-------------

As with regular 1 dimensional wham, each line of the metadata file should
either be blank, begin with a '#', or have the following format:

   /path/to/timeseries/file loc_win_x loc_win_y spring_x spring_y [correl_time] [temperature]

 /path/to/timeseries/file  name of one of the time series files.

 loc win x and loc win y   are the locations x_0 of the minimum of the
                           biasing terms in the first and second dimensions
                           of the reaction coordinate

 spring x spring y         spring constants k in 1/2*k*(x-x_0)^2

 correl_time               correlation time of data (only used for MC bootstrapping)

 temperature               temperature of simulation; if not supplied the value from the
                           commandline is assumed


time series file
----------------

  # comments
  time  x(t) y(t)  [E_pot]


 time            is ignored
 x(t)            value of observable x
 y(t)            value of observable y

 E_pot           potential energy; must be supplied if the temperature was supplied in
                 the metadata file


free energy output
------------------

#X              Y               Free            +/-             Pro             +/-
30.500000       80.500000       inf     0.000000        0.000000        0.000000


"""

import os, errno
import re
try:
    # pysqlite2 was the one one used for development
    from pysqlite2 import dbapi2 as sqlite
except ImportError:
    # newer versions
    import sqlite3 as sqlite
import numpy
#import gzip
import angles,config
try:
    config.angles['reference_state']
except AttributeError:
    raise ImportError("Your config.py file is missing the 'angles'. Please setup your AdK project again:\n\n"
                      "  ./bin/install.sh\n\n"
                      "This should get you the new version. Sorry for the inconvenience. :-*")

RMolarGasConstant = 8.314472 # J mol-1 K-1   # http://physics.nist.gov/cgi-bin/cuu/Value?r
RMolarGasConstant_kcal = 0.239005736 * 1e-3 * RMolarGasConstant  # 1 joule = 0.239005736 calories

class Project(object):
    """This class collects the functionality to prepare and run wham-2d.

    Typically one sets up the project with a angle database

       P = Project(db)

    Then write files (required because wham-2d is run externally and cannot
    access the db directly):

       P.write_meta_file()
       P.write_data_files()


    Run wham-2d

       P.run_wham()

    and display results

       P.F.plot()

    All methods take arguments for tweaking, The final result is a free energy
    object, P.F.
    
    """

    # find equilibration steps in the header of a wham time series
    # '# nequilibration = 2000'
    #
    # '# length = 1000'
    # '# length = inf'     <--- if no length limit either inf or
    #  ---                 <--- absent (legacy)
    ts_nequilibration_re = re.compile(r'^\s*#\s*('
                                      '(nequilibration\s*=\s*(?P<NEQUIL>\d+))'
                                      '|'
                                      '(length\s*=\s*(?P<LENGTH>\d+|inf))'
                                      ')')
    # '# firstframe = 2000'
    # '# lastframe = 5000' --- NOT USED
    ts_framelimits_re = re.compile(r'^\s*#\s*('
                                   '(firstframe\s*=\s*(?P<FIRST>\d+))'
                                   '|'
                                   '(lastframe\s*=\s*(?P<LAST>\d+))'
                                   ')')
    
    def __init__(self,db,name=None,
                 topdir=os.path.join(config.basedir,'pmf','wham'),
                 wham_exe='wham-2d',temperature=300,
                 nequilibration=2000, length=None,
                 **limits):
        """Set up a wham-2d project from a AngleDB database.

        db            pmf-only database (must have trajectories table)
        name          name of the project; default contains nequilibration and temperature
        topdir        store all generated files under <topdir>/<name>
        wham_exe      path to Alan Grossfield's wham-2d
        temperature   temperature in Kelvin of the umbrella simulations
        nequilibration       number of steps to discard as equilibration; eg
                             look at the decay of the MMFP GEO energy
        length               use a stretch of length frames; None takes all [None]
        
        limits: 
           minNMP,maxNMP        any limit set will replace the ones derived from the maximum
           minLID,maxLID        ...extent of the data limits (**limits)

        CAVEATS:
          - If the users sets the name/path wrong and data filese are overwritten then wham
            will probably crash because it will run out of data.
          - nequilibration is really the frame number of the first frame used in wham
            (frame numbers start at 1 so really equilibration is nequilibration-1; this was
            an oversight which we are keeping to stay in line with the previous 
        """
        
        self.db = db              # AngleDB
        if name is None:
            if length is None:
                # standard old name
                name = 'wham_eq%(nequilibration)d_T%(temperature)gK'
            else:
                # new name with length
                name = 'wham_eq%(nequilibration)d_L%(length)d_T%(temperature)gK'
            self.name = str(name) % vars()
        self.topdir = topdir
        self.workdir = os.path.join(self.topdir, self.name)
        self.datadir = os.path.join(self.workdir,'timeseries')
        self.wham = wham_exe      # should check that wham-2d can be found
        self.temperature = temperature
        self.nequilibration = nequilibration
        self.length = length
        if self.length is None:
            self.length = float('inf')  # fudge so that we can always use first <= frame <= last
        self.firstframe = self.nequilibration # should be nequil-1 but we keep in line w/legacy 
        self.lastframe = self.firstframe + self.length  # frames are 1-based (so no -1)

        for d in (self.workdir,self.datadir):
            try:
                os.makedirs(d)
            except OSError,err:
                if err.errno != errno.EEXIST:
                    raise
        # check db
        try:
            db.numrows
        except:
            # try as filename and load
            import angles
            try:
                dbnew = angles.AngleDB(db)                
                dbnew.numrows
            except:
                raise ValueError('db '+str(db)+' is not a functional database')
            self.db = dbnew
        if self.db.numwindows == 0:
            raise RuntimeError("There are no umbrella windows defined in the db. "
                               "Run db.add_metadata() first.")

        # limits
        # use angleranges wherever input contains None
        # write a metafile  that excludes windows that are completely outside the limits
        (minNMP,maxNMP),(minLID,maxLID) = self.db.angleranges(**limits)
        self.limits = {'minNMP':minNMP, 'maxNMP':maxNMP, 'minLID':minLID, 'maxLID':maxLID}

        # add new auxiliary table to db (holds filenames of data files that contain
        # data within the desired limits AND first/lastframe)
        # (make name a valid SQL identifier by replacing decimal dots with _)
        self.whamfilesname = "first_%(firstframe)r__last_%(lastframe)r__" % vars(self) + \
                             "___".join([str(k)+"__"+str(v).replace('.','_')
                                         for k,v in self.limits.items()])
        # now SQL translates __whamfiles__
        self.SQL_transformation_rules = {'__whamfiles__': 'whamfilesname',}

        # meta file and free energies are stored in their own
        # directory; in this way the directory names and filenames
        # ensure that file are only overwritten if the input
        # parameters are the same
        self.metadir = os.path.join(self.workdir,'wham2d___'+self.whamfilesname)
        self.metafile = os.path.join(self.metadir,'wham2d.meta')
        try:
            os.makedirs(self.metadir)
        except OSError,err:
            if err.errno != errno.EEXIST:
                raise


        # List of data files that will be used as input to the wham-2d code:
        # (1) make sure that the table is clean:
        try:
            self.SQL("""DROP TABLE __whamfiles__""")
        except sqlite.OperationalError:
            pass     
        # (2) fill tables with names of data files which
        # - are umbrella windows (NOT __meta__.forceconstant ISNULL)
        # - have more than nequilibration data points
        # - have at least one data point within the limits (excluding equilibration!)
        #     check if at least one of the four corners is inside the limits
        #     lower left (0): minNMP,minLID; lower right (1): maxNMP,minLID
        #     upper left (2): minNMP,maxLID; upper right (3): maxNMP,maxLID
        #     point A inside limits: minlimit_i < A_i < maxlimit_i (AND over i=NMP,LID)
        # However, with this 'rectangular comparison' algorithm we can still include trajectories
        # that actually don't have points; these are filtered out in a second step.
        #        
        # (2a) Rectangular filter
        #
        # TODO: do I need __data__ here of can I get away with __self__ (what it was before)
        #       (athough may not make a huge difference, perhaps 5% longer)
        #
        self.SQL("""CREATE TEMP TABLE __whamfiles__ AS
                    SELECT TID, filename, COUNT(filename) AS numdata,
                           min(NMP) AS minNMP, max(NMP) AS maxNMP,
                           min(LID) AS minLID, max(LID) AS maxLID 
                      FROM __data__
         NATURAL LEFT JOIN __meta__
                     WHERE frame BETWEEN ? AND ?
                           AND forceconstant NOTNULL
                  GROUP BY filename
                    HAVING numdata > 0
                       AND ( 
                            ( (minNMP BETWEEN %(minNMP)g AND %(maxNMP)g) AND
                              (minLID BETWEEN %(minLID)g AND %(maxLID)g) )
                        OR  ( (maxNMP BETWEEN %(minNMP)g AND %(maxNMP)g) AND
                              (minLID BETWEEN %(minLID)g AND %(maxLID)g) )
                        OR  ( (minNMP BETWEEN %(minNMP)g AND %(maxNMP)g) AND
                              (maxLID BETWEEN %(minLID)g AND %(maxLID)g) )
                        OR  ( (maxNMP BETWEEN %(minNMP)g AND %(maxNMP)g) AND
                              (maxLID BETWEEN %(minLID)g AND %(maxLID)g) )
                           )
                    """ % self.limits, self.firstframe, self.lastframe)

        # (2b) Fine grained frame filter
        # Now check trajectories individually
        # and remove trajectories any from __whamfiles__ that have no
        # data points in the limits
        # TODO: THIS IS VERY SLOW WITH BIG DBs.
        for filename in [x[0] for x in self.SQL("""SELECT filename FROM __whamfiles__""")]:
            # check whole trajectory in limits
            result = self.SQL("""SELECT COUNT(*)
                          FROM __data__
                         WHERE filename = ?
                               AND frame BETWEEN ? AND ?
                               AND NMP BETWEEN ? AND ?
                               AND LID BETWEEN ? AND ?
                           """,
                              filename, self.firstframe, self.lastframe,
                              self.limits['minNMP'],self.limits['maxNMP'],
                              self.limits['minLID'],self.limits['maxLID'])
            N_datainlimits, = result[0]
            if N_datainlimits == 0:
                self.SQL("""DELETE FROM __whamfiles__ WHERE filename = ?""", filename)
                # print "DEBUG  deleted %r from __whamfiles__" % filename

        # (3) remember filenames (for convenience)
        self.filenames = self.SQL("""SELECT filename FROM __whamfiles__ ORDER BY filename ASC""",
                                     asrecarray=True).filename


    def _transform_SQL(self,SQL):
        """Replaces special strings such as __self__ with table or view name;
           rules are in self.SQL_transformation_rules"""
        # ugly kludge: copied, pasted & modified
        # angles.BaseSelection's _transform_SQL and SQL in order to
        # have transparent replacement of __whamfiles__ on a
        # per-Project instance basis (if we only add it to
        # db.SQL_transformation_rules then instances will overwrite
        # each other) (Maybe make Project also a BaseSelection?)
        _SQL = SQL
        for target,replacement in self.SQL_transformation_rules.items():
            try:
                _SQL = _SQL.replace(target, str(self.__getattribute__(replacement)))
            except (AttributeError,TypeError):
                pass
        return _SQL
    
    def SQL(self,SQL,*sqlvalues,**kwargs):
        """Execute SQL statement and return results.

        SQL(SQL_statement, *values)

        SQL_statement is a legal SQL statement (see
        http://www.sqlite.org/lang.html for the SQL syntax understood by the
        sqlite3 database engine). It can contain the '?' place holder. In this
        case, *values will be interpolated into the SQL string using the
        dbapi.

        Note: __self__ is replaced by the table name.
              __meta__ is replaced by the name of the auxiliary trajectories table.
              __whamfiles__ is replaced by the filenames of the data files which contain
                            data within the limits; this depends on each Project instance

        :Arguments:
        *values            one or more values that are referenced with '?' in the SQL
        verbose            show sql statement [False]
        asrecarray         return a numpy rec array [False]
        """
        # db function replaces __self__ and __meta__, we only do __whamfiles__
        return self.db.SQL(self._transform_SQL(SQL),*sqlvalues,**kwargs)

    def write(self,**kwargs):
        """Write meta and data files with default settings."""
        self.write_meta_file(**kwargs)
        self.write_data_files()

    def write_meta_file(self,absolute_path=False):
        """Write the metadata file, using the trajectories table in the database.

        absolute_path         True: write full path to file
                              False: only include datadir name in path

        NOTE: 
        - The force constants in Charmm are in units 'kcal*mol**-1*radian**-2' but
          the wham code needs 'kcal*mol**-1*degree**-2'.
        """
        from numpy import pi
        k_rad2deg = 2 * (pi/180.)**2      # conversion factor 1/rad**2 --> 1/deg**2, 2x for Charmm

        meta = open(self.metafile,'w')
        # Comment lines can't be too long or wham-2d crashes with weird error messages!
        meta.write("# Input file for Alan Grossfield's wham-2d\n"
                   "# written by $Id: wham.py 4107 2010-04-27 12:15:54Z root $\n"
                   "# AngleDB: %(filename)r\n" % vars(self.db) +
                   "# nequilibration: %(nequilibration)d\n"
                   "# length:         %(length)r\n"
                   "# firstframe:     %(firstframe)r\n"
                   "# lastframe:      %(lastframe)r\n" % vars(self) + 
                   "# limits: NMP %(minNMP)g %(maxNMP)g   LID %(minLID)g %(maxLID)g\n" % self.limits + 
                   "# timeseries\t\t\t\tNMP_ref  LID_ref\t NMP_k  LID_k\n")

        SQL = """SELECT filename,NMP_ref,LID_ref,
                        forceconstant * %(k_rad2deg)f AS NMP_k,
                        forceconstant * %(k_rad2deg)f AS LID_k
                   FROM __meta__ AS M
                   NATURAL LEFT JOIN __whamfiles__ AS W  WHERE W.filename NOT NULL
                   ORDER BY NMP_ref,LID_ref ASC""" % vars()
        try:
            for (p,NMP_ref,LID_ref,NMP_k,LID_k) in self.SQL(SQL):
                fn = self.create_timeseries_path(p,absolute_path=absolute_path)
                meta.write("%(fn)s\t %(NMP_ref)g  %(LID_ref)g\t  %(NMP_k)g  %(LID_k)g\n" % vars())
        finally:
            meta.close()
        print "Wrote meta file '%(metafile)s'." % vars(self)
        
    def write_data_files(self,gzip=False):
        """Write out all time series that appear in the meta file.

        If a time series already exists AND has been produced for the same
        nequilibration then the file is not written but the existing one is
        kept. This speads up re-analysis tremendously.
        """

        nfiles = len(self.filenames)    # only take files that contain data within limits
        nequilibration = self.nequilibration   # aliased for vars()
        length = self.length
        for ifile,filename in enumerate(self.filenames):
            fn = self.create_timeseries_path(filename)
            ifile += 1
            progress = 100.0*ifile/nfiles            
            if self.timeseries_exists(fn):
                print "[%(progress)5.2f%% %(ifile)4d/%(nfiles)d] Kept timeseries file '%(fn)s'." % vars()
                continue
            timeseries = open(fn,'w')
            try:
                timeseries.write("""# NMP-LID timeseries, written by $Id: wham.py 4107 2010-04-27 12:15:54Z root $\n"""
                                 """# pickle = "%(filename)s"\n"""
                                 """# nequilibration = %(nequilibration)d\n"""
                                 """# length = %(length)r\n""" % vars())
                SQL = """SELECT frame,NMP,LID FROM __self__
                          WHERE filename = ?
                                AND frame BETWEEN ? AND ?
                       ORDER BY frame ASC"""
                for t,NMP,LID in self.db.SQL(SQL,filename,self.firstframe,self.lastframe):
                    # SQL does dbAPI interpolation and it is MUCH faster to slurp a whole
                    # data set instead of using the db iterator over single rows (x100?)
                    timeseries.write("%(t).1f \t%(NMP)f \t%(LID)f\n" % vars())
            finally:
                timeseries.close()
            t -= nequilibration
            print "[%(progress)5.2f%% %(ifile)4d/%(nfiles)d] Wrote timeseries file '%(fn)s' with %(t)d entries." % vars()

    def timeseries_exists(self,fn):
        """Return True if a time series of the desired length exists."""
        try:
            ts = open(fn)
        except IOError:
            return False
        # find '# nequilibration = 2000'
        nequilibration = None
        length = float('inf')  # default/legacy is length = inf
        lineno = 0
        maxlines = 100  # only read up to that many header lines
        try:
            for line in ts:
                lineno += 1
                m = self.ts_nequilibration_re.match(line)
                if m:
                    # set nequilibration and length if possible
                    try:
                        nequilibration = int(m.group('NEQUIL'))
                    except TypeError:
                        pass
                    try:
                        length = float(m.group('LENGTH'))  # only float can get inf
                        try:
                            length = int(length)
                        except OverflowError:
                            pass   # we got inf, alright
                    except TypeError:
                        pass                        
                if lineno > maxlines:
                    break
        finally:
            ts.close()
        return nequilibration == self.nequilibration and length == self.length

    def create_timeseries_name(self,p):
        return os.path.splitext(os.path.basename(p))[0] + '.dat'

    def create_timeseries_path(self,p,absolute_path=True):
        fn = self.create_timeseries_name(p)
        if absolute_path:
            return os.path.join(self.datadir, fn)
        else:
            return os.path.join(os.path.basename(self.datadir), fn)

    def run_wham(self,binwidth=1,tolerance=0.001,freefile=None,
                 temperature=None,metafile=None,padding=0):
        """Run the external wham-2d program and set up a FreeEnergy object.

        With default arguments the umbrella sampled data is unbiased across the
        whole domain. Set arguments such as binwidth and tolerance to influence
        the run.

        wham-2d must be in your shell's PATH; if wham-2d is not installed
        download it from Alan Grossfield's web page at
        http://membrane.urmc.rochester.edu/Software/WHAM/WHAM.html

        #NOTE: By default we are writing gzipped data files. Either
        #unzip files manually or use a patched version of wham-2d that
        #can read gzipped files. (Ask Oliver <orbeckst@gmail.com>.)
        # (NOT YET)

        Keyword arguments:
        
           freefile             result file with free energies in kcal/mol
           binwidth             width of a bin; the bin numbers in each directions
                                are calculated from the domain and the binwidth
           tolerance            iterate to self consistency until the free energy weights
                                change by less than tolerance kT
           temperature          temperature in K at which simualtions where run [300 K]

        """
        from subprocess import call

        if freefile is None:
            _name = 'free_bw%(binwidth)g_tol%(tolerance)g.dat' % vars()
            freefile = os.path.join(self.metadir, _name)
        if metafile is None:
            metafile = self.metafile
        if temperature is None:
            temperature = self.temperature

        # limits from self
        minNMP,maxNMP,minLID,maxLID = [self.limits[k] for k in ('minNMP','maxNMP','minLID','maxLID')]

        nbins_NMP = numpy.ceil((maxNMP-minNMP)/binwidth)
        nbins_LID = numpy.ceil((maxLID-minLID)/binwidth)

        self.whamargs = {'freefile':freefile,'metafile':metafile,
                         'temperature':temperature,
                         'nbins_NMP': nbins_NMP, 'nbins_LID':nbins_LID,
                         'tolerance': tolerance,
                         'limits':self.limits,
                         }                         

        args = "Px=0 %(minNMP)g %(maxNMP)g %(nbins_NMP)d  "\
               "Py=0 %(minLID)g %(maxLID)g %(nbins_LID)d  "\
               "%(tolerance)g %(temperature)g %(padding)d %(metafile)s %(freefile)s" % vars()
        wham_command = self.wham + " " + args
        print "run_wham(): "+ wham_command
        oldir = os.getcwd()
        try:
            os.chdir(self.workdir)
            retcode = call(wham_command, shell=True)
            if retcode == 0:
                pass
            elif retcode < 0:
                raise RuntimeError("run_wham(): wham was terminated by signal %d" % -retcode)
            else:
                raise RuntimeError("run_wham(): FAILURE, wham returned %d" % retcode)
#         except OSError, e:
#             raise OSError("Execution of wham failed:", e)
        finally:
            os.chdir(oldir)

        print "run_wham(): generated output free energy file '%(freefile)s'." %vars()

        whamargs = dict([(k,v) for k,v in self.whamargs.items() if k != 'freefile'])
        self.F = FreeEnergy(freefile,**whamargs)
        print "run_wham(): added free energy object self.F"

    def __repr__(self):
        return "<PMF.wham.Project with "+str(len(self.filenames))+ \
               " angle timeseries within limits "+str(self.limits)+">"


class FreeEnergy(object):
    """Load and analyze the output of wham-2d.

      F = FreeEnergy('pmf/wham/free.dat')
      F.plot()

    The 2D energy landscape can be accessed as ``F.free_energy`` and the X
    (NMP) and Y (LID) values in ``F.X`` and ``F.Y``.  (The points are stored in
    C-order, i.e. the last dimension varying fastest.)

    Free_energy, probability, X, and Y can be restricted to a subdomain by
    using the set_datalimits() method. The original values are always stored
    in attributes with underscore prepended, eg _X, _Y, _free_energy.

    Interpolation functions ``W(NMP,LID)`` and ``P(NMP,LID)`` provided cubic
    splines on the defined data range.
    """
    
    def __init__(self,freefile,reference=config.angles['reference_state'],**whamargs):
        """Load wham-2d output file and shift energy to reference state.        

        freefile    free energy file produced by wham-2d
        reference   (NMP,LID) coordinates of state that is used as the reference state;
        **whamargs  the arguments that were passed to wham-2d (or rather run_wham)
                    (eg temperature=300)

        """
        self.filename = freefile
        self.whamargs = whamargs   # parameters for the wham-2d run
        self.temperature = whamargs.setdefault('temperature',300.0)
        
        self._init_read_free()    # part of __init__, do not call on its own
        self._init_datalimits()   # call set_datalimits() with args to change
        
        # shift energy relative to reference state
        try:
            self.reference_state = reference
            self.reference = self.get_reference_energy(*reference)
        except IndexError:
            import warnings
            self.reference_state = None
            self.reference = 0
            warnings.warn("Cannot find reference state in the map. "
                          "Free energy will not be shifted at all.")
        self._free_energy -= self.reference
        self._init_interpolation_functions()  # have everything to set up self.W and self.P        

    def _init_read_free(self):
        free = open(self.filename,'r')
        _X = []
        _Y = []
        _FE = []
        try:
            for line in free:
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                _x,_y,_F,DF,P,DP = map(float, line.split())
                _X.append(_x)
                _Y.append(_y)
                _FE.append(_F)
        finally:
            free.close()

        # midpoints of bins from (x,y) coordinates
        self._midpoints = (numpy.unique(_X),   # original NMP midpoints
                           numpy.unique(_Y))   # original LID midpoints
        self._X = self._midpoints[0]
        self._Y = self._midpoints[1]
                           
        nx,ny = len(self._X), len(self._Y)        
        F = numpy.array(_FE,order='C').reshape((nx, ny))
        # original free energy
        self._free_energy = numpy.ma.array(F, mask=numpy.logical_not(numpy.isfinite(F)), fill_value=1000);

        # reconstruct what input histogram2d would need
        self._edges = (self._mid2edges(self._X),    # NMP bin edges
                       self._mid2edges(self._Y))    # LID bin edges


    def _mid2edges(self,m):
        # construct 2 additional points before first and and after last midpoint
        # so that m[0] and m[-1] will end up at the centre of the first and last bin
        # (simply subtract/add first/last binwidth: m' = m[0] - (m[1]-m[0])
        x = numpy.concatenate(([2*m[0]-m[1]], m, [2*m[-1]-m[-2]]))
        return 0.5*(x[1:] + x[:-1])  # edges are the midpoints of the extended midpoints

    def _init_interpolation_functions(self):
        # need to be called when datalimits change
        try:
            self.W = self.interpolationFunction(self.free_energy)  # interpolation function for PMF
            #self.P = self.interpolationFunction(self.probability,cval=0)  # ... and for probability
            #self.P.__doc__ = "B-spline function over the probability density, p(NMP,LID)"
            self.P = self.interpolatedProbability()
        except AttributeError:
            raise AttributeError('Still missing data attributes. Run __init__ properly.')

    def _init_datalimits(self):
        """Setup basic datalimits structures but not functions etc."""
        try:
            self.set_datalimits()
        except AttributeError:
            pass

    def set_datalimits(self,minNMP=None,maxNMP=None,minLID=None,maxLID=None):
        """Restrict the range of data that are actually used when self.free_energy and probability are used"""
        # unelegant copy&paste crap:
        if minNMP is None:
            minNMP = self._X[0]
        if maxNMP is None:
            maxNMP = self._X[-1]
        if minLID is None:
            minLID = self._Y[0]
        if maxLID is None:
            maxLID = self._Y[-1]

        limits = numpy.array([[minNMP, maxNMP], [minLID, maxLID]])
        datalimits = []
        # find index
        for d,(dmin,dmax) in enumerate(limits):
            e = self._edges[d]                   # |...|...|
            amin = e.searchsorted(dmin,'left')  - 1   # different treatment of outliers left/right
            amax = e.searchsorted(dmax,'right') - 1   #
            assert amin >= 0
            assert amax < len(self._midpoints[d])
            datalimits.append(slice(amin,amax+1))
        self.datalimits = tuple(datalimits)

        # complicated things that need updating when limits change;        
        # they can fail and raise an AttributeError (which is ignored
        # when run from _init_datalimits)
        self._init_interpolation_functions()


    def get_reference_energy(self,refNMP,refLID):
        """Returns the value of the free energy at the reference point."""

        # find bin of ref: the one with a single count
        ref,e1,e2 = numpy.histogram2d([refNMP],[refLID], bins=self.edges, normed=False)        
        return self.free_energy[ref == 1][0]

    def plot(self,title=r'Free energy in kcal/mol',NMPlimits=None,LIDlimits=None,
             min_F=None, max_F=12, use_spline=True, resolution=0.5, with_colorbar=True,
             **kwargs):
        """Make a 2D contour plot of the free energy landscape.

        title         plot title
        NMPlimits     (xmin,xmax) tuple in degrees, changing the view in the x-dimension
        LIDlimits     (ymin,ymax)
        min_F         lower cutoff for the plot at this height value in kcal/mol (unit of wham-2d)
        max_F         upper cutoff for the plot at this height value in kcal/mol (unit of wham-2d)
        use_spline    use spline-interpolated data instead of original unbiased histogram [True]
        resolution    when using splines, resolution in degrees along each axis [0.5]
        with_colorbar [True]
        **kwargs      arguments for pylab.plot
        """
        import pylab

        if use_spline:
            _NMP,_LID = numpy.mgrid[self.NMP.min():self.NMP.max():resolution,
                                  self.LID.min():self.LID.max():resolution]
            NMP,LID = numpy.ogrid[self.NMP.min():self.NMP.max():resolution,
                                  self.LID.min():self.LID.max():resolution]
            NMP,LID = NMP.ravel(), LID.ravel()
            W = self.W(_NMP,_LID)
        else:
            NMP = self.NMP
            LID = self.LID
            W = self.free_energy

        kwargs.setdefault('cmap',pylab.cm.jet)
        if min_F is None:
            min_F = -self.reference
        kwargs.setdefault('vmin',min_F)
        kwargs.setdefault('vmax',max_F)
        kwargs['origin'] = 'lower'
        kwargs['aspect'] = 'equal'
        kwargs['extent'] = (NMP.min(),NMP.max(), LID.min(),LID.max())
        ncontours = kwargs.pop('ncountours',100)

        pylab.imshow(W.T, **kwargs)
        # imshow is faster/less complicated. contourf not fully working
        #kwargs.setdefault('extend','both')
        #pylab.contourf(NMP,LID, W.T, numpy.linspace(min_F,max_F,ncontours),**kwargs)
        if with_colorbar:
            pylab.colorbar(extend='max')  # ticks=numpy.arange(0,max_F)

        ax = pylab.gca()
        #ax.set_aspect(kwargs['aspect'])  # when using contourf
        degreeFormatter = pylab.matplotlib.ticker.FormatStrFormatter(r'%d$^\circ$')
        ax.xaxis.set_major_formatter(degreeFormatter)
        ax.yaxis.set_major_formatter(degreeFormatter)
        
        pylab.contour(NMP,LID,W.T,numpy.arange(min_F,max_F,0.5),
                      colors='k',linewidths=0.2,alpha=0.5)        
        pylab.contour(NMP,LID,W.T,numpy.arange(min_F,max_F,1),
                      colors='k',linewidths=0.5)
        
        pylab.xlabel(r'angle NMP-CORE')
        pylab.ylabel(r'angle LID-CORE')

        if NMPlimits:
            pylab.xlim(NMPlimits)
        if LIDlimits:
            pylab.xlim(LIDlimits)

        if title:
            pylab.title(title)

    def interpolationFunction(self,data=None,spline_order=3,cval=None):
        """Returns a function W(x,y) that interpolates any values on the PMF.

        W = interpolationFunction(self,data=None,spline_order=3,cval=None)

        cval is set to data.max() by default (works for the PMF) but for the
        probability this should be 0 or data.min(). cval cannot be chosen too large
        or too small or NaN because otherwise the spline interpolation breaks down
        near the region and produces wild oscillations.
        """
        # see http://www.scipy.org/Cookbook/Interpolation
        from scipy import ndimage

        if data is None:
            data = self.free_energy

        if cval is None:
            cval = data.max()
        _data = data.filled(cval)   # fill with max, not 1000 to keep spline happy

        coeffs = ndimage.spline_filter(_data,order=spline_order)
        x0,y0 = self.X[0], self.Y[0]
        dx,dy = self.X[1] - x0, self.Y[1] - y0
        def _transform(cnew, c0, dc):
            return (numpy.atleast_1d(cnew) - c0)/dc
        def interpolatedF(NMP,LID):
            """B-spline function over the PMF, W(NMP,LID).

            Example usage:
            >>> XX,YY = numpy.mgrid[40:75:0.5,96:150:0.5]
            >>> FF = FreeEnergy.W(XX,YY)
            >>> contourf(XX,YY,FF, numpy.linspace(-3,11,100),extend='both')
            """
            return ndimage.map_coordinates(coeffs,
                                           numpy.array([_transform(NMP,x0,dx),_transform(LID,y0,dy)]),
                                           prefilter=False,mode='constant',cval=cval)
        return interpolatedF            

    def X():
        doc = """NMP coordinates of bins (within data limits)"""
        def fget(self):
            return self._X[self.datalimits[0]]
        return locals()
    X = property(**X())
    NMP = X

    def Y():
        doc = """LID coordinates of bins (within data limits)"""
        def fget(self):
            return self._Y[self.datalimits[1]]
        return locals()
    Y = property(**Y())
    LID = Y

    def edges():
        doc = "edges that can be used for histogram2d (note: orig in _edges); assigning sets limits"
        def fget(self):
            return (self._mid2edges(self.X), self._mid2edges(self.Y))
        def fset(self,edges):
            raise NotImplemented("Use set_datalimits() to change the limits for now")
        return locals()
    edges = property(**edges())

    def limits():
        doc = "lower and upper limits of the data, array([[min0,min1], [max0,max1]]), from edges"
        def fget(self):
            e = self.edges
            return numpy.array([ map(numpy.min, e), map(numpy.max, e) ])
        return locals()
    limits = property(**limits())

    def limits_dict():
        doc = "lower and upper limits of the data, oldstyle dict, from edges"
        def fget(self):
            lim = self.limits
            return {'minNMP':lim[0,0],'maxNMP':lim[1,0],
                    'minLID':lim[0,1],'maxLID':lim[1,1]}
        return locals()
    limits_dict = property(**limits_dict())

    def free_energy():
        doc = "Free energy array, restricted to defined data rectangle (see set_datalimits)"
        def fget(self):
            return self._free_energy[self.datalimits]
        return locals()
    free_energy = property(**free_energy())


    def probability():
        doc = """probability density (in degree**(-2), derived from PMF (in kcal/mol)

        P(NMP,LID) = exp(-W(NMP/LID)/RT) / dNMPdLID

        where dNMPdLID is the area element.
        """
        def fget(self):
            X,Y = self.edges
            dXdY = numpy.outer(X[1:] - X[:-1],Y[1:] - Y[:-1])
            RT = RMolarGasConstant_kcal * self.temperature
            return numpy.exp(-self.free_energy/RT)/dXdY
        return locals()
    probability = property(**probability())

    def interpolatedProbability(self):
        eX,eY = self.edges
        # calculating of density is a bit iffy: I effectively assume
        # that all original bins have the same area
        dXdY = numpy.outer(eX[1:] - eX[:-1], eY[1:] - eY[:-1]).mean()
        def P(X,Y,temperature=self.temperature):
            """Probability derived from the spline-interpolated PMF."""
            RT = RMolarGasConstant_kcal * temperature        
            return numpy.exp(-self.W(X,Y)/RT)/dXdY
        return P
            

class FreeEnergyContainer(object):
    """Wrap free energy array and edges into this class and provide it
    as the 'filename' argument to the DerivedFreeEnergy class."""
    def __init__(self,F,X,Y):
        self.F = F   # array
        self.X = X   # NMP midpoints
        self.Y = Y   # LID midpoints
        
class DerivedFreeEnergy(FreeEnergy):
    def _init_read_free(self):
        # subverted... interpret self.filename as a numpy array
        try:
            FEC = self.filename  # FreeEnergyContainer?
            F = FEC.F
        except AttributeError:
            super(DerivedFreeEnergy,self)._init_read_free()
            return
        self._midpoints = (FEC.X, FEC.Y)
        self._X = FEC.X
        self._Y = FEC.Y
        self._free_energy = numpy.ma.array(F, mask=numpy.logical_not(numpy.isfinite(F)),
                                           fill_value=1000);
        # reconstruct what input histogram2d would need
        self._edges = (self._mid2edges(self._X),    # NMP bin edges
                       self._mid2edges(self._Y))    # LID bin edges


class BlockScanner(object):
    """Create WHAM Projects with frames in blocks.

    Note that you need to know how many blocks fit into the windows;
    input values are used without check against the data.
    """
    def __init__(self,db,**kwargs):
        """BS = BlockScanner(db, nblocks=5, lenblock=1000, stride=lenblock, offset=0)

        firstframes = stride * numpy.arange(nblocks) + offset
        """
        self.db = db
        self.nblocks = kwargs.pop('nblocks',5)
        self.lenblock = kwargs.pop('lenblock',1000)
        self.stride = kwargs.pop('stride',self.lenblock)
        self.offset = kwargs.pop('offset',0)

        self.projects = {}      # index by starting frame

        kwargs['length'] = self.lenblock
        firstframes = self.stride * numpy.arange(self.nblocks) + self.offset
        for nequil in firstframes:
            print "-- init(): starting frame = %d + %d" % (nequil, self.lenblock)
            kwargs['nequilibration'] = nequil
            self.projects[nequil] = Project(db,**kwargs)

        self.firstframes = self.projects.keys()
        self.firstframes.sort()
        
    def sorted_projects(self,msg=None,firstframes=None):
        if firstframes is None:
            firstframes = self.firstframes
        for i in firstframes:
            if msg is not None:
                print "-- %s(): starting frame = %d + %d" % (msg, i, self.lenblock)
            yield i,self.projects[i]

    def sorted_F(self,**kwargs):
        for nequil,P in self.sorted_projects(**kwargs):
            yield nequil,P.F
    
    def write(self):
        for nequil,P in self.sorted_projects(msg='write'):
            P.write()

    def run_wham(self,**kwargs):
        kwargs.pop('freefile',None)
        for nequil,P in self.sorted_projects(msg='run_wham'):
            P.run_wham(**kwargs)
        self.aggregate(**kwargs)

    def aggregate(self,start=None,stop=None,step=None, **whamargs):
        """Combine blocks to get free_energy and free_energy_std.

        Use the slice argument to select blocks.
        """
        firstframes = self.firstframes[slice(start,stop,step)]
        self.currentblocks = firstframes
        Fmean = numpy.mean([F.free_energy
                            for nequil,F in self.sorted_F(firstframes=firstframes)], axis=0)
        Fstd = numpy.std([F.free_energy
                          for nequil,F in self.sorted_F(firstframes=firstframes)], axis=0)

        # hackish... fake the free file to use the FreeEnergy class:
        nequil,F = self.sorted_F().next()  # just need one for the X/Y points; all the same
        FEC = FreeEnergyContainer(Fmean, F.X, F.Y)            # hack (1)
        self.free_energy = DerivedFreeEnergy(FEC,**whamargs)  # hack (2)
        FEC = FreeEnergyContainer(Fstd, F.X, F.Y)
        self.free_energy_std = DerivedFreeEnergy(FEC,**whamargs)
        
    def plot(self,**kwargs):
        import pylab
        ncols = 2 
        nrows = numpy.ceil(len(self.projects)/2.0)   # UGLY!! if we have many plots....
        kwargs['with_colorbar'] = False
        for iplot,(nequil,F) in enumerate(self.sorted_F(msg='plot')):
            pylab.subplot(nrows,ncols,iplot+1)
            kwargs['title'] = "frames %d --> %d" % (nequil, nequil + self.lenblock-1)
            F.plot(**kwargs)

    def plot_mean(self,**kwargs):
        # make it use FreeEnergy ...
        self.free_energy.plot(**kwargs)

    def plot_std(self,**kwargs):
        self.free_energy_std.plot(**kwargs)
        
