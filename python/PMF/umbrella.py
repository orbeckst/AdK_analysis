# $Id: umbrella.py 3173 2009-03-26 12:20:44Z oliver $
"""This module automates the process of setting up umbrella sampling
windows. It solves the following problem:

  Given a target windows position in (NMP,LID)-angle space, find a
  suitable starting frame from existing trajectories, extract the
  frame, prepare the input scripts to be submitted to a SGE queuing
  system, and update the appropriate meta data files (cons_values.txt
  and create a *.txt file for the reference values).

It requires an angle database (angles for all frames of trajectories
--- see PMF.angles) that includes paths to the dcds. This can be set
up with the

  db.add_dcds(...,include_all=True)

or with the add_dcds_to(db) function.

Typical usage:

  db = PMF.angles.setup()          # build/load database with all simulations
  PMF.umbrella.add_dcds_to(db)     # add all necessary metadata

  U = PMF.umbrella.UmbrellaSetup(db, NAME, forceconstant=3000, ...)

  U.extract_frame(55,120)          # extract pdb for NMP=55,LID=120 to pdb
  
  # extract a rectangular grid
  U.extract_grid(35,80,5,  95,140,5)

  # write required files
  U.write_autometa(XXX)   # FIXME: filename XXX?
  U.write_parameters()

"""

import os, errno
import numpy
from MDAnalysis import Universe

import config

# filenames hardcoded into the SGE submission script
DEFAULT_WINDOWS = 'windows.autometa'
DEFAULT_PARAMETERS = 'umbrella_parameters.pickle'

try:
    DCD_DEFAULT_PATTERNS = config.dcd_default_patterns
except AttributeError:
    raise ImportError("Sorry, config.py changed. Please do './bin/install.sh' again. :-*")

class UmbrellaSetup(object):
    """A class to automatically generate AdK umbrella windows from an AngleDB and dcds.

    The basic idea is to find starting configurations via the
    database, extract the closest frame from the original dcd, and
    write an entry into a 'autometa' file that can be used in turn to
    add the results to the database.

    The output is a directory that contains everything to submit an SGE array job.
    """

    def __init__(self,db,name,forceconstant=3000,force=False,
                 autoumbrella_dir=None,libdir=None,jobdir=None,
                 psf=None,sge=None,inp=None,forcefield=None):
        """Set up a class to automatically generate AdK umbrella windows from a AngleDB.

        U = UmbrellaSetup(db,name,autoumbrella_dir=DIR,psf=PSF,sge=[SGE1, SGE2, ...],
                          forcefield=[RTF, PRM, INP, ...])

        Arguments:

          db             AngleDB database (run add_dcds_to(db) first!)
          name           name of the project; all files will be stored under
                         <autoumbrella_dir>/<name>
        
        Keyword args:

          forceconstant       in kcal/mol / radians**2 (same for NMP and LID)
          force               False: fail if overwriting a project [default]
                              True:  clobber

          The following keywords should normally be left at the defaults:
          
          autoumbrella_dir    top directory
          libdir              directory where forcefields and other input is stored
          jobdir              directory under which the project directory 'name' is created
          psf                 Charmm psf file for the dcds
          sge                 LIST of sge scripts
          forcefield          LIST of all required force field files

        Bugs:
        * input is not checked; if it says LIST it must be a list!
        """
        
        self.db = db
        self.name = name
        self.forceconstant = forceconstant  # same for one setup
        self.autoumbrella_dir = autoumbrella_dir or os.path.join(config.basedir,'pmf','autoumbrella')
        self.libdir = libdir or self._top('lib')
        self.jobdir = jobdir or self._top('jobs')
        self.psf = psf or self._lib('open_final.psf')
        self.inp = inp or self._lib('pmf_angles.inp')
        self.sge = sge or [ #self._lib('pmf_autoumbrella.sge'),
            self._lib('pmf_angles_darthtater.sge'),
            self._lib('pmf_angles_deathspud.sge'),self._lib('pmf_angles_timberwulf.sge'),]
        self.forcefield = forcefield or \
            [self._lib('top_all22_prot.inp'), self._lib('par_all22_prot.inp'),
             self._lib('acepar22_prot.inp'),]

        self.copyfiles = [self.psf, self.inp] + self.sge + self.forcefield
        self.autometa = Autometa(self.name)   # each generated window's meta data is added

        self.project_dir = os.path.join(self.jobdir, name)
        # maybe force --> unlink dir first?
        try:
            os.makedirs(self.project_dir)
        except OSError,err:
            if err.errno == errno.EEXIST:
                if not force:
                    raise OSError('Project directory %(project_dir)r exists. Use force=True to overwrite files in it.' % vars(self))
            else:
                raise
        print "Project directory: %(project_dir)r" % vars(self)
        self.copy_libfiles()


    def _top(self,*args):
        return os.path.join(self.autoumbrella_dir,*args)

    def _lib(self,*args):
        return os.path.join(self.libdir,*args)

    def copy_libfiles(self):
        import shutil
        for f in self.copyfiles:
            shutil.copy2(f, self.project_dir)
        print "Added standard files (Charmm script, force field, ...) to the project directory."
        
    def write(self):
        """Write all files. Old files will be overwritten."""
        self.write_autometa()
        self.write_parameters()

    def write_autometa(self,filename=None):
        """Write all meta data about the generated frames to filename (overwriting previous content)

        Do NOT CHANGE the filename; the SGE script expects to find 'windows.autometa'.
        """
        if filename is None:
            filename = os.path.join(self.project_dir, DEFAULT_WINDOWS)
        am_filename = self.autometa.write(filename,with_linenumbers=True)
        print "Wrote autometa file %(am_filename)r." % vars()
        return am_filename


    def write_parameters(self,filename=None):
        """Unpickle this to get important project parameters.

        Do NOT CHANGE the filename; the SGE script expects to find 'umbrella_parameters.pickle'.
        """
        import cPickle
        if filename is None:
            filename = os.path.join(self.project_dir, DEFAULT_PARAMETERS)        
        # generate 'anonymous' keys
        def to_anon_dict(l,key):
            return dict( [(str(key)+'_'+str(i), val) for i,val in enumerate(l)] )
        # structure mimics what staging.JOB expects
        parameters = {'inputfiles':
                      {'psf': self.psf,
                       'inp': self.inp,
                       },
                      'outputfiles':
                      {},
                      'variables':
                      {'project_dir': self.project_dir,},
                      }
        parameters['inputfiles'].update(to_anon_dict(self.forcefield, 'FF'))
        cPickle.dump(parameters, open(filename,'w'))
        print "Wrote parameter pickle file %(filename)r." % vars()
        return filename
                      

    # copied & modified from PMF.angles
    def qrange(self,NMP=None,LID=None,asrecarray=False,verbose=False):
        """Find all frames for which the angles are within the given range.

        qrange(NMP,LID) --> list of dcd paths and frames

        Example:
          qrange(NMP=(None,80),LID=(90,95.5))

        Angle ranges are given as tuples; if one boundary is None then it is ignored.
        Boundaries are inclusive, ie LID_min <= angle <= LID_max.

        Results are ordered by distance from the centre of the desired
        range (but note that unset ranges are set to 0).
        """
        _operators = (">=", "<=")
        where_clauses = []
        centers = {'NMP': 0, 'LID': 0}
        def process_limits(column,limits):
            if limits is not None:
                where_clauses.extend(["%(column)s %(op)s %(val)r" % vars()
                                      for val,op in zip(limits,_operators) if val is not None])
                centers[column] = numpy.mean( [val for val in limits if val is not None] )
        process_limits('LID',LID)
        process_limits('NMP',NMP)
        where = " AND ".join(where_clauses)

        return self.db.SQL("""SELECT dcdpath,frame,NMP,LID,filename,
                                     distance(NMP,LID,?,?) AS distance
                             FROM __self__
                NATURAL LEFT JOIN __meta__
                            WHERE NOT dcdpath ISNULL AND """+\
                               where + " " +\
                      """ORDER BY distance ASC""",
                           centers['NMP'], centers['LID'],
                           verbose=verbose,asrecarray=asrecarray)
        
    def qaround(self,NMP=None,LID=None,delta=1.0,asrecarray=False,verbose=False):
        """Find all frames in the range [angle-delta, angle+delta].

        qaround(NMP=<angle>,LID=<angle>,delta=1.0) --> list of dcds and frames

        Results are ordered by distance from the desired coordinates
        (This only works meaningfully when both NMP and LID
        coordinates are supplied because unset coordinates are set to
        0.)
        """
        
        if LID is None and NMP is None:
            raise ValueError("At least one of LID or NMP should be provided.")

        DELTA = numpy.array([-delta,delta])
        LIDrange = NMPrange = None
        if LID:
            LIDrange = LID + DELTA
        if NMP:
            NMPrange = NMP + DELTA
        return self.qrange(LID=LIDrange,NMP=NMPrange,asrecarray=asrecarray,verbose=verbose)

    def qbest(self,NMP,LID,delta=1.0,**kwargs):
        """Return closest match (out of all matches within delta)."""
        x = self.qaround(NMP=NMP,LID=LID,delta=delta,**kwargs)
        if len(x) == 0:
            raise RuntimeError('qbest(): No match within delta=%(delta)g around NMP=%(NMP)g,LID=%(LID)g. Increase delta.' % vars())
        return x[0]

    def extract_frame(self,NMP,LID,**kwargs):
        """Extract best-matching frame NMP,LID from corresponding dcd.

        pdbpath = extract_frame(NMP,LID, forceconstant=self.forceconstant,**kwargs)

        The force constant (in kcal/mol/deg**2) can be changed from
        the default; this simply results in a new entry in the
        autometa file if the pdb file already exists.

        If no match can be found a RuntimeError is raised; try again
        with higher delta (see qbest()).

        Keyword args:

        forceconstant       angular force constant (for Charmm's MMFP GEO)
        delta               radius around NMP,LID in which the best match is searched
        
        """
        kwargs['asrecarray'] = True
        forceconstant = kwargs.pop('forceconstant',self.forceconstant)
        
        x = self.qbest(NMP,LID,**kwargs)
        # check for some strange crashes where x is a tuple: problem with generating
        # a numpy.recarray from a large record list
        assert type(x) == numpy.core.records.record
        print "Best match for NMP=%g LID=%g at dist=%.2f: frame %d of '%s'" \
              % (NMP,LID,x.distance,x.frame,x.dcdpath)
        universe = Universe(self.psf, x.dcdpath)
        pdbpath = extract_pdb_from_universe(universe,x.frame,directory=self.project_dir)
        pdbname = os.path.basename(pdbpath)
        print "  --> NMP=%g LID=%g  '%s'" % (x.NMP, x.LID, pdbname)
        self.autometa.append( (pdbname, forceconstant, NMP, LID) )
        print "  --> appended to meta list with forceconstant=%g" % forceconstant
        return pdbpath

    def extract_list(self,angles,**kwargs):
        """Extract all frames in list angles = [(NMP1,LID1), (NMP2,LID2), ...].

        pdbpaths = extract_list(angles,**kwargs)

        See extract_frame() for details.
        """
        pdbpaths = []
        nwindows = len(angles)
        print "\n========== Extracting %(nwindows)d frames ===========" % vars()
        print "(This can take a while. Do not forget to run U.write() afterwards.)"
        for i,(NMP,LID) in enumerate(angles):
            print "[%5.1f%%] window %4d out of %5d" % (100.0*(i+1)/nwindows, i+1, nwindows)
            p = self.extract_frame(NMP,LID,**kwargs)
            pdbpaths.append(p)
        print "============ extraction complete =================\n"\
              "(extract more frames and then run U.write() to write the metafile."
        return pdbpaths

    def extract_grid(self, NMPmin,NMPmax,NMPstep, LIDmin,LIDmax,LIDstep, **kwargs):
        """Extract all frames on the grid (end points included).

        pdbpaths = extract_grid(NMPmin,NMPmax,NMPstep, LIDmin,LIDmax,LIDstep)

        Arguments:
        delta       maximum radius in angle space for frames to be included (in degree)

        See extract_frame() and grid2d() for details. Note that you might need a
        larger delta value if you do not have many starting structure in the
        database. Smaller delta are better for performance.
        """
        return self.extract_list(grid2d(NMPmin,NMPmax,NMPstep, LIDmin,LIDmax,LIDstep),**kwargs)


def grid2d(xmin,xmax,xstep, ymin,ymax,ystep):
    """Returns a list of sample point on the specified grid."""
    # no idea how to do this nicely in numpy
    def generate(start,stop,step):
        assert start < stop
        assert step > 0
        num = int(numpy.round(float((stop-start))/step) + 1)
        return numpy.linspace(start,stop,num, endpoint=True)

    X = generate(xmin,xmax,xstep)
    Y = generate(ymin,ymax,ystep)
    return numpy.array( [(x,y) for x in X for y in Y] )

# use alias to make it more look like wham.Project    
Project = UmbrellaSetup
    
class Autometa(object):
    """Hold meta data for generated frames.

    AM = Autometa()
    AM.append( (pdbname, forceconstant,NMP_ref,LID_ref) )

    The input is checked:
    - A record must be a 5-tuple (or a ValueError is raised)
    - It must not have been inserted into the table previously
      (or a IntegrityError is raised)

    (The underlying data structure is an in-memory SQL table. It only
    allows appending and iterating.)
    """

    def __init__(self,name,**kwargs):
        """Initialize the autometa table. name=RUN_NAME is required"""        
        from pysqlite2 import dbapi2 as sqlite
        self.name = str(name)
        self.connection = sqlite.connect(":memory:")
        self.cursor = self.connection.cursor()
        self.cursor.execute("""CREATE TABLE autometa (idx INTEGER PRIMARY KEY, pdbname,forceconstant,NMP_ref,LID_ref)""")
        self.cursor.execute("""CREATE UNIQUE INDEX window ON autometa (pdbname,forceconstant,NMP_ref,LID_ref)""")
        self.connection.commit()

    def _checked(self,x):
        """Return the input tuple if it is valid input for a Autometa record."""
        try:
            pdbname,forceconstant,NMP_ref,LID_ref = x
        except (ValueError,TypeError):
            raise ValueError('Argument MUST be (pdbname,forceconstant,NMP_ref,LID_ref)')
        return x

    def append(self,x):
        self.cursor.execute(\
            """INSERT INTO autometa (pdbname,forceconstant,NMP_ref,LID_ref) VALUES (?,?,?,?)""",
            self._checked(x))
        self.connection.commit()
        
    def extend(self,iterable):
        self.cursor.executemany(\
            """INSERT INTO autometa (pdbname,forceconstant,NMP_ref,LID_ref) VALUES (?,?,?,?)""",
            [self._checked(x) for x in iterable])
        self.connection.commit()

    def tolist(self):
        return list(self)

    def write(self,filename,comment=None,with_linenumbers=True):
        """Write meta file; with_linenumbers=True enumerates the lines (SHOULD BE SET!)."""
        import util
        # definition of the format of an autometa file
        fmt = "%(pdbname)-30s %(forceconstant)10g   %(NMP_ref)7g  %(LID_ref)7g   %(dcdname)s\n"
        if with_linenumbers:
            fmt = "%(linenumber)4d " + fmt
        if comment is None:
            comment = ""
        else:
            comment = '# ' + str(comment) + '\n'

        filename = util.filename(filename,ext='autometa')
        autometa = open(filename,'w')
        try:
            autometa.write("# autometa file for autogenerated umbrella windows, written by\n"
                           "# $Id: umbrella.py 3173 2009-03-26 12:20:44Z oliver $\n" + \
                           comment + \
                           "#    NAME=" + str(self.name) + "\n"\
                           "# id pdbfilename                   forceconstant   NMP_ref  LID_ref   dcdname\n"
                           "#                                [kcal/mol/rad**2] [degree] [degree]         \n")
            for linenumber,(pdbname,forceconstant,NMP_ref,LID_ref) in enumerate(self):
                linenumber += 1
                # format must match the SGE script: run_number=$(printf "%03d" ${SGE_TASK_ID})
                dcdname = "%s_%03d.dcd" % (self.name, linenumber)                
                autometa.write(fmt % vars())
        finally:
            autometa.close()
        return filename

    def __iter__(self):
        for row in self.connection.execute("""SELECT pdbname,forceconstant,NMP_ref,LID_ref FROM autometa ORDER BY idx"""):
            yield row
    def __len__(self):
        self.cursor.execute('SELECT COUNT(*) FROM autometa')
        length, = self.cursor.fetchone()
        return length
    def __repr__(self):
        return str(self.tolist())
    def __del__(self):
        self.cursor.execute('DROP TABLE autometa')

def extract_pdb_from_universe(universe,frame,directory=os.path.curdir):
    """Extract frame from dcd in universe and write to a pdb.

    pdbpath = extract_pdb_from_universe(universe, framenumber)

    The filename of the pdb file is auto generated:

       DCDNAME_Fnnnn.pdb  where nnnnn is the frame number

    and saved under directory.

    Note that we assume Charmm numbering for trajectories (starting at 1).
    """

    trjname = universe.trajectory.filename
    pdbname = '%s_f%05d.pdb' % (os.path.splitext(os.path.basename(trjname))[0], frame)
    pdbpath = os.path.join(directory,pdbname)

    print "Extracting frame %(frame)d from dcd %(trjname)r into pdb %(pdbname)r." % vars()
    
    all = universe.selectAtoms("all")
    # advance to frame
    # IMPORTANT: MDAnalysis frames start at 0, Charmm and frame numbering here starts at 1
    assert frame > 0
    universe.trajectory[frame-1]
    try:
        all.write(pdbpath)
    except AttributeError:
        # Needs latest version of MDAnalysis devel >= r143 for the pdb frame writer
        raise ImportError("Please upgrade to MDAnalysis with revision >= 143.")
    return pdbpath


def add_dcds_to(db,*filepatterns):
    """Convenience function to add dcds in known places to a AngleDB

    add_dcds_to(db,'/path/to/other/*.dcd', '../yet/another/path/*.dcd')

    The default patterns are specific to our setup. The dcd repository
    is on greenwulf so it really makes only sense to run all this on
    greenwulf.

    Note that the blacklist is taken into account to avoid starting
    windows from perhaps distored conformations.

    See txt/trajectories.txt for details.

    $Id: umbrella.py 3173 2009-03-26 12:20:44Z oliver $
    """

    filepatterns = DCD_DEFAULT_PATTERNS + list(filepatterns)
    
    db.add_metadata(include_all=True)
    db.add_dcds(*filepatterns)
        
def show_missing(db):
    """Return all rows from __meta__ containing angle trajectories without matching dcdpath."""
    #return db.SQL("SELECT * FROM __meta__ WHERE dcdpath ISNULL",asrecarray=True)
    return db.SQL("""SELECT * FROM __meta__
          NATURAL LEFT JOIN (SELECT filename,trajectory FROM __self__ GROUP BY filename)
          WHERE dcdpath ISNULL""",
                  asrecarray=True)

def show_assigned(db):
    """Return all rows from __meta__ containing angle trajectories with matching dcdpath."""
    return db.SQL("""SELECT * FROM __meta__
          NATURAL LEFT JOIN (SELECT filename,trajectory FROM __self__ GROUP BY filename)
          WHERE NOT dcdpath ISNULL""",
                  asrecarray=True)


def show_notmeta(db):
    """Show filenames that are not in meta (blacklist already removed), so only due to bugs."""
    import angles, config
    blacklist = angles.BlacklistMap()
    
    records =  db.SQL("""SELECT filename FROM __self__  GROUP BY filename
                     EXCEPT
                     SELECT filename FROM __meta__""")
    result = [fn for fn, in records if not blacklist.match(fn)]
    return result


def enable_test_environment():
    """Switch DCD_DEFAULT_PATTERNS to DCD_FAKE_PATTERNS
    
    Testing on deathspud:
    ====================

       import PMF.umbrella
       from PMF.umbrella import enable_test_environment, add_dcds_to, UmbrellaSetup

       enable_test_environment()
       add_dcds_to(db, PMF.umbrella.fake_patterns)

       U = UmbrellaSetup(db,'XXTEST')

       # only have oc330 ...
       U.qbest(68.64,146.8)   # testing...

       # extract two frames
       U.extract_frame(55,120)
       U.extract_frame(68.64,146.8)

    Extracting NMP=55 LID=120: dist=0.04: frame 39 of u'/home/oliver/Projects/DIMS/AdK_external_fake/greenwulf/project/oc330.dcd'
    """

    global DCD_DEFAULT_PATTERNS, DCD_FAKE_PATTERNS, disabled_DCD_DEFAULT_PATTERNS
    
    DCD_FAKE_PATTERNS = config.dcd_testing_patterns
    disabled_DCD_DEFAULT_PATTERNS = DCD_DEFAULT_PATTERNS[:]
    DCD_DEFAULT_PATTERNS = DCD_FAKE_PATTERNS[:]

    print "Switched DCD_DEFAULT_PATTERNS to the ones for deathspud testing."
    

