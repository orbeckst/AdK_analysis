# $Id: observables.py 4107 2010-04-27 12:15:54Z root $
# -*- coding: latin-1 -*-
"""Functions to analyze observables in the angle db.

Modelled after the energy functions in wham.py.

Observables are classes.

They are typically instantiated with a database or a selection as the first
argument. In some cases (for 'derived variables') the appropriate tables in the
database have to be computed first. See the ``analysis`` module for the
appropriate Creator classes.

All observables are held in the observables registry
``OBSERVABLES``. Typically, one simply instantiates a ObservablesGroup and uses
the methods of the group to analyse and plot the results.

Example
-------

 db = PMF.angles.AnglesDB('db/pmf.db')

 FRET = PMF.observables.ObservableGroup('FRET',db,'pmf',binwidth=0.5)
 FRET.plot()

 F = PMF.wham.FreeEnergy('.../free.dat')
 ProbFRET = PMF.observables.ObservableGroup('FRET_probability',db,'pmf',PMF=F)
  
"""

import os
import numpy
import config
from analysis import Observable,AngleprojectedObservable,PMFProjectedObservable

# use a class factory to create observables
# column:    the SQL column name that is used to read data;
#            the class name is constructed as 'column' with an optional prefix
# qualifier: string to uniquely define this observable from other of the
#            same observable_type
# legend:    string used in matplotlib legends
_saltbridges = numpy.rec.fromrecords(
    [('SB1', '1_D33-R156', r'D33-R156'),
     ('SB2', '2_K97-N190', r'K97-N190'),
     ('SB3', '3_D54-K157', r'D54-K157'),
     ('SB4', '4_D54-R156', r'D54-R156'),
     ('SB5', '5_D54-R167', r'D54-R167'),
     ('SB6', '6_R36-D158', r'R36-D158'),
     ('SB7', '7_R36-E170', r'R36-E170'),
     ('SB8', '8_K57-D158', r'K57-D158'),
     ('SB9', '9_K57-E170', r'K57-E170'),
     ('SB10', '10_D118-K136', r'D118-K136'),
     ('SB11', '11_D158-R36', r'D158-R36'),
     ('SB12', '12_D84-K13', r'D84-K13'),
     ('SB13', '13_E185-K97', r'E185-K97'),
     ('SB14', '14_E44-K47', r'E44-K47'),
     ], names='column,qualifier,legend')

_saltbridge_energies = numpy.rec.fromrecords(
    [('ESB1', '1_D33-R156', r'D33-R156'),
     ('ESB2', '2_K97-N190', r'K97-N190'),
     ('ESB3', '3_D54-K157', r'D54-K157'),
     ('ESB4', '4_D54-R156', r'D54-R156'),
     ('ESB5', '5_D54-R167', r'D54-R167'),
     ('ESB6', '6_R36-D158', r'R36-D158'),
     ('ESB7', '7_R36-E170', r'R36-E170'),
     ('ESB8', '8_K57-D158', r'K57-D158'),
     ('ESB9', '9_K57-E170', r'K57-E170'),
     ('ESB10', '10_D118-K136', r'D118-K136'),
     ('ESB11', '11_D158-R36', r'D158-R36'),
     ('ESB12', '12_D84-K13', r'D84-K13'),
     ('ESB13', '13_E185-K97', r'E185-K97'),
     ('ESB14', '14_E44-K47', r'E44-K47'),
     ], names='column,qualifier,legend')

_saltbridge_indicators = numpy.rec.fromrecords(
    [('ISB1', '1_D33-R156', r'D33-R156'),
     ('ISB2', '2_K97-N190', r'K97-N190'),
     ('ISB3', '3_D54-K157', r'D54-K157'),
     ('ISB4', '4_D54-R156', r'D54-R156'),
     ('ISB5', '5_D54-R167', r'D54-R167'),
     ('ISB6', '6_R36-D158', r'R36-D158'),
     ('ISB7', '7_R36-E170', r'R36-E170'),
     ('ISB8', '8_K57-D158', r'K57-D158'),
     ('ISB9', '9_K57-E170', r'K57-E170'),
     ('ISB10', '10_D118-K136', r'D118-K136'),
     ('ISB11', '11_D158-R36', r'D158-R36'),
     ('ISB12', '12_D84-K13', r'D84-K13'),
     ('ISB13', '13_E185-K97', r'E185-K97'),
     ('ISB14', '14_E44-K47', r'E44-K47'),
     ], names='column,qualifier,legend')

_FRETs = numpy.rec.fromrecords(
    [('FRET_Kern', 'Kern',r'I52 -K145 [Henzler-Wildman]',
      r'I52 -K145 [Henzler-Wildman et al, Nature 450 (2007), 838]'),
     ('FRET_Hanson', 'Hanson', r'A127-A194 [Hanson]',
      r'A127-A194 [Hanson et al, PNAS 104 (2007)]'),
     ('FRET_Sinev', 'Sinev', r'A55-V169 [Sinev]',
      r'A55-V169 [Sinev et al, Biochemistry 35 (1996)]'),
     ], names='column,qualifier,legend,doc')


_RMSDs = numpy.rec.fromrecords(
    [('cRMSD',"closed", r"RMSD($X(t)$, 1AKE/closed)",
      """RMSD(X(t), 1AKE) (distance from standard CLOSED state)."""),
     ('oRMSD',"open", r"RMSD($X(t)$, 4AKE/open)",
       """RMSD(X(t), 4AKE) (distance from standard OPEN state)."""),
     ('DeltaRMSD', "DeltaRMSD", r"$\Delta$RMSD = RMSD($X(t)$, 1AKE) - RMSD($X(t)$, 4AKE)",
      """DeltaRMSD = RMSD(X(t), 1AKE) - RMSD(X(t), 4AKE)"""),
     ], names='column,qualifier,legend,doc')

_contacts = numpy.rec.fromrecords(
    [('q1', 'q1', r"$q_1$",
      """Fraction of contacts relative to reference state 1 (1AKE minimized)."""),
     ('q2', 'q2', r"$q_2$",
      """Fraction of contacts relative to reference state 2 (4AKE minimized)."""),
     ], names='column,qualifier,legend,doc')

#------------------------------------------------------------
# observables registry
#
# add lists containing observables (use observable_type as key)

class ObservablesRegistry(object):
    def __init__(self,**kwargs):
        super(ObservablesRegistry,self).__init__(**kwargs)
        self.registry = {}

    def register(self,observable_type, cls, cls_name, cls_dict=None):
        """Create the observable class from a factory function."""
        if not observable_type in self.registry:
            self.registry[observable_type] = {}
        if cls_dict is None:
            # use class directly
            self.registry[observable_type][cls_name] = cls
        else:
            # modify class with cls_dict
            self.registry[observable_type][cls_name] = type(cls_name, (cls,), cls_dict)

    def register_records(self,observable_type,cls,records,name_prefix=''):
        """Register all observables in record array.

        observable_type   identifies the observable (used as part of filenames!)
        cls               Observable class
        records           numpy record array, holding entries describing observables

        The record array requires at a minimum the columns

          column       name of the SQL column
          qualifier    name that distinguish the entry from other observables
                       of the same type (is used as parts of filenames!)
          legend       legend for plots
          doc          optional addition to the __doc__ string
        """
        for r in records:
            cls_name = name_prefix + str(r.column)
            doc = str(cls.__doc__) + ': ' + str(r.legend)
            try:
                doc += '\n' + str(r.doc)
            except:
                pass
            self.register(observable_type, cls, cls_name,
                          {'observable_column': r.column,
                           'qualifier': r.qualifier,
                           'legend': r.legend,
                           '__doc__': doc,
                           })

    def update_dict(self,vardict):
        """Add all observables to a module dictionary.

        update_dict(globals())

        This must be executed at the module level.
        """
        for o in self.values():
            vardict.update(o)
        
    def keys(self):
        return self.registry.keys()
    def values(self):
        return self.registry.values()
    def items(self):
        return self.registry.items()
    def __getitem__(self,x):
        return self.registry[x]
    def __iter__(self):
        return self.registry.__iter__()
    def __len__(self):
        return len(self.registry)


OBSERVABLES = ObservablesRegistry()

# observables
# by convention 1D observables prefix the class name with 'DeltaRMSD_'
class DeltaRMSD_Observable(Observable):
    x_column = 'DeltaRMSD'
    default_binwidth = 0.2
    xlabel = ur'$\Delta\rho = \rho_c - \rho_o$ (Å)'

class DeltaRMSD_Saltbridge(DeltaRMSD_Observable):
    """Observable: saltbridge vs DeltaRMSD"""
    observable_type = 'DeltaRMSD_saltbridge'

OBSERVABLES.register_records('DeltaRMSD_saltbridge',
                             DeltaRMSD_Saltbridge,
                             _saltbridges,
                             name_prefix='DRMSD_')

class DeltaRMSD_SaltbridgeIndicator(DeltaRMSD_Observable):
    """Observable: probability of saltbridge vs DeltaRMSD"""
    observable_type = 'DeltaRMSD_saltbridge_indicator'

OBSERVABLES.register_records('DeltaRMSD_saltbridge_indicator',
                             DeltaRMSD_SaltbridgeIndicator,
                             _saltbridge_indicators,
                             name_prefix='DRMSD_')


class DeltaRMSD_FRET(DeltaRMSD_Observable):
    """Observable: FRET residue CA-CA distance vs DeltaRMSD"""    
    observable_type = 'DeltaRMSD_FRET'

OBSERVABLES.register_records('DeltaRMSD_FRET',
                             DeltaRMSD_FRET,
                             _FRETs,
                             name_prefix='DRMSD_'
                             )

class DeltaRMSD_contacts(DeltaRMSD_Observable):
    """Observable: fractional contacts q vs DeltaRMSD"""    
    observable_type = 'DeltaRMSD_contacts'

OBSERVABLES.register_records('DeltaRMSD_contacts',
                             DeltaRMSD_contacts,
                             _contacts,
                             name_prefix='DRMSD_'
                             )

#------------------------------------------------------------
# definition of PMF/Probability projections on an observable
# does not work because one needs to project...
class FRET_probability(PMFProjectedObservable):
    """Observable: probability from PMF projected on FRET distance"""    
    observable_type = 'FRET_probability'    
    default_binwidth = 0.2
    xlabel = ur'FRET CA-CA distance (Å)'
    ylabel = ur'probability density (Å$^-1$)'
    _title = {'observable': r'$P(d)\ (\AA^-1)$',
              'stdev':      r'$P(d)/\sqrt(N)\ (\AA^-1)$',
              'counts':     r'counts $N$',
              }

OBSERVABLES.register_records('FRET_probability',
                             FRET_probability,
                             _FRETs,
                             name_prefix='probability_'
                             )
### XXX: ObservableGroup 'FRET_probability' is broken ####
###      but individual classes are ok
###      (Problem: need individual Angleprojected observables for init)


#------------------------------------------------------------
# definition of Angle-projected observables
#

class PotentialEnergy(AngleprojectedObservable):
    """Potential energy (without umbrella restraints) projected on NMP and LID angles."""
    observable_type = "energy"
    qualifier = "Upot"
    observable_column = "energy" # could leave empty and use the class magic + case insensitive SQL
    major_contour = 50.0
    default_shift = True
    default_Nmin = 100
    _title = {'observable': r'$\Delta U$ (kcal/mol)',
              'stdev':      r'$\sigma_{\Delta U}$ (kcal/mol)',
              'counts':     r'histogram counts $N$',
              }

OBSERVABLES.register('Energy',PotentialEnergy,'U')


class RMSD(AngleprojectedObservable):
    """RMSD projected on NMP and LID angles."""
    observable_type = "RMSD"
    major_contour = 1.0

OBSERVABLES.register_records('RMSD', RMSD, _RMSDs)


# FRET distances (CA-CA) from the literature
class FRET(AngleprojectedObservable):
    """FRET CA-CA distances projected on NMP and LID angles."""
    observable_type = 'FRET'
    major_contour = 1.0
    
OBSERVABLES.register_records('FRET', FRET, _FRETs)


# salt bridges (see txt/salt-bridge_info.txt)
# 1:D33-R156  2:K97-N190  3:D54-K157  4:D54-R156  5:D54-R167  
# 6:R36-D158  7:R36-E170  8:K57-D158  9:K57-E170 10:D118-K136
class Saltbridge(AngleprojectedObservable):
    """Salt bridge distances projected on NMP and LID angles."""    
    observable_type = 'saltbridge'
    major_contour = 1.0    # Angstroem

OBSERVABLES.register_records('saltbridges', Saltbridge, _saltbridges)


class SaltbridgeEnergy(AngleprojectedObservable):
    """Salt bridge interaction energy projected on NMP and LID angles."""        
    observable_type = 'saltbridge_energy'
    major_contour = 10.0    # kcal/mol
    _title = {'observable': r'$E$ (kcal/mol)',
              'stdev':      r'$\sigma_E$ (kcal/mol)',
              'counts':     r'counts $N$',
              }

OBSERVABLES.register_records('saltbridge_energy', SaltbridgeEnergy, _saltbridge_energies)


class SaltbridgeIndicator(AngleprojectedObservable):
    """Salt bridge existence projected on NMP and LID angles."""        
    observable_type = 'saltbridge_indicator'
    major_contour = 0.2    # probability
    _title = {'observable': r'$I_{sb}$',
              'stdev':      r'$\sigma_{I_{sb}}$',
              'counts':     r'counts $N$',
              }

OBSERVABLES.register_records('saltbridge_indicator', SaltbridgeIndicator, _saltbridge_indicators)

# fractional contacts relative to reference states
class Contacts(AngleprojectedObservable):
    """q1 and q2 projected on NMP and LID angles."""
    observable_type = 'Contacts'
    major_contour = 0.01
    
OBSERVABLES.register_records('Contacts', Contacts, _contacts)


# done with defining Observable classes
del _saltbridges
del _saltbridge_energies
del _saltbridge_indicators
del _FRETs
del _RMSDs
del _contacts

# add them all to the module scope
OBSERVABLES.update_dict(globals())


class ObservableGroup(object):
    """A ObservableGroup bundles similar observables and allows analysis and
    plotting of all of them together. Observables are stored in the observables
    registry OBSERVABLES. New observables can be added; see PMF.observables for
    how to do it.

    Instantiation creates the observables and can take a while, depending on
    the size of the database.

    >>>  og = ObservableGroup(observable_type,db,runtype)

    One can access the individual observables as attributes for introspection
    or simply use the objects plot() and plot_all() methods to create figures
    with standard settings in the standard location under figs/.
    """ + """
    Currently the following ``observable_type``s are set up:
    %r""" % OBSERVABLES.keys()

    
    def __init__(self,observable_type,db,runtype,**kwargs):
        """Calculate all observables in the group observable_type from db.

        observable_type
        db               AngleDB (with observables added via db.insertmany_XXXs(...))
                         or an appropriate selection which contains the required columns
        runtype          string identifier, eg 'dims', 'pmf', ...

        **kwargs  passed to Observable (eg binwidth=2)

        Note that for a big db building the observables can take a long time.     
        """

        self.observable_type = observable_type
        self.runtype = runtype
        try:
            self._classes = OBSERVABLES[observable_type].values()
        except KeyError:
            raise ValueError('observable_type must be one of %r' % OBSERVABLES.keys())
        print "-- Analyzing %s for %r (can take a while...)" % (observable_type,
                                                                [cls.qualifier for cls in self._classes])
        self.observables = dict([(cls.__name__, cls(db,**kwargs)) for  cls in self._classes])

        # add observables as attributes XXX/but prefix with 'O' to avoid potential clashes/
        self.__dict__.update(self.observables)
        

    def plot(self,**kwargs):
        """Create graphs for all defined observables.

        plot(mode=MODE,**kwargs)

        mode      MODE = 'observable', 'stdev', 'counts'
        overlay   function that is called with no arguments, eg
                  overlay = lambda : plot(....)        
        **kwargs  passed to Observable._auto_plots()
                  title is ignored

        Figures are saved in the default location.
        """
        mode = kwargs.pop('mode','observable')
        title = kwargs.pop('title','None')  # ignored
        figs = os.path.join(config.basedir,'figs',self.observable_type)
        for o in self.observables.values():
            print ("-- %(observable_type)s: "+str(o.qualifier)+" for runtype=%(runtype)r") % vars(self)
            figdir = os.path.join(figs,o.qualifier)
            filebasename = o.filebasename(self.runtype)
            o._auto_plots(mode,filebasename,figdir,kwargs)
            print  "--- directory: %s" % figdir


    def plot_all(self,**kwargs):
        """Create observable, standard deviation and count graphs for all defined observables.

        plot_all(**kwargs)


        overlay   function that is called with no arguments, eg
                  overlay = lambda : plot(....)
        **kwargs  passed to to Observable.plot_all()
                  (although only title and cmap can be set for plot)

        Figures are saved in the default location.
        """        
        figs = os.path.join(config.basedir,'figs',self.observable_type)
        plotargs = {}   # pull out the ones used for plotting
        plotargs['title'] = kwargs.pop('title',None)
        plotargs['cmap'] = kwargs.pop('cmap',None)
        plotargs['runtype'] = self.runtype
        plotargs['overlay'] = kwargs.pop('overlay',None)
        # ... and erase the ones that don't make sense (because they will not apply
        # to all three plot modes, hence we're better off with the automatic values
        # TODO: 'do the right thing' and apply where it makes sense or select mode?
        kwargs.pop('min_contour',None)
        kwargs.pop('max_contour',None)
        for o in self.observables.values():
            print ("-- %(observable_type)s: "+str(o.qualifier)+" for runtype=%(runtype)r") % vars(self)
            plotargs['figdir'] = os.path.join(figs,o.qualifier)
            o.plot_all(**plotargs)
    

# quick hacks

import pylab
from analysis import Projector
class FRETprojector_group(object):
    def __init__(self,PMF,FRET,**kwargs):
        self.P = {'Hanson':            Projector(PMF=PMF,X=FRET.FRET_Hanson,**kwargs),
                  'Henzler-Wildman':   Projector(PMF=PMF,X=FRET.FRET_Kern,**kwargs),
                  'Sinev':             Projector(PMF=PMF,X=FRET.FRET_Sinev,**kwargs),
                  }
        for name,obj in self.P.items():
            attr_name = 'FRET_'+str(name).replace('-','_')
            self.__dict__[attr_name] = obj

    def plot_probabilities(self,**kwargs):
        self._plot('observable',**kwargs)
        pylab.ylabel(ur'probability (Å$^{-1}$)')
        print "Save to 'figs/FRET/projection/probabilities_splines.pdf'"
        return pylab.gca()        

    def plot_PMF(self,**kwargs):
        def prob2pmf(x):
            y = numpy.ma.array(x, mask=(x<=0))
            return -numpy.log(y/y.max())
        kwargs['func'] = prob2pmf
        self._plot('observable',**kwargs)
        pylab.ylabel(ur'$\mathcal{W}(d)/kT$')
        print "Save to 'figs/FRET/projection/pmf_splines.pdf'"
        return pylab.gca()        

    def plot_areas(self,**kwargs):
        self._plot('bincounts',**kwargs)
        pylab.ylabel(ur'contour area ($^{\circ2}$)')
        print "Save to 'figs/FRET/projection/areas_splines.pdf'"        
        return pylab.gca()

    def _plot(self,observable,**kwargs):
        from pylab import figure,xlabel,plot,legend
        func = kwargs.pop('func',lambda x: x)
        figure()
        for o in ('Hanson','Henzler-Wildman','Sinev'):
            p = self.P[o]
            data = func(p.__getattribute__(observable))
            plot(p.xvalues, data, lw=2, label=str(o), **kwargs)
        legend()
        xlabel(ur'FRET CA-CA distance (Å)')
        return pylab.gca()
    
    
    
