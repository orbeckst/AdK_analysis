#!/usr/bin/env python
# $Id: analysis.py 4102 2010-04-26 21:10:32Z oliver $
"""Process DIMS AdK transitions: extract the angle between LID and CORE domain
and NMP and CORE domain and save the angle pairs as pickled numpy arrays.

Format of array
---------------

 shape (2,N)               N = number of saved steps in trajectory
 [LID_angle, NMP_angle]    angles in degree

Note on angle definitions
-------------------------

The definitions of the angles are described in
txt/observables.txt. They *must* match between here and the MMFP GEO
restraints in pmf/pmf_angles.inp and pmf/pmf_temp.inp.

The corresponding CHARMM code:

   MMFP
   GEO  MAXGEO 5000 -
        sphere RCM angle -
        harmonic symmetric force  @FORCECONSTANT tref @NMP dtoff 0.0 -
        select resid 115:125 .and. backbone end     select resid 90:100 .and. backbone end     select resid 35:55 .and. backbone end 
   GEO  sphere RCM angle -
        harmonic symmetric force  @FORCECONSTANT tref @LID dtoff 0.0 -
        select resid 179:185 .and. backbone end     select resid 115:125 .and. backbone end     select resid 125:153 .and. backbone end 
   END

"""

import config    # configuration for this installation, AdK.config

import os
import numpy
from MDAnalysis import Universe
from Scientific.Geometry import Vector


try:
	from numpy import degrees
except ImportError:
	_rad2deg = 180./numpy.pi
	def degrees(x):
		return _rad2deg * x

# not working properly yet / not used
import fcntl, errno
class LockedFile(object):
	def __new__(cls, filename, mode='r'):
		"""Return a locked file object or raise IOError if no lock can be acquired
		
		f = LockedFile(filename)
		if f is None: raise IOError('Cannot get a lock')

		This implementation uses os.O_EXLOCK; this is not sufficient on nfs file
		systems. See the discussion in open(2) [http://linux.die.net/man/2/open].
		"""
		
		try:
			fd = os.open(filename,os.O_RDONLY)
		except IOError,err:
			# still a race between try and next os.open
			if err.errno == errno.ENOENT and 'w' in mode:
				fd = os.open(filename,os.O_WRONLY)
			else:
				raise
		try:
			fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
		except IOError,err:
			if err.errno in (errno.EACCES, errno.EAGAIN):
				raise IOError("Cannot acquire lock for %r" % filename)
			else:
				raise
		return os.fdopen(fd,mode)

	def __init__(self, filename, mode='r'):
		self.filename = name

	def lock(self):
		"""Raises IOError if no lock can be acquired"""
		try:
			fcntl.flock(self, fcntl.LOCK_EX | fcntl.LOCK_NB)
		except IOError,err:
			if err.errno in (errno.EACCES, errno.EAGAIN):
				raise IOError("Cannot acquire lock for %r" % self.name)
			else:
				raise
				
	def unlock(self):
		fcntl.flock(self, fcntl.LOCK_UN)

	def __repr__(self):
		return "LockedFile(%(filename)r,mode=%(mode)r)" % vars(self)
			

def angles_from_dcd(psfname,dcdname,targetdir=os.path.curdir,force=False,infix="_pmf"):
        """Extract LID and NMP angles from dcd and write a pickled numpy array.

        Arguments:
        psfname         psf for the trajectory
        dcdname         trajectory	
        targetdir       write pickle files to this directory [.]
	infix           pickle filename is '<dcdname><infix>_angles.pickle' ["_pmf"]
        force           True: always compute the trajectory and overwrite existing
                        False: only compute if the pickle file does not exist yet [False]
        """
        
	anglename = os.path.splitext(os.path.basename(dcdname))[0] + infix + "_angles.pickle"
	# normalize name (because some dcds already have a _pmf appended)
	savefile = os.path.join(targetdir, anglename.replace('_pmf_pmf','_pmf'))
	if os.path.isfile(savefile) and not force:
		print "File %(savefile)s already exists, skipping %(dcdname)s." % vars()
		return True
	try:
		universe = Universe(psfname,dcdname)
	except:
		print "Failed to load Universe(%s,%s)" % (psfname,dcdname)
		return False
		
	print "Analyzing: dcdfile = %(dcdname)s ---> %(savefile)s" % vars()
	try:
		tot_d = []
		for ts in universe.dcd:
			# angles: (top,core,bot)  
			# v0 = top - core, v1 = bot - core
			# cos alpha = v0*v1/v0v1
                        
                        # These defintions match the ones in the MMFP GEO restraints
                        # of pmf_angles.inp or pmf_bat.inp
			nbot = universe.selectAtoms("resid 115:125 and backbone").centerOfGeometry()
			ntop = universe.selectAtoms("resid 35:55 and backbone").centerOfGeometry()
			ncore = universe.selectAtoms("resid 90:100 and backbone").centerOfGeometry()

			lidc = universe.selectAtoms("resid 115:125 and backbone").centerOfGeometry()
			lidt = universe.selectAtoms("resid 125:153 and backbone").centerOfGeometry()
			lidb = universe.selectAtoms("resid 179:185 and backbone").centerOfGeometry()
			
			l_v0 = Vector(lidb - lidc)
			l_v1 = Vector(lidt - lidc)
			n_v0 = Vector(ntop - ncore)
			n_v1 = Vector(nbot - ncore)

			LID_angle = degrees( l_v1.angle(l_v0) )
			NMP_angle = degrees( n_v1.angle(n_v0) )

			tot_d.append([LID_angle,NMP_angle])
                        if  ts.frame % 10 == 0 or ts.frame == universe.dcd.numframes:
                            print "%s  frame %6d/%6d\r" % (dcdname, ts.frame, universe.dcd.numframes),
		dist_profile = numpy.array(tot_d)
		dist_profile.dump(savefile)
		print "Saved %d angle-pairs to '%s'        " % (len(tot_d),savefile)
		del dist_profile
		del tot_d
	except:
		print "Failed to process"
		return False
        return True



from subprocess import Popen,call
import tempfile

def CHARMM_LIB(*args):
	return os.path.join(config.basedir,'lib','charmm',*args)
ENERGY_ANGLES = CHARMM_LIB('energy_angles.inp')
FRET_DISTANCES = CHARMM_LIB('fret_distances.inp')
SALTBRIDGES = CHARMM_LIB('saltbridge_distances.inp')
INTER_SB = CHARMM_LIB('inter_sb.inp')
MG_DIST = CHARMM_LIB('mg_distances.inp')

def energy_from_dcd(**kwargs):
	return _run_charmm_script(ENERGY_ANGLES,'_energy.dat',**kwargs)

def inter_sb_from_dcd(**kwargs):
	return _run_charmm_script(INTER_SB,'_inter_sb.dat',**kwargs)

def fret_from_dcd(**kwargs):
	return _run_charmm_script(FRET_DISTANCES,'_fret.dat',**kwargs)

def saltbridges_from_dcd(**kwargs):
	return _run_charmm_script(SALTBRIDGES,'_saltbridges.dat',**kwargs)

def mg_distance_from_dcd(**kwargs):
        return _run_charmm_script(MG_DIST,'_distance.dat',**kwargs)


def _run_charmm_script(charmm_script,datasuffix,
		       psf=None,dcdname=None,targetdir=os.path.curdir, force=False,charmm='c35a1-dims',
		       toppar=os.path.join(config.basedir,'lib','charmm'),keeplog=False):
	if psf is None or dcdname is None:
		raise ValueError("psf and/or dcdname are missing")
	output = os.path.splitext(os.path.basename(dcdname))[0] + datasuffix
	output = os.path.join(targetdir, output)
	output_bz = output+'.bz2'    
	if (os.path.isfile(output) or os.path.isfile(output_bz)) and not force:
		print "File %(output)r or %(output_bz)r already exists, skipping %(dcdname)r." % vars()
		return None
	print "Charmm: %(dcdname)r --> %(output)r" % vars()
	# double quotes are important to keep case in Charmm and must not be removed so
	# we avoid using /bin/sh implicitly and build the command as a list
	cmd = [charmm, 'PSF="'+psf+'"', 'DCD="'+dcdname+'"', 'OUTPUT="'+output+'"',
	       'TOPPAR="'+toppar+'"']
	script = open(charmm_script,'r')
	log,logname = tempfile.mkstemp(suffix='.log',prefix='charmm_')  # is deleted if all went ok
	try:
		rc = Popen(cmd, stdin=script, stdout=log, shell=False).wait()
	finally:
		script.close()
		os.close(log)

	if rc > 0 or not os.path.isfile(output):
		print "WARNING: Charmm failed, returncode = %(rc)d for %(dcdname)r (maybe look at %(logname)r?)" % vars()
		return None
	if not keeplog:
		os.unlink(logname)    # remove log because all went well
	else:
		print "DEBUG mode: keeping logfile %(logname)r." % vars()
    
	print "bzip2:   %(output)r --> " % vars(),
	cmd = ['bzip2']
	if force:
		cmd.append('-f')
	cmd.append(output)
	rc = call(cmd)
	if rc == 0:
		print "%r" % output_bz
		return output_bz
	else:
		print "FAILED"
		return output

RMS_CO_ANGLES = os.path.join(config.basedir,'lib','charmm','rms_cl.inp')
RMS_OC_ANGLES = os.path.join(config.basedir,'lib','charmm','rms_calc_cl.inp')

def rms_from_dcd(psf,dcdname,targetdir=os.path.curdir, force=False,charmm='c35a1-dims',
                    toppar=os.path.join(config.basedir,'lib','charmm'),keeplog=False):
        import arrayio
        import subprocess
        charmm_script = RMS_CO_ANGLES
        charmm_script2 = RMS_OC_ANGLES
        output = os.path.splitext(os.path.basename(dcdname))[0] + 'co_rms.dat'
        output = os.path.join(targetdir, output)
        output2 = os.path.splitext(os.path.basename(dcdname))[0] + 'oc_rms.dat'
        output2 = os.path.join(targetdir, output2)
        output3 = os.path.splitext(os.path.basename(dcdname))[0]
        output3 = os.path.join(targetdir, output3)
        if (os.path.isfile(output)) and not force:
                print "File %(output)r already exists, skipping %(dcdname)r." % vars()
                return None
        print "Charmm: %(dcdname)r --> %(output)r" % vars()
        if (os.path.isfile(output2)) and not force:
                print "File %(output2)r already exists, skipping %(dcdname)r." % vars()
                return None
        print "Charmm: %(dcdname)r --> %(output2)r" % vars()
        # double quotes are important to keep case in Charmm and must not be removed so
        # we avoid using /bin/sh implicitly and build the command as a list
        cmd = [charmm, 'PSF="'+psf+'"', 'DCD="'+dcdname+'"', 'OUTPUT="'+output+'"',
               'TOPPAR="'+toppar+'"']
        cmd2 = [charmm, 'PSF="'+psf+'"', 'DCD="'+dcdname+'"', 'OUTPUT="'+output2+'"',
               'TOPPAR="'+toppar+'"']
        script = open(charmm_script,'r')
        script2 = open(charmm_script2,'r')
        log = open(output,'w')
        log2 = open(output2,'w')
        try:
                rc_co = Popen(cmd,stdin=script, stdout=subprocess.PIPE, shell=False) #.wait()
                rc_1 = Popen(["grep","THUS RMS DIFF IS"],stdin=rc_co.stdout,stdout=log)
                out1 = rc_1.communicate()[0]
                rc_oc = Popen(cmd2, stdin=script2, stdout=subprocess.PIPE, shell=False) #.wait()
                rc_2 = Popen(["grep","THUS RMS DIFF IS"],stdin=rc_oc.stdout,stdout=log2)
                out2 = rc_2.communicate()[0]
                #print "hi"
        finally:
                script.close()
                script2.close()
                log.close()
                log2.close()
	tot_d = []
	tot_rms = arrayio.readFloatArrayWithSkip(output,[4])
	tot_rms2 = arrayio.readFloatArrayWithSkip(output2,[4])
	for man in range(len(tot_rms)-1):
		rms = tot_rms[man]
		rms2 = tot_rms2[man]
		deltaRMSD = rms2 - rms
		#print deltaRMSD
		tot_d.append([rms,rms2,deltaRMSD]) ##,angles[:,0][man],angles[:,1][man]])
	dist_profile = numpy.array([tot_d])
	dist_profile.dump(output3+'_deltaRMSD.pickle')

        if not os.path.isfile(output):
                print "WARNING: Charmm failed, returncode = %(rc_co)r for %(dcdname)r (maybe look at %(log)r?)" % vars()
                return None
        #if not keeplog:
        #        os.unlink(log)    # remove log because all went well
        else:
                print "DEBUG mode: keeping logfile %(log)r." % vars()

	return output, output2, output3


def contacts_from_dcd(psfname, dcdname, refpdb1=None, refpdb2=None, radius=8.0, 
		      targetdir=os.path.curdir, force=False, infix=""):
	"""Contact (q1-q2) analysis.

	Supply the reference pdbs; if refpdb1 is None then the first
	frme of the dcd is used; if refpdb2 is None then the last
	frame is used. These pdbs are extracted to the targetdir.
	"""
	from contacts import ContactAnalysis
	C = ContactAnalysis(psfname, dcdname, ref1=refpdb1, ref2=refpdb2, radius=radius,
			    targetdir=targetdir, infix=infix, force=force)
	output = C.run(store=False)
	if output:
		print "Analyzed: dcdfile = %(dcdname)r ---> %(output)r" % vars()
	return output
