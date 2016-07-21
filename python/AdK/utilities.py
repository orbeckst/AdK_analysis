#!/usr/bin/env python
# $Id: utilities.py 2175 2008-08-13 06:10:14Z oliver $

__doc__ = """
Read a pdb file and extract the chain data in Modeller format:

>P1;1AKE_B
structure:1ake: 1 :B: 214 :B::::

Uses Bio.PDB from http://biopython.org
See http://biopython.org/docs/api/public/
and http://www.binf.ku.dk/~thamelry/biopdb_faq.pdf
"""

from Bio.PDB.PDBParser import PDBParser
import Bio.PDB
import os,sys
import sre

#------------------------------------------------------------
# Global variables
#------------------------------------------------------------

verbose = 3


#------------------------------------------------------------
# functions
#------------------------------------------------------------


def msg(level,m):
    """
       Print a message if the verbosity level  >= the level of the message.
       The verbosity level 'verbose' must be declared globally
    """
    if (verbose >= level):
        if (verbose > 3):
            """If in debug mode then show debug level with message"""
            m = "[dbglvl%3d] %s" % (level,m)
        sys.stderr.write("%s\n" % m)


def choose_altLoc (structure, altLoc="A"):
    """
    Chose one set of disorder specifiers on a per atom basis.

    Parameters: structure=<structure>
                    structure from the PDBParser
                altLoc=<string> ("A")
                    string describing the disordered set
    """

    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                if residue.is_disordered():
                    for atom in residue.get_list():
                        if atom.is_disordered():
                            if atom.disordered_has_id(altLoc):
                                atom.disordered_select(altLoc)

class PDB2Modeller(object):
    def __init__(self,pdb_id,pdbfile):
        import os.path
        regex_chain = sre.compile('[A-Z]')
        p = PDBParser(PERMISSIVE=1)
        self.structure = p.get_structure(pdb_id, pdbfile)
        msg(3,"Read pdb file '%s'" % pdbfile)
    
        # Not needed, just as a reminder how to select one set: 
        ### choose_altLoc(structure, altLoc)

        pdbfile_id = os.path.splitext(pdbfile)[0]

        self.data = {}

        # Only look at chains with an identifier
        chains = [c for c in 
                  Bio.PDB.Selection.unfold_entities(self.structure,'C') 
                  if regex_chain.match(c.id) is not None]
        for c in chains:
            chain_id = c.id
            residues = [r.id for r in c.get_iterator() if r.id[0] == ' ']
            if len(residues) == 0:
                continue
            first,last = residues[0], residues[-1]
            first_resid = first[1]
            last_resid = last[1]
            modeller_id = "%s_%s" % (pdb_id.upper(), chain_id)
            self.data[modeller_id] = ModellerHeader(modeller_id,pdb_id,pdbfile_id,
                                                    chain_id,first_resid,last_resid)
            msg(3,str(self.data[modeller_id]))

    def show(self):
        for modeller_id in sorted(self.data):
            print self.data[modeller_id].header()

    def write(self,filename,mode='w'):
        f = open(filename,mode)
        for modeller_id in sorted(self.data):
            f.write(self.data[modeller_id].header())
        f.close()


class ModellerHeader(object):
    def __init__(self,modeller_id,pdb_id,pdbfile_id,chain_id,first_resid,last_resid):
        self.modeller_id = modeller_id
        self.pdb_id = pdb_id
        self.pdbfile_id = pdbfile_id
        self.chain_id = chain_id
        self.first_resid = first_resid
        self.last_resid = last_resid

    def header(self):
        return """>P1;%(modeller_id)s\n"""\
            """structure:%(pdbfile_id)s: %(first_resid)d :%(chain_id)s: %(last_resid)d :%(chain_id)s::::\n""" % self.__dict__

    def __str__(self):
        return "%(modeller_id)s %(pdbfile_id)s: %(first_resid)d :%(chain_id)s: %(last_resid)d :%(chain_id)s" % self.__dict__

    def __repr__(self):
        return """ModellerHeader(modeller_id='%(modeller_id)s',pdb_id='%(pdb_id)s',"""\
            """pdbfile_id='%(pdbfile_id)s',chain_id='%(chain_id)s',"""\
            """first_resid=%(first_resid)d,last_resid=%(last_resid)d)""" % self.__dict__


class  all_pdb2modeller(object):
    """Convenience class to extract the Modeller headers from a number
    of pdb files and optionally patch these headers into an existing
    pir alignment file:

    >>> all = AdK.utilities.all_pdb2modeller('*.pdb')
    >>> all.modify_alignment('../../seq/AdK_withlid.pir')
    """

    def __init__(self,globpattern="*.pdb"):
        """Extract the Modeller line from all files matching globpattern."""
        import glob
        pdb_pattern = sre.compile('(?P<pdb_id>\d\w\w\w).*\.(pdb|ent)')
        self.pdbfiles = glob.glob(globpattern)
        self.p2ms = []

        for pdbfile in self.pdbfiles:
            m = pdb_pattern.match(pdbfile)
            if m is None:
                msg(0,"%(pdbfile)s does not include a valid PDB ID in filename "
                "so we skip it." % locals())
                continue
            pdb_id = m.group('pdb_id')
            msg(0,"%(pdb_id)s" % locals())  
            self.p2ms.append(PDB2Modeller(pdb_id,pdbfile))

        self.headers = {}   # lookup headers by modeller_id
        for x in self.p2ms:
            self.headers.update(x.data)

    def show(self):
        for x in self.p2ms:
            x.show()
            
    def write(self,filename='pdb2modeller.txt'):
        try:
            os.unlink(filename)
        except:
            pass
        for x in self.p2ms:
            x.write(filename,'a')

    def modify_alignment(self,pirfile):
        """Insert/replace Modeller lines in pir file.

        This only works if the PIR file already has the proper
        Modeller ID in the header line (tcoffee produces this
        automatically when ran on pdb files). See the example below.

           Modeller ID
            vvvvvv
        >P1;1ANK_A|PDBID|CHAIN|SEQUENCE
        (empty line)

        >P1;1ANK_A
        structure:1ank: 1 :A: 214 :A::::
        """

        modellerfile = pirfile + '_modeller'
        pir = open(pirfile,'r')
        modeller = open(modellerfile,'w')

        header_pat = sre.compile('>P1;(?P<modeller_id>\d\w\w\w_[A-Z])')
        
        for line in pir:
            m = header_pat.match(line)
            if m is not None:
                modeller_id = m.group('modeller_id')   # header line
                pir.next() # discard next line (empty in original PIR file)
                # substitute with new header
                line = self.headers[modeller_id].header()
                msg(3,line)
            modeller.write(line)
                
        pir.close()
        modeller.close()

        os.rename(modellerfile,pirfile)
        msg(3,"Done: patched %(pirfile)s" % locals())
