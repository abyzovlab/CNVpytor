""" cnvpytor.trio

class Trio: CNVpytor class for trio analysis
"""
from __future__ import absolute_import, print_function, division

from .io import *
from .utils import *
from .genome import Genome
import numpy as np
import logging


_logger = logging.getLogger("cnvpytor.trio")

class Trio:
    """
    Class for trio analysis
    """
    def __init__(self, files):
        """
        Initialize a Trio object

        Parameters
        ----------
        filenames : list of str
            CNVpytor filenames: child, parent1, parent2
        """
        _logger.debug("Trio class init")
        self.filenames = files
        if len(files) != 3:
            _logger.debug("Trio class init: not a trio")
            raise Exception("Trio object requires 3 filenames")

        self.io = [IO(f, ro=False) for f in files]
        self.io_gc = self.io[0]
        self.io_mask = self.io[0]
        if self.io[0].signal_exists(None, None, "reference genome") \
                and self.io[0].signal_exists(None, None, "use reference"):
            rg_name = np.array(self.io[0].get_signal(None, None, "reference genome")).astype("str")[0]
            if rg_name in Genome.reference_genomes:
                Genome.detected_genome = rg_name
                rg_use = self.io[0].get_signal(None, None, "use reference")
                if "gc_file" in Genome.reference_genomes[rg_name] and rg_use[0] == 1:
                    _logger.debug("Using GC content from database for reference genome '%s'." % rg_name)
                    self.io_gc = IO(Genome.reference_genomes[rg_name]["gc_file"], ro=True, buffer=True)
                if "mask_file" in Genome.reference_genomes[rg_name] and rg_use[1] == 1:
                    _logger.debug("Using strict mask from database for reference genome '%s'." % rg_name)
                    self.io_mask = IO(Genome.reference_genomes[rg_name]["mask_file"], ro=True, buffer=True)

    def trio_phase(self, parents=False, callset=None):
        """ Phase SNPs in a trio.

        Parameters
        ----------
        parents : bool
            If True, phase parents.
        callset : str or None
            It will assume SNP data if None. Otherwise it will assume SNV data
            stored under the name provided by callset variable.

        Returns
        -------
        None

        """

        for c in self.io[0].snp_chromosomes():
            if c in self.io[1].snp_chromosomes() and c in self.io[2].snp_chromosomes():
                _logger.info("Phasing SNP data for chromosome '%s'." % c)
                ch, fa, mo = {}, {}, {}
                pos, ref, alt, nref, nalt, gt, flag, qual = self.io[0].read_snp(c, callset=callset)
                for i in range(len(pos)):
                    k = str(pos[i])+":"+ref[i]+">"+alt[i]
                    ch[k] = (gt[i] % 4,i)
                pos1, ref1, alt1, nref1, nalt1, gt1, flag1, qual1 = self.io[1].read_snp(c, callset=callset)
                for i in range(len(pos1)):
                    k = str(pos1[i])+":"+ref1[i]+">"+alt1[i]
                    fa[k] = (gt1[i] % 4,i)
                pos2, ref2, alt2, nref2, nalt2, gt2, flag2, qual2 = self.io[2].read_snp(c, callset=callset)
                for i in range(len(pos2)):
                    k = str(pos2[i])+":"+ref2[i]+">"+alt2[i]
                    mo[k] = (gt2[i] % 4,i)

                _logger.info("Phasing SNP data for child '%s'." % c)
                rpos, rref, ralt, rnref, rnalt, rgt, rflag, rqual = [], [], [], [], [], [], [], []
                for i in range(len(pos)):
                    k = str(pos[i])+":"+ref[i]+">"+alt[i]
                    resgt = ch[k][0]
                    if ch[k][0]==0 or ch[k][0]==3:
                        resgt = ch[k][0] + 4
                    else:
                        if k in fa and k in mo and fa[k][0]!=0 and mo[k][0]!=0:
                            if fa[k][0]==3 and mo[k][0] in [1,2]:
                                resgt = 5
                                #_logger.debug("Phasing 5 (1/1 0/1). SNP: %s:%s" % (c, k))
                            elif mo[k][0]==3 and fa[k][0] in [1,2]:
                                resgt = 6
                                #_logger.debug("Phasing 6 (0/1 1/1). SNP: %s:%s" % (c, k))
                            else:
                                _logger.debug("Unable to phase. Not unique - skipping. SNP: %s:%s" % (c,k))
                        elif k in fa and fa[k][0]!=0:
                            resgt = 5
                            #_logger.debug("Phasing 5 (?/1 0/0). SNP: %s:%s" % (c, k))
                        elif k in mo and mo[k][0]!=0:
                            resgt = 6
                            #_logger.debug("Phasing 6 (0/0 ?/1). SNP: %s:%s" % (c, k))
                        else:
                            _logger.debug("Unable to phase. Not present in parents - skipping. SNP: %s:%s" % (c,k))
                    if resgt != -1:
                        rpos.append(pos[i])
                        rref.append(ref[i])
                        ralt.append(alt[i])
                        rnref.append(nref[i])
                        rnalt.append(nalt[i])
                        rgt.append(gt[i])
                        rflag.append(flag[i])
                        rqual.append(qual[i])
                        rgt[-1] = resgt
                _logger.info("Writing phased SNP data for child '%s'." % c)
                self.io[0].save_snp(c, rpos, rref, ralt, rnref, rnalt, rgt, rflag, rqual, update=True, callset=callset)
                if parents:
                    _logger.info("Phasing SNP data for parent 1 '%s'." % c)
                    rpos, rref, ralt, rnref, rnalt, rgt, rflag, rqual = [], [], [], [], [], [], [], []
                    for i in range(len(pos1)):
                        k = str(pos1[i])+":"+ref1[i]+">"+alt1[i]
                        resgt = fa[k][0]
                        if fa[k][0]==0 or fa[k][0]==3:
                            resgt = fa[k][0] + 4
                        else:
                            if k in ch and ch[k][0]!=0 and (k not in mo or mo[k][0]==0):
                                resgt = 5
                            elif k in ch and ch[k][0]==3:
                                resgt = 5
                            elif (k not in ch or ch[k][0]==0):
                                resgt = 6
                            else:
                                _logger.debug("Unable to phase. Not unique. SNP: %s:%s" % (c,k))
                        if resgt != -1:
                            rpos.append(pos1[i])
                            rref.append(ref1[i])
                            ralt.append(alt1[i])
                            rnref.append(nref1[i])
                            rnalt.append(nalt1[i])
                            rgt.append(gt1[i])
                            rflag.append(flag1[i])
                            rqual.append(qual1[i])
                            rgt[-1] = resgt
                    _logger.info("Writing phased SNP data for parent 1 '%s'." % c)
                    self.io[1].save_snp(c, rpos, rref, ralt, rnref, rnalt, rgt, rflag, rqual, update=True,
                                        callset=callset)

                    _logger.info("Phasing SNP data for parent 2 '%s'." % c)
                    rpos, rref, ralt, rnref, rnalt, rgt, rflag, rqual = [], [], [], [], [], [], [], []
                    for i in range(len(pos2)):
                        k = str(pos2[i])+":"+ref2[i]+">"+alt2[i]
                        resgt = mo[k][0]
                        if mo[k][0]==0 or mo[k][0]==3:
                            resgt = mo[k][0] + 4
                        else:
                            if k in ch and ch[k][0]!=0 and (k not in fa or fa[k][0]==0):
                                resgt = 6
                            elif k in ch and ch[k][0]==3:
                                resgt = 6
                            elif (k not in ch or ch[k][0]==0):
                                resgt = 5
                            else:
                                _logger.debug("Unable to phase. Not unique. SNP: %s:%s" % (c,k))
                        if resgt != -1:
                            rpos.append(pos2[i])
                            rref.append(ref2[i])
                            ralt.append(alt2[i])
                            rnref.append(nref2[i])
                            rnalt.append(nalt2[i])
                            rgt.append(gt2[i])
                            rflag.append(flag2[i])
                            rqual.append(qual2[i])
                            rgt[-1] = resgt

                    _logger.info("Writing phased SNP data for parent 2 '%s'." % c)
                    self.io[2].save_snp(c, rpos, rref, ralt, rnref, rnalt, rgt, rflag, rqual, update=True,
                                        callset=callset)
            else:
                _logger.info("Chromosome '%s' not present in parents data." % c)
