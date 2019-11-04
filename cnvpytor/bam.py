"""
class Bam: BAM/CRAM/SAM reading class

"""
from __future__ import absolute_import, print_function, division
from .genome import Genome
import pysam
import numpy as np
import logging
import os
import random

_logger = logging.getLogger("cnvpytor.bam")


class Bam:

    def __init__(self, filename, max_fragment_len=5000, max_read_len=300, reference_filename=False):
        """
        Opens BAM/CRAM/SAM file, reads chromosome names/lengths from header file and detects reference genome

        Parameters
        ----------
        filename : str
            Name of the BAM/CRAM/SAM file.
        max_fragment_len : int
            Maximal fragment length used in distribution calculation (default: 5000)
        max_read_len : int
            Maximal read length used in distribution calculation (default: 300).

        """
        self.max_read_len = max_read_len
        self.max_frg_len = max_fragment_len
        self.reference_genome = None
        self.filename = filename
        if filename[-4:] == ".bam":
            self.file = pysam.AlignmentFile(filename, "rb")
        elif filename[-4:] == ".sam":
            self.file = pysam.AlignmentFile(filename, "r")
        elif filename[-5:] == ".cram":
            if reference_filename:
                self.file = pysam.AlignmentFile(filename, "rc",reference_filename=reference_filename)
            else:
                self.file = pysam.AlignmentFile(filename, "rc")
        else:
            _logger.warning("Unsuported file type: " + filename)
        self.len = {}
        if self.file:
            _logger.info("File: " + filename + " successfully open")
            self.reference_genome = Genome.detect_genome(self.file.header.references, self.file.header.lengths)
            if self.reference_genome:
                _logger.info("Detected reference genome: " + self.reference_genome)
            for c, l in zip(self.file.header.references, self.file.header.lengths):
                self.len[c] = l

    def get_chr_len(self):
        """ Get chromosome names and lengths.

        Returns
        -------
        chrs : list of str
            Chromosome names from BAM/CRAM/SAM header.
        len : list of str
            Chromosome lengths from BAM/CRAM/SAM header.

        """
        return self.file.header.references, self.file.header.lengths

    def read_chromosome(self, chr_name):
        """
        Reads chromosome RD data and calculates read vs template length distribution

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.

        Returns
        -------
        rd_p : numpy.ndarray
            RD parity array (100bp bins).
        rd_u : numpy.ndarray
            RD unique array (100bp bins).
        his_read_frg : numpy.ndarray
            2D distribution of read and template lengths.

        """
        if not (chr_name in self.len):
            _logger.warning("Can not find chromosome '%s' in file '%s'." % (chr_name, self.filename))
            return None, None, None
        _logger.debug("Reading chromosome %s from filename %s" % (chr_name, self.filename))
        n = self.len[chr_name] // 100 + 1
        rd_p = np.zeros(n)
        rd_u = np.zeros(n)
        his_read_frg = np.zeros((self.max_read_len, self.max_frg_len))
        try:
            for r in self.file.fetch(chr_name, multiple_iterators=True):
                assert isinstance(r, pysam.libcalignedsegment.AlignedSegment)
                if r.template_length and r.reference_length:
                    fl = abs(r.template_length)
                    rl = r.reference_length
                    if (rl < self.max_read_len) and (fl < self.max_frg_len):
                        his_read_frg[rl][fl] += 1
                    if r.is_unmapped or r.is_secondary or r.is_duplicate:
                        continue
                if r.reference_start and r.reference_end and not (r.mapping_quality is None):
                    mid = (r.reference_end + r.reference_start) // 200
                    rd_p[mid] += 1
                    if r.mapping_quality > 0:
                        rd_u[mid] += 1
        except IOError:
            _logger.error("Error while reading file '%s'" % self.filename)
            exit(0)
        except ValueError:
            _logger.error("Error while reading file '%s'. Probably index is missing." % self.filename)
            exit(0)
        return rd_p, rd_u, his_read_frg

    def pileup(self, chr_name, pos, ref, alt, tmp_file=".cnvpytor"):
        if not (chr_name in self.len):
            _logger.warning("Can not find chromosome '%s' in file '%s'." % (chr_name, self.filename))
            return
        _logger.debug("Pileup chromosome %s from filename %s" % (chr_name, self.filename))
        tmp_file += "_" + str(random.randint(0, 1e10)) + "_" + chr_name
        f = open(tmp_file, "w")
        for i in pos:
            print(chr_name, i, file=f)
        f.close()
        mpile = pysam.mpileup("-r", chr_name, "-l", tmp_file, self.filename)
        os.remove(tmp_file)
        pos_seq = dict([(int(x.split("\t")[1]), x.split("\t")[4].upper()) for x in mpile.split("\n") if x != ""])
        nref = [0] * len(pos)
        nalt = [0] * len(pos)
        for ix in range(len(pos)):
            if pos[ix] in pos_seq:
                nref[ix] = pos_seq[pos[ix]].count(ref[ix])
                nalt[ix] = pos_seq[pos[ix]].count(alt[ix])
        return nref, nalt
