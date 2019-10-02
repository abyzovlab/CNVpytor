""" cnvpytor.fasta

class Fasta: reading fa.gz file

"""
from .genome import Genome
import pysam
import logging
import re

_logger = logging.getLogger("cnvpytor.fasta")


class Fasta:

    def __init__(self, filename):
        """
        Opens FASTA file, reads chromosome names/lengths and detects reference genome

        Parameters
        ----------
        filename : str
            Name of the BAM/CRAM/SAM file.

        """
        self.reference_genome = None
        self.filename = filename
        self.file = None
        try:
            self.file = pysam.FastaFile(filename)
        except IOError:
            _logger.error("Problem opening file '%s' " % filename)
            exit(0)
        except ValueError:
            _logger.error("Index for filename '%s' is missing!" % filename)
            exit(0)

        self.len = {}
        if self.file:
            _logger.info("File: " + filename + " successfully open")
            self.reference_genome = Genome.detect_genome(self.file.references, self.file.lengths)
            if self.reference_genome:
                _logger.info("Detected reference genome: " + self.reference_genome)
            for c, l in zip(self.file.references, self.file.lengths):
                self.len[c] = l

    def get_chr_len(self):
        """
        Get chromosome names and lengths.

        Returns
        -------
        chrs : list of str
            Chromosome names from BAM/CRAM/SAM header.
        len : list of str
            Chromosome lengths from BAM/CRAM/SAM header.

        """
        return self.file.references, self.file.lengths

    def read_chromosome_gc(self, chr_name):
        """
        Reads chromosome GC/AT content

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.

        Returns
        -------
        gc : list of int
            Binned GC content (100bp bins).
        at : list of int
            Binned AT content (100bp bins).

        """
        if not (chr_name in self.len):
            _logger.warning("Can not find chromosome '%s' in fasta file '%s'." % (chr, self.filename))
            return None, None
        _logger.debug("Reading chromosome: %s" % chr_name)
        seq = self.file.fetch(chr_name).upper()

        gc = [seq.count("G", i, i + 100) + seq.count("C", i, i + 100)
              for i in range(0, len(seq), 100)]
        at = [seq.count("A", i, i + 100) + seq.count("T", i, i + 100)
              for i in range(0, len(seq), 100)]
        n = self.len[chr_name] // 100 + 1
        if len(gc) < n:
            gc.append(0)
            at.append(0)
        tot = len(seq)
        sgc = sum(gc)
        sat = sum(at)
        snn = tot - sgc - sat
        _logger.info(
            "GC/AT/N content: %.1f%% / %.1f%% / %.1f%%" % (100. * sgc / tot, 100. * sat / tot, 100. * snn / tot))
        return gc, at

    def read_chromosome_mask_p_regions(self, chr_name):
        """
        Reads chromosome strict mask P regions.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.

        Returns
        -------
        p_mask : list of (int, int)
            List of strict mask P regions: [(start_1, end_1), (start_2, end_2),...]

        """
        if not (chr_name in self.len):
            _logger.warning("Can not find chromosome '%s' in fasta file '%s'." % (chr_name, self.filename))
            return None
        _logger.debug("Reading chromosome: %s" % chr_name)
        seq = self.file.fetch(chr_name).upper()
        p = re.compile("([P]+)")
        re_iterator = p.finditer(seq)
        return [i.span() for i in re_iterator]
