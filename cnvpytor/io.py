""" cnvpytor.io

class IO: Reading/writing CNVpytor files (extension .pytor) using h5py library.

"""
from __future__ import absolute_import, print_function, division
from .genome import Genome
from .utils import *
import logging
import numpy as np
import h5py
import re

_logger = logging.getLogger("cnvpytor.io")

FLAG_AUTO = 0x0001
FLAG_SEX = 0x0002
FLAG_MT = 0x0004
FLAG_GC_CORR = 0x0010
FLAG_AT_CORR = 0x0020
FLAG_USEMASK = 0x0100
FLAG_USEID = 0x0200
FLAG_USEHAP = 0x0400


class IO:
    signals = {
        "RD p": "%(chr)s_rd_p",
        "RD u": "%(chr)s_rd_u",
        "GC/AT": "%(chr)s_gc",
        "mask": "%(chr)s_mask",
        "GC corr": "gc_corr%(flag)s",
        "RD p dist": "dist_rd_p%(flag)s",
        "RD u dist": "dist_rd_u%(flag)s",
        "RD GC dist": "dist_rd_gc%(flag)s",
        "RD stat": "rd_stat%(flag)s",
        "RD": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s",
        "RD unique": "his_rd_u_%(chr)s_%(bin_size)d%(rd_flag)s",
        "RD raw": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_raw",
        "RD l1": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l1",
        "RD l2": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l2",
        "RD l3": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l3",
        "RD partition": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s",
        "RD call": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s_merge",
        "GC": "%(chr)s_gc_%(bin_size)",
        "SNP pos": "%(chr)s_snp_pos",
        "SNP desc": "%(chr)s_snp_desc",
        "SNP counts": "%(chr)s_snp_counts",
        "SNP qual": "%(chr)s_snp_qual",
        "SNP bin_count": "snp_bafc_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP baf": "snp_baf_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP maf": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP likelihood": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i1": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i2": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i3": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i4": "snp_i4_bafc_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP likelihood partition": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP maf partition": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i1 partition": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i2 partition": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i3 partition": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i4 partition": "snp_i4_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP likelihood call": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP maf call": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i1 call": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i2 call": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i3 call": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i4 call": "snp_i4_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "RD chromosomes": "rd_chromosomes",
        "SNP chromosomes": "snp_chromosomes",
        "read frg dist": "read_frg_len",
        "reference genome": "reference_genome",
        "use reference": "use_reference"
    }

    def __init__(self, filename):
        """
        Opens CNVpytor file for reading/writing

        Parameters
        ----------
        filename : str
            Name of the file.

        """
        self.filename = filename
        _logger.debug("Opening h5 file '%s'" % self.filename)
        try:
            self.file = h5py.File(filename, "r+")
        except IOError:
            self.file = h5py.File(filename, "w")

    def __del__(self):
        _logger.debug("Closing h5 file '%s'" % self.filename)
        self.file.close()

    @staticmethod
    def suffix_rd_flag(flags):
        """
        Converts binary flags into suffix used in RD signal names.

        Parameters
        ----------
        flags : int
            Binary flag (FLAG_GC_CORR = 0x0010, FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100).

        Returns
        -------
        s : str
            Suffix string used in RD signal names.

        """
        s = ""
        if flags & FLAG_AT_CORR:
            s += "_AT"
        if flags & FLAG_GC_CORR:
            s += "_GC"
        if flags & FLAG_USEMASK:
            s += "_mask"
        return s

    @staticmethod
    def suffix_snp_flag(flags):
        """
        Converts binary flags into suffix used in SNP signal names.

        Parameters
        ----------
        flags : int
            Binary flag (FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        s : str
            Suffix string used in SNP signal names.

        """
        s = ""
        if flags & FLAG_USEMASK:
            s += "_mask"
        if flags & FLAG_USEID:
            s += "_id"
        if flags & FLAG_USEHAP:
            s += "_hap"
        return s

    @staticmethod
    def suffix_flag(flags):
        """
        Converts binary flags into suffix used in distribution signal names.

        Parameters
        ----------
        flags : int
            Binary flag (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004).

        Returns
        -------
        s : str
            Suffix string used in distribution signal names.

        """
        s = ""
        if flags & FLAG_AUTO:
            s += "_auto"
        if flags & FLAG_SEX:
            s += "_sex"
        if flags & FLAG_MT:
            s += "_mt"
        return s

    def chromosomes_with_signal(self, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """
        Returns list of chromosomes with signal stored in CNVpytor file

        Parameters
        ----------
        bin_size : int or None
            Bin size or None.
        signal : str
            Signal name.
        flags : int
            Binary flag
            (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004, FLAG_GC_CORR = 0x0010
            FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        chrs : list of str
            List of chromosome names.

        """
        search_string = "^" + self.signal_name("(.*)", bin_size, signal, flags) + "$"
        chrs = []
        for key in self.file.keys():
            res = re.findall(search_string, key)
            if len(res) > 0:
                chrs.append(res[0])
        return chrs

    def signal_name(self, chr_name, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """
        Returns h5py variable name for a given signal.

        Parameters
        ----------
        chr_name : str or None
            Name of the chromosome or None.
        bin_size : int or None
            Bin size or None.
        signal : str
            Signal name.
        flags : int
            Binary flag
            (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004, FLAG_GC_CORR = 0x0010
            FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        sig_name : str
            Signal h5py variable name

        """
        if signal in self.signals:
            return self.signals[signal] % {"chr": chr_name, "bin_size": bin_size, "rd_flag": self.suffix_rd_flag(flags),
                                           "snp_flag": self.suffix_snp_flag(flags), "flag": self.suffix_flag(flags)}
        else:
            logging.warning("Signal '%s' does not exists!" % signal)
            return None

    def signal_exists(self, chr_name, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """
        Checks does signal exist.

        Parameters
        ----------
        chr_name : str or None
            Name of the chromosome or None.
        bin_size : int or None
            Bin size or None.
        signal : str
            Signal name.
        flags : int
            Binary flag
            (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004, FLAG_GC_CORR = 0x0010
            FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        exists : bool
            True if signal exists in CNVpytor file

        """
        signame = self.signal_name(chr_name, bin_size, signal, flags)
        if not signame:
            return False
        return signame in self.file

    def create_signal(self, chr_name, bin_size, signal, data, flags=0):
        """
        Stores signal data into CNVpytor file and returns data set instance.

        Parameters
        ----------
        chr_name : str or None
            Name of the chromosome or None.
        bin_size : int or None
            Bin size or None.
        signal : str
            Signal name.
        data : numpy.ndarray
            Array contains data.
        flags : int
            Binary flag
            (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004, FLAG_GC_CORR = 0x0010
            FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        data_set : h5py._hl.dataset.Dataset
            Data set instance.

        """
        signame = self.signal_name(chr_name, bin_size, signal, flags)
        if not signame:
            return None
        if signame in self.file:
            del self.file[signame]
        ds = self.file.create_dataset(signame, data.shape, dtype=str(data.dtype), compression="gzip",
                                      compression_opts=9, data=data)
        self._flush()
        return ds

    def update_signal(self, chr_name, bin_size, signal, data, flags=0):
        """
        Updates signal data in CNVpytor file and returns data set instance.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome or None.
        bin_size : int
            Bin size or None.
        signal : str
            Signal name.
        data : numpy.ndarray
            Array contains data.
        flags : int
            Binary flag
            (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004, FLAG_GC_CORR = 0x0010
            FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        data_set : h5py._hl.dataset.Dataset
            Data set instance.

        """
        signame = self.signal_name(chr_name, bin_size, signal, flags)
        if not signame:
            return None
        if not (signame in self.file):
            _logger.warning("Signal %s does not exist in file %s!" % (signame, self.filename))
            return None
        self.file[signame] = data
        self._flush()
        return self.file[signame]

    def get_signal(self, chr_name, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """
        Reads signal data from CNVpytor file and returns pointer to data set.

        Parameters
        ----------
        chr_name : str or None
            Name of the chromosome or None.
        bin_size : int or None
            Bin size or None.
        signal : str
            Signal name.
        flags : int
            Binary flag
            (FLAG_AUTO = 0x0001, FLAG_SEX = 0x0002, FLAG_MT = 0x0004, FLAG_GC_CORR = 0x0010
            FLAG_AT_CORR = 0x0020, FLAG_USEMASK = 0x0100, FLAG_USEID = 0x0200, FLAG_USEHAP = 0x0400).

        Returns
        -------
        array : numpy.nparray
            Array contains data.

        """
        signame = self.signal_name(chr_name, bin_size, signal, flags)
        if not signame:
            return None
        if not (signame in self.file):
            logging.debug("Signal '%s' does not exist in file '%s'!" % (signame, self.filename))
            return []
        return np.array(self.file[signame])

    def _flush(self):
        """
        Flush pyh5 file.

        Returns
        -------
        None

        """
        self.file.flush()

    def ls(self):
        """
        Prints content of CNVpytor file.

        Returns
        -------
        None

        """
        print("Filename '%s':" % self.filename)
        print(", ".join(self.file.keys()))
        print()

    @staticmethod
    def save_root_trees(root_filename):
        """
        Save RD and VCF data into root file. Requires ROOT installed.

        Parameters
        ----------
        root_filename : sr
            Name of the root file.

        Returns
        -------
        None

        Not implemented yet!

        """
        try:
            import ROOT
        except ImportError:
            logging.warning("ROOT package is not installed - root file not saved!")
        else:
            # Save RD and SNP data into root file - TODO
            _logger.debug(root_filename, ROOT.__version__)

    def save_rd(self, chr_name, rd_p, rd_u):
        """
        Compress and stores RD data into CNVpytor file and returns data set instances.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        rd_p : numpy.ndarray
            Array with RD parity data.
        rd_u : numpy.ndarray
            Array with RD unique data.

        Returns
        -------
        ds_p : h5py._hl.dataset.Dataset
            Data set instance with RD parity signal.
        ds_u : h5py._hl.dataset.Dataset
            Data set instance with RD unique signal.

        """
        data_type = "uint32" if Genome.is_mt_chrom(chr_name) else "uint16"
        crd_p, crd_u = rd_compress(rd_p, rd_u, data_type)
        ds_p = self.create_signal(chr_name, None, "RD p", crd_p)
        ds_u = self.create_signal(chr_name, None, "RD u", crd_u)
        if not (chr_name in self.rd_chromosomes()):
            rd_chroms = self.rd_chromosomes()
            rd_chroms.append(chr_name)
            self.create_signal(None, None, "RD chromosomes", np.array([np.string_(x) for x in rd_chroms]))
        return ds_p, ds_u

    def save_vcf(self, chr_name, pos, ref, alt, nref, nalt, gt, flag, qual):
        """
        Compress and stores SNP data into CNVpytor file.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        pos : list of int
            List of SNP positions.
        ref : list of str
            List of SNP reference base (A, T, G, C or .).
        alt : list of str
            List of SNP alternative base (A, T, G, C or .).
        nref : list of int
            Count of reads contains reference SNP.
        nalt : list of int
            Count of reads contains alternative SNP.
        gt : list of int
            List of genotypes (0 - "0/0", 1 - "0/1", 3- "1/1", 4 - "0|0" , 5 - "0|1", 6 - "1|0", 7 - "1|1").
        flag : list of int
            Binary flag: first bit 1 - SNP exists in database, second bit 1 - SNP in P region of strict mask.
        qual : list of int
            SNP quality (scale 0 - 255).

        Returns
        -------
            None
        """
        snp_pos, snp_desc, snp_counts, snp_qual = snp_compress(pos, ref, alt, nref, nalt, gt, flag, qual)
        self.create_signal(chr_name, None, "SNP pos", snp_pos)
        self.create_signal(chr_name, None, "SNP desc", snp_desc)
        self.create_signal(chr_name, None, "SNP counts", snp_counts)
        self.create_signal(chr_name, None, "SNP qual", snp_qual)
        if not (chr_name in self.snp_chromosomes()):
            snp_chroms = self.snp_chromosomes()
            snp_chroms.append(chr_name)
            self.create_signal(None, None, "SNP chromosomes", np.array([np.string_(x) for x in snp_chroms]))

    def read_rd(self, chr_name):
        """
        Reads RD signals

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.

        Returns
        -------
        rd_p : numpy.ndarray
            Array with RD parity data.
        rd_u : numpy.ndarray
            Array with RD unique data.

        """
        crd_p = self.get_signal(chr_name, None, "RD p")
        crd_u = self.get_signal(chr_name, None, "RD u")
        rd_p, rd_u = rd_decompress(crd_p, crd_u)
        return rd_p, rd_u

    def rd_chromosomes(self):
        """
        Lists all chromosomes with RD signal stored in CNVpytor file.

        Returns
        -------
        chrs : list of str
            List of chromosome names with RD signal.

        """
        return list(np.array(self.get_signal(None, None, "RD chromosomes")).astype("str"))

    def gc_chromosomes(self):
        """
        Lists all chromosomes with GC/AT content data stored in CNVpytor file.

        Returns
        -------
        chrs : list of str
            List of chromosome names with GC/AT content data.

        """
        return self.chromosomes_with_signal(None, "GC/AT")

    def snp_chromosomes(self):
        """
        Lists all chromosomes with SNP signal stored in CNVpytor file.

        Returns
        -------
        chrs : list of str
            List of chromosome names with SNP signal.

        """
        return list(np.array(self.get_signal(None, None, "SNP chromosomes")).astype("str"))

    def mask_chromosomes(self):
        """
        Lists all chromosomes with strict P mask stored in CNVpytor file.

        Returns
        -------
        chrs : list of str
            List of chromosome names with strict P mask.

        """
        return self.chromosomes_with_signal(None, "mask")

    def rd_chromosome_name(self, name):
        """
        Finds name of the chromosome used for RD signal variable name.

        Parameters
        ----------
        name : str
            Chromosome name

        Returns
        -------
        chr : str or None
            Synonym for provided chromosome name used for RD signal.
            If such chromosome does not exist returns None.

        """
        rdcs = self.rd_chromosomes()
        if name in rdcs:
            return name
        if Genome.extended_chrom_name(name) in rdcs:
            return Genome.extended_chrom_name(name)
        if Genome.canonical_chrom_name(name) in rdcs:
            return Genome.canonical_chrom_name(name)
        for rdc in rdcs:
            if Genome.canonical_chrom_name(name) == Genome.canonical_chrom_name(rdc):
                return rdc
        return None

    def snp_chromosome_name(self, name):
        """
        Finds name of the chromosome used for SNP signal variable name.

        Parameters
        ----------
        name : str
            Chromosome name

        Returns
        -------
        chr : str or None
            Synonym for provided chromosome name used for SNP signal.
            If such chromosome does not exist returns None.

        """
        snpcs = self.snp_chromosomes()
        if name in snpcs:
            return name
        if Genome.extended_chrom_name(name) in snpcs:
            return Genome.extended_chrom_name(name)
        if Genome.canonical_chrom_name(name) in snpcs:
            return Genome.canonical_chrom_name(name)
        for snpc in snpcs:
            if Genome.canonical_chrom_name(name) == Genome.canonical_chrom_name(snpc):
                return snpc
        return None
