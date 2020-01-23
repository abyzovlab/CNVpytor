""" cnvpytor.io

class IO: Reading/writing CNVpytor files (extension .pytor) using h5py library.

"""
from __future__ import absolute_import, print_function, division
from .genome import Genome
from .utils import *
from .version import __version__
import datetime
import logging
import os.path
import io
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


class Signals(object):
    signals = {
        "RD p": "%(chr)s_rd_p",
        "RD u": "%(chr)s_rd_u",
        "GC/AT": "%(chr)s_gc",
        "mask": "%(chr)s_mask",
        "GC corr": "gc_corr_%(bin_size)d%(flag)s",
        "RD p dist": "dist_rd_p_%(bin_size)d%(flag)s",
        "RD u dist": "dist_rd_u_%(bin_size)d%(flag)s",
        "RD GC dist": "dist_rd_gc_%(bin_size)d%(flag)s",
        "RD stat": "rd_stat_%(bin_size)d%(flag)s",
        "RD": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s",
        "RD unique": "his_rd_u_%(chr)s_%(bin_size)d%(rd_flag)s",
        "RD raw": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_raw",
        "RD l1": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l1",
        "RD l2": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l2",
        "RD l3": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l3",
        "RD partition": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s",
        "RD call": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s_merge",
        "RD mosaic segments": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s_mosaic_segments",
        "RD mosaic call": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s_mosaic_call",
        "RD level": "rd_level_%(bin_size)d%(flag)s",
        "GC": "%(chr)s_gc_%(bin_size)",
        "SNP pos": "%(chr)s_snp_pos",
        "SNP desc": "%(chr)s_snp_desc",
        "SNP counts": "%(chr)s_snp_counts",
        "SNP qual": "%(chr)s_snp_qual",
        "SNP bin count 0|0": "snp_bafc_00_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin count 0|1": "snp_bafc_01_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin count 1|0": "snp_bafc_10_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin count 1|1": "snp_bafc_11_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin reads 0|0": "snp_readc_00_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin reads 0|1": "snp_readc_01_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin reads 1|0": "snp_readc_10_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP bin reads 1|1": "snp_readc_11_%(chr)s_%(bin_size)d%(snp_flag)s",
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
        "SNP likelihood segments": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s_segments",
        "SNP likelihood call": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP maf call": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i1 call": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i2 call": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i3 call": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i4 call": "snp_i4_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "somatic SNP pos": "somatic_%(name)s_%(chr)s_snp_pos",
        "somatic SNP desc": "somatic_%(name)s_%(chr)s_snp_desc",
        "somatic SNP counts": "somatic_%(name)s_%(chr)s_snp_counts",
        "somatic SNP qual": "somatic_%(name)s_%(chr)s_snp_qual",
        "RD chromosomes": "rd_chromosomes",
        "SNP chromosomes": "snp_chromosomes",
        "chromosome lengths": "chr_len",
        "read frg dist": "read_frg_len",
        "reference genome": "reference_genome",
        "use reference": "use_reference"
    }

    def __init__(self):
        pass

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
        if flags & FLAG_GC_CORR:
            s += "_GC"
        return s

    def signal_name(self, chr_name, bin_size, signal, flags=0, name=''):
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
            try:
                return self.signals[signal] % {"chr": chr_name, "bin_size": bin_size,
                                               "rd_flag": self.suffix_rd_flag(flags),
                                               "snp_flag": self.suffix_snp_flag(flags), "flag": self.suffix_flag(flags),
                                               "name": name}
            except TypeError:
                return None
        else:
            return None


class IO(Signals):

    def __init__(self, filename, ro=False, buffer=False, create=True):
        """
        Opens CNVpytor file for reading/writing

        Parameters
        ----------
        filename : str
            Name of the file.
        ro : bool
            Opens file in read-only mode. Default: False.
        buffer : bool
            It will copy hdf5 file in RAM buffer before opening. Works with read-only mode. Default: False.
        create : bool
            It will create file when set and when file does not exist. Otherwise, if file does not exist
            it will log an error. Default: True.

        """
        Signals.__init__(self)
        self.filename = filename
        self.file = None
        _logger.debug("Opening h5 file '%s'" % self.filename)
        if ro:
            try:
                if buffer:
                    with open(filename, 'rb') as f:
                        self.bytesio = io.BytesIO(f.read())
                    self.file = h5py.File(self.bytesio, "r")
                else:
                    self.file = h5py.File(filename, "r")
                _logger.debug("File '%s' successfully opened in read-only mode." % self.filename)
            except IOError:
                _logger.error("Unable to open file %s!" % filename)
                exit(0)
        elif os.path.exists(filename):
            try:
                self.file = h5py.File(filename, "r+")
                _logger.debug("File '%s' successfully opened." % self.filename)
            except IOError:
                _logger.error("Unable to open file %s!" % filename)
                exit(0)
        elif create:
            try:
                self.file = h5py.File(filename, "w")
                now = datetime.datetime.now()
                # Meta data
                self.add_meta_attribute('Version', __version__)
                self.add_meta_attribute('Date', now.strftime("%Y-%m-%d %H:%M"))
                _logger.debug("File '%s' successfully created." % self.filename)
            except IOError:
                _logger.error("Unable to create file %s!" % filename)
                exit(0)
        else:
            _logger.error("File %s is missing!" % filename)
            exit(0)

    def __del__(self):
        _logger.debug("Closing h5 file '%s'" % self.filename)
        if self.file:
            self.file.close()

    def chromosomes_with_signal(self, bin_size, signal, flags=0, name=''):
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
        search_string = "^" + self.signal_name("(.[^_]*)", bin_size, signal, flags, name) + "$"
        chrs = []
        for key in self.file.keys():
            res = re.findall(search_string, key)
            if len(res) > 0:
                chrs.append(res[0])
        return chrs

    def chromosomes_bin_sizes_with_signal(self, signal, flags=0, name=''):
        """
        Returns list of chromosome bin_size pairs with signal stored in CNVpytor file

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
        chrs_bss : list of (str, int)
            List of tuples (chromosome name, bin size).

        """
        search_string = "^" + self.signal_name("(.[^_]*)", 17110806, signal, flags, name) + "$"
        search_string = search_string.replace("17110806", "(.[0-9]*)")
        chrs_bss = []
        for key in self.file.keys():
            res = re.findall(search_string, key)
            if len(res) > 0:
                chrs_bss.append(res[0])
        return chrs_bss

    def signal_exists(self, chr_name, bin_size, signal, flags=0, name=''):
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
        signame = self.signal_name(chr_name, bin_size, signal, flags, name)
        if not signame:
            return False
        return signame in self.file

    def create_signal(self, chr_name, bin_size, signal, data, flags=0, name=''):
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
        signame = self.signal_name(chr_name, bin_size, signal, flags, name)
        if not signame:
            return None
        if signame in self.file:
            del self.file[signame]
        ds = self.file.create_dataset(signame, data.shape, dtype=str(data.dtype), compression="gzip",
                                      compression_opts=9, data=data)
        self._flush()
        return ds

    def update_signal(self, chr_name, bin_size, signal, data, flags=0, name=''):
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
        signame = self.signal_name(chr_name, bin_size, signal, flags, name)
        if not signame:
            return None
        if not (signame in self.file):
            _logger.warning("Signal %s does not exist in file %s!" % (signame, self.filename))
            return None
        self.file[signame] = data
        self._flush()
        return self.file[signame]

    def get_signal(self, chr_name, bin_size, signal, flags=0, name=''):
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
        signame = self.signal_name(chr_name, bin_size, signal, flags, name)
        if not signame:
            return None
        if not (signame in self.file):
            _logger.debug("Signal '%s' does not exist in file '%s'!" % (signame, self.filename))
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
        print()
        print("Filename '%s'" % self.filename)
        print("-----------" + "-" * len(self.filename))

        parameter_list = ['Date', 'Version']
        if parameter_list and self.file.attrs.keys():
            date = self.file.attrs['Date']
            version = self.file.attrs['Version']
            print("File created: {} using CNVpytor ver {}\n".format(date, version))

        print("Chromosomes with RD signal: " + ", ".join(self.rd_chromosomes()))
        print()
        print("Chromosomes with SNP signal: " + ", ".join(self.snp_chromosomes()))
        print()
        if self.signal_exists(None, None, "reference genome") and self.signal_exists(None, None, "use reference"):
            rg_name = np.array(self.get_signal(None, None, "reference genome")).astype("str")[0]
            print("Using reference genome: " + rg_name + " [ ", end='')
            if rg_name in Genome.reference_genomes:
                rg_use = self.get_signal(None, None, "use reference")
                if "gc_file" in Genome.reference_genomes[rg_name] and rg_use[0] == 1:
                    print("GC: yes, ", end='')
                else:
                    print("GC: no, ", end='')
                if "mask_file" in Genome.reference_genomes[rg_name] and rg_use[1] == 1:
                    print("mask: yes ]")
                else:
                    print("mask: no ]")
        else:
            print("Reference genome is not set.")
        print()

        chr_bs = self.chromosomes_bin_sizes_with_signal("RD")
        chrs = {}
        bss = []
        for c, b in chr_bs:
            if c not in chrs:
                chrs[c] = []
            chrs[c].append(int(b))
            if int(b) not in bss:
                bss.append(int(b))

        print("Chromosomes with RD histograms [bin sizes]: " + ", ".join(chrs.keys()) + " " + str(sorted(bss)))
        print()
        chr_bs = self.chromosomes_bin_sizes_with_signal("SNP likelihood", FLAG_USEMASK)
        chrs = {}
        bss = []
        for c, b in chr_bs:
            if c not in chrs:
                chrs[c] = []
            chrs[c].append(int(b))
            if int(b) not in bss:
                bss.append(int(b))

        print("Chromosomes with SNP histograms [bin sizes]: " + ", ".join(chrs.keys()) + " " + str(sorted(bss)))
        print()
        chr_len = list(np.array(self.get_signal(None, None, "chromosome lengths")).astype("str"))
        chr_len = dict(zip(chr_len[::2], chr_len[1::2]))
        print("Chromosome lengths: " + str(chr_len))

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

    def save_rd(self, chr_name, rd_p, rd_u, chromosome_length=None):
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
        _logger.info("Saving chromosome RD data for chromosome '%s'." % chr_name)
        data_type = "uint32" if Genome.is_mt_chrom(chr_name) else "uint16"
        crd_p, crd_u = rd_compress(rd_p, rd_u, data_type)
        snp_name = self.snp_chromosome_name(chr_name)
        if not (snp_name is None):
            if snp_name == chr_name:
                _logger.info("Detecting SNP data in file '%s' for the same chromosome." % self.filename)
            else:
                _logger.info(
                    "Detecting RD data in file '%s' for the same chromosome with different name '%s'. SNP name will be used." % (
                        self.filename, snp_name))
                chr_name = snp_name
        if chromosome_length is not None:
            self.set_chromosome_length(chr_name, chromosome_length)
        ds_p = self.create_signal(chr_name, None, "RD p", crd_p)
        ds_u = self.create_signal(chr_name, None, "RD u", crd_u)
        if not (chr_name in self.rd_chromosomes()):
            rd_chroms = self.rd_chromosomes()
            rd_chroms.append(chr_name)
            self.create_signal(None, None, "RD chromosomes", np.array([np.string_(x) for x in rd_chroms]))
        return ds_p, ds_u

    def save_snp(self, chr_name, pos, ref, alt, nref, nalt, gt, flag, qual, update=False, callset=None,
                 chromosome_length=None):
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
        if callset is None:
            _logger.info("Saving SNP data for chromosome '%s'." % chr_name)
        else:
            _logger.info("Saving somatic '%s' SNV data for chromosome '%s'." % (callset, chr_name))
        snp_pos, snp_desc, snp_counts, snp_qual = snp_compress(pos, ref, alt, nref, nalt, gt, flag, qual)
        rd_name = self.rd_chromosome_name(chr_name)
        if not update and not (rd_name is None):
            if rd_name == chr_name:
                _logger.info("Detecting RD data in file '%s' for the same chromosome." % self.filename)
            else:
                _logger.info(
                    "Detecting RD data in file '%s' for the same chromosome with different name '%s'. RD name will be used." % (
                        self.filename, rd_name))
                chr_name = rd_name
        if not self.is_chromosome_length_set(chr_name):
            if chromosome_length is not None:
                self.set_chromosome_length(chr_name, chromosome_length)
            elif chromosome_length is None:
                self.set_chromosome_length(chr_name, pos[-1] + 1)

        if callset is None:
            self.create_signal(chr_name, None, "SNP pos", snp_pos)
            self.create_signal(chr_name, None, "SNP desc", snp_desc)
            self.create_signal(chr_name, None, "SNP counts", snp_counts)
            self.create_signal(chr_name, None, "SNP qual", snp_qual)
        else:
            self.create_signal(chr_name, None, "somatic SNP pos", snp_pos, name=callset)
            self.create_signal(chr_name, None, "somatic SNP desc", snp_desc, name=callset)
            self.create_signal(chr_name, None, "somatic SNP counts", snp_counts, name=callset)
            self.create_signal(chr_name, None, "somatic SNP qual", snp_qual, name=callset)
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

    def read_snp(self, chr_name, callset=None):
        """
        Reads SNP signals

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.

        Returns
        -------
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

        """
        if callset is None:
            snp_pos = self.get_signal(chr_name, None, "SNP pos")
            snp_desc = self.get_signal(chr_name, None, "SNP desc")
            snp_counts = self.get_signal(chr_name, None, "SNP counts")
            snp_qual = self.get_signal(chr_name, None, "SNP qual")
        else:
            snp_pos = self.get_signal(chr_name, None, "somatic SNP pos", name=callset)
            snp_desc = self.get_signal(chr_name, None, "somatic SNP desc", name=callset)
            snp_counts = self.get_signal(chr_name, None, "somatic SNP counts", name=callset)
            snp_qual = self.get_signal(chr_name, None, "somatic SNP qual", name=callset)
        pos, ref, alt, nref, nalt, gt, flag, qual = snp_decompress(snp_pos, snp_desc, snp_counts, snp_qual)
        return pos, ref, alt, nref, nalt, gt, flag, qual

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

    def set_chromosome_length(self, chromosome, length):
        """

        Parameters
        ----------
        chromosome : str
            Chromosome name
        length: int
            Chromosome length

        Returns
        -------
            None

        """
        chr_len = list(np.array(self.get_signal(None, None, "chromosome lengths")).astype("str"))
        chr_len = dict(zip(chr_len[::2], chr_len[1::2]))
        chr_len[chromosome] = str(length)
        self.create_signal(None, None, "chromosome lengths",
                           np.array([np.string_(x) for s in chr_len.items() for x in s]))

    def get_chromosome_length(self, chromosome):
        """

        Parameters
        ----------
        chromosome : str
            Chromosome name

        Returns
        -------
        len : int or None
            Chromosome length

        """
        chr_len = list(np.array(self.get_signal(None, None, "chromosome lengths")).astype("str"))
        chr_len = dict(zip(chr_len[::2], chr_len[1::2]))
        if chromosome in chr_len:
            return int(chr_len[chromosome])
        return None

    def is_chromosome_length_set(self, chromosome):
        """

        Parameters
        ----------
        chromosome : str
            Chromosome name

        Returns
        -------
        bool
            True if chromosome length is set

        """
        chr_len = list(np.array(self.get_signal(None, None, "chromosome lengths")).astype("str"))
        chr_len = dict(zip(chr_len[::2], chr_len[1::2]))
        return chromosome in chr_len

    def add_meta_attribute(self, attribute, value):
        self.file.attrs[attribute] = str(value)

    def read_meta_attribute(self):
        for attribute, value in self.file.attrs.items():
            print("{}: {}".format(attribute, value))
