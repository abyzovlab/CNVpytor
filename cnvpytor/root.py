""" cnvpytor.root

class Root: main CNVpytor class
"""
from __future__ import absolute_import, print_function, division

from .io import *
from .bam import Bam
from .vcf import Vcf
from .fasta import Fasta
from .genome import Genome
from .viewer import anim_plot_likelihood, anim_plot_rd
import numpy as np
import logging
import os.path

_logger = logging.getLogger("cnvpytor.root")


class Root:

    def __init__(self, filename, create=False, max_cores=8):
        """
        Class constructor opens CNVpytor file. The class implements all core CNVpytor calculations.

        Parameters
        ----------
        filename : str
            CNVpytor filename
        create : bool
            Creates cnvpytor file if not exists
        max_cores : int
            Maximal number of cores used in parallel calculations

        """
        _logger.debug("App class init: filename '%s'; max cores %d." % (filename, max_cores))
        self.io = IO(filename, create=create)
        self.max_cores = max_cores
        self.io_gc = self.io
        self.io_mask = self.io
        if self.io.signal_exists(None, None, "reference genome") and self.io.signal_exists(None, None, "use reference"):
            rg_name = np.array(self.io.get_signal(None, None, "reference genome")).astype("str")[0]
            if rg_name in Genome.reference_genomes:
                rg_use = self.io.get_signal(None, None, "use reference")
                if "gc_file" in Genome.reference_genomes[rg_name] and rg_use[0] == 1:
                    _logger.debug("Using GC content from database for reference genome '%s'." % rg_name)
                    self.io_gc = IO(Genome.reference_genomes[rg_name]["gc_file"], ro=True, buffer=True)
                if "mask_file" in Genome.reference_genomes[rg_name] and rg_use[1] == 1:
                    _logger.debug("Using strict mask from database for reference genome '%s'." % rg_name)
                    self.io_mask = IO(Genome.reference_genomes[rg_name]["mask_file"], ro=True, buffer=True)

    def _read_bam(self, bf, chroms, reference_filename=False):
        bamf = Bam(bf, reference_filename=reference_filename)
        if bamf.reference_genome:
            self.io.create_signal(None, None, "reference genome", np.array([np.string_(bamf.reference_genome)]))
            self.io.create_signal(None, None, "use reference", np.array([1, 1]).astype("uint8"))

        def read_chromosome(x):
            chrs, lens = x
            _logger.info("Reading data for chromosome %s with length %d" % (chrs, lens))
            l_rd_p, l_rd_u, l_his_read_frg = bamf.read_chromosome(chrs)
            return l_rd_p, l_rd_u, l_his_read_frg

        chrname, chrlen = bamf.get_chr_len()
        chr_len = [(c, l) for (c, l) in zip(chrname, chrlen) if len(chroms) == 0 or c in chroms]
        if self.max_cores == 1:
            count = 0
            cum_his_read_frg = None
            for cl in chr_len:
                rd_p, rd_u, his_read_frg = read_chromosome(cl)
                if not rd_p is None:
                    if cum_his_read_frg is None:
                        cum_his_read_frg = his_read_frg
                    else:
                        cum_his_read_frg += his_read_frg
                    self.io.save_rd(cl[0], rd_p, rd_u, chromosome_length=cl[1])
                    count += 1
            if not cum_his_read_frg is None:
                self.io.create_signal(None, None, "read frg dist", cum_his_read_frg)
            return count
        else:
            from .pool import parmap
            res = parmap(read_chromosome, chr_len, cores=self.max_cores)
            count = 0
            cum_his_read_frg = None
            for c, r in zip(chr_len, res):
                if not r[0] is None:
                    if cum_his_read_frg is None:
                        cum_his_read_frg = r[2]
                    else:
                        cum_his_read_frg += r[2]
                    self.io.save_rd(c[0], r[0], r[1], chromosome_length=c[1])
                    count += 1
            if not cum_his_read_frg is None:
                self.io.create_signal(None, None, "read frg dist", cum_his_read_frg)
            return count

    def rd_stat(self, chroms=[]):
        """ Calculate RD statistics for 100 bp bins.

        Parameters
        ----------
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.

        """
        rd_stat = []
        auto, sex, mt = False, False, False

        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_gc_chromosomes[rd_name] = c

        for c in rd_gc_chromosomes:
            if (len(chroms) == 0 or c in chroms):
                rd_p, rd_u = self.io.read_rd(c)
                max_bin = max(int(10 * np.mean(rd_u) + 1), int(10 * np.mean(rd_p) + 1))
                dist_p, bins = np.histogram(rd_p, bins=range(max_bin + 1))
                dist_u, bins = np.histogram(rd_u, bins=range(max_bin + 1))
                n_p, m_p, s_p = fit_normal(bins[1:-1], dist_p[1:])[0]
                n_u, m_u, s_u = fit_normal(bins[1:-1], dist_u[1:])[0]
                _logger.info(
                    "Chromosome '%s' stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (c, m_p, s_p))
                _logger.info(
                    "Chromosome '%s' stat - RD unique distribution gaussian fit:  %.2f +- %.2f" % (c, m_u, s_u))
                rd_stat.append((c, m_p, s_p, m_u, s_u))
                auto = auto or Genome.is_autosome(c)
                sex = sex or Genome.is_sex_chrom(c)
                mt = mt or Genome.is_mt_chrom(c)

        if auto:
            max_bin_auto = int(
                max(map(lambda x: 5 * x[1] + 5 * x[2], filter(lambda x: Genome.is_autosome(x[0]), rd_stat))))
            _logger.debug("Max RD for autosomes calculated: %d" % max_bin_auto)
            dist_p_auto = np.zeros(max_bin_auto)
            dist_u_auto = np.zeros(max_bin_auto)
            dist_p_gc_auto = np.zeros((max_bin_auto, 101))
            bins_auto = range(max_bin_auto + 1)
        if sex:
            max_bin_sex = int(
                max(map(lambda x: 5 * x[1] + 5 * x[2], filter(lambda x: Genome.is_sex_chrom(x[0]), rd_stat))))
            _logger.debug("Max RD for sex chromosomes calculated: %d" % max_bin_sex)
            dist_p_sex = np.zeros(max_bin_sex)
            dist_u_sex = np.zeros(max_bin_sex)
            dist_p_gc_sex = np.zeros((max_bin_sex, 101))
            bins_sex = range(max_bin_sex + 1)

        if mt:
            max_bin_mt = int(
                max(map(lambda x: 5 * x[1] + 5 * x[2], filter(lambda x: Genome.is_mt_chrom(x[0]), rd_stat))))
            _logger.debug("Max RD for mitochondria calculated: %d" % max_bin_mt)
            if max_bin_mt < 100:
                max_bin_mt = 100
            bin_size_mt = max_bin_mt // 100
            bins_mt = range(0, max_bin_mt // bin_size_mt * bin_size_mt + bin_size_mt, bin_size_mt)
            n_bins_mt = len(bins_mt) - 1
            _logger.debug("Using %d bin size, %d bins for mitochondria chromosome." % (bin_size_mt, n_bins_mt))
            dist_p_mt = np.zeros(n_bins_mt)
            dist_u_mt = np.zeros(n_bins_mt)
            dist_p_gc_mt = np.zeros((n_bins_mt, 101))

        for c in rd_gc_chromosomes:
            if len(chroms) == 0 or c in chroms:
                _logger.info("Chromosome '%s' stat " % c)
                rd_p, rd_u = self.io.read_rd(c)
                if Genome.is_autosome(c):
                    dist_p, bins = np.histogram(rd_p, bins=bins_auto)
                    dist_u, bins = np.histogram(rd_u, bins=bins_auto)
                    gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                    gcp = gcp_decompress(gcat)
                    dist_p_gc, xbins, ybins = np.histogram2d(rd_p, gcp, bins=(bins_auto, range(102)))
                    dist_p_auto += dist_p
                    dist_u_auto += dist_u
                    dist_p_gc_auto += dist_p_gc
                elif Genome.is_sex_chrom(c):
                    dist_p, bins = np.histogram(rd_p, bins=bins_sex)
                    dist_u, bins = np.histogram(rd_u, bins=bins_sex)
                    gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                    gcp = gcp_decompress(gcat)
                    dist_p_gc, xbins, ybins = np.histogram2d(rd_p, gcp, bins=(bins_sex, range(102)))
                    dist_p_sex += dist_p
                    dist_u_sex += dist_u
                    dist_p_gc_sex += dist_p_gc
                elif Genome.is_mt_chrom(c):
                    dist_p, bins = np.histogram(rd_p, bins=bins_mt)
                    dist_u, bins = np.histogram(rd_u, bins=bins_mt)
                    gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                    gcp = gcp_decompress(gcat)
                    dist_p_gc, xbins, ybins = np.histogram2d(rd_p, gcp, bins=(bins_mt, range(102)))
                    dist_p_mt += dist_p
                    dist_u_mt += dist_u
                    dist_p_gc_mt += dist_p_gc

        if auto:
            n_auto, m_auto, s_auto = fit_normal(np.array(bins_auto[1:-1]), dist_p_auto[1:])[0]
            _logger.info("Autosomes stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (m_auto, s_auto))
            stat_auto = np.array([max_bin_auto, 1, max_bin_auto, n_auto, m_auto, s_auto])
            self.io.create_signal(None, 100, "RD p dist", dist_p_auto, flags=FLAG_AUTO)
            self.io.create_signal(None, 100, "RD u dist", dist_u_auto, flags=FLAG_AUTO)
            self.io.create_signal(None, 100, "RD GC dist", dist_p_gc_auto, flags=FLAG_AUTO)
            self.io.create_signal(None, 100, "RD stat", stat_auto, flags=FLAG_AUTO)
            gc_corr_auto = calculate_gc_correction(dist_p_gc_auto, m_auto, s_auto)
            self.io.create_signal(None, 100, "GC corr", gc_corr_auto, flags=FLAG_AUTO)

        if sex:
            n_sex, m_sex, s_sex = fit_normal(np.array(bins_sex[1:-1]), dist_p_sex[1:])[0]
            _logger.info("Sex chromosomes stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (m_sex, s_sex))
            stat_sex = np.array([max_bin_sex, 1, max_bin_sex, n_sex, m_sex, s_sex])
            self.io.create_signal(None, 100, "RD p dist", dist_p_sex, flags=FLAG_SEX)
            self.io.create_signal(None, 100, "RD u dist", dist_u_sex, flags=FLAG_SEX)
            self.io.create_signal(None, 100, "RD GC dist", dist_p_gc_sex, flags=FLAG_SEX)
            self.io.create_signal(None, 100, "RD stat", stat_sex, flags=FLAG_SEX)
            gc_corr_sex = calculate_gc_correction(dist_p_gc_sex, m_sex, s_sex)
            self.io.create_signal(None, 100, "GC corr", gc_corr_sex, flags=FLAG_SEX)

        if mt:
            n_mt, m_mt, s_mt = fit_normal(np.array(bins_mt[1:-1]), dist_p_mt[1:])[0]
            _logger.info("Mitochondria stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (m_mt, s_mt))
            if auto:
                _logger.info("Mitochondria stat - number of mitochondria per cell:  %.2f +- %.2f" % (
                    2. * m_mt / m_auto, 2. * s_mt / m_auto + s_auto * m_mt / (m_auto * m_auto)))
            stat_mt = np.array([max_bin_mt, bin_size_mt, n_bins_mt, n_mt, m_mt, s_mt])
            self.io.create_signal(None, 100, "RD p dist", dist_p_mt, flags=FLAG_MT)
            self.io.create_signal(None, 100, "RD u dist", dist_u_mt, flags=FLAG_MT)
            self.io.create_signal(None, 100, "RD GC dist", dist_p_gc_mt, flags=FLAG_MT)
            self.io.create_signal(None, 100, "RD stat", stat_mt, flags=FLAG_MT)
            gc_corr_mt = calculate_gc_correction(dist_p_gc_mt, m_mt, s_mt, bin_size_mt)
            self.io.create_signal(None, 100, "GC corr", gc_corr_mt, flags=FLAG_MT)

    def _read_vcf(self, vcf_file, chroms, sample='', use_index=False, no_counts=False, ad_tag="AD", gt_tag="GT",
                  callset=None):
        vcff = Vcf(vcf_file)
        chrs = [c for c in vcff.get_chromosomes() if len(chroms) == 0 or c in chroms]

        use_index = use_index and os.path.exists(vcf_file + ".tbi")

        def save_data(chr, pos, ref, alt, nref, nalt, gt, flag, qual):
            if (len(chroms) == 0 or chr in chroms) and (not pos is None) and (len(pos) > 0):
                self.io.save_snp(chr, pos, ref, alt, nref, nalt, gt, flag, qual, callset=callset,
                                 chromosome_length=vcff.lengths[chr])
            # TODO: Stop reading if all from chrom list are read.

        def save_data_no_counts(chr, pos, ref, alt, gt, flag, qual):
            if (len(chroms) == 0 or chr in chroms) and (not pos is None) and (len(pos) > 0):
                self.io.save_snp(chr, pos, ref, alt, np.zeros_like(pos), np.zeros_like(pos), gt, flag, qual,
                                 callset=callset, chromosome_length=vcff.lengths[chr])

        if use_index:
            count = 0
            for c in chrs:
                _logger.info("Reading variant data for chromosome %s" % c)
                if no_counts:
                    pos, ref, alt, gt, flag, qual = vcff.read_chromosome_snp_no_counts(c, sample, gt_tag=gt_tag)
                    nref, nalt = np.zeros_like(pos), np.zeros_like(pos)
                else:
                    pos, ref, alt, nref, nalt, gt, flag, qual = vcff.read_chromosome_snp(c, sample, ad_tag=ad_tag,
                                                                                         gt_tag=gt_tag)

                if not pos is None and len(pos) > 0:
                    self.io.save_snp(c, pos, ref, alt, nref, nalt, gt, flag, qual, callset=callset,
                                     chromosome_length=vcff.lengths[c])
                    count += 1
            return count
        else:
            if no_counts:
                return vcff.read_all_snp_no_counts(save_data_no_counts, sample, gt_tag=gt_tag)
            else:
                return vcff.read_all_snp(save_data, sample, ad_tag=ad_tag, gt_tag=gt_tag)

    def rd(self, bamfiles, chroms=[], reference_filename=False):
        """ Read chromosomes from bam/sam/cram file(s) and store in cnvpytor file.

        Parameters
        ----------
        bamfiles : list of str
            List of bam/sam/cram files
        chroms : list of str
            List of chromosomes. Import all available if empty.
        reference_filename : str
            Only for CRAM files - reference genome filename

        """
        for bf in bamfiles:
            self._read_bam(bf, chroms, reference_filename=reference_filename)
            self.io.add_meta_attribute("BAM", bf)

        if self.io.signal_exists(None, None, "reference genome"):
            rg_name = np.array(self.io.get_signal(None, None, "reference genome")).astype("str")[0]
            if "mask_file" in Genome.reference_genomes[rg_name]:
                _logger.info("Strict mask for reference genome '%s' found in database." % rg_name)
                self.io_mask = IO(Genome.reference_genomes[rg_name]["mask_file"])
            if "gc_file" in Genome.reference_genomes[rg_name]:
                _logger.info("GC content for reference genome '%s' found in database." % rg_name)
                self.io_gc = IO(Genome.reference_genomes[rg_name]["gc_file"])
                self.rd_stat()

    def vcf(self, vcf_files, chroms=[], sample='', no_counts=False, ad_tag="AD", gt_tag="GT", callset=None,
            use_index=True):
        """ Read SNP data from variant file(s) and store in .cnvpytor file

        Parameters
        ----------
        vcf_files : list of str
            List of variant filenames.
        chroms : list of str
            List of chromosomes. Import all available if empty.
        sample : str
            Name of the sample. It will read first sample if empty string is provided.
        no_counts : bool
            Do not read AD counts if true.
        ad_tag : str
            AD tag (default 'AD').
        gt_tag : str
            GT tag (default 'GT').
        callset : str or None
            It will assume SNP data if None. Otherwise it will assume SNV data and
            store under the name provided by callset variable.
        use_index : bool
            Use index file for vcf parsing. Default is False.

        """
        for vcf_file in vcf_files:
            self._read_vcf(vcf_file, chroms, sample, no_counts=no_counts, ad_tag=ad_tag, gt_tag=gt_tag,
                           callset=callset, use_index=use_index)
            self.io.add_meta_attribute("VCF", vcf_file)

    def rd_from_vcf(self, vcf_file, chroms=[], sample='', ad_tag="AD", dp_tag="DP", use_index=False):
        """
        Read RD from variant file(s) and store in .cnvpytor file

        Parameters
        ----------
        vcf_file : str
            Variant filename
        chroms : list of str
            List of chromosomes. Imports all available if empty.
        sample : str
            Name of the sample. It will read first sample if empty string is provided.
        ad_tag : str
            AD tag (default 'AD')
        dp_tag : str
            DP tag (default 'DP')
        use_index : bool
            Use index file for vcf parsing. Default is False.

        """
        vcff = Vcf(vcf_file)
        chrs = [c for c in vcff.get_chromosomes() if len(chroms) == 0 or c in chroms]

        def save_data(chr, rd_p, rd_u):
            if (len(chroms) == 0 or chr in chroms) and (not rd_p is None) and (len(rd_p) > 0):
                self.io.save_rd(chr, rd_p, rd_u, chromosome_length=vcff.lengths[chr])
            # TODO: Stop reading if all from chroms list are already read.

        if use_index:
            count = 0
            for c in chrs:
                _logger.info("Reading variant data for chromosome %s" % c)

                rd_p, rd_u = vcff.read_chromosome_rd(c, sample, ad_tag=ad_tag, dp_tag=dp_tag)
                if not rd_p is None and len(rd_p) > 0:
                    self.io.save_rd(chr, rd_p, rd_u, chromosome_length=vcff.lengths[chr])
                    count += 1
            self.io.add_meta_attribute("RD from VCF", vcf_file)
            return count
        else:
            count = vcff.read_all_rd(save_data, sample, ad_tag=ad_tag, dp_tag=dp_tag)
            self.io.add_meta_attribute("RD from VCF", vcf_file)
            return count

    def _pileup_bam(self, bamfile, chroms, pos, ref, alt, nref, nalt, reference_filename):
        _logger.info("Calculating pileup from bam file '%s'." % bamfile)
        bamf = Bam(bamfile, reference_filename=reference_filename)

        def pileup_chromosome(x):
            c, l, snpc = x
            _logger.info("Pileup chromosome %s with length %d" % (c, l))
            return bamf.pileup(c, pos[snpc], ref[snpc], alt[snpc])

        chrname, chrlen = bamf.get_chr_len()
        chr_len = [(c, l, self.io.snp_chromosome_name(c)) for (c, l) in zip(chrname, chrlen) if
                   self.io.snp_chromosome_name(c) in chroms]

        if self.max_cores == 1:
            for cl in chr_len:
                r = pileup_chromosome(cl)
                nref[cl[2]] = [x + y for x, y in zip(nref[cl[2]], r[0])]
                nalt[cl[2]] = [x + y for x, y in zip(nalt[cl[2]], r[1])]

        else:
            from .pool import parmap
            res = parmap(pileup_chromosome, chr_len, cores=self.max_cores)
            for cl, r in zip(chr_len, res):
                nref[cl[2]] = [x + y for x, y in zip(nref[cl[2]], r[0])]
                nalt[cl[2]] = [x + y for x, y in zip(nalt[cl[2]], r[1])]

    def pileup(self, bamfiles, chroms=[], reference_filename=False):
        """
        Read bam/sam/cram file(s) and pile up SNP counts at positions imported with vcf(vcf_files, ...) method

        Parameters
        ----------
        bamfiles : list of str
            Bam files.
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        reference_filename : str
            Only for CRAM files - reference genome filename
        """
        chrs = [c for c in self.io.snp_chromosomes() if (len(chroms) == 0 or c in chroms)]

        pos, ref, alt, nref, nalt, gt, flag, qual = {}, {}, {}, {}, {}, {}, {}, {}
        for c in chrs:
            _logger.info("Decompressing and setting all SNP counts to zero for chromosome '%s'." % c)
            pos[c], ref[c], alt[c], nref[c], nalt[c], gt[c], flag[c], qual[c] = self.io.read_snp(c)
            nref[c] = [0] * len(nref[c])
            nalt[c] = [0] * len(nalt[c])
        for bf in bamfiles:
            self._pileup_bam(bf, chrs, pos, ref, alt, nref, nalt, reference_filename=reference_filename)

        for c in chrs:
            _logger.info("Saving SNP data for chromosome '%s' in file '%s'." % (c, self.io.filename))
            self.io.save_snp(c, pos[c], ref[c], alt[c], nref[c], nalt[c], gt[c], flag[c], qual[c])

    def rd_from_snp(self, chroms=[], callset=None):
        """ Create RD signal from already imported SNP signal

        Parameters
        ----------
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        callset : str or None
            It will assume SNP data if None. Otherwise it will assume SNV data
            stored under the name provided by callset variable.

        Returns
        -------

        """

        chromosomes = [c for c in self.io.snp_chromosomes() if (len(chroms) == 0 or c in chroms)]

        snp_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            snp_name = self.io.snp_chromosome_name(c)
            if snp_name is not None:
                snp_mask_chromosomes[snp_name] = c

        for c in chromosomes:
            _logger.info("Calculating RD signal for chromosome '%s' using SNP counts." % c)
            pos, ref, alt, nref, nalt, gt, flag, qual = self.io.read_snp(c, callset=callset)
            n = self.io.get_chromosome_length(c) // 100 + 1
            rd = np.zeros(n)
            ncg = self.io.get_chromosome_length(c) // 10000 + 1
            rdcg = np.zeros(ncg)
            rdc = np.zeros(ncg)
            for p, c1, c2 in zip(pos, nref, nalt):
                rdcg[(p - 1) // 10000] += c1 + c2
                rdc[(p - 1) // 10000] += 1
            rdcg = rdcg / rdc
            for i in range(n):
                rd[i] = rdcg[i // 100]
            self.io.save_rd(c, rd, rd)

    def gc(self, filename, chroms=[], make_gc_genome_file=False):
        """ Read GC content from reference fasta file and store in .cnvnator file
                Arguments:
                    * FASTA filename
                    * list of chromosome names
        """
        fasta = Fasta(filename)
        chrname, chrlen = fasta.get_chr_len()
        chr_len_rdn = []
        if make_gc_genome_file:
            _logger.debug("GC data for all reference genome contigs will be read and stored in cnvnator file.")
            for (c, l) in zip(chrname, chrlen):
                chr_len_rdn.append((c, l, c))
        else:
            for (c, l) in zip(chrname, chrlen):
                rd_name = self.io.rd_chromosome_name(c)
                if len(chroms) == 0 or (rd_name in chroms) or (c in chroms):
                    if rd_name is None:
                        _logger.debug("Not found RD signal for chromosome '%s'. Skipping..." % c)
                    else:
                        _logger.debug("Found RD signal for chromosome '%s' with name '%s'." % (c, rd_name))
                        chr_len_rdn.append((c, l, rd_name))
        count = 0
        for c, l, rdn in chr_len_rdn:
            _logger.info("Reading data for chromosome '%s' with length %d" % (c, l))
            gc, at = fasta.read_chromosome_gc(c)
            if not gc is None:
                self.io.create_signal(rdn, None, "GC/AT", gc_at_compress(gc, at))
                count += 1
                _logger.info("GC data for chromosome '%s' saved in file '%s'." % (rdn, self.io.filename))
        if count > 0:
            _logger.info("Running RD statistics on chromosomes with imported GC data.")
            self.rd_stat()
        return count

    def copy_gc(self, filename, chroms=[]):
        """ Copy GC content from another cnvnator file
                Arguments:
                    * cnvnator filename
                    * list of chromosome names
        """
        io_src = IO(filename)
        chr_rdn = []
        for c in io_src.chromosomes_with_signal(None, "GC/AT"):
            rd_name = self.io.rd_chromosome_name(c)
            if len(chroms) == 0 or (rd_name in chroms) or (c in chroms):
                if rd_name is None:
                    _logger.debug("Not found RD signal for chromosome '%s'. Skipping..." % c)
                else:
                    _logger.debug("Found RD signal for chromosome '%s' with name '%s'." % (c, rd_name))
                    chr_rdn.append((c, rd_name))
        count = 0
        for c, rdn in chr_rdn:
            _logger.info(
                "Copying GC data for chromosome '%s' from file '%s' in file '%s'." % (c, filename, self.io.filename))
            self.io.create_signal(rdn, None, "GC/AT", io_src.get_signal(c, None, "GC/AT").astype(dtype="uint8"))
        return count

    def mask(self, filename, chroms=[], make_mask_genome_file=False):
        """ Read strict mask fasta file and store in .cnvnator file
                Arguments:
                    * FASTA filename
                    * list of chromosome names
        """
        fasta = Fasta(filename)
        chrname, chrlen = fasta.get_chr_len()
        chr_len_rdn = []
        if make_mask_genome_file:
            _logger.debug("Strict mask data for all reference genome contigs will be read and stored in cnvnator file.")
            for (c, l) in zip(chrname, chrlen):
                chr_len_rdn.append((c, l, c))
        else:
            for (c, l) in zip(chrname, chrlen):
                rd_name = self.io.rd_chromosome_name(c)
                snp_name = self.io.snp_chromosome_name(c)
                if len(chroms) == 0 or (rd_name in chroms) or (snp_name in chroms) or (c in chroms):
                    if (rd_name is None) and (snp_name is None):
                        _logger.debug("Not found RD or SNP signal for chromosome '%s'. Skipping..." % c)
                    elif (rd_name is None):
                        _logger.debug("Found SNP signal for chromosome '%s' with name '%s'." % (c, snp_name))
                        chr_len_rdn.append((c, l, snp_name))
                    else:
                        _logger.debug("Found RD signal for chromosome '%s' with name '%s'." % (c, rd_name))
                        chr_len_rdn.append((c, l, rd_name))
        count = 0
        for c, l, rdn in chr_len_rdn:
            _logger.info("Reading data for chromosome '%s' with length %d" % (c, l))
            mask = fasta.read_chromosome_mask_p_regions(c)
            if not mask is None:
                self.io.create_signal(rdn, None, "mask", mask_compress(mask))
                count += 1
                _logger.info("Strict mask data for chromosome '%s' saved in file '%s'." % (rdn, self.io.filename))
        return count

    def copy_mask(self, filename, chroms=[]):
        """ Copy strict mask data from another cnvnator file
                Arguments:
                    * cnvnator filename
                    * list of chromosome names
        """
        io_src = IO(filename)
        chr_rdn = []
        for c in io_src.chromosomes_with_signal(None, "mask"):
            rd_name = self.io.rd_chromosome_name(c)
            snp_name = self.snp_chromosome_name(c)
            if len(chroms) == 0 or (rd_name in chroms) or (snp_name in chroms) or (c in chroms):
                if (rd_name is None) and (snp_name is None):
                    _logger.debug("Not found RD or SNP signal for chromosome '%s'. Skipping..." % c)
                elif (rd_name is None):
                    _logger.debug("Found SNP signal for chromosome '%s' with name '%s'." % (c, snp_name))
                    chr_rdn.append((c, snp_name))
                else:
                    _logger.debug("Found RD signal for chromosome '%s' with name '%s'." % (c, rd_name))
                    chr_rdn.append((c, rd_name))
        count = 0
        for c, rdn in chr_rdn:
            _logger.info(
                "Copying strict mask data for chromosome '%s' from file '%s' in file '%s'." % (
                    c, filename, self.io.filename))
            self.io.create_signal(rdn, None, "mask", io_src.get_signal(c, None, "mask").astype(dtype="uint32"))
        return count

    def set_reference_genome(self, rg):
        """ Manually sets reference genome.

        Parameters
        ----------
        rg: str
            Name of reference genome

        Returns
        -------
            True if genome exists in database

        """
        if rg in Genome.reference_genomes:
            _logger.info("Reference genome '%s' found in database." % rg)
            self.io.create_signal(None, None, "reference genome", np.array([np.string_(rg)]))
            self.io.create_signal(None, None, "use reference", np.array([1, 1]).astype("uint8"))
            if "mask_file" in Genome.reference_genomes[rg]:
                _logger.info("Strict mask for reference genome '%s' found in database." % rg)
            if "gc_file" in Genome.reference_genomes[rg]:
                _logger.info("GC content for reference genome '%s' found in database." % rg)
        else:
            _logger.warning("Reference genome '%s' not found in database." % rg)

    def calculate_histograms(self, bin_sizes, chroms=[]):
        """
        Calculates RD histograms and store data into cnvpytor file.

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.

        """
        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None and (len(chroms) == 0 or (rd_name in chroms) or (c in chroms)):
                rd_gc_chromosomes[rd_name] = c

        rd_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None and (len(chroms) == 0 or (rd_name in chroms) or (c in chroms)):
                rd_mask_chromosomes[rd_name] = c

        for bin_size in bin_sizes:
            his_stat = []
            auto, sex, mt = False, False, False
            bin_ratio = bin_size // 100

            for c in self.io.rd_chromosomes():
                if c in rd_gc_chromosomes:
                    _logger.info("Calculating histograms using bin size %d for chromosome '%s'." % (bin_size, c))
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    rd_p, rd_u = self.io.read_rd(c)
                    his_p = np.concatenate(
                        (rd_p, np.zeros(bin_ratio - len(rd_p) + (len(rd_p) // bin_ratio * bin_ratio))))
                    his_p.resize((len(his_p) // bin_ratio, bin_ratio))
                    his_p = his_p.sum(axis=1)
                    his_u = np.concatenate(
                        (rd_u, np.zeros(bin_ratio - len(rd_u) + (len(rd_u) // bin_ratio * bin_ratio))))
                    his_u.resize((len(his_u) // bin_ratio, bin_ratio))
                    his_u = his_u.sum(axis=1)
                    if bin_ratio == 1:
                        his_p = his_p[:-1]
                        his_u = his_u[:-1]
                    self.io.create_signal(c, bin_size, "RD", his_p)
                    self.io.create_signal(c, bin_size, "RD unique", his_u)

                    max_bin = max(int(10 * np.mean(his_u) + 1), int(10 * np.mean(his_p) + 1))
                    if max_bin < 10000:
                        max_bin = 10000
                    rd_bin_size = max_bin // 10000
                    rd_bins = range(0, max_bin // rd_bin_size * rd_bin_size + rd_bin_size, rd_bin_size)
                    dist_p, bins = np.histogram(his_p, bins=rd_bins)
                    dist_u, bins = np.histogram(his_u, bins=rd_bins)
                    n_p, m_p, s_p = fit_normal(bins[1:-1], dist_p[1:])[0]
                    n_u, m_u, s_u = fit_normal(bins[1:-1], dist_u[1:])[0]
                    _logger.info(
                        "Chromosome '%s' bin size %d stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (
                            c, bin_size, m_p, s_p))
                    _logger.info(
                        "Chromosome '%s' bin size %d stat - RD unique distribution gaussian fit:  %.2f +- %.2f" % (
                            c, bin_size, m_u, s_u))
                    his_stat.append((c, m_p, s_p, m_u, s_u))
                    auto = auto or Genome.is_autosome(c)
                    sex = sex or Genome.is_sex_chrom(c)
                    mt = mt or Genome.is_mt_chrom(c)

            mt = mt and (bin_size <= 1000)

            if auto:
                max_bin_auto = int(
                    max(map(lambda x: 5 * x[1] + 5 * x[2], filter(lambda x: Genome.is_autosome(x[0]), his_stat))))
                _logger.debug("Max RD for autosome histograms calculated: %d." % max_bin_auto)
                if max_bin_auto < 1000:
                    max_bin_auto = 1000
                bin_size_auto = max_bin_auto // 1000
                bins_auto = range(0, max_bin_auto // bin_size_auto * bin_size_auto + bin_size_auto, bin_size_auto)
                n_bins_auto = len(bins_auto) - 1
                _logger.debug(
                    "Using %d RD bin size, %d bins for autosome RD - GC distributions." % (bin_size_auto, n_bins_auto))
                dist_p_auto = np.zeros(n_bins_auto)
                dist_p_gccorr_auto = np.zeros(n_bins_auto)
                dist_u_auto = np.zeros(n_bins_auto)
                dist_p_gc_auto = np.zeros((n_bins_auto, 101))

            if sex:
                max_bin_sex = int(
                    max(map(lambda x: 5 * x[1] + 5 * x[2], filter(lambda x: Genome.is_sex_chrom(x[0]), his_stat))))
                _logger.debug("Max RD for sex chromosome histograms calculated: %d." % max_bin_sex)
                if max_bin_sex < 1000:
                    max_bin_sex = 1000
                bin_size_sex = max_bin_sex // 1000
                bins_sex = range(0, max_bin_sex // bin_size_sex * bin_size_sex + bin_size_sex, bin_size_sex)
                n_bins_sex = len(bins_sex) - 1
                _logger.debug(
                    "Using %d RD bin size, %d bins for sex chromosome RD - GC distributions." % (
                        bin_size_sex, n_bins_sex))
                dist_p_sex = np.zeros(n_bins_sex)
                dist_p_gccorr_sex = np.zeros(n_bins_sex)
                dist_u_sex = np.zeros(n_bins_sex)
                dist_p_gc_sex = np.zeros((n_bins_sex, 101))

            if mt:
                max_bin_mt = int(
                    max(map(lambda x: 5 * x[1] + 5 * x[2], filter(lambda x: Genome.is_mt_chrom(x[0]), his_stat))))
                _logger.debug("Max RD for mitochondria histogram calculated: %d." % max_bin_mt)
                if max_bin_mt < 1000:
                    max_bin_mt = 1000
                bin_size_mt = max_bin_mt // 1000
                bins_mt = range(0, max_bin_mt // bin_size_mt * bin_size_mt + bin_size_mt, bin_size_mt)
                n_bins_mt = len(bins_mt) - 1
                _logger.debug("Using %d bin size, %d bins for mitochondria chromosome." % (bin_size_mt, n_bins_mt))
                dist_p_mt = np.zeros(n_bins_mt)
                dist_p_gccorr_mt = np.zeros(n_bins_mt)
                dist_u_mt = np.zeros(n_bins_mt)
                dist_p_gc_mt = np.zeros((n_bins_mt, 101))

            _logger.info("Calculating global statistics.")
            for c in self.io.rd_chromosomes():
                if c in rd_gc_chromosomes:
                    _logger.debug("Chromosome '%s'." % c)
                    his_p = self.io.get_signal(c, bin_size, "RD")
                    his_u = self.io.get_signal(c, bin_size, "RD unique")
                    if Genome.is_autosome(c):
                        dist_p, bins = np.histogram(his_p, bins=bins_auto)
                        dist_u, bins = np.histogram(his_u, bins=bins_auto)
                        gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                        dist_p_gc, xbins, ybins = np.histogram2d(his_p, gcp_decompress(gcat, bin_ratio),
                                                                 bins=(bins_auto, range(102)))
                        dist_p_auto += dist_p
                        dist_u_auto += dist_u
                        dist_p_gc_auto += dist_p_gc
                    elif Genome.is_sex_chrom(c):
                        dist_p, bins = np.histogram(his_p, bins=bins_sex)
                        dist_u, bins = np.histogram(his_u, bins=bins_sex)
                        gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                        dist_p_gc, xbins, ybins = np.histogram2d(his_p, gcp_decompress(gcat, bin_ratio),
                                                                 bins=(bins_sex, range(102)))
                        dist_p_sex += dist_p
                        dist_u_sex += dist_u
                        dist_p_gc_sex += dist_p_gc
                    elif Genome.is_mt_chrom(c) and mt:
                        dist_p, bins = np.histogram(his_p, bins=bins_mt)
                        dist_u, bins = np.histogram(his_u, bins=bins_mt)
                        gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                        dist_p_gc, xbins, ybins = np.histogram2d(his_p, gcp_decompress(gcat, bin_ratio),
                                                                 bins=(bins_mt, range(102)))
                        dist_p_mt += dist_p
                        dist_u_mt += dist_u
                        dist_p_gc_mt += dist_p_gc

            if auto:
                n_auto, m_auto, s_auto = fit_normal(np.array(bins_auto[1:-1]), dist_p_auto[1:])[0]
                _logger.info("Autosomes stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (m_auto, s_auto))
                stat_auto = np.array([max_bin_auto, bin_size_auto, n_bins_auto, n_auto, m_auto, s_auto])
                self.io.create_signal(None, bin_size, "RD p dist", dist_p_auto, flags=FLAG_AUTO)
                self.io.create_signal(None, bin_size, "RD u dist", dist_u_auto, flags=FLAG_AUTO)
                self.io.create_signal(None, bin_size, "RD GC dist", dist_p_gc_auto, flags=FLAG_AUTO)
                self.io.create_signal(None, bin_size, "RD stat", stat_auto, flags=FLAG_AUTO)
                gc_corr_auto = calculate_gc_correction(dist_p_gc_auto, m_auto, s_auto, bin_size_auto)
                self.io.create_signal(None, bin_size, "GC corr", gc_corr_auto, flags=FLAG_AUTO)

            if sex:
                n_sex, m_sex, s_sex = fit_normal(np.array(bins_sex[1:-1]), dist_p_sex[1:])[0]
                _logger.info(
                    "Sex chromosomes stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (m_sex, s_sex))
                stat_sex = np.array([max_bin_sex, bin_size_sex, n_bins_sex, n_sex, m_sex, s_sex])
                self.io.create_signal(None, bin_size, "RD p dist", dist_p_sex, flags=FLAG_SEX)
                self.io.create_signal(None, bin_size, "RD u dist", dist_u_sex, flags=FLAG_SEX)
                self.io.create_signal(None, bin_size, "RD GC dist", dist_p_gc_sex, flags=FLAG_SEX)
                self.io.create_signal(None, bin_size, "RD stat", stat_sex, flags=FLAG_SEX)
                gc_corr_sex = calculate_gc_correction(dist_p_gc_sex, m_sex, s_sex, bin_size_sex)
                self.io.create_signal(None, bin_size, "GC corr", gc_corr_sex, flags=FLAG_SEX)

            if mt:
                n_mt, m_mt, s_mt = fit_normal(np.array(bins_mt[1:-1]), dist_p_mt[1:])[0]
                _logger.info("Mitochondria stat - RD parity distribution gaussian fit:  %.2f +- %.2f" % (m_mt, s_mt))
                if auto:
                    _logger.info("Mitochondria stat - number of mitochondria per cell:  %.2f +- %.2f" % (
                        2. * m_mt / m_auto, 2. * s_mt / m_auto + s_auto * m_mt / (m_auto * m_auto)))
                stat_mt = np.array([max_bin_mt, bin_size_mt, n_bins_mt, n_mt, m_mt, s_mt])
                self.io.create_signal(None, bin_size, "RD p dist", dist_p_mt, flags=FLAG_MT)
                self.io.create_signal(None, bin_size, "RD u dist", dist_u_mt, flags=FLAG_MT)
                self.io.create_signal(None, bin_size, "RD GC dist", dist_p_gc_mt, flags=FLAG_MT)
                self.io.create_signal(None, bin_size, "RD stat", stat_mt, flags=FLAG_MT)
                gc_corr_mt = calculate_gc_correction(dist_p_gc_mt, m_mt, s_mt, bin_size_mt)
                self.io.create_signal(None, bin_size, "GC corr", gc_corr_mt, flags=FLAG_MT)

            for c in self.io.rd_chromosomes():
                if c in rd_gc_chromosomes and (mt or not Genome.is_mt_chrom(c)):
                    _logger.info(
                        "Calculating GC corrected RD histogram using bin size %d for chromosome '%s'." % (bin_size, c))
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    his_p = self.io.get_signal(c, bin_size, "RD")
                    _logger.debug("Calculating GC corrected RD")
                    gc_corr = self.io.get_signal(None, bin_size, "GC corr", flag)
                    gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                    his_p_corr = his_p / np.array(list(map(lambda x: gc_corr[int(x)], gcp_decompress(gcat, bin_ratio))))
                    self.io.create_signal(c, bin_size, "RD", his_p_corr, flags=FLAG_GC_CORR)
                    if Genome.is_autosome(c):
                        dist_p_gc, bins = np.histogram(his_p_corr, bins=bins_auto)
                        dist_p_gccorr_auto += dist_p_gc
                    elif Genome.is_sex_chrom(c):
                        dist_p_gc, bins = np.histogram(his_p_corr, bins=bins_sex)
                        dist_p_gccorr_sex += dist_p_gc
                    elif Genome.is_mt_chrom(c) and mt:
                        dist_p_gc, bins = np.histogram(his_p_corr, bins=bins_mt)
                        dist_p_gccorr_mt += dist_p_gc

            if auto:
                n_auto, m_auto, s_auto = fit_normal(np.array(bins_auto[1:-1]), dist_p_gccorr_auto[1:])[0]
                _logger.info(
                    "Autosomes stat - RD parity distribution gaussian fit after GC correction:  %.2f +- %.2f" % (
                        m_auto, s_auto))
                stat_auto = np.array([max_bin_auto, bin_size_auto, n_bins_auto, n_auto, m_auto, s_auto])
                self.io.create_signal(None, bin_size, "RD p dist", dist_p_gccorr_auto, flags=(FLAG_AUTO | FLAG_GC_CORR))
                self.io.create_signal(None, bin_size, "RD stat", stat_auto, flags=(FLAG_AUTO | FLAG_GC_CORR))

            if sex:
                n_sex, m_sex, s_sex = fit_normal(np.array(bins_sex[1:-1]), dist_p_gccorr_sex[1:])[0]
                _logger.info(
                    "Sex chromosomes stat - RD parity distribution gaussian fit after GC correction:  %.2f +- %.2f" % (
                        m_sex, s_sex))
                stat_sex = np.array([max_bin_sex, bin_size_sex, n_bins_sex, n_sex, m_sex, s_sex])
                self.io.create_signal(None, bin_size, "RD p dist", dist_p_gccorr_sex, flags=(FLAG_SEX | FLAG_GC_CORR))
                self.io.create_signal(None, bin_size, "RD stat", stat_sex, flags=(FLAG_SEX | FLAG_GC_CORR))

            if mt:
                n_mt, m_mt, s_mt = fit_normal(np.array(bins_mt[1:-1]), dist_p_gccorr_mt[1:])[0]
                _logger.info(
                    "Mitochondria stat - RD parity distribution gaussian fit after GC correction:  %.2f +- %.2f" % (
                        m_mt, s_mt))
                if auto:
                    _logger.info(
                        "Mitochondria stat - number of mitochondria per cell after GC correction:  %.2f +- %.2f" % (
                            2. * m_mt / m_auto, 2. * s_mt / m_auto + s_auto * m_mt / (m_auto * m_auto)))
                stat_mt = np.array([max_bin_mt, bin_size_mt, n_bins_mt, n_mt, m_mt, s_mt])
                self.io.create_signal(None, bin_size, "RD p dist", dist_p_gccorr_mt, flags=(FLAG_MT | FLAG_GC_CORR))
                self.io.create_signal(None, bin_size, "RD stat", stat_mt, flags=(FLAG_MT | FLAG_GC_CORR))

            for c in self.io.rd_chromosomes():
                if (c in rd_gc_chromosomes) and (c in rd_mask_chromosomes):
                    _logger.info("Calculating P-mask histograms using bin size %d for chromosome '%s'." % (bin_size, c))
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    rd_p, rd_u = self.io.read_rd(c)
                    mask = mask_decompress(self.io_mask.get_signal(rd_mask_chromosomes[c], None, "mask"))
                    bin_ratio = bin_size // 100
                    n_bins = len(rd_p) // bin_ratio + 1
                    his_p = np.zeros(n_bins)
                    his_p_corr = np.zeros(n_bins)
                    his_u = np.zeros(n_bins)
                    his_count = np.zeros(n_bins)
                    for (b, e) in mask:
                        for p in range((b - 1) // 100 + 2, (e - 1) // 100):
                            his_p[p // bin_ratio] += rd_p[p]
                            his_u[p // bin_ratio] += rd_u[p]
                            his_count[p // bin_ratio] += 1
                    np.seterr(divide='ignore', invalid='ignore')
                    his_p = his_p / his_count * bin_ratio
                    his_u = his_u / his_count * bin_ratio
                    if bin_ratio == 1:
                        his_p = his_p[:-1]
                        his_u = his_u[:-1]
                    gc_corr = self.io.get_signal(None, bin_size, "GC corr", flag)
                    gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                    his_p_corr = his_p / np.array(list(map(lambda x: gc_corr[int(x)], gcp_decompress(gcat, bin_ratio))))
                    self.io.create_signal(c, bin_size, "RD", his_p, flags=FLAG_USEMASK)
                    self.io.create_signal(c, bin_size, "RD unique", his_u, flags=FLAG_USEMASK)
                    self.io.create_signal(c, bin_size, "RD", his_p_corr, flags=FLAG_GC_CORR | FLAG_USEMASK)

    def partition(self, bin_sizes, chroms=[], use_gc_corr=True, use_mask=False, repeats=3, genome_size=2.9e9):
        """
        Calculates segmentation of RD signal.

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        use_gc_corr : bool
            Use GC corrected signal if True. Default: True.
        use_mask : bool
            Use P-mask filter if True. Default: False.

        """
        bin_bands = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128]

        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_gc_chromosomes[rd_name] = c

        rd_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_mask_chromosomes[rd_name] = c

        for bin_size in bin_sizes:
            for c in self.io.rd_chromosomes():
                if (c in rd_gc_chromosomes or not use_gc_corr) and (c in rd_mask_chromosomes or not use_mask) and (
                        len(chroms) == 0 or (c in chroms)):
                    flag_stat = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    if use_gc_corr:
                        flag_stat |= FLAG_GC_CORR
                    flag_rd = (FLAG_GC_CORR if use_gc_corr else 0) | (FLAG_USEMASK if use_mask else 0)
                    if self.io.signal_exists(c, bin_size, "RD stat", flag_stat) and self.io.signal_exists(c, bin_size,
                                                                                                          "RD",
                                                                                                          flag_rd):
                        _logger.info("Calculating partition using bin size %d for chromosome '%s'." % (bin_size, c))
                        stat = self.io.get_signal(c, bin_size, "RD stat", flag_stat)
                        mean = stat[4]
                        std = stat[5]
                        rd = self.io.get_signal(c, bin_size, "RD", flag_rd)
                        rd = np.nan_to_num(rd)
                        masked = np.zeros_like(rd, dtype=bool)
                        levels = np.copy(rd)

                        for bin_band in bin_bands:
                            _logger.debug("Bin band is %d." % bin_band)
                            levels[np.logical_not(masked)] = rd[np.logical_not(masked)]
                            nm_levels = levels[np.logical_not(masked)]
                            mask_borders = [0]
                            count = 0
                            for i in range(len(masked)):
                                if masked[i]:
                                    if count > 0:
                                        mask_borders.append(mask_borders[-1] + count - 1)
                                        count = 0
                                else:
                                    count += 1

                            kk = np.arange(3 * bin_band + 1)
                            exp_kk = kk * np.exp(-0.5 * kk ** 2 / bin_band ** 2)
                            for step in range(repeats):
                                isig = np.ones_like(nm_levels) * 4. / std ** 2
                                isig[nm_levels >= (mean / 4)] = mean / std ** 2 / nm_levels[nm_levels >= (mean / 4)]

                                def calc_grad(k):
                                    if k < len(nm_levels):
                                        return exp_kk[k] * np.exp(
                                            np.concatenate(
                                                (-0.5 * ((nm_levels - np.roll(nm_levels, -k)) ** 2 * isig)[:-k],
                                                 [0] * k))) - np.concatenate(([0] * k, exp_kk[
                                            k] * np.exp(-0.5 * ((nm_levels - np.roll(nm_levels, k)) ** 2 * isig)[k:])))
                                    else:
                                        return np.zeros_like(nm_levels)

                                if True:
                                    grad = np.sum([calc_grad(k) for k in range(1, 3 * bin_band + 1)], axis=0)
                                else:
                                    import multiprocessing
                                    from multiprocessing import sharedctypes
                                    grad = np.ctypeslib.as_ctypes(np.zeros(len(nm_levels)))
                                    shared_array = sharedctypes.RawArray(grad._type_, grad)

                                    def calc_grad_par(ks):
                                        tmp = np.ctypeslib.as_array(shared_array)
                                        for k in ks:
                                            if k < len(nm_levels):
                                                tmp += exp_kk[k] * np.exp(
                                                    np.concatenate(
                                                        (-0.5 * ((nm_levels - np.roll(nm_levels, -k)) ** 2 * isig)[:-k],
                                                         [0] * k))) - np.concatenate(([0] * k, exp_kk[
                                                    k] * np.exp(
                                                    -0.5 * ((nm_levels - np.roll(nm_levels, k)) ** 2 * isig)[k:])))

                                    jobs = []
                                    nj = self.max_cores

                                    for i in range(nj):
                                        ks = [k for k in range(1, 3 * bin_band + 1) if k % nj == i]
                                        if len(ks) > 0:
                                            process = multiprocessing.Process(target=calc_grad_par, args=(ks,))
                                            jobs.append(process)
                                            process.start()
                                    for j in jobs:
                                        j.join()
                                    grad = np.ctypeslib.as_array(shared_array)

                                border = [i for i in range(grad.size - 1) if grad[i] < 0 and grad[i + 1] >= 0]
                                border.append(grad.size - 1)
                                border = sorted(list(set(border + mask_borders)))
                                pb = 0
                                for b in border:
                                    nm_levels[pb:b + 1] = np.mean(nm_levels[pb:b + 1])
                                    pb = b + 1

                            levels[np.logical_not(masked)] = nm_levels
                            border = [0] + list(np.argwhere(np.abs(np.diff(levels)) > 0.01)[:, 0] + 1) + [levels.size]
                            masked = np.zeros_like(rd, dtype=bool)

                            for i in range(1, len(border)):
                                seg = [border[i - 1], border[i]]
                                seg_left = [border[i - 1], border[i - 1]]
                                if i > 1:
                                    seg_left[0] = border[i - 2]
                                else:
                                    continue
                                seg_right = [border[i], border[i]]
                                if i < (len(border) - 1):
                                    seg_right[1] = border[i + 1]
                                else:
                                    continue
                                n = seg[1] - seg[0]
                                n_left = seg_left[1] - seg_left[0]
                                n_right = seg_right[1] - seg_right[0]
                                if n <= 1:
                                    continue
                                seg_mean = np.mean(levels[seg[0]:seg[1]])
                                seg_std = np.std(levels[seg[0]:seg[1]])
                                if (n_right <= 15) or (n_left <= 15) or (n <= 15):
                                    ns = 1.8 * np.sqrt(levels[seg_left[0]] / mean) * std
                                    if np.abs(levels[seg_left[0]] - levels[seg[0]]) < ns:
                                        continue
                                    ns = 1.8 * np.sqrt(levels[seg_right[0]] / mean) * std
                                    if np.abs(levels[seg_right[0]] - levels[seg[0]]) < ns:
                                        continue
                                else:
                                    seg_left_mean = np.mean(levels[seg_left[0]:seg_left[1]])
                                    seg_left_std = np.std(levels[seg_left[0]:seg_left[1]])
                                    seg_right_mean = np.mean(levels[seg_right[0]:seg_right[1]])
                                    seg_right_std = np.std(levels[seg_right[0]:seg_right[1]])
                                    if t_test_2_samples(seg_mean, seg_std, n, seg_left_mean, seg_left_std, n_left) > (
                                            0.01 / genome_size * bin_size * (n + n_left)):
                                        continue
                                    if t_test_2_samples(seg_mean, seg_std, n, seg_right_mean, seg_right_std,
                                                        n_right) > (
                                            0.01 / genome_size * bin_size * (n + n_right)):
                                        continue
                                if t_test_1_sample(mean, seg_mean, seg_std, n) > 0.05:
                                    continue
                                masked[seg[0]:seg[1]] = True
                                levels[seg[0]:seg[1]] = np.mean(rd[seg[0]:seg[1]])

                        self.io.create_signal(c, bin_size, "RD partition", levels, flags=flag_rd)

    def call(self, bin_sizes, chroms=[], print_calls=False, use_gc_corr=True, use_mask=False, genome_size=2.9e9,
             genome_cnv_fraction=0.01):
        """
        CNV caller based on the segmented RD signal.

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        print_calls : bool
            Print to stdout list of calls if true.
        use_gc_corr : bool
            Use GC corrected signal if True. Default: True.
        use_mask : bool
            Use P-mask filter if True. Default: False.

        Returns
        -------
        calls: dist
            Dictionary bin_size -> list of calls
            Each call is list with values:
                type ("deletion" or "duplication"),
                chromosome name, start, end,
                size, cnv level
                e1, e2, e3, e4 - estimations of p-Value
                q0 - percentage of uniquely mapped reads in cnv region

        """
        ret = {}
        normal_genome_size = genome_size * (1 - genome_cnv_fraction)
        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_gc_chromosomes[rd_name] = c

        rd_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_mask_chromosomes[rd_name] = c

        for bin_size in bin_sizes:
            ret[bin_size] = []
            for c in self.io.rd_chromosomes():
                if (c in rd_gc_chromosomes or not use_gc_corr) and (c in rd_mask_chromosomes or not use_mask) and (
                        len(chroms) == 0 or (c in chroms)):
                    flag_stat = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    flag_auto = FLAG_AUTO
                    if use_gc_corr:
                        flag_stat |= FLAG_GC_CORR
                        flag_auto |= FLAG_GC_CORR
                    flag_rd = (FLAG_GC_CORR if use_gc_corr else 0) | (FLAG_USEMASK if use_mask else 0)
                    if self.io.signal_exists(c, bin_size, "RD stat", flag_stat) and \
                            self.io.signal_exists(c, bin_size, "RD", flag_rd) and \
                            self.io.signal_exists(c, bin_size, "RD partition", flag_rd):
                        _logger.debug("Calculating CNV calls using bin size %d for chromosome '%s'." % (bin_size, c))
                        stat = self.io.get_signal(c, bin_size, "RD stat", flag_stat)
                        mean = stat[4]
                        std = stat[5]
                        rd = self.io.get_signal(c, bin_size, "RD", flag_rd)
                        rd = np.nan_to_num(rd)
                        levels = self.io.get_signal(c, bin_size, "RD partition", flag_rd)
                        delta = 0.25
                        if Genome.is_sex_chrom(c) and self.io.signal_exists(c, bin_size, "RD stat", flag_auto):
                            stat_auto = self.io.get_signal(c, bin_size, "RD stat", flag_auto)
                            if stat_auto[4] * 0.66 > mean:
                                _logger.debug("Assuming male individual!")
                                delta *= 2

                        _logger.debug("Merging levels with relative difference smaller than %f.", delta)
                        delta *= mean
                        done = False
                        while not done:
                            done = True
                            # border - list of borders between segments
                            border = [0] + list(np.argwhere(np.abs(np.diff(levels)) > 0.01)[:, 0] + 1) + [levels.size]
                            for ix in range(len(border) - 2):
                                if ix < len(border) - 2:
                                    v1 = np.abs(levels[border[ix]] - levels[border[ix + 1]])
                                    if v1 < delta:
                                        v2 = v1 + 1
                                        v3 = v1 + 1
                                        if ix > 0:
                                            v2 = np.abs(levels[border[ix]] - levels[border[ix - 1]])
                                        if ix < len(border) - 3:
                                            v3 = np.abs(levels[border[ix + 1]] - levels[border[ix + 2]])
                                        if v1 < v2 and v1 < v3:
                                            done = False
                                            levels[border[ix]:border[ix + 2]] = np.mean(
                                                levels[border[ix]:border[ix + 2]])
                                            del border[ix + 1]

                        _logger.debug("Calling segments")
                        min = mean - delta
                        max = mean + delta

                        flags = [""] * len(levels)
                        segments = []

                        b = 0
                        while b < len(levels):
                            b0 = b
                            bs = b
                            while b < len(levels) and levels[b] < min:
                                b += 1
                            be = b
                            if be > bs + 1:
                                adj = adjustToEvalue(mean, std, rd, bs, be, 0.05 * bin_size / normal_genome_size)
                                if adj is not None:
                                    bs, be = adj
                                    segments.append([bs, be, -1])
                                    flags[bs:be] = ["D"] * (be - bs)

                            bs = b
                            while b < len(levels) and levels[b] > max:
                                b += 1
                            be = b
                            if be > bs + 1:
                                adj = adjustToEvalue(mean, std, rd, bs, be, 0.05 * bin_size / normal_genome_size)
                                if adj is not None:
                                    bs, be = adj
                                    segments.append([bs, be, +1])
                                    flags[bs:be] = ["A"] * (be - bs)
                            if b == b0:
                                b += 1

                        _logger.debug("Calling additional deletions")
                        b = 0
                        while b < len(levels):
                            while b < len(levels) and flags[b] != "":
                                b += 1
                            bs = b
                            while b < len(levels) and levels[b] < min:
                                b += 1
                            be = b
                            if be > bs + 1:
                                if gaussianEValue(mean, std, rd, bs, be) < 0.05 / normal_genome_size:
                                    segments.append([bs, be, -1])
                                    flags[bs:be] = ["d"] * (be - bs)
                                b -= 1
                            b += 1

                        b = 0
                        if b < len(levels):
                            cf = flags[b]
                        bs = 0
                        merge = rd.copy()
                        while b < len(levels):
                            while flags[b] == cf:
                                b += 1
                                if b >= len(flags):
                                    break
                            if b > bs:
                                merge[bs:b] = np.mean(merge[bs:b])
                            if b < len(levels):
                                cf = flags[b]
                            bs = b

                        self.io.create_signal(c, bin_size, "RD call", merge, flags=flag_rd)

                        _logger.debug("Calculate/print calls")
                        b = 0
                        while b < len(levels):
                            cf = flags[b]
                            if cf == "":
                                b += 1
                                continue
                            bs = b
                            while b < len(levels) and cf == flags[b]:
                                b += 1
                            cnv = np.mean(rd[bs:b]) / mean
                            if cf.upper() == "D":
                                etype = "deletion"
                            else:
                                etype = "duplication"
                            start = bin_size * bs + 1
                            end = bin_size * b
                            size = end - start + 1
                            e1 = getEValue(mean, std, rd, bs, b) * normal_genome_size / bin_size
                            e2 = gaussianEValue(mean, std, rd, bs, b) * normal_genome_size
                            e3, e4 = 1, 1
                            tmp = int(1000. / bin_size + 0.5)
                            if bs + tmp < b - tmp:
                                e3 = getEValue(mean, std, rd, bs + tmp, b - tmp) * normal_genome_size / bin_size
                                e4 = gaussianEValue(mean, std, rd, bs + tmp, b - tmp) * normal_genome_size
                            rd_p = self.io.get_signal(c, bin_size, "RD")
                            rd_u = self.io.get_signal(c, bin_size, "RD unique")
                            q0 = -1
                            if sum(rd_p[bs:b]) > 0:
                                q0 = (sum(rd_p[bs:b]) - sum(rd_u[bs:b])) / sum(rd_p[bs:b])
                            if print_calls:
                                print("%s\t%s:%d-%d\t%d\t%.4f\t%e\t%e\t%e\t%e\t%.4f\t" % (
                                    etype, c, start, end, size, cnv, e1, e2, e3, e4, q0))
                            ret[bin_size].append([etype, c, start, end, size, cnv, e1, e2, e3, e4, q0])
        return ret

    def call_mosaic(self, bin_sizes, chroms=[], use_gc_corr=True, use_mask=False, omin=None,
                    max_distance=0.3, anim=""):
        """
        CNV caller based on likelihood merger.

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        use_gc_corr : bool
            Use GC corrected signal if True. Default: True.
        use_mask : bool
            Use P-mask filter if True. Default: False.

        """
        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_gc_chromosomes[rd_name] = c

        rd_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_mask_chromosomes[rd_name] = c

        for bin_size in bin_sizes:
            if omin is None:
                overlap_min = 0.05 * bin_size / 3e9
            else:
                overlap_min = omin

            for c in self.io.rd_chromosomes():
                if (c in rd_gc_chromosomes or not use_gc_corr) and (c in rd_mask_chromosomes or not use_mask) and (
                        len(chroms) == 0 or (c in chroms)):
                    flag_stat = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    flag_auto = FLAG_AUTO
                    if use_gc_corr:
                        flag_stat |= FLAG_GC_CORR
                        flag_auto |= FLAG_GC_CORR
                    flag_rd = (FLAG_GC_CORR if use_gc_corr else 0) | (FLAG_USEMASK if use_mask else 0)
                    if self.io.signal_exists(c, bin_size, "RD stat", flag_stat) and \
                            self.io.signal_exists(c, bin_size, "RD", flag_rd):
                        _logger.info("Calculating mosaic calls using bin size %d for chromosome '%s'." % (bin_size, c))
                        stat = self.io.get_signal(c, bin_size, "RD stat", flag_stat)
                        mean = stat[4]
                        std = stat[5]
                        rd = self.io.get_signal(c, bin_size, "RD", flag_rd)
                        bins = len(rd)
                        valid = np.isfinite(rd)
                        level = rd[valid]
                        error = np.sqrt(level) ** 2 + std ** 2
                        loc_fl = np.min(list(zip(np.abs(np.diff(level))[:-1], np.abs(np.diff(level))[1:])), axis=1)
                        loc_fl = np.concatenate(([0], loc_fl, [0]))
                        error += (loc_fl / 2) ** 2
                        error = np.sqrt(error)
                        level = list(level)
                        error = list(error)
                        segments = [[i] for i in range(bins) if np.isfinite(rd[i])]
                        overlaps = [normal_overlap(level[i], error[i], level[i + 1], error[i + 1]) for i in
                                    range(len(segments) - 1)]

                        iter = 0
                        if anim != "":
                            anim_plot_rd(level, error, segments, bins, iter, anim + c + "_0_" + str(bin_size), 0,
                                         0, mean)

                        while len(overlaps) > 0:
                            maxo = max(overlaps)
                            if maxo < overlap_min:
                                break
                            i = overlaps.index(maxo)
                            nl, ne = normal_merge(level[i], error[i], level[i + 1], error[i + 1])
                            level[i] = nl
                            error[i] = ne
                            segments[i] += segments[i + 1]
                            del level[i + 1]
                            del error[i + 1]
                            del segments[i + 1]
                            del overlaps[i]
                            if i < len(overlaps):
                                overlaps[i] = normal_overlap(level[i], error[i], level[i + 1], error[i + 1])
                            if i > 0:
                                overlaps[i - 1] = normal_overlap(level[i - 1], error[i - 1], level[i], error[i])
                            iter = iter + 1
                            if anim != "" and (iter % 100) == 0:
                                anim_plot_rd(level, error, segments, bins, iter, anim + c + "_0_" + str(bin_size), maxo,
                                             mino, mean)

                        iter = 0
                        ons = -1

                        _logger.info("Second stage. Number of segments: %d." % len(level))

                        while True:
                            overlaps = [normal_overlap(level[i], error[i], level[j], error[j]) for i in
                                        range(len(level)) for j in range(i + 1, len(level)) if
                                        (segments[j][0] - segments[i][-1]) < max_distance * (
                                                len(segments[i]) + len(segments[j]))]
                            if len(overlaps) == 0:
                                break

                            maxo = max(overlaps)
                            if maxo < overlap_min:
                                break
                            i, j = 0, 1
                            while i < len(segments) - 1:

                                if (segments[j][0] - segments[i][-1]) < max_distance * (
                                        len(segments[i]) + len(segments[j])) and normal_overlap(level[i], error[i],
                                                                                                level[j],
                                                                                                error[j]) == maxo:
                                    nl, ne = normal_merge(level[i], error[i], level[j], error[j])
                                    level[i] = nl
                                    error[i] = ne
                                    segments[i] += segments[j]
                                    segments[i] = sorted(segments[i])
                                    del level[j]
                                    del error[j]
                                    del segments[j]

                                    if j >= len(segments):
                                        i += 1
                                        j = i + 1
                                else:
                                    j += 1
                                    if j >= len(segments):
                                        i += 1
                                        j = i + 1
                            iter = iter + 1
                            if anim != "" and (iter % 100) == 0:
                                anim_plot_rd(level, error, segments, bins, iter, anim + c + "_1_" + str(bin_size), maxo,
                                             mino, mean)
                            _logger.info("Iteration: %d. Number of segments: %d." % (iter, len(level)))
                            if ons == len(segments):
                                break
                            ons = len(segments)

                        for i in range(len(segments)):
                            i1 = level[i] / mean
                            i2 = t_test_1_sample(mean, level[i], error[i], len(segments[i]))
                            if i2 is not None and i2 < (0.05 * bin_size / 3e9) and abs(i1 - 1.) > 0.01:
                                print(c + ":" + str(segments[i][0] * bin_size + 1) + "-" + str(
                                    segments[i][-1] * bin_size + bin_size),
                                      (segments[i][-1] - segments[i][0] + 1) * bin_size, bin_size, len(segments[i]),
                                      i1, i2)

                        self.io.create_signal(c, bin_size, "RD mosaic segments",
                                              data=segments_code(segments), flags=flag_rd)
                        self.io.create_signal(c, bin_size, "RD mosaic call",
                                              data=np.array([level, error], dtype="float32"), flags=flag_rd)

    def mask_snps(self, callset=None):
        """ Flags SNPs in P-region (sets second bit of the flag to 1 for SNP inside P region, or to 0 otherwise).
        Requires imported mask data or recognized reference genome with mask data.

        Parameters
        ----------
        callset : str or None
            It will assume SNP data if None. Otherwise it will assume SNV data
            stored under the name provided by callset variable.

        """
        snp_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            snp_name = self.io.snp_chromosome_name(c)
            if snp_name is not None:
                snp_mask_chromosomes[snp_name] = c

        for c in self.io.snp_chromosomes():
            if c in snp_mask_chromosomes:
                _logger.info("Masking SNP data for chromosome '%s'." % c)
                pos, ref, alt, nref, nalt, gt, flag, qual = self.io.read_snp(c, callset=callset)
                mask = mask_decompress(self.io_mask.get_signal(snp_mask_chromosomes[c], None, "mask"))
                mask_ix = 0
                for snp_ix in range(len(pos)):
                    while mask_ix < len(mask) and mask[mask_ix][1] < pos[snp_ix]:
                        mask_ix += 1
                    if mask_ix < len(mask) and mask[mask_ix][0] < pos[snp_ix] and pos[snp_ix] < (
                            mask[mask_ix][1] + 1):
                        flag[snp_ix] = flag[snp_ix] | 2
                    else:
                        flag[snp_ix] = flag[snp_ix] & 1
                if len(pos) > 0:
                    self.io.save_snp(c, pos, ref, alt, nref, nalt, gt, flag, qual, update=True, callset=callset)

    def calculate_baf(self, bin_sizes, chroms=[], use_mask=True, use_id=False, use_phase=False, res=200,
                      reduce_noise=False, blw=0.8):
        """
        Calculates BAF histograms and store data into cnvpytor file.

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        use_mask : bool
            Use P-mask filter if True. Default: True.
        use_id : bool
            Use id flag filter if True. Default: False.
        use_phase : bool
            Use phasing information if True and available. Default: False.
        res: int
            Likelihood function resolution. Default: 200.


        """
        snp_flag = (FLAG_USEMASK if use_mask else 0) | (FLAG_USEID if use_id else 0)
        for c in self.io.snp_chromosomes():
            if len(chroms) == 0 or c in chroms:
                _logger.info("Calculating BAF histograms for chromosome '%s'." % c)
                pos, ref, alt, nref, nalt, gt, flag, qual = self.io.read_snp(c)
                max_bin = {}
                count00 = {}
                count01 = {}
                count10 = {}
                count11 = {}
                reads00 = {}
                reads01 = {}
                reads10 = {}
                reads11 = {}
                baf = {}
                maf = {}
                likelihood = {}
                i1 = {}
                i2 = {}
                lh_x = np.arange(1.0 / res, 1., 1.0 / res)
                for bs in bin_sizes:
                    max_bin[bs] = (pos[-1] - 1) // bs + 1
                    count00[bs] = np.zeros(max_bin[bs])
                    count01[bs] = np.zeros(max_bin[bs])
                    count10[bs] = np.zeros(max_bin[bs])
                    count11[bs] = np.zeros(max_bin[bs])
                    reads00[bs] = np.zeros(max_bin[bs])
                    reads01[bs] = np.zeros(max_bin[bs])
                    reads10[bs] = np.zeros(max_bin[bs])
                    reads11[bs] = np.zeros(max_bin[bs])
                    baf[bs] = np.zeros(max_bin[bs])
                    maf[bs] = np.zeros(max_bin[bs])
                    likelihood[bs] = np.ones((max_bin[bs], res - 1)).astype("float") / (res - 1)
                    i1[bs] = np.zeros(max_bin[bs])
                    i2[bs] = np.zeros(max_bin[bs])

                for i in range(len(pos)):
                    if (nalt[i] + nref[i]) > 0 and (not use_id or (flag[i] & 1)) and (not use_mask or (flag[i] & 2)):
                        if gt[i] == 1 or gt[i] == 5 or gt[i] == 6:
                            for bs in bin_sizes:
                                b = (pos[i] - 1) // bs
                                if use_phase and (gt[i] == 5):
                                    reads01[bs][b] += nref[i]
                                    reads10[bs][b] += nalt[i]
                                    count10[bs][b] += 1
                                    snp_baf = 1.0 * nalt[i] / (nalt[i] + nref[i])
                                    likelihood[bs][b] *= beta(nalt[i], nref[i], lh_x, phased=True)
                                    s = np.sum(likelihood[bs][b])
                                    if s != 0.0:
                                        likelihood[bs][b] /= s
                                elif use_phase and (gt[i] == 6):
                                    reads01[bs][b] += nalt[i]
                                    reads10[bs][b] += nref[i]
                                    count01[bs][b] += 1
                                    snp_baf = 1.0 * nref[i] / (nalt[i] + nref[i])
                                    likelihood[bs][b] *= beta(nref[i], nalt[i], lh_x, phased=True)
                                    s = np.sum(likelihood[bs][b])
                                    if s != 0.0:
                                        likelihood[bs][b] /= s
                                else:
                                    snp_baf = 1.0 * nalt[i] / (nalt[i] + nref[i])
                                    reads01[bs][b] += nalt[i]
                                    reads10[bs][b] += nref[i]
                                    count01[bs][b] += 1
                                    if reduce_noise:
                                        likelihood[bs][b] *= beta(nalt[i] + (1 if nalt[i] < nref[i] else 0),
                                                                  nref[i] + (1 if nref[i] < nalt[i] else 0), lh_x)
                                    else:
                                        likelihood[bs][b] *= beta(nalt[i] * blw, nref[i] * blw, lh_x)
                                    s = np.sum(likelihood[bs][b])
                                    if s != 0.0:
                                        likelihood[bs][b] /= s

                                baf[bs][b] += snp_baf
                                maf[bs][b] += 1.0 - snp_baf if snp_baf > 0.5 else snp_baf
                        else:
                            for bs in bin_sizes:
                                b = (pos[i] - 1) // bs
                                if use_phase and (gt[i] == 7):
                                    count11[bs][b] += 1
                                    reads11[bs][b] += nalt[i]
                                    reads00[bs][b] += nref[i]
                                elif use_phase and (gt[i] == 4):
                                    count00[bs][b] += 1
                                else:
                                    count11[bs][b] += 1
                                    reads11[bs][b] += nalt[i]
                                    reads00[bs][b] += nref[i]

                for bs in bin_sizes:
                    for i in range(max_bin[bs]):
                        count = count01[bs][i] + count10[bs][i]
                        if count > 0:
                            baf[bs][i] /= count
                            maf[bs][i] /= count
                            max_lh = np.amax(likelihood[bs][i])
                            ix = np.where(likelihood[bs][i] == max_lh)[0][0]
                            i1[bs][i] = 1.0 * (res // 2 - 1 - ix) / res if ix <= (res // 2 - 1) else 1.0 * (
                                    ix - res // 2 + 1) / res
                            i2[bs][i] = likelihood[bs][i][res // 2 - 1] / max_lh
                    _logger.info("Saving BAF histograms with bin size %d for chromosome '%s'." % (bs, c))
                    self.io.create_signal(c, bs, "SNP bin count 0|0", count00[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin count 0|1", count01[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin count 1|0", count10[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin count 1|1", count11[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin reads 0|0", reads00[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin reads 0|1", reads01[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin reads 1|0", reads10[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP bin reads 1|1", reads11[bs].astype("uint16"), snp_flag)
                    self.io.create_signal(c, bs, "SNP baf", baf[bs].astype("float32"), snp_flag)
                    self.io.create_signal(c, bs, "SNP maf", maf[bs].astype("float32"), snp_flag)
                    self.io.create_signal(c, bs, "SNP likelihood", likelihood[bs].astype("float32"), snp_flag)
                    self.io.create_signal(c, bs, "SNP i1", i1[bs].astype("float32"), snp_flag)
                    self.io.create_signal(c, bs, "SNP i2", i2[bs].astype("float32"), snp_flag)

    def call_baf(self, bin_sizes, chroms=[], use_mask=True, use_id=False, odec=0.9, omin=None, mcount=None,
                 max_distance=0.1, anim=""):
        """ CNV caller based on BAF likelihood mearger

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        use_mask : bool
            Use P-mask filter if True. Default: False.
        use_id : bool
            Use ID filter if True. Default: False.
        """

        for bin_size in bin_sizes:
            if omin is None:
                overlap_min = 1e-7 * bin_size
            else:
                overlap_min = omin
            if mcount is None:
                min_count = bin_size // 10000
            else:
                min_count = mcount

            snp_flag = (FLAG_USEMASK if use_mask else 0) | (FLAG_USEID if use_id else 0)
            for c in self.io.snp_chromosomes():
                if len(chroms) == 0 or c in chroms and self.io.signal_exists(c, bin_size, "SNP likelihood", snp_flag):
                    # _logger.info(
                    #    "Caling CNV-s by merging BAF likelihood with bin size %d for chromosome '%s'." % (bin_size, c))
                    likelihood = list(self.io.get_signal(c, bin_size, "SNP likelihood", snp_flag).astype("float64"))
                    snp_hets = self.io.get_signal(c, bin_size, "SNP bin count 0|1", snp_flag)
                    snp_hets += self.io.get_signal(c, bin_size, "SNP bin count 1|0", snp_flag)

                    bins = len(likelihood)
                    res = likelihood[0].size
                    segments = [[i] for i in range(bins) if
                                snp_hets[i] >= min_count and np.sum(likelihood[i]) > 0.0]
                    likelihood = [likelihood[i] for i in range(bins) if
                                  snp_hets[i] >= min_count and np.sum(likelihood[i]) > 0.0]
                    overlaps = [likelihood_overlap(likelihood[i], likelihood[i + 1]) for i in range(len(segments) - 1)]

                    iter = 0
                    while len(overlaps) > 0:
                        maxo = max(overlaps)
                        mino = max(maxo * odec, overlap_min)
                        if maxo < overlap_min:
                            break
                        i = 0
                        while i < len(overlaps):
                            if overlaps[i] > mino:
                                nlh = likelihood[i] * likelihood[i + 1]
                                likelihood[i] = nlh / np.sum(nlh)
                                segments[i] += segments[i + 1]
                                del likelihood[i + 1]
                                del segments[i + 1]
                                del overlaps[i]
                                if i < len(overlaps):
                                    overlaps[i] = likelihood_overlap(likelihood[i], likelihood[i + 1])
                                if i > 0:
                                    overlaps[i - 1] = likelihood_overlap(likelihood[i - 1], likelihood[i])
                            else:
                                i = i + 1

                        iter = iter + 1
                        if anim != "":
                            anim_plot_likelihood(likelihood, segments, bins, res, iter,
                                                 anim + c + "_0_" + str(bin_size), maxo, mino)

                    overlaps = [[
                        likelihood_overlap(likelihood[i], likelihood[j])
                        if (segments[j][0] - segments[i][-1]) < max_distance * (
                                len(segments[i]) + len(segments[j])) and i < j
                        else 0 for j in range(len(segments))]
                        for i in range(len(segments))]

                    iter = 0
                    ons = -1

                    while True:
                        overlaps = [likelihood_overlap(likelihood[i], likelihood[j]) for i in range(len(likelihood))
                                    for j in range(i + 1, len(likelihood)) if
                                    (segments[j][0] - segments[i][-1]) < max_distance * (
                                            len(segments[i]) + len(segments[j]))]
                        if len(overlaps) == 0:
                            break
                        maxo = max(overlaps)
                        mino = max(maxo * odec, overlap_min)
                        if maxo < overlap_min:
                            break
                        i, j = 0, 1
                        while i < len(segments) - 1:
                            #                            if overlaps[i][j] > mino:
                            if likelihood_overlap(likelihood[i], likelihood[j]) > mino and (
                                    segments[j][0] - segments[i][-1]) < max_distance * (
                                    len(segments[i]) + len(segments[j])):
                                nlh = likelihood[i] * likelihood[j]
                                likelihood[i] = nlh / np.sum(nlh)
                                segments[i] += segments[j]
                                segments[i] = sorted(segments[i])
                                del likelihood[j]
                                del segments[j]

                                if j >= len(segments):
                                    i += 1
                                    j = i + 1
                            else:
                                j += 1
                                if j >= len(segments):
                                    i += 1
                                    j = i + 1
                        iter = iter + 1
                        if anim != "":
                            anim_plot_likelihood(likelihood, segments, bins, res, iter,
                                                 anim + c + "_1_" + str(bin_size), maxo, mino)

                        if ons == len(segments):
                            break
                        ons = len(segments)

                    for i in range(len(segments)):
                        i1, i2 = likelihood_baf_pval(likelihood[i])

                        print(c + ":" + str(segments[i][0] * bin_size + 1) + "-" + str(
                            segments[i][-1] * bin_size + bin_size),
                              (segments[i][-1] - segments[i][0] + 1) * bin_size, bin_size, len(segments[i]),
                              i1, i2)

                    self.io.create_signal(c, bin_size, "SNP likelihood segments",
                                          data=segments_code(segments), flags=snp_flag)
                    self.io.create_signal(c, bin_size, "SNP likelihood call",
                                          data=np.array(likelihood, dtype="float32"), flags=snp_flag)

    def call_2d(self, bin_sizes, chroms=[], use_gc_corr=True, rd_use_mask=False, snp_use_mask=True, snp_use_id=False,
                odec=0.9, omin=None, mcount=None, max_distance=0.1, anim=""):
        """
        CNV caller using combined RD and BAF sigal based on likelihood merger.

        Parameters
        ----------
        bin_sizes : list of int
            List of histogram bin sizes
        chroms : list of str
            List of chromosomes. Calculates for all available if empty.
        use_gc_corr : bool
            Use GC corrected signal if True. Default: True.
        rd_use_mask : bool
            Use P-mask filter for RD if True. Default: False.
        snp_use_mask : bool
            Use P-mask filter for SNP if True. Default: True.
        snp_use_id : bool
            Use ID filter for SNP if True. Default: False.

        """

        snp_flag = (FLAG_USEMASK if snp_use_mask else 0) | (FLAG_USEID if snp_use_id else 0)
        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_gc_chromosomes[rd_name] = c

        rd_mask_chromosomes = {}
        for c in self.io_mask.mask_chromosomes():
            rd_name = self.io.rd_chromosome_name(c)
            if not rd_name is None:
                rd_mask_chromosomes[rd_name] = c

        for bin_size in bin_sizes:
            if omin is None:
                overlap_min = 0.05 * bin_size / 3e9
            else:
                overlap_min = omin
            if mcount is None:
                min_count = bin_size // 10000
            else:
                min_count = mcount

            for c in self.io.rd_chromosomes():
                if (c in rd_gc_chromosomes or not use_gc_corr) and (c in rd_mask_chromosomes or not rd_use_mask) and (
                        self.io.signal_exists(c, bin_size, "SNP likelihood", snp_flag)) and (
                        len(chroms) == 0 or (c in chroms)):
                    flag_stat = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    flag_auto = FLAG_AUTO
                    if use_gc_corr:
                        flag_stat |= FLAG_GC_CORR
                        flag_auto |= FLAG_GC_CORR
                    flag_rd = (FLAG_GC_CORR if use_gc_corr else 0) | (FLAG_USEMASK if rd_use_mask else 0)
                    if self.io.signal_exists(c, bin_size, "RD stat", flag_stat) and \
                            self.io.signal_exists(c, bin_size, "RD", flag_rd):
                        _logger.info("Calculating 2d calls using bin size %d for chromosome '%s'." % (bin_size, c))
                        stat = self.io.get_signal(c, bin_size, "RD stat", flag_stat)
                        mean = stat[4]
                        std = stat[5]
                        rd = self.io.get_signal(c, bin_size, "RD", flag_rd)
                        rd_bins = len(rd)

                        snp_likelihood = list(
                            self.io.get_signal(c, bin_size, "SNP likelihood", snp_flag).astype("float64"))
                        snp_hets = self.io.get_signal(c, bin_size, "SNP bin count 0|1", snp_flag)
                        snp_hets += self.io.get_signal(c, bin_size, "SNP bin count 1|0", snp_flag)

                        snp_bins = len(snp_likelihood)
                        res = snp_likelihood[0].size
                        bins = min(rd_bins, snp_bins)

                        segments = [[i] for i in range(bins) if
                                    snp_hets[i] >= min_count and np.sum(snp_likelihood[i]) > 0.0 and np.isfinite(rd[i])]

                        level = [rd[i] for i in range(bins) if
                                 snp_hets[i] >= min_count and np.sum(snp_likelihood[i]) > 0.0 and np.isfinite(rd[i])]
                        level = np.array(level)
                        ilevel = level.copy()
                        error = np.sqrt(level) ** 2 + std ** 2
                        loc_fl = np.min(list(zip(np.abs(np.diff(level))[:-1], np.abs(np.diff(level))[1:])), axis=1)
                        loc_fl = np.concatenate(([0], loc_fl, [0]))
                        error += (loc_fl / 2) ** 2
                        error = np.sqrt(error)
                        level = list(level)
                        error = list(error)

                        likelihood = [snp_likelihood[i] for i in range(bins) if
                                      snp_hets[i] >= min_count and np.sum(snp_likelihood[i]) > 0.0 and np.isfinite(
                                          rd[i])]

                        overlaps = [normal_overlap(level[i], error[i], level[i + 1], error[i + 1]) * likelihood_overlap(
                            likelihood[i], likelihood[i + 1]) for i in range(len(segments) - 1)]

                        iter = 0
                        if anim != "":
                            anim_plot_rd(level, error, segments, bins, iter, anim + c + "_0_" + str(bin_size), 0,
                                         0, mean)

                        while len(overlaps) > 0:
                            maxo = max(overlaps)
                            if maxo < overlap_min:
                                break
                            i = overlaps.index(maxo)
                            nl, ne = normal_merge(level[i], error[i], level[i + 1], error[i + 1])
                            nlh = likelihood[i] * likelihood[i + 1]
                            level[i] = nl
                            error[i] = ne
                            likelihood[i] = nlh / np.sum(nlh)
                            segments[i] += segments[i + 1]
                            del level[i + 1]
                            del error[i + 1]
                            del segments[i + 1]
                            del likelihood[i + 1]
                            del overlaps[i]
                            if i < len(overlaps):
                                overlaps[i] = normal_overlap(level[i], error[i], level[i + 1],
                                                             error[i + 1]) * likelihood_overlap(likelihood[i],
                                                                                                likelihood[i + 1])
                            if i > 0:
                                overlaps[i - 1] = normal_overlap(level[i - 1], error[i - 1], level[i],
                                                                 error[i]) * likelihood_overlap(likelihood[i - 1],
                                                                                                likelihood[i])
                            iter = iter + 1
                            if anim != "" and (iter % 100) == 0:
                                anim_plot_rd(level, error, segments, bins, iter, anim + c + "_0_" + str(bin_size), maxo,
                                             mino, mean)

                        iter = 0
                        ons = -1

                        _logger.info("Second stage. Number of segments: %d." % len(level))

                        while True:
                            overlaps = [normal_overlap(level[i], error[i], level[j], error[j]) * likelihood_overlap(
                                likelihood[i], likelihood[j]) for i in range(len(level)) for j in
                                        range(i + 1, len(level)) if
                                        (segments[j][0] - segments[i][-1]) < max_distance * (
                                                len(segments[i]) + len(segments[j]))]
                            if len(overlaps) == 0:
                                break

                            maxo = max(overlaps)
                            if maxo < overlap_min:
                                break
                            i, j = 0, 1
                            while i < len(segments) - 1:

                                if (segments[j][0] - segments[i][-1]) < max_distance * (
                                        len(segments[i]) + len(segments[j])) and \
                                        normal_overlap(level[i], error[i], level[j], error[j]) * likelihood_overlap(
                                    likelihood[i], likelihood[j]) == maxo:
                                    nl, ne = normal_merge(level[i], error[i], level[j], error[j])
                                    nlh = likelihood[i] * likelihood[j]
                                    level[i] = nl
                                    error[i] = ne
                                    likelihood[i] = nlh / np.sum(nlh)
                                    segments[i] += segments[j]
                                    segments[i] = sorted(segments[i])
                                    del level[j]
                                    del error[j]
                                    del segments[j]

                                    if j >= len(segments):
                                        i += 1
                                        j = i + 1
                                else:
                                    j += 1
                                    if j >= len(segments):
                                        i += 1
                                        j = i + 1
                            iter = iter + 1
                            if anim != "" and (iter % 100) == 0:
                                anim_plot_rd(level, error, segments, bins, iter, anim + c + "_1_" + str(bin_size), maxo,
                                             mino, mean)
                            _logger.info("Iteration: %d. Number of segments: %d." % (iter, len(level)))
                            if ons == len(segments):
                                break
                            ons = len(segments)

                        # rd = [np.nan] * bins
                        # for i in range(len(segments)):
                        #     for b in segments[i]:
                        #         rd[b] = level[i]
                        # plt.plot(rd)
                        # plt.show()
                        # return

                        for i in range(len(segments)):
                            rd_mean = level[i] / mean
                            rd_p = t_test_1_sample(mean, level[i], error[i], len(segments[i]))
                            baf_mean, baf_p = likelihood_baf_pval(likelihood[i])
                            print(c + ":" + str(segments[i][0] * bin_size + 1) + "-" + str(
                                segments[i][-1] * bin_size + bin_size),
                                  (segments[i][-1] - segments[i][0] + 1) * bin_size, bin_size, len(segments[i]),
                                  rd_mean, rd_p, baf_mean, baf_p)

                        self.io.create_signal(c, bin_size, "RD mosaic segments",
                                              data=segments_code(segments), flags=flag_rd)
                        self.io.create_signal(c, bin_size, "RD mosaic call",
                                              data=np.array([level, error], dtype="float32"), flags=flag_rd)
                        self.io.create_signal(c, bin_size, "SNP likelihood segments",
                                              data=segments_code(segments), flags=snp_flag)
                        self.io.create_signal(c, bin_size, "SNP likelihood call",
                                              data=np.array(likelihood, dtype="float32"), flags=snp_flag)

        return

    def ls(self):
        """ Print to stdout content of cnvpytor file.

        """
        self.io.ls()
