""" cnvpytor.vcf

class Vcf: reads SNP data from VCF file

"""
import pysam
import logging
import numpy as np
from .utils import int1, gt_from_list, gt_from_str
_logger = logging.getLogger("cnvpytor.vcf")


class Vcf:

    def __init__(self, filename):
        """
        Opens VCF file, reads chromosome/sample names

        Parameters
        ----------
        filename : str
            Name of the VCF file.

        """
        self.filename = filename
        try:
            self.file = pysam.VariantFile(filename)
        except IOError:
            _logger.error("Problem opening file '%s'!" % filename)
            exit(0)
        except ValueError:
            _logger.error("Problem opening file '%s'!" % filename)
            exit(0)

        self.chrs = []
        self.samples = []
        self.lengths = {}
        if self.file:
            _logger.info("File: " + filename + " successfully open")
            self.chrs = list(self.file.header.contigs)
            self.samples = list(self.file.header.samples)
            for c in self.chrs:
                self.lengths[c] = self.file.header.contigs[c].length
            _logger.debug("Header contigs: %s" % ", ".join(self.chrs))
            _logger.debug("Header samples: %s" % ", ".join(self.samples))

    def get_chromosomes(self):
        """
        Get chromosome names.

        Returns
        -------
        chrs : list of str
            List of chromosome names from VCF header

        """
        return self.chrs

    def get_samples(self):
        """
        Get sample names.

        Returns
        -------
        samples : list of str
            List of sample names from VCF header

        """
        return self.samples

    def read_chromosome_snp(self, chr_name, sample='', ad_tag='AD', gt_tag='GT', filter=True):
        """
        Read SNP/indel data for given chromosome and sample.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        sample : str
            Name of the sample (first one if empty - default)
        ad_tag : str
            AD tag in VCF file (default: 'AD')
        gt_tag : str
            GT tag in VCF file (default: 'GT')
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

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
        if sample == '':
            sample = self.samples[0]
        pos = []
        ref = []
        alt = []
        nref = []
        nalt = []
        gt = []
        flag = []
        qual = []
        filter_stat = {}
        alphabet = ['A', 'T', 'G', 'C', '.']
        try:
            for rec in self.file.fetch(chr_name):
                if len(rec.filter.keys()) == 0:
                    if "No filter (.)" in filter_stat:
                        filter_stat["No filter (.)"] += 1
                    else:
                        filter_stat["No filter (.)"] = 1
                for f in rec.filter.keys():
                    if f in filter_stat:
                        filter_stat[f] += 1
                    else:
                        filter_stat[f] = 1
                if ("PASS" in rec.filter.keys() or not filter) and rec.alts and len(rec.alts) >= 1 and (
                        gt_tag in rec.samples[sample].keys()) and (
                        ad_tag in rec.samples[sample].keys()) and len(rec.samples[sample][gt_tag]) > 1 and len(
                    rec.samples[sample][ad_tag]) > 1:
                    if (len(rec.samples[sample][ad_tag]) > 1) and (rec.ref in alphabet) and (rec.alts[0] in alphabet):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        if rec.qual == "." or rec.qual is None:
                            qual.append(0)
                        else:
                            qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        try:
                            nref.append(int(rec.samples[sample][ad_tag][0]))
                            nalt.append(int(rec.samples[sample][ad_tag][1]))
                        except ValueError:
                            nref.append(0)
                            nalt.append(0)
                        try:
                            UNICODE_EXISTS = bool(type(unicode))
                        except NameError:
                            unicode = str
                        if isinstance(rec.samples[sample][gt_tag], str) or isinstance(rec.samples[sample][gt_tag],
                                                                                      unicode):
                            gt.append(gt_from_str(rec.samples[sample][gt_tag]))
                        else:
                            gt.append(gt_from_list(rec.samples[sample][gt_tag], rec.samples[sample].phased))
        except ValueError:
            _logger.error("Variant file reading problem. Probably index file is missing or corrupted.")
            exit(0)
        _logger.debug("Chromosome '%s' read. Number of variants to store: %d." % (chr_name, len(pos)))
        _logger.debug("Variants filter field statistics:")
        for f in filter_stat:
            _logger.debug(" * '%s' : %d" % (f, filter_stat[f]))
        return pos, ref, alt, nref, nalt, gt, flag, qual

    def read_chromosome_snp_no_counts(self, chr_name, sample='', gt_tag='GT', filter=True):
        """
        Read SNP/indel data without counts (AD tag) for given chromosome and sample.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        sample : str
            Name of the sample (first one if empty - default)
        gt_tag : str
            GT tag in VCF file (default: 'GT')
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

        Returns
        -------
        pos : list of int
            List of SNP positions.
        ref : list of str
            List of SNP reference base (A, T, G, C or .).
        alt : list of str
            List of SNP alternative base (A, T, G, C or .).
        gt : list of int
            List of genotypes (0 - "0/0", 1 - "0/1", 3- "1/1", 4 - "0|0" , 5 - "0|1", 6 - "1|0", 7 - "1|1").
        flag : list of int
            Binary flag: first bit 1 - SNP exists in database, second bit 1 - SNP in P region of strict mask.
        qual : list of int
            SNP quality (scale 0 - 255).

        """
        if sample == '':
            sample = self.samples[0]
        pos = []
        ref = []
        alt = []
        gt = []
        flag = []
        qual = []
        alphabet = ['A', 'T', 'G', 'C', '.']
        try:
            for rec in self.file.fetch(chr_name):
                if ("PASS" in rec.filter.keys() or not filter) and rec.alts and len(rec.alts) >= 1 and (
                        gt_tag in rec.samples[sample].keys()) and len(rec.samples[sample][gt_tag]) > 1:
                    if (rec.ref in alphabet) and (rec.alts[0] in alphabet):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        if rec.qual == "." or rec.qual is None:
                            qual.append(0)
                        else:
                            qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        try:
                            UNICODE_EXISTS = bool(type(unicode))
                        except NameError:
                            unicode = str
                        if isinstance(rec.samples[sample][gt_tag], str) or isinstance(rec.samples[sample][gt_tag],
                                                                                      unicode):
                            gt.append(gt_from_str(rec.samples[sample][gt_tag]))
                        else:
                            gt.append(gt_from_list(rec.samples[sample][gt_tag], rec.samples[sample].phased))

        except ValueError:
            _logger.error("Variant file reading problem. Probably index file is missing or corrupted.")
            exit(0)
        return pos, ref, alt, gt, flag, qual

    def read_chromosome_rd(self, chr_name, sample='', ad_tag='AD', dp_tag='DP', filter=False):
        """
        Read RD data for given chromosome and sample.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        sample : str
            Name of the sample (first one if empty - default)
        ad_tag : str
            AD tag in VCF file (default: 'AD')
        dp_tag : str
            DP tag in VCF file (default: 'DP')
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

        Returns
        -------
        rd_p : list of int
            binned RD signal from DP tag (100 bins)
        rd_u : list of int
            binned RD signal from AD tag (100 bins)

        """
        if sample == '':
            sample = self.samples[0]
        pos = []
        rd_ad = []
        rd_dp = []
        try:
            for rec in self.file.fetch(chr_name):
                if ("PASS" in rec.filter.keys() or not filter) and (ad_tag in rec.samples[sample].keys()) and len(
                        rec.samples[sample][ad_tag]) > 1 and (dp_tag in rec.samples[sample].keys()) and (
                        rec.samples[sample][dp_tag] is not None):
                    try:
                        rd1 = sum(map(int, rec.samples[sample][ad_tag]))
                        rd2 = int(rec.samples[sample][dp_tag])
                        pos.append(rec.pos)
                        rd_ad.append(rd1)
                        rd_dp.append(rd2)
                    except ValueError:
                        pass
        except ValueError:
            _logger.error("Variant file reading problem. Probably index file is missing or corrupted.")
            exit(0)
        n = 0
        if len(pos) == 0:
            return None, None
        n = pos[-1] + 1
        if self.lengths[chr_name] is not None:
            n = self.lengths[chr_name]
        n = n // 100 + 1
        rd_p = np.zeros(n)
        rd_u = np.zeros(n)
        count = np.zeros(n)
        for p, rd1, rd2 in zip(pos, rd_dp, rd_ad):
            cgb = (p - 1) // 10000 * 100
            rd_p[cgb] += rd1
            rd_u[cgb] += rd2
            count[cgb] += 1
        for i in range(n // 100 + 1):
            if (100 * i) < n:
                if count[100 * i] != 0:
                    rd_p[i * 100:(i + 1) * 100] = rd_p[i * 100] / count[i * 100]
                    rd_u[i * 100:(i + 1) * 100] = rd_u[i * 100] / count[i * 100]
                else:
                    rd_p[i * 100:(i + 1) * 100] = np.nan
                    rd_u[i * 100:(i + 1) * 100] = np.nan
        return rd_p, rd_u

    def read_all_snp(self, callback, sample='', ad_tag='AD', gt_tag='GT', filter=True):
        """
        Read SNP/indel data for given sample name.

        Parameters
        ----------
        callback : callable
            Function to call after read a chromosome:
            callback(chrom, pos, ref, alt, nref, nalt, gt, flag, qual)
        sample : str
            Name of the sample (first one if empty - default)
        ad_tag : str
            AD tag in VCF file (default: 'AD')
        gt_tag : str
            GT tag in VCF file (default: 'GT')
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

        Returns
        -------
        count : int
            Number of read chromosomes

        """
        if sample == '':
            sample = self.samples[0]
        pos = []
        ref = []
        alt = []
        nref = []
        nalt = []
        gt = []
        flag = []
        qual = []
        filter_stat = {}
        last_chrom = None
        count = 0
        alphabet = ['A', 'T', 'G', 'C', '.']

        try:
            for rec in self.file.fetch():
                if last_chrom is None:
                    last_chrom = rec.chrom
                if last_chrom != rec.chrom:
                    _logger.info("Chromosome '%s' read. Number of variants to store: %d." % (last_chrom, len(pos)))
                    _logger.debug("Variants filter field statistics:")
                    for f in filter_stat:
                        _logger.debug("    * '%s' : %d" % (f, filter_stat[f]))
                    callback(last_chrom, pos, ref, alt, nref, nalt, gt, flag, qual)
                    pos = []
                    ref = []
                    alt = []
                    nref = []
                    nalt = []
                    gt = []
                    flag = []
                    qual = []
                    filter_stat = {}
                    count += 1
                if len(rec.filter.keys()) == 0:
                    if "No filter (.)" in filter_stat:
                        filter_stat["No filter (.)"] += 1
                    else:
                        filter_stat["No filter (.)"] = 1
                for f in rec.filter.keys():
                    if f in filter_stat:
                        filter_stat[f] += 1
                    else:
                        filter_stat[f] = 1

                if ("PASS" in rec.filter.keys() or not filter) and rec.alts and len(rec.alts) >= 1 and (
                        gt_tag in rec.samples[sample].keys()) and (
                        ad_tag in rec.samples[sample].keys()) and len(rec.samples[sample][gt_tag]) > 1 and len(
                    rec.samples[sample][ad_tag]) > 1:
                    if (len(rec.samples[sample][ad_tag]) > 1) and (rec.ref in alphabet) and (
                            rec.alts[0] in alphabet) and (rec.samples[sample][gt_tag][0] is not None) and (
                            rec.samples[sample][gt_tag][1] is not None):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        if rec.qual == "." or rec.qual is None:
                            qual.append(0)
                        else:
                            qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        try:
                            nref.append(int(rec.samples[sample][ad_tag][0]))
                            nalt.append(int(rec.samples[sample][ad_tag][1]))
                        except ValueError:
                            nref.append(0)
                            nalt.append(0)
                        try:
                            UNICODE_EXISTS = bool(type(unicode))
                        except NameError:
                            unicode = str
                        if isinstance(rec.samples[sample][gt_tag], str) or isinstance(rec.samples[sample][gt_tag],
                                                                                      unicode):
                            gt.append(gt_from_str(rec.samples[sample][gt_tag]))
                        else:
                            gt.append(gt_from_list(rec.samples[sample][gt_tag], rec.samples[sample].phased))

                last_chrom = rec.chrom
            _logger.debug("Chromosome '%s' read. Number of variants to store: %d." % (last_chrom, len(pos)))
            _logger.debug("Variants filter field statistics:")
            for f in filter_stat:
                _logger.debug(" * '%s' : %d" % (f, filter_stat[f]))
            callback(last_chrom, pos, ref, alt, nref, nalt, gt, flag, qual)
            count += 1
            return count
        except ValueError:
            _logger.error("Variant file reading problem.")
            exit(0)

    def read_all_snp_no_counts(self, callback, sample='', gt_tag='GT', filter=True):
        """
        Read SNP/indel data without counts (AD tag) for given sample name.

        Parameters
        ----------
        callback : callable
            Function to call after read a chromosome:
            callback(chrom, pos, ref, alt, gt, flag, qual)
        sample : str
            Name of the sample (first one if empty - default)
        gt_tag : str
            GT tag in VCF file (default: 'GT')
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

        Returns
        -------
        count : int
            Number of read chromosomes

        """
        if sample == '':
            sample = self.samples[0]
        pos = []
        ref = []
        alt = []
        gt = []
        flag = []
        qual = []
        last_chrom = None
        count = 0
        alphabet = ['A', 'T', 'G', 'C', '.']

        try:
            for rec in self.file.fetch():
                if last_chrom is None:
                    last_chrom = rec.chrom
                if last_chrom != rec.chrom:
                    callback(last_chrom, pos, ref, alt, gt, flag, qual)
                    pos = []
                    ref = []
                    alt = []
                    gt = []
                    flag = []
                    qual = []
                    count += 1

                if ("PASS" in rec.filter.keys() or not filter) and rec.alts and len(rec.alts) >= 1 and (
                        gt_tag in rec.samples[sample].keys()) and len(rec.samples[sample][gt_tag]) > 1:
                    if (rec.ref in alphabet) and (rec.alts[0] in alphabet) and (
                            rec.samples[sample][gt_tag][0] is not None) and (
                            rec.samples[sample][gt_tag][1] is not None):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        if rec.qual == "." or rec.qual is None:
                            qual.append(0)
                        else:
                            qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        try:
                            UNICODE_EXISTS = bool(type(unicode))
                        except NameError:
                            unicode = str
                        if isinstance(rec.samples[sample][gt_tag], str) or isinstance(rec.samples[sample][gt_tag],
                                                                                      unicode):
                            gt.append(gt_from_str(rec.samples[sample][gt_tag]))
                        else:
                            gt.append(gt_from_list(rec.samples[sample][gt_tag], rec.samples[sample].phased))

                last_chrom = rec.chrom
            if len(pos) > 0:
                callback(last_chrom, pos, ref, alt, gt, flag, qual)
                count += 1
            return count
        except ValueError:
            _logger.error("Variant file reading problem.")
            exit(0)

    def read_all_snp_positions(self, callback, filter=True):
        """
        Read SNP/indel data positions.

        Parameters
        ----------
        callback : callable
            Function to call after read a chromosome:
            callback(chrom, pos, ref, alt)
        sample : str
            Name of the sample (first one if empty - default)
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

        Returns
        -------
        count : int
            Number of read chromosomes

        """
        pos = []
        ref = []
        alt = []
        last_chrom = None
        count = 0
        alphabet = ['A', 'T', 'G', 'C', '.']

        try:
            for rec in self.file.fetch():
                if last_chrom is None:
                    last_chrom = rec.chrom
                if last_chrom != rec.chrom:
                    callback(last_chrom, pos, ref, alt)
                    pos = []
                    ref = []
                    alt = []
                    count += 1
                if ("PASS" in rec.filter.keys() or not filter) and rec.alts and len(rec.alts) >= 1:
                    if (rec.ref in alphabet) and (rec.alts[0] in alphabet):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                last_chrom = rec.chrom
            if len(pos) > 0:
                callback(last_chrom, pos, ref, alt)
                count += 1
            return count
        except ValueError:
            _logger.error("Variant file reading problem.")
            exit(0)

    def read_all_rd(self, callback, sample='', ad_tag='AD', dp_tag='DP', filter=False):
        """
        Read RD data for given chromosome and sample.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        sample : str
            Name of the sample (first one if empty - default)
        ad_tag : str
            AD tag in VCF file (default: 'AD')
        dp_tag : str
            DP tag in VCF file (default: 'DP')
        filter : bool
            If True it will read only variants with PASS filter, otherwise all

        Returns
        -------
        rd_p : list of int
            binned RD signal from DP tag (100 bins)
        rd_u : list of int
            binned RD signal from AD tag (100 bins)

        """
        if sample == '':
            sample = self.samples[0]
        pos = []
        rd_ad = []
        rd_dp = []
        last_chrom = None
        chr_count = 0

        try:
            for rec in self.file.fetch():
                if last_chrom is None:
                    last_chrom = rec.chrom
                if last_chrom != rec.chrom and len(pos) > 0:
                    n = pos[-1] + 1
                    if self.lengths[last_chrom] is not None:
                        n = self.lengths[last_chrom]
                    n = n // 100 + 1
                    rd_p = np.zeros(n)
                    rd_u = np.zeros(n)
                    count = np.zeros(n)
                    for p, rd1, rd2 in zip(pos, rd_dp, rd_ad):
                        cgb = (p - 1) // 10000 * 100
                        rd_p[cgb] += rd1
                        rd_u[cgb] += rd2
                        count[cgb] += 1
                    for i in range(n // 100 + 1):
                        if (100 * i) < n:
                            if count[100 * i] != 0:
                                rd_p[i * 100:(i + 1) * 100] = rd_p[i * 100] / count[i * 100]
                                rd_u[i * 100:(i + 1) * 100] = rd_u[i * 100] / count[i * 100]
                            else:
                                rd_p[i * 100:(i + 1) * 100] = np.nan
                                rd_u[i * 100:(i + 1) * 100] = np.nan
                    callback(last_chrom, rd_p, rd_u)
                    pos = []
                    rd_ad = []
                    rd_dp = []
                    chr_count += 1

                if ("PASS" in rec.filter.keys() or not filter) and (ad_tag in rec.samples[sample].keys()) and len(
                        rec.samples[sample][ad_tag]) > 1 and (
                        dp_tag in rec.samples[sample].keys()) and (rec.samples[sample][dp_tag] is not None):
                    try:
                        rd1 = sum(map(int, rec.samples[sample][ad_tag]))
                        rd2 = int(rec.samples[sample][dp_tag])
                        pos.append(rec.pos)
                        rd_ad.append(rd1)
                        rd_dp.append(rd2)
                    except ValueError:
                        pass
                last_chrom = rec.chrom
            if len(pos) > 0:
                n = pos[-1] + 1
                if self.lengths[last_chrom] is not None:
                    n = self.lengths[last_chrom]
                n = n // 100 + 1
                rd_p = np.zeros(n)
                rd_u = np.zeros(n)
                count = np.zeros(n)

                for p, rd1, rd2 in zip(pos, rd_dp, rd_ad):
                    cgb = (p - 1) // 10000 * 100
                    rd_p[cgb] += rd1
                    rd_u[cgb] += rd2
                    count[cgb] += 1
                for i in range(n // 100 + 1):
                    if (100 * i) < n:
                        if count[100 * i] != 0:
                            rd_p[i * 100:(i + 1) * 100] = rd_p[i * 100] / count[i * 100]
                            rd_u[i * 100:(i + 1) * 100] = rd_u[i * 100] / count[i * 100]
                        else:
                            rd_p[i * 100:(i + 1) * 100] = np.nan
                            rd_u[i * 100:(i + 1) * 100] = np.nan
                callback(last_chrom, rd_p, rd_u)
                chr_count += 1
            return chr_count
        except ValueError:
            _logger.error("Variant file reading problem. Probably index file is missing or corrupted.")
            exit(0)


class CreateVCF:
    def __init__(self, filename, reference_genome, chromosome_length_dct, c_date):
        """
        Create a VCF file

        Parameters
        ----------
        filename : str
            Name of the VCF file.
        reference_genome: dict
            reference genome dct
        chromosome_length_dct: dict
            chromosome length dct; format name: length
        c_date : date
            current date
        """
        self.filename = filename
        self.reference_genome = reference_genome
        self.chromosome_length_dct = chromosome_length_dct
        try:
            self.vcf = pysam.VariantFile(filename, "w",  header=self.vcf_header(current_date=c_date))
        except IOError:
            _logger.error("Problem opening file '%s'!" % filename)
            exit(0)
        except ValueError:
            _logger.error("Problem opening file '%s'!" % filename)
            exit(0)

    def vcf_header(self, current_date):
        if self.reference_genome:
            rg = self.reference_genome["name"]
        else:
            rg = "unknown"

        # Create a VCF header
        vcfh = pysam.VariantHeader()

        # add source
        vcfh.add_line(f"##fileDate={current_date}")
        vcfh.add_line(f"##reference={rg}")
        vcfh.add_line("##source=CNVpytor")
        # INFO
        vcfh.add_meta('INFO', items=[('ID', "END"), ('Number', 1), ('Type', 'Integer'),
                                     ('Description', 'End position of the variant described in this record')])
        vcfh.add_meta('INFO', items=[('ID', "IMPRECISE"), ('Number', 0), ('Type', 'Flag'),
                                     ('Description', 'Imprecise structural variation')])
        vcfh.add_meta('INFO', items=[('ID', "SVLEN"), ('Number', 1), ('Type', 'Integer'),
                                     ('Description', 'Difference in length between REF and ALT alleles')])
        vcfh.add_meta('INFO', items=[('ID', "SVTYPE"), ('Number', 1), ('Type', 'String'),
                                     ('Description', 'Type of structural variant')])
        vcfh.add_meta('INFO', items=[('ID', "pytorRD"), ('Number', 1), ('Type', 'Float'),
                                     ('Description', 'Normalized RD')])
        vcfh.add_meta('INFO', items=[('ID', "pytorP1"), ('Number', 1), ('Type', 'Float'),
                                     ('Description', 'e-val by t-test')])
        vcfh.add_meta('INFO', items=[('ID', "pytorP2"), ('Number', 1), ('Type', 'Float'),
                                     ('Description', 'e-val by Gaussian tail')])
        vcfh.add_meta('INFO', items=[('ID', "pytorP3"), ('Number', 1), ('Type', 'Float'),
                                     ('Description', 'e-val by t-test (middle)')])
        vcfh.add_meta('INFO', items=[('ID', "pytorP4"), ('Number', 1), ('Type', 'Float'),
                                     ('Description', 'Type of structural variant')])
        vcfh.add_meta('INFO', items=[('ID', "pytorQ0"), ('Number', 1), ('Type', 'Float'),
                                     ('Description', 'Fraction of reads with 0 mapping quality')])
        vcfh.add_meta('INFO', items=[('ID', "pytorPN"), ('Number', 1), ('Type', 'Integer'),
                                     ('Description', 'Fraction of N bases')])
        vcfh.add_meta('INFO', items=[('ID', "pytorDG"), ('Number', 1), ('Type', 'Integer'),
                                     ('Description', 'Distance to nearest gap in reference genome')])
        vcfh.add_meta('INFO', items=[('ID', "pytorCL"), ('Number', 1), ('Type', 'String'),
                                     ('Description', 'Caller method')])
        vcfh.add_meta('INFO', items=[('ID', "SAMPLES"), ('Number', 1), ('Type', 'String'),
                                      ('Description', 'Sample genotyped to have the variant')])

        # vcfh.alts({'ID': 'Del', 'Description': 'Deletion'})
        # ALT
        vcfh.add_meta('ALT', items=[('ID', "DEL"), ('Description', 'Deletion')])

        vcfh.add_meta('ALT', items=[('ID', "DUP"), ('Description', 'Duplication')])
        vcfh.add_meta('ALT', items=[('ID', "LOH"), ('Description', 'Copy number neutral loss of heterozygosity')])
        # FORMAT
        vcfh.add_meta('FORMAT', items=[('ID', "GT"), ('Number', 1), ('Type', 'String'), ('Description', 'Genotype')])
        vcfh.add_meta('FORMAT', items=[('ID', "CN"), ('Number', 1), ('Type', 'Float'),
                                       ('Description', 'Copy number genotype for imprecise events')])

        return vcfh

    def insert_all_contigs(self):
        for c, c_length in self.chr_len_dict.items():
            self.vcf.header.contigs.add(c, length=c_length)

    def insert_records(self, calls):
        samples = []
        chr_list = []
        for call in calls:
            call_chr = str(call[3])
            sample = call[0]
            if sample not in samples:
                self.vcf.header.add_sample(sample)
                samples.append(sample)

            if call_chr not in chr_list:
                chr_list.append(call_chr)
                if call_chr in self.chromosome_length_dct:
                    self.vcf.header.contigs.add(call_chr, length=self.chromosome_length_dct[call_chr])
                else:
                    self.vcf.header.contigs.add(call_chr)

        for idx, call in enumerate(calls):
            if len(call) == 0:
                continue

            record_id = "CNVpytor_" + {"deletion": "del", "duplication": "dup", "cnnloh": "loh"}[call[2]] + str(idx)
            alt = {"deletion": "<DEL>", "duplication": "<DUP>", "cnnloh": "<LOH>"}[call[2]]

            call_chr = str(call[3])
            start_loci = int(call[4])
            stop_loci = int(call[5])
            if start_loci != 1:
                start_loci = start_loci - 1

            info_dct = {"END": stop_loci, "IMPRECISE": True, "SVLEN": int(call[6]), "SVTYPE": alt[1:4],
                        "pytorRD": call[7], "pytorP1": call[8], "pytorP2": call[9], "pytorP3": call[10],
                        "pytorP4": call[11], "pytorQ0": call[12], "pytorPN": int(call[13]), "pytorDG": int(call[14]),
                        "pytorCL": call[1]}

            r = self.vcf.new_record(contig=call_chr, start=start_loci, stop=stop_loci, alleles=('N', alt), id=record_id,
                                    filter="PASS", info=info_dct)

            # add sample information
            for sample in samples:
                if sample == call[0]:
                    if call[2] == "deletion" and call[7] < 0.25:
                        r.samples[sample]['GT'] = (1, 1)
                        r.samples[sample]['CN'] = 0
                    elif call[2] == "deletion" and call[7] > 0.25:
                        r.samples[sample]['GT'] = (0, 1)
                        r.samples[sample]['CN'] = 1
                    elif call[2] == "duplication" and call[7] <= 1.75:
                        r.samples[sample]['GT'] = (0, 1)
                        r.samples[sample]['CN'] = 3
                    elif call[2] == "duplication" and 1.75 < call[7] <= 2.25:
                        r.samples[sample]['GT'] = (None, 1)
                        r.samples[sample]['CN'] = 4
                    elif call[2] == "duplication" and call[7] > 2.25:
                        r.samples[sample]['GT'] = (None, 1)
                        r.samples[sample]['CN'] = round(call[7] * 2)
                    else:
                        r.samples[sample]['GT'] = (None, None)
                        r.samples[sample]['CN'] = None

            self.vcf.write(r)

        # Close the VCF file
        self.vcf.close()

