""" cnvpytor.vcf

class Vcf: reads SNP data from VCF file

"""
import pysam
import logging

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
        if self.file:
            _logger.info("File: " + filename + " successfully open")
            self.chrs = list(self.file.header.contigs)
            self.samples = list(self.file.header.samples)
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

    def read_chromosome_snp(self, chr_name, sample=''):
        """
        Read SNP/indel data for given chromosome and sample.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        sample : str
            Name of the sample (first one if empty - default)

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
        alphabet = ['A', 'T', 'G', 'C', '.']
        try:
            for rec in self.file.fetch(chr_name):
                if "PASS" in rec.filter.keys() and rec.alts and len(rec.alts) == 1 and (
                        "GT" in rec.samples[sample].keys()) and (
                        "AD" in rec.samples[sample].keys()) and len(rec.samples[sample]["GT"]) > 1 and len(
                    rec.samples[sample]["AD"]) > 1:
                    if (len(rec.samples[sample]["AD"]) > 1) and (rec.ref in alphabet) and (rec.alts[0] in alphabet):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        nref.append(rec.samples[sample]["AD"][0])
                        nalt.append(rec.samples[sample]["AD"][1])
                        gt.append(rec.samples[sample]["GT"][0] * 2 + rec.samples[sample]["GT"][1])
                        if rec.samples[sample].phased:
                            gt[-1] += 4
        except ValueError:
            _logger.error("Variant file reading problem. Probably index file is missing or corrupted.")
            exit(0)
        return pos, ref, alt, nref, nalt, gt, flag, qual

    def read_chromosome_snp_no_counts(self, chr_name, sample=''):
        """
        Read SNP/indel data without counts (AD tag) for given chromosome and sample.

        Parameters
        ----------
        chr_name : str
            Name of the chromosome.
        sample : str
            Name of the sample (first one if empty - default)

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
                if "PASS" in rec.filter.keys() and rec.alts and len(rec.alts) == 1 and (
                        "GT" in rec.samples[sample].keys()) and len(rec.samples[sample]["GT"]) > 1:
                    if (rec.ref in alphabet) and (rec.alts[0] in alphabet):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        gt.append(rec.samples[sample]["GT"][0] * 2 + rec.samples[sample]["GT"][1])
                        if rec.samples[sample].phased:
                            gt[-1] += 4
        except ValueError:
            _logger.error("Variant file reading problem. Probably index file is missing or corrupted.")
            exit(0)
        return pos, ref, alt, gt, flag, qual

    def read_all_snp(self, callback, sample=''):
        """
        Read SNP/indel data for given sample name.

        Parameters
        ----------
        callback : callable
            Function to call after read a chromosome:
            callback(chrom, pos, ref, alt, nref, nalt, gt, flag, qual)

        sample : str
            Name of the sample (first one if empty - default)

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
        last_chrom = None
        count = 0
        alphabet = ['A', 'T', 'G', 'C', '.']

        try:
            for rec in self.file.fetch():
                if last_chrom is None:
                    last_chrom = rec.chrom
                if last_chrom != rec.chrom:
                    callback(last_chrom, pos, ref, alt, nref, nalt, gt, flag, qual)
                    pos = []
                    ref = []
                    alt = []
                    nref = []
                    nalt = []
                    gt = []
                    flag = []
                    qual = []
                    count += 1

                if "PASS" in rec.filter.keys() and rec.alts and len(rec.alts) == 1 and (
                        "GT" in rec.samples[sample].keys()) and (
                        "AD" in rec.samples[sample].keys()) and len(rec.samples[sample]["GT"]) > 1 and len(
                    rec.samples[sample]["AD"]) > 1:
                    if (len(rec.samples[sample]["AD"]) > 1) and (rec.ref in alphabet) and (
                            rec.alts[0] in alphabet) and (rec.samples[sample]["GT"][0] is not None) and (
                            rec.samples[sample]["GT"][1] is not None):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)   # Assign P region as default (-mask_snps updates this flag)
                        qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        nref.append(rec.samples[sample]["AD"][0])
                        nalt.append(rec.samples[sample]["AD"][1])
                        gt.append(rec.samples[sample]["GT"][0] * 2 + rec.samples[sample]["GT"][1])
                        if rec.samples[sample].phased:
                            gt[-1] += 4
                last_chrom = rec.chrom
            if len(pos) > 0:
                callback(last_chrom, pos, ref, alt, nref, nalt, gt, flag, qual)
                count += 1
            return count
        except ValueError:
            _logger.error("Variant file reading problem.")
            exit(0)

    def read_all_snp_no_counts(self, callback, sample=''):
        """
        Read SNP/indel data without counts (AD tag) for given sample name.

        Parameters
        ----------
        callback : callable
            Function to call after read a chromosome:
            callback(chrom, pos, ref, alt, gt, flag, qual)

        sample : str
            Name of the sample (first one if empty - default)

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

                if "PASS" in rec.filter.keys() and rec.alts and len(rec.alts) == 1 and (
                        "GT" in rec.samples[sample].keys()) and len(rec.samples[sample]["GT"]) > 1:
                    if (rec.ref in alphabet) and (rec.alts[0] in alphabet) and (
                            rec.samples[sample]["GT"][0] is not None) and (rec.samples[sample]["GT"][1] is not None):
                        pos.append(rec.pos)
                        ref.append(rec.ref)
                        alt.append(rec.alts[0])
                        flag.append(2)  # Assign P region as default (-mask_snps updates this flag)
                        qual.append(int(rec.qual / 10))  # divide QUAL by factor 10 and truncate to one byte
                        if qual[-1] > 255:
                            qual[-1] = 255
                        gt.append(rec.samples[sample]["GT"][0] * 2 + rec.samples[sample]["GT"][1])
                        if rec.samples[sample].phased:
                            gt[-1] += 4
                last_chrom = rec.chrom
            if len(pos) > 0:
                callback(last_chrom, pos, ref, alt, gt, flag, qual)
                count += 1
            return count
        except ValueError:
            _logger.error("Variant file reading problem.")
            exit(0)
