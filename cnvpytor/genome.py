""" cnvpytor.genome

class Genome: detect reference / parity / naming functions / reference genome data files

"""
from __future__ import absolute_import, print_function, division
from .utils import *
from .utils import get_install_path
from collections import OrderedDict
import logging
import os

_logger = logging.getLogger("cnvpytor.genome")


class Genome:
    pytor_datapath = get_install_path().joinpath("data")
    reference_genomes = {
        "hg19": {
            "name": "GRCh37",
            "species": "human",
            "chromosomes": OrderedDict(
                [("chr1", (249250621, "A")), ("chr2", (243199373, "A")), ("chr3", (198022430, "A")),
                 ("chr4", (191154276, "A")), ("chr5", (180915260, "A")), ("chr6", (171115067, "A")),
                 ("chr7", (159138663, "A")), ("chr8", (146364022, "A")), ("chr9", (141213431, "A")),
                 ("chr10", (135534747, "A")), ("chr11", (135006516, "A")), ("chr12", (133851895, "A")),
                 ("chr13", (115169878, "A")), ("chr14", (107349540, "A")), ("chr15", (102531392, "A")),
                 ("chr16", (90354753, "A")), ("chr17", (81195210, "A")), ("chr18", (78077248, "A")),
                 ("chr19", (59128983, "A")), ("chr20", (63025520, "A")), ("chr21", (48129895, "A")),
                 ("chr22", (51304566, "A")), ("chrX", (155270560, "S")), ("chrY", (59373566, "S")),
                 ("chrM", (16571, "M"))]),
            "gc_file": pytor_datapath.joinpath("gc_hg19.pytor"),
            "mask_file": pytor_datapath.joinpath("mask_hg19.pytor"),
            "ensembl_api_region": "https://grch37.rest.ensembl.org/overlap/region/human/{region}?content-type=application/json;feature=gene;"
        },
        "hg38": {
            "name": "GRCh38",
            "species": "human",
            "chromosomes": OrderedDict(
                [("chr1", (248956422, "A")), ("chr2", (242193529, "A")), ("chr3", (198295559, "A")),
                 ("chr4", (190214555, "A")), ("chr5", (181538259, "A")), ("chr6", (170805979, "A")),
                 ("chr7", (159345973, "A")), ("chr8", (145138636, "A")), ("chr9", (138394717, "A")),
                 ("chr10", (133797422, "A")), ("chr11", (135086622, "A")), ("chr12", (133275309, "A")),
                 ("chr13", (114364328, "A")), ("chr14", (107043718, "A")), ("chr15", (101991189, "A")),
                 ("chr16", (90338345, "A")), ("chr17", (83257441, "A")), ("chr18", (80373285, "A")),
                 ("chr19", (58617616, "A")), ("chr20", (64444167, "A")), ("chr21", (46709983, "A")),
                 ("chr22", (50818468, "A")), ("chrX", (156040895, "S")), ("chrY", (57227415, "S")),
                 ("chrM", (16569, "M"))]),
            "gc_file": pytor_datapath.joinpath("gc_hg38.pytor"),
            "mask_file": pytor_datapath.joinpath("mask_hg38.pytor"),
            "ensembl_api_region": "https://rest.ensembl.org/overlap/region/human/{region}?content-type=application/json;feature=gene;"
        },
        "chm13": {
            "name": "CHM13 v2.0",
            "species": "human",
            "chromosomes": OrderedDict([
                ("chr1", (248387328, "A")), ("chr2", (242696752, "A")), ("chr3", (201105948, "A")),
                ("chr4", (193574945, "A")), ("chr5", (182045439, "A")), ("chr6", (172126628, "A")),
                ("chr7", (160567428, "A")), ("chr8", (146259331, "A")), ("chr9", (150617247, "A")),
                ("chr10", (134758134, "A")), ("chr11", (135127769, "A")), ("chr12", (133324548, "A")),
                ("chr13", (113566686, "A")), ("chr14", (101161492, "A")), ("chr15", (99753195, "A")),
                ("chr16", (96330374, "A")), ("chr17", (84276897, "A")), ("chr18", (80542538, "A")),
                ("chr19", (61707364, "A")), ("chr20", (66210255, "A")), ("chr21", (45090682, "A")),
                ("chr22", (51324926, "A")), ("chrX", (154259566, "S")), ("chrY", (62460029, "S")),
                ("chrMT", (16569, "M"))
            ]),
            "gc_file": pytor_datapath.joinpath("gc_chm13v2.0.pytor"),
        },
        "chm13v1.1": {
            "name": "CHM13 v1.1 + hg38 chrY",
            "species": "human",
            "chromosomes": OrderedDict([
                ("chr1", (248387328, "A")), ("chr2", (242696752, "A")), ("chr3", (201105948, "A")),
                ("chr4", (193574945, "A")), ("chr5", (182045439, "A")), ("chr6", (172126628, "A")),
                ("chr7", (160567428, "A")), ("chr8", (146259331, "A")), ("chr9", (150617247, "A")),
                ("chr10", (134758134, "A")), ("chr11", (135127769, "A")), ("chr12", (133324548, "A")),
                ("chr13", (113566686, "A")), ("chr14", (101161492, "A")), ("chr15", (99753195, "A")),
                ("chr16", (96330374, "A")), ("chr17", (84276897, "A")), ("chr18", (80542538, "A")),
                ("chr19", (61707364, "A")), ("chr20", (66210255, "A")), ("chr21", (45090682, "A")),
                ("chr22", (51324926, "A")), ("chrX", (154259566, "S")), ("chrY", (57227415, "S")),
                ("chrMT", (16569, "M"))
            ]),
            "gc_file": pytor_datapath.joinpath("gc_chm13v1.1.pytor"),
        },
        "kn99": {
            "name": "Cryptococcus neoformans var. grubii KN99 GCA_002216725.1",
            "species": "cryptococcus neoformans",
            "chromosomes": OrderedDict([
                ("CP022321.1", (2291500, "A")), ("CP022322.1", (1621676, "A")), ("CP022323.1", (1574972, "A")),
                ("CP022324.1", (1084805, "A")), ("CP022325.1", (1814975, "A")), ("CP022326.1", (1422463, "A")),
                ("CP022327.1", (1399209, "A")), ("CP022328.1", (1398693, "A")), ("CP022329.1", (1186813, "A")),
                ("CP022330.1", (1059962, "A")), ("CP022331.1", (1562107, "A")), ("CP022332.1", (774060, "A")),
                ("CP022333.1", (756017, "A")), ("CP022334.1", (942472, "A")), ("CP022335.1", (24923, "M"))
            ]),
            "gc_file": pytor_datapath.joinpath("gc_kn99.pytor"),
        }
    }

    detected_genome = None

    @staticmethod
    def canonical_chrom_name(name):
        """
        Removes prefix chr, chrom or chromosome

        Parameters
        ----------
        name : str
            Name of the chromosome

        Returns
        -------
        cname : str
            Canonical chromosome name.

        """
        cname = name.upper().replace("CHROMOSOME", "").replace("CHROM", "").replace("CHR", "")
        if cname == "MT":
            cname = "M"
        if cname.startswith('NC_'):
            try:
                chrom_num = int(float(cname[3:]))
                if chrom_num <= 22:
                    cname = str(chrom_num)
                elif chrom_num == 23:
                    cname = 'X'
                elif chrom_num == 24:
                    cname = 'Y'
                elif chrom_num == 12920:
                    cname = 'M'
            except ValueError:
                cname = cname
        return cname

    @staticmethod
    def extended_chrom_name(name):
        """
        Add 'chr' prefix to the chromosome name

        Parameters
        ----------
        name : str
            Name of the chromosome

        Returns
        -------
        ename : str
            Extended chromosome name.

        """
        return "chr" + Genome.canonical_chrom_name(name)

    @classmethod
    def check_resources(cls):
        """
        Check do resource files exist.

        Returns
        -------
        ok : bool
            Returns True if all resource files exist.

        """
        _logger.debug("Checking reference genome resource files.")
        for i in cls.reference_genomes:
            if "gc_file" in cls.reference_genomes[i] and not os.path.exists(cls.reference_genomes[i]["gc_file"]):
                return False
            if "mask_file" in cls.reference_genomes[i] and not os.path.exists(cls.reference_genomes[i]["mask_file"]):
                return False
        return True

    @classmethod
    def download_resources(cls):
        """
        Download missing resource files from gitHub.

        Returns
        -------
        None

        """
        _logger.info("Updating reference genome resource files...")
        for i in cls.reference_genomes:
            if "gc_file" in cls.reference_genomes[i] and not os.path.exists(cls.reference_genomes[i]["gc_file"]):
                _logger.info("Detecting missing GC resource file for reference genome '%s'" % i)
                res = cls.reference_genomes[i]["gc_file"]
                fn = res.split("/")[-1]
                url = "https://github.com/abyzovlab/CNVpytor/raw/master/cnvpytor/data/" + fn
                if is_downloadable(url):
                    _logger.info("Downloading GC resource file: %s", fn)
                    try:
                        download(url, res)
                        _logger.info("File downloaded.")
                    except Exception as e:
                        _logger.error("Problem with downloading/saving resource files.")
                        _logger.error("Exception details: " + str(e))

                else:
                    _logger.warning("GC resource file is not downloadable!")
            if "mask_file" in cls.reference_genomes[i] and not os.path.exists(cls.reference_genomes[i]["mask_file"]):
                _logger.info("Detecting missing MASK resource file for reference genome '%s'" % i)
                res = cls.reference_genomes[i]["mask_file"]
                fn = res.split("/")[-1]
                url = "https://github.com/abyzovlab/CNVpytor/raw/master/cnvpytor/data/" + fn
                if is_downloadable(url):
                    _logger.info("Downloading MASK resource file: %s", fn)
                    try:
                        download(url, res)
                        _logger.info("File downloaded.")
                    except Exception as e:
                        _logger.error("Problem with downloading/saving resource files.")
                        _logger.error("Exception details: " + str(e))
        _logger.info("Done.")

    @classmethod
    def is_autosome(cls, name):
        """
        Checks is chromosome with given name listed as autosome in the reference genome.
        If reference genome is not detected, returns True if name is not equal to
        'M', 'X', 'Y' or 'SEX' and if does not contain 'GL' or 'NC'.

        Parameters
        ----------
        name : str
            Name of the chromosome

        Returns
        -------
        bool
            Return True if chromosome is autosome

        """
        if cls.detected_genome is None:
            return (not cls.is_sex_chrom(name)) and (not cls.is_mt_chrom(name)) and (not ("GL" in name.upper())) and (
                not ("NC" in name.upper()))
        if cls.extended_chrom_name(name) in cls.reference_genomes[cls.detected_genome]["chromosomes"]:
            return cls.reference_genomes[cls.detected_genome]["chromosomes"][cls.extended_chrom_name(name)][1] == "A"
        elif name in cls.reference_genomes[cls.detected_genome]["chromosomes"]:
            return cls.reference_genomes[cls.detected_genome]["chromosomes"][name][1] == "A"
        return False

    @classmethod
    def is_sex_chrom(cls, name):
        """
        Checks is chromosome with given name listed as sex chromosome in the reference genome.
        If reference genome is not detected, returns True if name is equal to 'X', 'Y' or 'SEX'.

        Parameters
        ----------
        name : str
            Name of the chromosome

        Returns
        -------
        bool
            Return True if chromosome is sex chromosome

        """
        if cls.detected_genome is None:
            return cls.canonical_chrom_name(name) in {"X", "Y", "SEX"}
        if cls.extended_chrom_name(name) in cls.reference_genomes[cls.detected_genome]["chromosomes"]:
            return cls.reference_genomes[cls.detected_genome]["chromosomes"][cls.extended_chrom_name(name)][1] == "S"
        elif name in cls.reference_genomes[cls.detected_genome]["chromosomes"]:
            return cls.reference_genomes[cls.detected_genome]["chromosomes"][name][1] == "S"
        return False

    @classmethod
    def is_mt_chrom(cls, name):
        """
        Checks is chromosome with given name listed as mitochondrial chromosome in the reference genome.
        If reference genome is not detected, returns True if name is equal to 'M' or 'MT'.

        Parameters
        ----------
        name : str
            Name of the chromosome

        Returns
        -------
        bool
            Return True if chromosome is mitochondrial chromosome

        """
        if cls.detected_genome is None:
            return cls.canonical_chrom_name(name) in {"M", "MT"}
        if cls.extended_chrom_name(name) in cls.reference_genomes[cls.detected_genome]["chromosomes"]:
            return cls.reference_genomes[cls.detected_genome]["chromosomes"][cls.extended_chrom_name(name)][1] == "M"
        elif name in cls.reference_genomes[cls.detected_genome]["chromosomes"]:
            return cls.reference_genomes[cls.detected_genome]["chromosomes"][name][1] == "M"

        return False

    @classmethod
    def detect_genome(cls, names, lengths):
        """
        Detects reference genome for given list od chromosome names and lengths.

        Parameters
        ----------
        names : list of str
            List of chromosome names.
        lengths : list of int
            List of chromosome lengths.

        Returns
        -------
        g : str or None
            Name of the reference genome if detected, otherwise None.

        """
        for g in cls.reference_genomes:
            found = True
            checked = False
            for c, l in zip(names, lengths):
                if ((c in cls.reference_genomes[g]["chromosomes"]) or (
                        cls.extended_chrom_name(c) in cls.reference_genomes[g]["chromosomes"])) and (
                        not cls.is_mt_chrom(c)):
                    checked = True
                    if c in cls.reference_genomes[g]["chromosomes"]:
                        found = found and (cls.reference_genomes[g]["chromosomes"][c][0] == l)
                    else:
                        found = found and (cls.reference_genomes[g]["chromosomes"][cls.extended_chrom_name(c)][0] == l)
            if checked and found:
                cls.detected_genome = g
                return g
        return None

    @classmethod
    def load_reference_genomes(cls, filename):
        """
        Load reference genomes from configuration file. File should be writen in format:

            |#File: example_ref_genome_conf.py
            |
            |import_reference_genomes = {
            |    "hg19": {
            |    "name": "GRCh37",
            |    "species": "human",
            |    "chromosomes": OrderedDict(
            |        [("chr1", (249250621, "A")), ("chr2", (243199373, "A")), ("chr3", (198022430, "A")),
            |        ("chr4", (191154276, "A")), ("chr5", (180915260, "A")), ("chr6", (171115067, "A")),
            |        ("chr7", (159138663, "A")), ("chr8", (146364022, "A")), ("chr9", (141213431, "A")),
            |        ("chr10", (135534747, "A")), ("chr11", (135006516, "A")), ("chr12", (133851895, "A")),
            |        ("chr13", (115169878, "A")), ("chr14", (107349540, "A")), ("chr15", (102531392, "A")),
            |        ("chr16", (90354753, "A")), ("chr17", (81195210, "A")), ("chr18", (78077248, "A")),
            |        ("chr19", (59128983, "A")), ("chr20", (63025520, "A")), ("chr21", (48129895, "A")),
            |        ("chr22", (51304566, "A")), ("chrX", (155270560, "S")), ("chrY", (59373566, "S")),
            |        ("chrM", (16571, "M"))]),
            |    "gc_file": "/path/gc_file.pytor",
            |    "mask_file": "/path/mask_file.pytor"
            |    }
            |}


        Parameters
        ----------
        filename : str
            Name of the configuration file

        Returns
        -------
        None

        """
        _logger.info("Reading configuration file '%s'." % filename)
        exec(open(filename).read(), globals())
        for g in import_reference_genomes:
            _logger.info("Importing reference genome data: '%s'." % g)
            cls.reference_genomes[g] = import_reference_genomes[g]

    @classmethod
    def available_reference_genomes(cls):
        """
        Returns available reference genome keys

        Returns
        -------
        None

        """
        return [rg for rg in cls.reference_genomes]

    @classmethod
    def print_reference_genomes(cls):
        """
        Prints loaded reference genomes

        Returns
        -------
        None

        """
        print("Available reference genomes:")
        for rg in cls.reference_genomes:
            print(f"* {rg}")
            print(f"\tName: {cls.reference_genomes[rg]['name']}")
            print(f"\tSpecies: {cls.reference_genomes[rg]['species']}")
            print(f"\tChromosomes: {len(cls.reference_genomes[rg]['chromosomes'])}")
            print(f"\tGC file specified: {'Yes' if 'gc_file' in cls.reference_genomes[rg] else 'Not'}")
            print(f"\tMask file specified: {'Yes' if 'mask_file' in cls.reference_genomes[rg] else 'Not'}")
            print(f"\tEnsembl API link specified: {'Yes' if 'ensembl_api_region' in cls.reference_genomes[rg] else 'Not'}")





