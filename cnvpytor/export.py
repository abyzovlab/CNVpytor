
import json
from pathlib import Path
from .io import *
from .genome import *


class Wiggle:
    def __init__(self, filename):
        self.filename = filename
        self.file = None
        import pyBigWig

        if not Path(filename).exists():
            try:
                self.file = pyBigWig.open(filename, 'w')
            except IOError as e:
                print("Unable to open file {}! Error: {}".format(filename, e))
            except RuntimeError as e:
                print("Unable to open file {}! Error: {}".format(filename, e))
        else:
            self.file = pyBigWig.open(filename)

    def add_header_list(self, chr_len_list):
        self.file.addHeader(chr_len_list)

    def add_fixedstep(self, chrom, position_int, value_list, span=1, step=1):
        self.file.addEntries(chrom, position_int, values=value_list, span=span, step=step)

    def create_wig(self, md5, chrom_list, bin_size, signal, flag):
        header_list = []

        for chrom in chrom_list:
            signal_details = md5.get_signal(chrom, bin_size, signal, flag)
            header_data = (chrom, signal_details.size * bin_size)
            header_list.append(header_data)
        self.add_header_list(header_list)

        for chrom in chrom_list:
            signal_details = md5.get_signal(chrom, bin_size, signal, flag)
            signal_value_list = signal_details[()]
            self.add_fixedstep(chrom, 0, signal_value_list, span=bin_size, step=bin_size)

    def __del__(self):

        if self.file:
            self.file.close()


class ExportJbrowse:

    rd_signal_dct = {
        "RD": {
            "FLAG": [0, 0x0010],
            "color": ["gray", "black"]
        },
        "RD partition": {
            "FLAG": [0x0010],
            "color": ["green"]
        },
        "RD call": {
            "FLAG": [0x0010],
            "color": ["red"]
        }
    }
    snp_signal_dct = {
        "SNP maf": {
            "FLAG": [0x0100]
        },
        "SNP baf": {
            "FLAG": [0x0100]
        },
        "SNP i1": {
            "FLAG": [0x0100]
        }
    }

    signal_dct = {
        "RD": "his_rd_p_%(bin_size)d%(rd_flag)s",
        "RD partition": "his_rd_p_%(bin_size)d_partition%(rd_flag)s",
        "RD call": "his_rd_p_%(bin_size)d_partition%(rd_flag)s_merge",
        "SNP baf": "snp_baf_%(bin_size)d%(snp_flag)s",
        "SNP maf": "snp_maf_%(bin_size)d%(snp_flag)s",
        "SNP i1": "snp_i1_%(bin_size)d%(snp_flag)s",
        "SNP i1 partition": "snp_i1_%(bin_size)d%(snp_flag)s_partition",

    }

    def __init__(self, file, dir_name):
        self.filename = file
        self.dir = dir_name
        self.io = IO(file, ro=True)

    @property
    def pytor_name(self):
        root_filename = Path(self.filename).resolve().stem
        return root_filename

    @property
    def export_main_dir(self):
        main_dir_name = "jbrowse_{}".format(self.pytor_name)
        main_dir = Path(self.dir).joinpath(main_dir_name)
        main_dir.mkdir(parents=True, exist_ok=True)
        return main_dir

    @property
    def export_data_dir(self):
        data_dir = self.export_main_dir.joinpath("bw")
        data_dir.mkdir(parents=True, exist_ok=True)
        return data_dir

    @property
    def export_seq_dir(self):
        seq_dir = self.export_main_dir.joinpath("seq")
        seq_dir.mkdir(parents=True, exist_ok=True)
        return seq_dir

    @property
    def export_tracklist_file(self):
        track_list = self.export_main_dir.joinpath("trackList.json")
        return track_list

    @property
    def export_tracks_file(self):
        track_list = self.export_main_dir.joinpath("tracks.conf")
        return track_list

    @property
    def export_ref_file(self):
        ref_file = self.export_seq_dir.joinpath("refSeqs.json")
        return ref_file

    def signal_name(self, bin_size, signal, flags=0):
        if signal in self.signal_dct:
            try:
                return self.signal_dct[signal] % {"bin_size": bin_size, "rd_flag": Signals().suffix_rd_flag(flags),
                                                  "snp_flag": Signals().suffix_snp_flag(flags),
                                                  "flag": Signals().suffix_flag(flags)}
            except TypeError:
                return None
        else:
            return None

    def rd_chr_bin(self):
        chr_bs = self.io.chromosomes_bin_sizes_with_signal("RD")
        chrs = {}
        bss = []
        for c, b in chr_bs:
            if c not in chrs:
                chrs[c] = []
            chrs[c].append(int(b))
            if int(b) not in bss:
                bss.append(int(b))
        return chrs, bss

    def snp_chr_bin(self):
        chr_bs = self.io.chromosomes_bin_sizes_with_signal("SNP likelihood", FLAG_USEMASK)
        chrs = {}
        bss = []
        for c, b in chr_bs:
            if c not in chrs:
                chrs[c] = []
            chrs[c].append(int(b))
            if int(b) not in bss:
                bss.append(int(b))
        return chrs, bss

    def rd_signal(self):
        rd_chr, rd_bin = self.rd_chr_bin()
        chrom_list = list(rd_chr)

        for signal_name, signal_dct in self.rd_signal_dct.items():
            for index, flag in enumerate(signal_dct['FLAG']):
                for bin_size in rd_bin:
                    signal = self.signal_name(bin_size, signal_name, flag)
                    bigwig_filename = "{}.bw".format(signal)
                    bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                    bigwig_file = str(bigwig_file)

                    wig = Wiggle(bigwig_file)
                    wig.create_wig(self.io, chrom_list, bin_size, signal_name, flag)

    def snp_signal(self):
        snp_chr, snp_bin = self.snp_chr_bin()
        chrom_list = list(snp_chr)

        for signal_name, signal_dct in self.snp_signal_dct.items():
            for index, flag in enumerate(signal_dct['FLAG']):
                for bin_size in snp_bin:
                    signal = self.signal_name(bin_size, signal_name, flag)
                    bigwig_filename = "{}.bw".format(signal)
                    bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                    bigwig_file = str(bigwig_file)

                    wig = Wiggle(bigwig_file)
                    wig.create_wig(self.io, chrom_list, bin_size, signal_name, flag)

    @staticmethod
    def add_config_reference():
        track_dct = {'formatVersion': 1, 'tracks': []}
        track_dct['tracks'].append({
            "category": "Reference sequence",
            "chunkSize": 20000,
            "key": "Reference sequence",
            "label": "DNA",
            "seqType": "dna",
            "storeClass": "JBrowse/Store/Sequence/StaticChunked",
            "type": "SequenceTrack",
            "urlTemplates": "seq/{refseq_dirpath}/{refseq}-"
        })
        return track_dct

    def add_config_track(self):
        rd_chr, rd_bin = self.rd_chr_bin()
        url_template_dct = []
        for signal_name, signal_dct in self.rd_signal_dct.items():
            if 'FLAG' in signal_dct:
                for index, flag in enumerate(signal_dct['FLAG']):
                    suffix_rd_flag = Signals.suffix_rd_flag(flag)
                    signal_id = "{}{}".format(signal_name, suffix_rd_flag)
                    scales = {}
                    for bin_size in rd_bin:
                        signal = self.signal_name(bin_size, signal_name, flag)
                        bigwig_filename = "{}.bw".format(signal)
                        bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                        bigwig_current_path = Path(bigwig_file.parent.name).joinpath(bigwig_file.name).as_posix()
                        scales[bin_size] = bigwig_current_path

                    url_template_dct.append({
                        "storeClass": "MultiScaleBigWig/Store/SeqFeature/MultiScaleBigWig",
                        "scales": scales,
                        "name": signal_id,
                        "color": signal_dct['color'][index]
                    })

        track_dct = self.add_config_reference()

        track_dct['tracks'].append({
            "category": self.pytor_name,
            'autoscale': 'local',
            "storeClass": "MultiBigWig/Store/SeqFeature/MultiBigWig",
            "showTooltips": True,
            "showLabels": True,
            "clickTooltips": True,
            "label": "RD",
            "type": "MultiBigWig/View/Track/MultiWiggle/MultiXYPlot",
            'urlTemplates': url_template_dct
        })
        return track_dct

    def add_snp_config_track(self):
        snp_dct_list = []
        snp_chr, snp_bin = self.snp_chr_bin()
        for signal_name, signal_dct in self.snp_signal_dct.items():
            for index, flag in enumerate(signal_dct['FLAG']):
                scales = {}
                for bin_size in snp_bin:
                    signal = self.signal_name(bin_size, signal_name, flag)
                    bigwig_filename = "{}.bw".format(signal)
                    bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                    bigwig_current_path = Path(bigwig_file.parent.name).joinpath(bigwig_file.name).as_posix()
                    scales[bin_size] = bigwig_current_path

                snp_dct = {
                    "category": self.pytor_name,
                    "autoscale": "global",
                    "storeClass": "MultiScaleBigWig/Store/SeqFeature/MultiScaleBigWig",
                    "scales": scales,
                    "showTooltips": True,
                    "label": signal_name,
                    "type": "JBrowse/View/Track/Wiggle/XYPlot",
                    "noFill": True
                }
                snp_dct_list.append(snp_dct)
        return snp_dct_list

    def create_tracklist_json(self):

        track_dct = self.add_config_track()
        for snp_dct in self.add_snp_config_track():
            track_dct['tracks'].append(snp_dct)

        with open(self.export_tracklist_file, 'w') as f:
            json.dump(track_dct, f, indent=2)
        return track_dct

    def create_reference_json(self):
        rg_name = np.array(self.io.get_signal(None, None, "reference genome")).astype("str")[0]
        if rg_name in Genome.reference_genomes:
            genome_dct = Genome.reference_genomes[rg_name]['chromosomes']
            chr_dct_list = []
            rdcs = self.io.rd_chromosomes()
            for key, value in genome_dct.items():
                if Genome.extended_chrom_name(key) in rdcs:
                    chr = Genome.extended_chrom_name(key)
                elif Genome.canonical_chrom_name(key) in rdcs:
                    chr = Genome.canonical_chrom_name(key)
                else:
                    chr = key
                dct_dct = {"end": value[0], "length": value[0], "name": chr, "start": 0}
                chr_dct_list.append(dct_dct)

            with open(self.export_ref_file, 'w') as f:
                json.dump(chr_dct_list, f, indent=2)

