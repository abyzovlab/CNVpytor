import json
from pathlib import Path
from .io import *


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
            # chrom = chrom.decode('UTF-8')
            #signal_name = md5.signal_name(chrom, bin_size, signal, flag)
            #signal_details = md5.get_signal(signal_name)
            signal_details = md5.get_signal(chrom, bin_size, signal, flag)
            header_data = (chrom, signal_details.size * bin_size)
            header_list.append(header_data)
        self.add_header_list(header_list)

        for chrom in chrom_list:
            # chrom = chrom.decode('UTF-8')

            # signal_name = md5.signal_name(chrom, bin_size, signal, flag)
            # signal_details = md5.get_signal(signal_name)
            signal_details = md5.get_signal(chrom, bin_size, signal, flag)
            signal_value_list = signal_details[()]
            self.add_fixedstep(chrom, 0, signal_value_list, span=bin_size, step=bin_size)

    def __del__(self):

        if self.file:
            self.file.close()


class ExportJbrowse:
    rd_signal_list = ['RD', 'RD partition', 'RD call']
    rd_flag_list = [[0, 0x0010], [0x0010], [0x0010]]
    snp_signal_list = ['SNP maf', 'SNP baf', 'SNP i1']

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
    def export_main_dir(self):
        root_filename = Path(self.filename).name
        main_dir_name = "jbrowse_{}".format(root_filename)
        main_dir = Path(self.dir).joinpath(main_dir_name)
        main_dir.mkdir(parents=True, exist_ok=True)
        return main_dir

    @property
    def export_data_dir(self):
        data_dir = self.export_main_dir.joinpath("bw")
        data_dir.mkdir(parents=True, exist_ok=True)
        return data_dir

    @property
    def export_config_file(self):
        track_list = self.export_main_dir.joinpath("trackList.json")
        return track_list

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

        for index, signal in enumerate(self.rd_signal_list):
            for flag in self.rd_flag_list[index]:
                for bin_size in rd_bin:
                    signal_name = self.signal_name(bin_size, signal, flag)
                    bigwig_filename = "{}.bw".format(signal_name)
                    bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                    bigwig_file = str(bigwig_file)

                    # print(bigwig_file)
                    wig = Wiggle(bigwig_file)
                    wig.create_wig(self.io, chrom_list, bin_size, signal, flag)

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

        for index, signal in enumerate(self.rd_signal_list):
            for flag in self.rd_flag_list[index]:
                suffix_rd_flag = Signals.suffix_rd_flag(flag)
                signal_id = "{}{}".format(signal, suffix_rd_flag)
                scales = {}
                for bin_size in rd_bin:
                    signal_name = self.signal_name(bin_size, signal, flag)
                    bigwig_filename = "{}.bw".format(signal_name)
                    bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                    bigwig_current_path = Path(bigwig_file.parent.name).joinpath(bigwig_file.name).as_posix()
                    # print(bigwig_current_path.)
                    # bigwig_file = str(bigwig_file)
                    scales[bin_size] = bigwig_current_path

                url_template_dct.append({
                    "storeClass": "MultiScaleBigWig/Store/SeqFeature/MultiScaleBigWig",
                    "scales": scales,
                    "name": signal_id
                })

        track_dct = self.add_config_reference()

        track_dct['tracks'].append({
            'autoscale': 'global',
            "storeClass": "MultiBigWig/Store/SeqFeature/MultiBigWig",
            "showTooltips": True,
            "showLabels": True,
            "clickTooltips": True,
            "label": "sample1 global",
            "type": "MultiBigWig/View/Track/MultiWiggle/MultiXYPlot",
            'urlTemplates': url_template_dct
        })
        return track_dct

    def create_tracklist_json(self):

        track_dct = self.add_config_track()
        with open(self.export_config_file, 'w') as f:
            json.dump(track_dct, f, indent=2)

    def create_reference_json(self):
        # signal_name = "reference_genome"
        # signal_name = "RD chromosomes"
        signal_name = "use reference"

        signal_data = np.array(self.io.get_signal(None, None, signal_name))
        print(signal_data)
        # signal_name = self.io.signal_name(chrom, bin_size, signal, flag)
        # signal_details = md5.get_signal(signal_name)



