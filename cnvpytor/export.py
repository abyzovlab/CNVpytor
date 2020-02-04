
import json
from pathlib import Path
from .io import *
from .genome import *

_logger = logging.getLogger("cnvpytor.export")


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

    def get_cnvpytor_signal(self, md5, chrom, bin_size, signal, flag):
        signal_details = md5.get_signal(chrom, bin_size, signal, flag)
        return signal_details

    def get_chrom_list(self, md5):
        chr_len = md5.get_signal(None, None, "chromosome lengths")
        chr_len_list = list(zip(chr_len[::2].astype(str), chr_len[1::2].astype(int)))
        return chr_len_list

    def create_wig_offset_transform(self, md5, chr_list, bin_size, signal, flag, offset):
        # add chr_list to add wig header
        self.add_header_list(chr_list)

        # add the data
        for (chrom, length) in chr_list:
            signal_details = md5.get_signal(chrom, bin_size, signal, flag)
            if isinstance(signal_details, np.ndarray):
                signal_value_list = signal_details[()]
                signal_value_list[signal_value_list != 0] += offset
                signal_value_list = np.absolute(signal_value_list)
                self.add_fixedstep(chrom, 0, signal_value_list, span=bin_size, step=bin_size)

    def create_wig(self, md5, chr_list, bin_size, signal, flag):
        # add chr_list to add wig header
        self.add_header_list(chr_list)

        # add the data
        for (chrom, length) in chr_list:
            signal_details = md5.get_signal(chrom, bin_size, signal, flag)
            if isinstance(signal_details, np.ndarray):
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
        "SNP baf": {
            "FLAG": [0x0100],
            "color": ["black"],
            "nonCont": [True],
        },
        "SNP i1": {
            "FLAG": [0x0100, 0x0100],
            "color": ["red", "red"],
            "nonCont": [True, True],
            "offset": [0.5, -0.5]
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
    def export_directory(self):
        main_dir_name = "jbrowse_{}".format(self.pytor_name)
        return main_dir_name

    @property
    def export_main_dir(self):
        main_dir = Path(self.dir).joinpath(self.export_directory)
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
        _logger.debug("Create Read depth related signals")
        rd_chr, rd_bin = self.rd_chr_bin()

        # get chr list
        chr_len = self.io.get_signal(None, None, "chromosome lengths")
        chr_list = list(zip(chr_len[::2].astype(str), chr_len[1::2].astype(int)))

        for signal_name, signal_dct in self.rd_signal_dct.items():
            for index, flag in enumerate(signal_dct['FLAG']):
                for bin_size in rd_bin:
                    signal = self.signal_name(bin_size, signal_name, flag)
                    bigwig_filename = "{}.bw".format(signal)
                    bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                    bigwig_file = str(bigwig_file)

                    wig = Wiggle(bigwig_file)
                    wig.create_wig(self.io, chr_list, bin_size, signal_name, flag)

    def snp_signal(self):
        _logger.debug("Create SNP related signals")
        snp_chr, snp_bin = self.snp_chr_bin()

        # get chr list
        chr_len = self.io.get_signal(None, None, "chromosome lengths")
        chr_list = list(zip(chr_len[::2].astype(str), chr_len[1::2].astype(int)))

        for signal_name, signal_dct in self.snp_signal_dct.items():
            for index, flag in enumerate(signal_dct['FLAG']):
                if "offset" in signal_dct:
                    offset = signal_dct['offset'][index]
                    for bin_size in snp_bin:
                        signal = self.signal_name(bin_size, signal_name, flag)

                        bigwig_filename = "{}_offset{}.bw".format(signal, offset)
                        bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                        bigwig_file = str(bigwig_file)
                        wig = Wiggle(bigwig_file)
                        wig.create_wig_offset_transform(self.io, chr_list, bin_size, signal_name, flag, offset)

                else:
                    for bin_size in snp_bin:
                        signal = self.signal_name(bin_size, signal_name, flag)
                        bigwig_filename = "{}.bw".format(signal)
                        bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                        bigwig_file = str(bigwig_file)

                        wig = Wiggle(bigwig_file)
                        wig.create_wig(self.io, chr_list, bin_size, signal_name, flag)

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

    def add_rd_config_track(self):
        _logger.debug("Get RD config track")
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
        _logger.debug("Get SNP config track info")
        snp_url_dct_list = []
        snp_chr, snp_bin = self.snp_chr_bin()
        for signal_name, signal_dct in self.snp_signal_dct.items():
            for index, flag in enumerate(signal_dct['FLAG']):
                suffix_flag = Signals.suffix_snp_flag(flag)
                scales = {}
                if "offset" in signal_dct:
                    offset = signal_dct['offset'][index]
                    signal_id = "{}{}{}".format(signal_name, suffix_flag, offset)
                    for bin_size in snp_bin:
                        signal = self.signal_name(bin_size, signal_name, flag)
                        bigwig_filename = "{}_offset{}.bw".format(signal, offset)
                        bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                        bigwig_current_path = Path(bigwig_file.parent.name).joinpath(bigwig_file.name).as_posix()
                        scales[bin_size] = bigwig_current_path
                else:
                    signal_id = "{}{}".format(signal_name, suffix_flag)
                    for bin_size in snp_bin:
                        signal = self.signal_name(bin_size, signal_name, flag)
                        bigwig_filename = "{}.bw".format(signal)
                        bigwig_file = self.export_data_dir.joinpath(bigwig_filename)
                        bigwig_current_path = Path(bigwig_file.parent.name).joinpath(bigwig_file.name).as_posix()
                        scales[bin_size] = bigwig_current_path

                snp_url_dct_list.append({
                    "storeClass": "MultiScaleBigWig/Store/SeqFeature/MultiScaleBigWig",
                    "scales": scales,
                    "name": signal_id,
                    "color": signal_dct['color'][index],
                    "nonCont": signal_dct['nonCont'][index]
                })
        track_dct = {
            "category": self.pytor_name,
            'autoscale': 'local',
            "storeClass": "MultiBigWig/Store/SeqFeature/MultiBigWig",
            "showTooltips": True,
            "showLabels": True,
            "clickTooltips": True,
            "max_score": 1,
            "label": "SNP",
            "type": "MultiBigWig/View/Track/MultiWiggle/MultiXYPlot",
            'urlTemplates': snp_url_dct_list,
        }
        return track_dct

    def create_tracklist_json(self):
        _logger.debug("Creates config file: {}".format(self.export_tracklist_file))

        # create rd config
        track_dct = self.add_rd_config_track()

        # create SNP config
        snp_track = self.add_snp_config_track()
        track_dct['tracks'].append(snp_track)

        with open(self.export_tracklist_file, 'w') as f:
            json.dump(track_dct, f, indent=2)
        return track_dct

    def create_reference_json(self):
        _logger.debug("Exporting reference details")
        # get signal details
        chr_len = list(np.array(self.io.get_signal(None, None, "chromosome lengths")).astype("str"))
        chr_dct = dict(zip(chr_len[::2], chr_len[1::2]))

        # create signal list in proper format
        chr_dct_list = []
        for chr, length in chr_dct.items():
            tmp_dct = {"end": length, "length": length, "name": chr, "start": 0}
            chr_dct_list.append(tmp_dct)

        # save it to file
        with open(self.export_ref_file, 'w') as f:
            json.dump(chr_dct_list, f, indent=2)

    def __del__(self):
        _logger.info("Export to: {}".format(self.export_main_dir))
        _logger.info("Copy this directory to jbrowse directory if export path is not set to jbrowse path")
        _logger.info("To access this via localhost: http://localhost/jbrowse/?data={}".format(self.export_directory))
