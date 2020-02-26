import json
from pathlib import Path
from .io import *
from .genome import *


_logger = logging.getLogger("cnvpytor.export")


class Wiggle:
    def __init__(self, filename):
        """
        creates bigwig file
        Parameters
        ----------
        filename : path
            Path for the bigwig filename
        """
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
        """
        Add header to the bigwig file
        Parameters
        ----------
        chr_len_list : list of tuple
            chromosome name and length list.

        Returns
        -------

        """
        self.file.addHeader(chr_len_list)

    def add_fixedstep(self, chrom, position_int, value_list, span=1, step=1):
        """
        Add fixed step formatted data
        Parameters
        ----------
        chrom : str
            chromosome name
        position_int : int
            start position
        value_list : list of values
            input values
        span : int
        step : int

        Returns
        -------

        """
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


class ExportJBrowse:

    rd_signal_dct = {
        "RD": {
            "FLAG": [0, 0x0010],
            "color": ["gray", "black"]
        },
        "RD partition": {
            "FLAG": [0x0010],
            "color": ["red"]
        },
        "RD call": {
            "FLAG": [0x0010],
            "color": ["green"]
        }
    }
    snp_signal_dct = {
        "SNP baf": {
            "FLAG": [0x0100],
            "color": ["gray"],
            "nonCont": [True],
        },
        "SNP i1": {
            "FLAG": [0x0100, 0x0100],
            "color": ["red", "red"],
            "nonCont": [True, True],
            "offset": [0.5, -0.5]
        },
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

    def __init__(self, files, dir_name):
        """
        Exports CNVpytor data
        Parameters
        ----------
        files : path
            CNVpytor files path
        dir_name: path
            Export directory path
        """
        self.files = files
        self.dir = Path(dir_name)
        self.io = [IO(f, ro=True) for f in files]
        self.export_dir = self.export_create_dir()

    @property
    def pytor_names(self):
        name_list = []
        for filename in self.files:
            name_list.append(Path(filename).resolve().stem)
        return name_list

    @property
    def export_directory(self):
        if self.dir.is_dir():
            if len(self.files) > 1:
                # for multiple input file
                default_name = self.dir.joinpath("cnvpytor_jbrowse_export")

            else:
                # for single_input_file
                default_name = self.dir.joinpath("jbrowse_{}".format(self.pytor_names[0]))

            if default_name.exists():
                tmp_name = default_name
                i = 1
                while default_name.exists():
                    update_name = "{}({})".format(tmp_name.name, i)
                    default_name = default_name.with_name(update_name)
                    i = i+1
            return default_name
        else:
            if self.dir.parent.exists():
                return self.dir
            else:
                _logger.error("Error: incorrect export path: {}".format(self.dir))
                exit(0)

    def export_create_dir(self):
        main_dir = self.export_directory
        main_dir.mkdir(parents=True, exist_ok=True)
        _logger.info("CNVpytor data exporting for JBrowse view in {}".format(main_dir))
        return main_dir

    @property
    def export_data_dir_list(self):
        data_dir = self.export_dir.joinpath("bw")
        data_dir.mkdir(parents=True, exist_ok=True)
        data_dir_list = []
        for root_name in self.pytor_names:
            root_data = data_dir.joinpath(root_name)
            root_data.mkdir(parents=True, exist_ok=True)
            data_dir_list.append(root_data)
        return data_dir_list

    @property
    def export_seq_dir(self):
        seq_dir = self.export_dir.joinpath("seq")
        seq_dir.mkdir(parents=True, exist_ok=True)
        return seq_dir

    @property
    def export_tracklist_file(self):
        track_list = self.export_dir.joinpath("trackList.json")
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

    def rd_chr_bin(self, root_io):
        chr_bs = root_io.chromosomes_bin_sizes_with_signal("RD")
        chrs = {}
        bss = []
        for c, b in chr_bs:
            if c not in chrs:
                chrs[c] = []
            chrs[c].append(int(b))
            if int(b) not in bss:
                bss.append(int(b))
        return chrs, bss

    def snp_chr_bin(self, root_io):
        chr_bs = root_io.chromosomes_bin_sizes_with_signal("SNP likelihood", FLAG_USEMASK)
        chrs = {}
        bss = []
        for c, b in chr_bs:
            if c not in chrs:
                chrs[c] = []
            chrs[c].append(int(b))
            if int(b) not in bss:
                bss.append(int(b))
        return chrs, bss

    @staticmethod
    def create_bigwig(root_io, bigwig_file, chr_list, bin_size, signal_name, flag, offset=None):
        wig = None
        for (chrom, length) in chr_list:
            signal_details = root_io.get_signal(chrom, bin_size, signal_name, flag)
            if isinstance(signal_details, np.ndarray):
                signal_value_list = signal_details[()]
                if offset is not None:
                    signal_value_list[signal_value_list != 0] += offset
                    signal_value_list = np.absolute(signal_value_list)

                if not isinstance(wig, Wiggle):
                    wig = Wiggle(bigwig_file)
                    wig.add_header_list(chr_list)

                wig.add_fixedstep(chrom, 0, signal_value_list, span=bin_size, step=bin_size)

    def rd_signal(self):
        _logger.debug("Create Read depth related signals")
        for root_index, root_io in enumerate(self.io):
            _logger.info("JBrowse export: RD related data for {}".format(self.pytor_names[root_index]))
            rd_chr, rd_bin = self.rd_chr_bin(root_io)

            # get chr list
            chr_len = root_io.get_signal(None, None, "chromosome lengths")
            chr_list = list(zip(chr_len[::2].astype(str), chr_len[1::2].astype(int)))

            for signal_name, signal_dct in self.rd_signal_dct.items():
                _logger.info("JBrowse export: RD signal {}".format(signal_name))
                for index, flag in enumerate(signal_dct['FLAG']):
                    for bin_size in rd_bin:
                        signal = self.signal_name(bin_size, signal_name, flag)
                        bigwig_filename = "{}.bw".format(signal)
                        bigwig_file = self.export_data_dir_list[root_index].joinpath(bigwig_filename)
                        bigwig_file = str(bigwig_file)

                        self.create_bigwig(root_io, bigwig_file, chr_list, bin_size, signal_name, flag)

    def snp_signal(self):
        _logger.debug("Create SNP related signals")
        for root_index, root_io in enumerate(self.io):
            _logger.info("JBrowse export: SNP related data for {}".format(self.pytor_names[root_index]))
            snp_chr, snp_bin = self.snp_chr_bin(root_io)

            # get chr list
            chr_len = root_io.get_signal(None, None, "chromosome lengths")
            chr_list = list(zip(chr_len[::2].astype(str), chr_len[1::2].astype(int)))

            for signal_name, signal_dct in self.snp_signal_dct.items():
                _logger.info("JBrowse export: SNP signal {}".format(signal_name))
                for index, flag in enumerate(signal_dct['FLAG']):
                    if "offset" in signal_dct:
                        offset = signal_dct['offset'][index]
                        for bin_size in snp_bin:
                            signal = self.signal_name(bin_size, signal_name, flag)
                            bigwig_filename = "{}_offset{}.bw".format(signal, offset)
                            bigwig_file = self.export_data_dir_list[root_index].joinpath(bigwig_filename)
                            bigwig_file = str(bigwig_file)

                            self.create_bigwig(root_io, bigwig_file, chr_list, bin_size, signal_name, flag,
                                               offset=offset)

                    else:
                        for bin_size in snp_bin:
                            signal = self.signal_name(bin_size, signal_name, flag)
                            bigwig_filename = "{}.bw".format(signal)
                            bigwig_file = self.export_data_dir_list[root_index].joinpath(bigwig_filename)
                            bigwig_file = str(bigwig_file)

                            self.create_bigwig(root_io, bigwig_file, chr_list, bin_size, signal_name, flag)

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
        track_dct_list = []
        for root_index, root_io in enumerate(self.io):
            rd_chr, rd_bin = self.rd_chr_bin(root_io)
            url_template_dct = []
            for signal_name, signal_dct in self.rd_signal_dct.items():
                if 'FLAG' in signal_dct:
                    for index, flag in enumerate(signal_dct['FLAG']):
                        suffix_rd_flag = Signals.suffix_rd_flag(flag)
                        signal_id = "{}_{}{}".format(self.pytor_names[root_index], signal_name, suffix_rd_flag)
                        scales = {}
                        for bin_size in rd_bin:
                            signal = self.signal_name(bin_size, signal_name, flag)
                            bigwig_filename = "{}.bw".format(signal)
                            bigwig_file = self.export_data_dir_list[root_index].joinpath(bigwig_filename)
                            bigwig_current_path = Path(bigwig_file.parent.parent.name).joinpath(bigwig_file.parent.name, bigwig_file.name).as_posix()
                            if bigwig_file.exists():
                                scales[bin_size] = bigwig_current_path

                        if len(scales) > 0:
                            url_template_dct.append({
                                "storeClass": "MultiScaleBigWig/Store/SeqFeature/MultiScaleBigWig",
                                "scales": scales,
                                "name": signal_id,
                                "color": signal_dct['color'][index],

                            })
            if len(url_template_dct) > 0:

                track_dct = {
                    "category": self.pytor_names[root_index],
                    'autoscale': 'local',
                    "storeClass": "MultiBigWig/Store/SeqFeature/MultiBigWig",
                    "showTooltips": True,
                    "showLabels": True,
                    "clickTooltips": True,
                    "key": "RD",
                    "label": "RD {}".format(self.pytor_names[root_index]),
                    "type": "MultiBigWig/View/Track/MultiWiggle/MultiXYPlot",
                    'useStdDev': True,
                    'urlTemplates': url_template_dct

                }
                track_dct_list.append(track_dct)
        return track_dct_list

    def add_snp_config_track(self):
        _logger.debug("Get SNP config track info")
        track_dct_list = []
        for root_index, root_io in enumerate(self.io):
            snp_url_dct_list = []
            snp_chr, snp_bin = self.snp_chr_bin(root_io)
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
                            bigwig_file = self.export_data_dir_list[root_index].joinpath(bigwig_filename)
                            bigwig_current_path = Path(bigwig_file.parent.parent.name).joinpath(bigwig_file.parent.name, bigwig_file.name).as_posix()
                            if bigwig_file.exists():
                                scales[bin_size] = bigwig_current_path
                    else:
                        signal_id = "{}{}".format(signal_name, suffix_flag)
                        for bin_size in snp_bin:
                            signal = self.signal_name(bin_size, signal_name, flag)
                            bigwig_filename = "{}.bw".format(signal)
                            bigwig_file = self.export_data_dir_list[root_index].joinpath(bigwig_filename)
                            bigwig_current_path = Path(bigwig_file.parent.parent.name).joinpath(bigwig_file.parent.name, bigwig_file.name).as_posix()
                            if bigwig_file.exists():
                                scales[bin_size] = bigwig_current_path
                    if len(scales) > 0:
                        snp_url_dct_list.append({
                            "storeClass": "MultiScaleBigWig/Store/SeqFeature/MultiScaleBigWig",
                            "scales": scales,
                            "name": signal_id,
                            "color": signal_dct['color'][index],
                            "nonCont": signal_dct['nonCont'][index]
                        })
            if len(snp_url_dct_list) > 0:
                track_dct = {
                    "category": self.pytor_names[root_index],
                    'autoscale': 'local',
                    "storeClass": "MultiBigWig/Store/SeqFeature/MultiBigWig",
                    "showTooltips": True,
                    "showLabels": True,
                    "clickTooltips": True,
                    "max_score": 1,
                    "key": "SNP",
                    "label": "SNP {}".format(self.pytor_names[root_index]),
                    "type": "MultiBigWig/View/Track/MultiWiggle/MultiXYPlot",
                    'urlTemplates': snp_url_dct_list,
                }
                track_dct_list.append(track_dct)
        return track_dct_list

    def create_tracklist_json(self):
        _logger.debug("Creates config file: {}".format(self.export_tracklist_file))

        # reference config
        track_dct = self.add_config_reference()

        # create rd config
        rd_track_list = self.add_rd_config_track()
        for rd_track in rd_track_list:
            track_dct['tracks'].append(rd_track)

        # create SNP config
        snp_track_list = self.add_snp_config_track()
        for snp_track in snp_track_list:
            track_dct['tracks'].append(snp_track)

        with open(self.export_tracklist_file, 'w') as f:
            json.dump(track_dct, f, indent=2)
        return track_dct

    def create_reference_json(self):
        _logger.debug("Exporting reference details")
        # get signal details
        chr_len = list(np.array(self.io[0].get_signal(None, None, "chromosome lengths")).astype("str"))
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

        _logger.info("JBrowse export: complete")
        _logger.info("Copy this directory to jbrowse directory if export path is not set to JBrowse path, "
                     "To access this via localhost: http://localhost/jbrowse/?data={}"
                     .format(self.export_directory.parent.name))
