""" cnvpytor.viewparams

Class ViewParams: parameters and description for Viewer class
"""
from __future__ import absolute_import, print_function, division

from .genome import *
import matplotlib.pyplot as plt
import logging

_logger = logging.getLogger("cnvpytor.viewparams")


class ViewParams(object):
    default = {
        "bin_size": None,
        "panels": ["rd"],
        "rd_partition": True,
        "rd_call": True,
        "rd_call_mosaic": False,
        "rd_use_mask": False,
        "rd_use_gc_corr": True,
        "rd_manhattan_range": [0, 2],
        "rd_manhattan_call": False,
        "snp_use_mask": True,
        "snp_use_id": False,
        "snp_use_phase": False,
        "rd_colors": "",
        "snp_colors": ["yellow", "orange", "cyan", "blue", "lime", "green", "yellow", "orange"],
        "baf_colors": "",
        "lh_colors": "",
        "plot_files": [],
        "plot_file": 0,
        "chrom": [],
        "style": None,
        "grid": "auto",
        "xkcd": False,
        "output_filename": "",
        "palette1": ["#555555", "#aaaaaa"],
        "palette2": ["#00ff00", "#0000ff"],
        "contrast": 20,
        "min_segment_size": 0
    }

    command_tree = {
        "set": {},
        "unset": {},
        "help": {},
        "save": None,
        "show": None,
        "quit": None,
        "rd": None,
        "likelihood": None,
        "baf": None,
        "snp": None,
        "info": None,
        "stat": None,
        "rdstat": None,
        "circular": None,
        "manhattan": None,
        "calls": None,
        "ls": None,
        "compare": None
    }

    param_help = {
        "help": help_format(
            topic="help",
            p_desc="Print help for a topic, command or parameter.",
            p_usage="help <topic>\n\nList of params: " +
                    ", ".join(default.keys()) +
                    "\n\nList of commands: " +
                    ", ".join(command_tree.keys()) +
                    "\n\nList of topics: plotting, signals",
            p_example="help bin_size\nhelp plotting"
        ),
        "plotting": help_format(
            topic="PLOTTING",
            p_desc="There are several types of plots:\n" +
                   "    * genome wide plots: rdstat, stat, rd, baf, snp, likelihood\n" +
                   "    * regions: genomic regions separated by space (new subplot) and comma (same subplot)\n" +
                   "    * manhattan: manhattan style whole genome plot\n" +
                   "    * calls: plot calls\n"
                   "    * circular: circular plots of read-depth and MAF",
            p_example="chr1:1M-20M\nchr21\nchr1:1M-20M,chr2:30M-45M chr21",
            p_see="regions, rdstat, stat, rd, baf, snp, likelihood, manhattan, calls, circular"
        ),
        "signals": help_format(
            topic="SIGNALS",
            p_desc="...",
            p_see="plotting"
        ),
        "set": help_format(
            topic="set",
            p_desc="Set parameter value. " +
                   "Requires argument or list of arguments except in the case when parameter type is bool.",
            p_usage="set param [value]",
            p_example="set bin_size 100000\nset rd_call_mosaic\nset panels rd likelihood",
            p_see="unset"
        ),
        "unset": help_format(
            topic="unset",
            p_desc="Set bool parameter to False.",
            p_usage="unset param",
            p_example="unset rd_call_mosaic",
            p_see="unset"
        ),
        "save": help_format(
            topic="save",
            p_desc="Save current figure to file. There is a shortcut (see examples): " +
                   "redirection '>' after plotting command.\n" +
                   "Available formats:\n" +
                   key_val_str(plt.gcf().canvas.get_supported_filetypes())[:-1],
            p_usage="save filename",
            p_example="save image.png\n manhattan > image.png\n chr1:20M-50M > image.png",
            p_see="output_filename"
        ),
        "show": help_format(
            topic="show",
            p_desc="Lists all parameters with values",
            p_usage="show",
            p_see="set, unset"
        ),
        "quit": help_format(
            topic="quit",
            p_desc="Exists interactive mode.",
            p_usage="quit",
            p_see="help"
        ),
        "rd": "",
        "likelihood": "",
        "snp": "",
        "baf": "",
        "stat": "",
        "rdstat": "",
        "circular": "",
        "manhattan": "",
        "calls": "",
        "ls": help_format(
            topic="ls",
            p_desc="Print content of pytor files",
            p_usage="ls",
            p_see="show"
        ),
        "compare": "",
        "bin_size": help_format(
            topic="bin_size",
            p_desc="Size of bins used for plotting",
            p_type="binsize_type (int divisible by 100)",
            p_default=str(default["bin_size"]),
            p_affects="all",
            p_example="set bin_size 100000",
            p_see="output_filename, xkcd"
        ),
        "panels": help_format(
            topic="panels",
            p_desc="List of panels to plot. Possible options are:\n" +
                   "    rd         - read depth plot,\n" +
                   "    likelihood - baf likelihood 2D plot,\n" +
                   "    baf        - binned baf, maf and likelihood peak position,\n" +
                   "    snp        - plot baf for each particular SNP",
            p_type="list of strings",
            p_default=str(default["panels"]),
            p_affects="region plot, manhattan",
            p_example="set panels rd likelihood",
            p_see="grid"
        ),
        "rd_partition": help_format(
            topic="rd_partition",
            p_desc="Enables plotting partition signal in rd plots",
            p_type="bool",
            p_default=str(default["rd_partition"]),
            p_affects="region plot, rd",
            p_example="set rd_partition\nunset rd_partition",
            p_see="rd_call, rd_call_mosaic"
        ),
        "rd_call": help_format(
            topic="rd_call",
            p_desc="Enables plotting call signal in rd plots",
            p_type="bool",
            p_default=str(default["rd_call"]),
            p_affects="region plot, rd",
            p_example="set rd_call\nunset rd_call",
            p_see="rd_partition, rd_call_mosaic"
        ),
        "rd_call_mosaic": help_format(
            topic="rd_call_mosaic",
            p_desc="Enables plotting mosaic call signal in rd plots",
            p_type="bool",
            p_default=str(default["rd_call_mosaic"]),
            p_affects="region plot, rd",
            p_example="set rd_call_mosaic\nunset rd_call_mosaic",
            p_see="rd_partition, rd_call"
        ),
        "rd_raw": "",
        "rd_use_mask": "",
        "rd_use_gc_corr": "",
        "rd_range": "",
        "rd_manhattan_range": "",
        "rd_manhattan_call": "",
        "snp_use_mask": "",
        "snp_use_id": "",
        "snp_use_phase": "",
        "rd_colors": "",
        "snp_colors": "",
        "baf_colors": "",
        "lh_colors": "",
        "plot_files": "",
        "plot_file": "",
        "chrom": "",
        "style": "",
        "grid": "",
        "xkcd": "",
        "output_filename": "",
        "palette1": "",
        "palette2": "",
        "contrast": "",
        "min_segment_size": "",
    }

    def __init__(self, params):
        for key in self.default:
            if key in params:
                setattr(self, key, params[key])
            else:
                setattr(self, key, self.default[key])

    def help(self, param):
        if param in self.param_help:
            print(self.param_help[param])
        else:
            print("\nUnknown parameter !\n")

    def set(self, param, args):
        if param in self.params and self.params[param] is False:
            self.__setattr__(param, True)
        elif param == "bin_size" and len(args) > 0:
            self.__setattr__(param, args[0])
        elif param == "contrast" and len(args) > 0:
            self.__setattr__(param, float(args[0]))
        elif param == "grid" and len(args) > 0:
            if args[0] == "auto":
                self.__setattr__(param, "auto")
            else:
                self.__setattr__(param, list(map(int, args)))
        elif param == "manhattan_range" and len(args) > 0:
            self.__setattr__(param, list(map(float, args)))
        elif param == "min_segment_size" and len(args) > 0:
            self.__setattr__(param, int(args[0]))
        elif param == "output_filename" and len(args) > 0:
            self.__setattr__(param, args[0])
        elif param == "plot_file" and len(args) > 0:
            self.__setattr__(param, int(args[0]))
        elif param == "plot_files":
            self.__setattr__(param, list(map(int, args)))
        elif param == "style":
            self.__setattr__(param, args[0])
        elif param in self.params and self.params[param] is not True:
            self.__setattr__(param, args)

    def unset(self, param):
        if param in self.params and self.params[param] is True:
            self.__setattr__(param, False)

    @property
    def params(self):
        dct = {}
        for key in self.default:
            dct[key] = getattr(self, key)

        return dct

    @property
    def bin_size_f(self):
        """
        Formatted bin_size (e.g. 1000 -> "1K", 10000000 -> "10M")

        Returns
        -------
        bin_size_f : str

        """
        return binsize_format(self.params["bin_size"])

    def __setattr__(self, name, value):

        if name == 'bin_size' and value is not None:
            try:
                value = binsize_type(value)
            except (ArgumentTypeError, ValueError):
                _logger.warning("bin_size should be integer divisible by 100")
        if name == 'style':
            if value in plt.style.available:
                plt.style.use("default")
                plt.style.use(value)

        if name == 'xkcd':
            if value:
                self.style = "classic"
                from matplotlib import patheffects
                plt.xkcd()
                plt.rcParams["path.effects"] = [patheffects.withStroke(linewidth=0.5, foreground="w")]

            elif hasattr(self, name) and self.xkcd:
                plt.rcdefaults()
                self.style = 'classic'

        super(ViewParams, self).__setattr__(name, value)
