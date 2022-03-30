""" cnvpytor.viewer

Class Viewer: ploting CNVpytor data
"""
from __future__ import absolute_import, print_function, division

from .io import *
from .utils import *
from .genome import *
from .viewparams import ViewParams, HelpDescription
from .annotator import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from scipy.cluster import hierarchy
from scipy.stats import beta

import numpy as np
import logging
import readline
import traceback
import os
import sys
import datetime

_logger = logging.getLogger("cnvpytor.viewer")


class Reader:
    def __init__(self, files):
        """ Class constructor opens cnvpytor files.

        Parameters
        ----------
        files : list of str
            List of cnvpytor filenames.

        """
        self.io = [IO(f, ro=True) for f in files]


class Show(Reader):
    def ls(self):
        """ Prints to stdout content of all cnvpytor files.

        """
        for i in self.io:
            i.ls()

    def gc_info(self):
        """ Prints to stdout gc content info of all cnvpytor files.

        """
        for i in self.io:
            i.gc_info(stdout=True)

    def meta(self):
        """ Prints to stdout meta tags of all cnvpytor files.

        """
        for i in self.io:
            i.print_meta_attribute()

    def info(self, bin_sizes):
        """ Prints to stdout RD info for all cnvpytor files.
        Columns are following:
            filename
            mean read length, stdev of read length in %
            mean template length, stdev of template length in %
            for each bin_size (including 100 always):
                rd level and corresponding stdev for each chromosome type (autosomes, sex chromosomes and mitochondria)

        """
        if 100 not in bin_sizes:
            bin_sizes = [100] + bin_sizes
        labels = ["FILE", "RL", "dRL[%]", "FL", "dFL[%]"]
        for bs in bin_sizes:
            labels.append("RD_AUTO_" + binsize_format(bs))
            labels.append("dRD_AUTO_" + binsize_format(bs) + "[%]")
            labels.append("RD_GC_AUTO_" + binsize_format(bs))
            labels.append("dRD_GC_AUTO_" + binsize_format(bs) + "[%]")
            labels.append("RD_XY_" + binsize_format(bs))
            labels.append("dRD_XY_" + binsize_format(bs) + "[%]")
            labels.append("RD_GC_XY_" + binsize_format(bs))
            labels.append("dRD_GC_XY_" + binsize_format(bs) + "[%]")
            if bs <= 500:
                labels.append("RD_MT_" + binsize_format(bs))
                labels.append("dRD_MT_" + binsize_format(bs) + "[%]")
                labels.append("RD_GC_MT_" + binsize_format(bs))
                labels.append("dRD_CG_MT_" + binsize_format(bs) + "[%]")
        print(("{:}\t{:}\t{:}\t{:}\t{:}\t" + "{:}\t" * (len(labels) - 5)).format(*tuple(labels)))
        for i in self.io:
            rfd = i.get_signal(None, None, "read frg dist")
            rd = np.sum(rfd, axis=1)
            fd = np.sum(rfd, axis=0)
            mrl = np.sum(rd * np.arange(rd.size)) / np.sum(rd)
            mfl = np.sum(fd * np.arange(fd.size)) / np.sum(fd)
            mrl2 = np.sum(rd * np.arange(rd.size) * np.arange(rd.size)) / np.sum(rd)
            mfl2 = np.sum(fd * np.arange(fd.size) * np.arange(fd.size)) / np.sum(fd)
            sdr = 100. * np.sqrt(mrl2 - mrl * mrl) / mrl
            sdf = 100. * np.sqrt(mfl2 - mfl * mfl) / mfl
            print("{:}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t".format(i.filename, mrl, sdr, mfl, sdf), end="")
            for bs in bin_sizes:
                for flag in [FLAG_AUTO, FLAG_SEX, FLAG_MT]:
                    if bs <= 500 or not flag == FLAG_MT:
                        if i.signal_exists(None, bs, "RD stat", flags=flag):
                            stat = i.get_signal(None, bs, "RD stat", flags=flag)
                            if stat[4] > 0:
                                stat[5] /= stat[4] / 100.
                            print("{:.2f}\t{:.2f}\t".format(stat[4], stat[5]), end="")
                        else:
                            print("{:}\t{:}\t".format("-", "-"), end="")
                        if i.signal_exists(None, bs, "RD stat", flags=(flag | FLAG_GC_CORR)):
                            stat = i.get_signal(None, bs, "RD stat", flags=(flag | FLAG_GC_CORR))
                            if stat[4] > 0:
                                stat[5] /= stat[4] / 100.
                            print("{:.2f}\t{:.2f}\t".format(stat[4], stat[5]), end="")
                        else:
                            print("{:}\t{:}\t".format("-", "-"), end="")
            print()


class Figure(ViewParams):
    def __init__(self, params, force_agg=False):
        """ Class implements visualisations and exports

        Parameters
        ----------
        params : dict
            Params to be passed to ViewParam class

        """
        if force_agg:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        ViewParams.__init__(self, params)
        self.fig = None
        self.fig_grid = None
        self.fig_sub_grid = None
        self.count = 0
        self.current = -1
        self.sg_current = -1

    def new_figure(self, panel_count, grid="auto", panel_size=None, title=None):
        """ Clear figure and create new plot layout.

        Parameters
        ----------
        panel_count : int
            Number of panels
        grid : str or (int, int)
            number of columns and rows (sx, sy) or "auto"
        panel_size : (float, float)
            size of a single panel (only when plots in file)

        """
        if panel_size is None:
            panel_size = self.panel_size
        if grid == "auto":
            grid = self.grid
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, dpi=self.dpi, facecolor='w', edgecolor='k')
        if title is not None:
            self.fig.suptitle(title, fontsize=16)
        sx, sy = self._get_grid(grid, panel_count)
        if self.output_filename != "":
            self.fig.set_figheight(panel_size[1] * sy)
            self.fig.set_figwidth(panel_size[0] * sx)
        self.fig_grid = gridspec.GridSpec(sy, sx, hspace=self.margins[5], wspace=self.margins[4])
        self.current = -1
        self.sg_current = -1

    def new_subgrid(self, panel_count, grid="auto", hspace=0, wspace=0):
        if grid == "auto":
            grid = self.subgrid
        sx, sy = self._get_grid(grid, panel_count)
        self.current += 1
        self.fig_sub_grid = gridspec.GridSpecFromSubplotSpec(sy, sx, subplot_spec=self.fig_grid[self.current],
                                                             wspace=wspace, hspace=hspace)
        self.sg_current = -1
        self.sg_current_ax = None

    def next_panel(self):
        """ Return axes of next panel

        Returns
        -------
        ax : matplotlib.axes.Axes
            Axes for a given panel
        """
        self.current += 1
        return self.fig.add_subplot(self.fig_grid[self.current])

    def next_subpanel(self, sharex=False):
        """ Return axes of next sub panel

        Returns
        -------
        ax : matplotlib.axes.Axes
            Axes for a given panel
        """
        self.sg_current += 1
        if self.sg_current == 0 or not sharex:
            self.sg_current_ax = self.fig.add_subplot(self.fig_sub_grid[self.sg_current])
        else:
            self.sg_current_ax = self.fig.add_subplot(self.fig_sub_grid[self.sg_current], sharex=self.sg_current_ax)
        return self.sg_current_ax

    def next_polar_panel(self):
        """ Return axes of next panel

        Returns
        -------
        ax : matplotlib.axes.Axes
            Axes for a given panel
        """
        self.current += 1
        return self.fig.add_subplot(self.fig_grid[self.current], projection="polar")

    def get_panel(self, i):
        """ Returns axes of a i-th panel

        Parameters
        ----------
        i : int
            Panel number

        Returns
        -------
        ax : matplotlib.axes.Axes
            Axes for a given panel

        """
        return self.fig.get_axes()[i]

    def _get_grid(self, grid, panel_count):
        if grid == "auto":
            sx, sy = self._panels_shape(panel_count)
        elif grid == "vertical":
            sx, sy = 1, panel_count
        elif grid == "horizontal":
            sx, sy = panel_count, 1
        else:
            sx, sy = tuple(grid)
        return sx, sy

    def fig_show(self, add_sufix=True, suffix=""):
        """ Plot figure. If output_filename is specified it will plot only into a file.

        Parameters
        ----------
        add_sufix : bool
            If true it will add sufix to output_filename in format prefix.sufix.count.extension
            where count is auto-incremented integer starting from 0 and
            prefix.extension is parsed from output_filename parameter.

        suffix : str
            Sufix used in filename.

        """
        bottom, top, left, right, wspace, hspace = self.margins
        plt.subplots_adjust(bottom=bottom, top=top, wspace=wspace, hspace=hspace, left=left, right=right)
        if self.output_filename != "":
            image_filename = self.output_filename
            if add_sufix:
                image_filename = self._image_filename(suffix)
            if image_filename is not None:
                try:
                    plt.savefig(image_filename, dpi=self.dpi)
                except:
                    _logger.warning("Figure is not saved due to an error!")
                plt.close(self.fig)
            else:
                _logger.warning("Figure is not saved!")
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def _image_filename(self, suffix):
        parts = self.output_filename.split(".")
        if parts[-1] not in ["png", "pdf", "jpg", "eps", "svg"]:
            _logger.warning("File extension should be: .jpg, .png, .svg, .eps or .pdf")
            return None
        if suffix == "":
            suffix = str(self.count).zfill(4)
        else:
            suffix += "." + str(self.count).zfill(4)
        self.count += 1
        parts[-1] = suffix + "." + parts[-1]
        return ".".join(parts)

    @staticmethod
    def _panels_shape(n):
        sx, sy = 1, 1
        if n == 2:
            sy = 2
        elif n in [3, 4]:
            sx, sy = 2, 2
        elif n in [5, 6]:
            sx, sy = 2, 3
        elif n in [7, 8, 9]:
            sx, sy = 3, 3
        elif n in [10, 11, 12]:
            sx, sy = 3, 4
        elif n in [13, 14, 15, 16]:
            sx, sy = 4, 4
        elif n in [17, 18, 19, 20]:
            sx, sy = 4, 5
        elif n in [21, 22, 23, 24]:
            sx, sy = 4, 6
        else:
            while sx * sy < n:
                sy += 1
                sx = int(2. * sy / 3 + 1.)
        return sx, sy


class Viewer(Show, Figure, HelpDescription):

    def __init__(self, files, params={}, force_agg=False, history_file_size=1000):
        """

        Parameters
        ----------
        files : list of str
            List of cnvpytor filenames
        params : dict
            List of parameters different than default to be passed to ViewParams class.

        """
        _logger.debug("Viewer class init: files [%s], params %s." % (", ".join(files), str(params)))
        Figure.__init__(self, params, force_agg=force_agg)
        Show.__init__(self, files)
        self.history_file_size = history_file_size
        self.cnvpytor_dir = os.path.expanduser('~/.cnvpytor')
        self.save_history = False
        if os.path.exists(self.cnvpytor_dir):
            if os.access(self.cnvpytor_dir, os.W_OK):
                self.save_history = True
            if os.path.exists(self.cnvpytor_dir + "/viewer.conf"):
                conf = eval(open(self.cnvpytor_dir + "/viewer.conf").read())
                for key in conf:
                    setattr(self, key, conf[key])

        self.io_gc = self.io[0]
        self.io_mask = self.io[0]
        self.reference_genome = None
        self.plot_files = list(range(len(files)))
        self.default["plot_files"] = list(range(len(files)))
        if self.io[0].signal_exists(None, None, "reference genome"):
            rg_name = np.array(self.io[0].get_signal(None, None, "reference genome")).astype("str")[0]
            Genome.detected_genome = rg_name
            self.reference_genome = Genome.reference_genomes[rg_name]
            if "mask_file" in Genome.reference_genomes[rg_name]:
                self.io_mask = IO(Genome.reference_genomes[rg_name]["mask_file"], ro=True, buffer=True)
            if "gc_file" in Genome.reference_genomes[rg_name]:
                self.io_gc = IO(Genome.reference_genomes[rg_name]["gc_file"], ro=True, buffer=True)

    def parse(self, command):
        current = "regions"
        regions = []

        for p in command:
            if p.isdigit() and (int(p) % 100) == 0:
                self.bin_size = int(p)
                if current == "rd":
                    self.rd()
                if current == "baf":
                    self.baf()
                if current == "likelihood":
                    self.likelihood()
                elif current == "manhattan":
                    self.global_plot()
                elif current == "calls":
                    if len(self.callers) > 0:
                        self.manhattan(plot_type=self.callers[0])
                elif current == "stat":
                    self.stat(int(p))
                elif current == "circular":
                    self.circular()
                elif current == "regions":
                    self.multiple_regions(regions)
                    regions = []
            elif p == "rdstat":
                self.stat()
            elif p == "snp":
                self.snp()
            elif p in ["rd", "baf", "manhattan", "calls", "stat", "regions", "likelihood", "circular"]:
                current = p
            elif current == "regions":
                regions.append(p)
            else:
                current = p

    def plot_command(self, command):
        self.interactive = False
        self.parse(command)

    def prompt(self):
        self.interactive = True

        chromosomes = set({})
        for f in self.io:
            chromosomes = chromosomes.union(set(f.rd_chromosomes()))
            chromosomes = chromosomes.union(set(f.snp_chromosomes()))
        for c in chromosomes:
            self.command_tree[c] = None
        self.command_tree["set"]["style"] = dict(zip(plt.style.available, [None] * len(plt.style.available)))
        if os.path.exists(self.cnvpytor_dir + "/history"):
            readline.read_history_file(self.cnvpytor_dir + "/history")

        readline.parse_and_bind("tab: complete")
        completer = PromptCompleter(self.command_tree)
        readline.set_completer(completer.complete)
        quit = False
        try:
            while not quit:
                prompt_str = ""
                if os.isatty(sys.stdin.fileno()):
                    prompt_str = "cnvpytor> "
                else:
                    self.interactive = False
                try:
                    line = raw_input(prompt_str)
                except NameError:
                    line = input(prompt_str)

                if line == "" or line[0] == "#":
                    continue

                if self.save_history and self.interactive:
                    readline.set_history_length(self.history_file_size)
                    readline.write_history_file(self.cnvpytor_dir + "/history")

                pre = line.split(">")
                f = pre[0].strip().split(" ")
                n = len(f)
                if len(line) == 0:
                    continue
                elif f[0] == "quit" or f[0] == "exit":
                    quit = True
                elif line[0] == "|":
                    try:
                        eval(compile(line[1:], '<string>', 'single'))
                    except Exception as e:
                        print(traceback.format_exc())
                elif f[0] == "save":
                    if n > 1:
                        try:
                            plt.savefig(f[1])
                        except ValueError:
                            _logger.warning("File extension should be: .jpg, .png, .svg, .eps or .pdf")
                        except:
                            _logger.warning("Figure is not saved due to an error!")

                elif f[0] in ["draw", "repaint", "update"]:
                    if n == 1:
                        self.fig.canvas.draw()
                elif f[0] == "ls":
                    self.ls()
                elif f[0] == "meta":
                    self.meta()
                elif f[0] == "show":
                    if n == 1:
                        self.show()
                elif f[0] == "set":
                    if n > 1:
                        self.set_param(f[1], f[2:])
                elif f[0] == "help" and n > 1:
                    self.help(f[1])
                elif f[0] == "help" and n == 1:
                    self.help("help")
                elif f[0] == "unset":
                    if n > 1:
                        self.unset(f[1])
                elif f[0] == "genotype":
                    if n > 1:
                        self.genotype_all([self.bin_size], f[1:], interactive=True)
                elif f[0] == "snv":
                    if n == 2:
                        self.snp(callset=f[1])
                    elif n == 1:
                        self.snp(callset="default")
                elif f[0] == "compare":
                    if n == 3:
                        self.compare(f[1], f[2], plot=self.plot)
                    elif n == 4:
                        self.compare(f[1], f[2], n_bins=int(f[3]), plot=self.plot)
                elif f[0] == "info":
                    if n > 1:
                        self.info(list(map(binsize_type, f[1:])))
                elif f[0] == "print":
                    if f[1] == "calls":
                        if self.print_filename == "":
                            self.print_calls()
                        else:
                            self.print_calls_file()
                    elif f[1] == "joint_calls" or f[1] == "merged_calls":
                        self.print_simple_merged_calls()

                else:
                    try:
                        if f[0] not in ["rdstat", "snp"]:
                            self.parse(f + [str(self.bin_size)])
                        else:
                            self.parse(f)
                        if len(pre) > 1:
                            fns = pre[1].strip().split(" ")
                            if fns[0] != "":
                                plt.savefig(fns[0], dpi=200)
                    except Exception as e:
                        print(traceback.format_exc())
        except (EOFError, KeyboardInterrupt):
            print()
            return

    def help(self, param):
        if param in self.param_help:
            print(self.param_help[param])
        else:
            print("\nUnknown parameter !\n")

    @staticmethod
    def set_style(style):
        if style in plt.style.available:
            plt.style.use("default")
            plt.style.use(style)

    def file_title(self, ix):
        if ix < len(self.file_titles):
            return self.file_titles[ix]
        else:
            return self.io[ix].filename.split("/")[-1].replace(".pytor", "")

    def show(self):
        print("\nParameters")
        for key in sorted(self.params.keys()):
            print("    * %s: %s" % (key, str(self.params[key])))
            if key == "plot_files":
                for i in range(len(self.io)):
                    print("            %d: %s" % (i, self.io[i].filename))
        print()

    def stat(self, his_bin_size=100, return_image=False):
        plt.clf()
        auto = self.io[self.plot_file].signal_exists(None, his_bin_size, "RD stat", FLAG_AUTO)
        sex = self.io[self.plot_file].signal_exists(None, his_bin_size, "RD stat", FLAG_SEX)
        mt = self.io[self.plot_file].signal_exists(None, his_bin_size, "RD stat", FLAG_MT) and (his_bin_size < 1001)
        if not (auto or sex or mt):
            return
        cond = [auto, sex, mt]
        stat_list = []
        n_cols = sum(map(int, cond))
        ix = 1
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(4 * n_cols, 8), dpi=90, facecolor='w', edgecolor='k')
        for t, c, flag in zip(["Autosomes", "X/Y", "Mitochondria"], cond, [FLAG_AUTO, FLAG_SEX, FLAG_MT]):
            if c:
                stat = self.io[self.plot_file].get_signal(None, his_bin_size, "RD stat", flag)
                stat_list.append(stat)
                max_rd = int(stat[0])
                bin_size = int(stat[1])
                n_bins = int(stat[2])
                lim_rd = int(max(2 * stat[4], stat[4] + 3 * stat[5]))
                _logger.info("RD stat for %s: %.2f +- %.2f" % (t, stat[4], stat[5]))
                if t == "Mitochondria" and auto:
                    _logger.info("RD stat for %s - number of mitochondria per cell: %.2f +- %.2f" % (
                        t, 2 * stat[4] / stat_list[0][4],
                        2 * stat[5] / stat_list[0][4] + stat_list[0][5] * stat[4] / (
                                stat_list[0][4] * stat_list[0][4])))
                his_p = self.io[self.plot_file].get_signal(None, his_bin_size, "RD p dist", flag)
                his_u = self.io[self.plot_file].get_signal(None, his_bin_size, "RD u dist", flag)
                his_rd_gc = self.io[self.plot_file].get_signal(None, his_bin_size, "RD GC dist", flag)
                gc_corr = self.io[self.plot_file].get_signal(None, his_bin_size, "GC corr", flag)
                ax = plt.subplot(2, n_cols, ix)
                ax.set_xlabel("RD")
                ax.set_ylabel("GC [%]")
                ax.xaxis.set_ticklabels([])
                ax.set_title(t)
                his_rd_gc[0][0] = 0
                ax.imshow(his_rd_gc[:lim_rd // bin_size, :].T, aspect="auto", interpolation='nearest', origin='lower')
                ax.plot(gc_corr * stat[4] / bin_size, range(101), "w-")

                ax = plt.subplot(2, n_cols, ix + n_cols)
                ax.set_ylabel("Normalised distribution")
                ax.set_xlabel("RD")
                ax.set_xlim([0, lim_rd])
                # ax.set_ylim([0, 1.1])
                bins = range(0, max_rd, bin_size)
                x = np.arange(0, max_rd // bin_size * bin_size, 0.1 * bin_size)
                plt.plot(x, normal(x, 1, stat[4], stat[5]), "g-")
                x = np.array(bins)
                plt.plot(x[:len(his_u)], his_u / stat[3], "y*")
                plt.plot(x[:len(his_p)], his_p / stat[3], "b*")
                ix += 1
        plt.subplots_adjust(bottom=0.08, top=0.95, wspace=0.25, hspace=0, left=0.05 * 3 / n_cols, right=0.95)
        if return_image:
            self.fig.canvas.draw()
            import PIL
            pil_image = PIL.Image.frombytes('RGB', self.fig.canvas.get_width_height(),
                                            self.fig.canvas.tostring_rgb())
            return pil_image
        elif self.output_filename != "":
            plt.savefig(self._image_filename("stat"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def rd(self):
        bin_size = self.bin_size
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_file].rd_chromosome_name(c)
            if self.io[self.plot_file].signal_exists(rd_chr, bin_size, "RD", 0) and \
                    self.io[self.plot_file].signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                chroms.append((rd_chr, l))
        self.new_figure(panel_count=len(chroms))
        for c, l in chroms:
            flag_rd = FLAG_USEMASK if self.rd_use_mask else 0
            mean, stdev = self.io[self.plot_file].rd_normal_level(bin_size, flag_rd | FLAG_GC_CORR)
            his_p = self.io[self.plot_file].get_signal(c, bin_size, "RD", flag_rd)
            his_p_corr = self.io[self.plot_file].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            his_p_seg = self.io[self.plot_file].get_signal(c, bin_size, "RD partition", flag_rd | FLAG_GC_CORR)
            his_p_call = self.io[self.plot_file].get_signal(c, bin_size, "RD call", flag_rd | FLAG_GC_CORR)
            his_p_mosaic_seg = self.io[self.plot_file].get_signal(c, bin_size, "RD mosaic segments",
                                                                  flag_rd | FLAG_GC_CORR)
            his_p_mosaic_seg = segments_decode(his_p_mosaic_seg)
            his_p_mosaic_call = self.io[self.plot_file].get_signal(c, bin_size, "RD mosaic call",
                                                                   flag_rd | FLAG_GC_CORR)
            his_p_mosaic_seg_2d = self.io[self.plot_file].get_signal(c, bin_size, "RD mosaic segments 2d",
                                                                     flag_rd | FLAG_GC_CORR)
            his_p_mosaic_seg_2d = segments_decode(his_p_mosaic_seg_2d)
            his_p_mosaic_call_2d = self.io[self.plot_file].get_signal(c, bin_size, "RD mosaic call 2d",
                                                                      flag_rd | FLAG_GC_CORR)
            his_p_mosaic = np.zeros_like(his_p) * np.nan
            if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call and (
                    "rd_mosaic" in self.callers):
                for seg, lev in zip(list(his_p_mosaic_seg), list(his_p_mosaic_call[0])):
                    for segi in seg:
                        his_p_mosaic[segi] = lev
            his_p_mosaic_2d = np.zeros_like(his_p) * np.nan
            if his_p_mosaic_call_2d is not None and len(his_p_mosaic_call_2d) > 0 and self.rd_call and (
                    "combined_mosaic" in self.callers):
                for seg, lev in zip(list(his_p_mosaic_seg_2d), list(his_p_mosaic_call_2d[0])):
                    for segi in seg:
                        his_p_mosaic_2d[segi] = lev
            ax = self.next_panel()
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // bin_size, 10e6 // bin_size), minor=[])
            if (self.rd_range[1] - self.rd_range[0]) < 30:
                ax.yaxis.set_ticks(np.arange(int(self.rd_range[0]), int(self.rd_range[1] + 1), 1) * mean / 2,
                                   minor=[])
            ax.set_ylim([self.rd_range[0] * mean / 2, self.rd_range[1] * mean / 2])
            n_bins = l // bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()
            if self.rd_raw:
                plt.step(his_p, self.rd_colors[0])
            if self.rd_corrected:
                plt.step(his_p_corr, self.rd_colors[1])
            if his_p_seg is not None and len(his_p_seg) > 0 and self.rd_partition:
                plt.step(his_p_seg, self.rd_colors[2])
            if his_p_call is not None and len(his_p_call) > 0 and self.rd_call:
                plt.step(his_p_call, self.rd_colors[3])
            if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call and (
                    "rd_mosaic" in self.callers):
                plt.step(his_p_mosaic, self.rd_colors[4])
            if his_p_mosaic_call_2d is not None and len(his_p_mosaic_call_2d) > 0 and self.rd_call and (
                    "combined_mosaic" in self.callers):
                plt.step(his_p_mosaic_2d, self.rd_colors[5])
        self.fig_show(suffix="rd")

    def rd_diff(self, file1, file2):
        bin_size = self.bin_size
        flag_rd = (FLAG_USEMASK if self.rd_use_mask else 0)
        mean, stdev = self.io[file1].rd_normal_level(bin_size, flag_rd | FLAG_GC_CORR)
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[file1].rd_chromosome_name(c)
            if self.io[file1].signal_exists(rd_chr, bin_size, "RD", 0) and \
                    self.io[file1].signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)) and (
                    len(self.chrom) == 0 or rd_chr in self.chrom):
                chroms.append((rd_chr, l))
        self.new_figure(panel_count=len(chroms))
        for c, l in chroms:
            flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
            stat1 = self.io[file1].get_signal(None, bin_size, "RD stat", flag)
            stat2 = self.io[file2].get_signal(None, bin_size, "RD stat", flag)
            if stat1 is None:
                _logger.error(
                    "Data for bin size %d is missing in file '%s'!" % (bin_size, self.io[file1].filename))
                return
            if stat2 is None:
                _logger.error(
                    "Data for bin size %d is missing in file '%s'!" % (bin_size, self.io[file2].filename))
                return

            his_p_corr1 = self.io[file1].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            his_p_corr2 = self.io[file2].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            ax = self.next_panel()
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            if (self.rd_range[1] - self.rd_range[0]) < 30:
                ax.yaxis.set_ticks(np.arange(int(self.rd_range[0]), int(self.rd_range[1] + 1), 1) * mean / 2,
                                   minor=[])
            ax.yaxis.set_ticks(np.arange(0, 2, 0.25), minor=[])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // bin_size, 10e6 // bin_size), minor=[])
            ax.set_ylim([0, 1])
            n_bins = l // bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()
            if self.markersize == "auto":
                plt.scatter(np.arange(len(his_p_corr1)), his_p_corr1 / stat1[4] - his_p_corr2 / stat2[4], marker='.',
                            s=10, alpha=0.7)
            else:
                plt.scatter(np.arange(len(his_p_corr1)), his_p_corr1 / stat1[4] - his_p_corr2 / stat2[4], marker='.',
                            s=self.markersize, alpha=0.7)
            # plt.step(np.abs(his_p_corr1 / stat1[4] - his_p_corr2 / stat2[4]), "k")
        self.fig_show(suffix="rd_diff")

    def likelihood(self):
        bin_size = self.bin_size
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
            FLAG_USEHAP if self.snp_use_phase else 0)
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        chroms = []
        if self.reference_genome is None:
            chroms = self.io[self.plot_file].snp_chromosomes()
        else:
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = self.io[self.plot_file].snp_chromosome_name(c)
                if self.io[self.plot_file].signal_exists(snp_chr, bin_size, "SNP likelihood", snp_flag) and (
                        Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                    chroms.append(snp_chr)
        self.new_figure(panel_count=len(chroms))
        for c in chroms:
            likelihood = self.io[self.plot_file].get_signal(c, bin_size, "SNP likelihood", snp_flag)
            img = np.array(likelihood).transpose()
            ax = self.next_panel()
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.imshow(img, aspect='auto')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticks(np.arange(0, likelihood.shape[0], 50), minor=[])
            ax.set_xlim([0, likelihood.shape[0]])
            if self.snp_call and ("baf_mosaic" in self.callers):
                likelihood = self.io[self.plot_file].get_signal(c, bin_size, "SNP likelihood call", snp_flag)
                segments = segments_decode(
                    self.io[self.plot_file].get_signal(c, bin_size, "SNP likelihood segments", snp_flag))
                call_pos = []
                call_i1 = []
                call_i2 = []
                call_c = []
                for s, lh in zip(segments, likelihood):
                    i1, i2, p = likelihood_pixels_pval(lh)
                    if i1 != i2 and len(s) > self.min_segment_size:
                        alpha = -np.log(p + 1e-40) / self.contrast
                        if alpha > 1:
                            alpha = 1
                        for pos in s:
                            call_pos.append(pos)
                            call_i1.append(min(i1, i2))
                            call_i2.append(max(i1, i2))
                            color = colors.to_rgb(self.lh_colors[0]) + (alpha,)
                            call_c.append(color)
                plt.scatter(call_pos, call_i1, s=self.lh_markersize, color=np.array(call_c), edgecolors='face',
                            marker=self.lh_marker)
                plt.scatter(call_pos, call_i2, s=self.lh_markersize, color=np.array(call_c), edgecolors='face',
                            marker=self.lh_marker)
            if self.snp_call and ("combined_mosaic" in self.callers):
                likelihood = self.io[self.plot_file].get_signal(c, bin_size, "SNP likelihood call 2d", snp_flag)
                segments = segments_decode(
                    self.io[self.plot_file].get_signal(c, bin_size, "SNP likelihood segments 2d", snp_flag))
                call_pos = []
                call_i1 = []
                call_i2 = []
                call_c = []
                for s, lh in zip(segments, likelihood):
                    i1, i2, p = likelihood_pixels_pval(lh)
                    if i1 != i2 and len(s) > self.min_segment_size:
                        alpha = -np.log(p + 1e-40) / self.contrast
                        if alpha > 1:
                            alpha = 1
                        for pos in s:
                            call_pos.append(pos)
                            call_i1.append(min(i1, i2))
                            call_i2.append(max(i1, i2))
                            color = colors.to_rgb(self.lh_colors[1]) + (alpha,)
                            call_c.append(color)
                plt.scatter(call_pos, call_i1, s=self.lh_markersize, color=np.array(call_c), edgecolors='face',
                            marker=self.lh_marker)
                plt.scatter(call_pos, call_i2, s=self.lh_markersize, color=np.array(call_c), edgecolors='face',
                            marker=self.lh_marker)
        self.fig_show(suffix="likelihood")

    def baf(self):
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
            FLAG_USEHAP if self.snp_use_phase else 0)
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            snp_chr = self.io[self.plot_file].snp_chromosome_name(c)
            if self.io[self.plot_file].signal_exists(snp_chr, self.bin_size, "SNP baf", snp_flag) and \
                    self.io[self.plot_file].signal_exists(snp_chr, self.bin_size, "SNP maf", snp_flag) and \
                    self.io[self.plot_file].signal_exists(snp_chr, self.bin_size, "SNP i1", snp_flag) and \
                    self.io[self.plot_file].signal_exists(snp_chr, self.bin_size, "SNP i2", snp_flag) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                chroms.append((snp_chr, l))

        self.new_figure(panel_count=len(chroms))
        for c, l in chroms:
            baf = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP baf", snp_flag)
            maf = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP maf", snp_flag)
            i1 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP i1", snp_flag)
            i2 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP i2", snp_flag)

            ax = self.next_panel()
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], minor=[])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // self.bin_size, 10e6 // self.bin_size), minor=[])
            ax.set_ylim([0, 1])
            n_bins = l // self.bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()
            ax.step(baf, self.baf_colors[0])
            ax.step(maf, self.baf_colors[1])
            ax.step(i1, self.baf_colors[2])
        self.fig_show(suffix="baf")

    def snp(self, plot_gt=None, plot_pmask=None, callset=None):
        if plot_pmask is None:
            plot_pmask = [0, 1]
        if plot_gt is None:
            plot_gt = [0, 1, 2, 3]
        chroms = []
        if self.reference_genome is None:
            chroms = self.io[self.plot_file].snp_chromosomes()
        else:
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = self.io[self.plot_file].snp_chromosome_name(c)
                if callset is None:
                    if self.io[self.plot_file].signal_exists(snp_chr, None, "SNP pos", 0) and \
                            self.io[self.plot_file].signal_exists(snp_chr, None, "SNP desc", 0) and \
                            self.io[self.plot_file].signal_exists(snp_chr, None, "SNP counts", 0) and \
                            self.io[self.plot_file].signal_exists(snp_chr, None, "SNP qual", 0) and \
                            (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                        chroms.append(snp_chr)
                else:
                    if self.io[self.plot_file].signal_exists(snp_chr, None, "somatic SNP pos", 0, name=callset) and \
                            self.io[self.plot_file].signal_exists(snp_chr, None, "somatic SNP desc", 0,
                                                                  name=callset) and \
                            self.io[self.plot_file].signal_exists(snp_chr, None, "somatic SNP counts", 0,
                                                                  name=callset) and \
                            self.io[self.plot_file].signal_exists(snp_chr, None, "somatic SNP qual", 0,
                                                                  name=callset) and \
                            (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                        chroms.append(snp_chr)
        self.new_figure(panel_count=len(chroms))
        for c in chroms:
            pos, ref, alt, nref, nalt, gt, flag, qual = self.io[self.plot_file].read_snp(c, callset=callset)
            hpos = []
            baf = []
            color = []
            qlpha = 0.7
            for i in range(len(pos)):
                if (nref[i] + nalt[i]) != 0:
                    if (gt[i] % 4 in plot_gt) and ((flag[i] >> 1) in plot_pmask):
                        hpos.append(pos[i])
                        if gt[i] % 4 != 2:
                            baf.append(1.0 * nalt[i] / (nref[i] + nalt[i]))
                        else:
                            baf.append(1.0 * nref[i] / (nref[i] + nalt[i]))
                        if self.snp_alpha_P:
                            alpha = None
                            color.append(colors.to_rgba(self.snp_colors[(gt[i] % 4) * 2 + 1], (flag[i] >> 1) * 0.4))
                        else:
                            color.append(self.snp_colors[(gt[i] % 4) * 2 + (flag[i] >> 1)])

            ax = self.next_panel()
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], minor=[])
            # l = max(pos)
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6), 10e6), minor=[])
            ax.set_ylim([0., 1.])
            ax.set_xlim([-0.05 * l, 1.05 * l])
            ax.grid()
            if self.markersize == "auto":
                ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=0.7)
            else:
                ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=0.7)
        self.fig_show(suffix="snp")

    def get_calls(self):
        bin_size = self.bin_size
        n = len(self.plot_files)
        ix = self.plot_files
        if self.annotate:
            annotator = Annotator(self.reference_genome)
        ret = []
        for caller in self.callers:
            if caller == "rd_mean_shift":
                for i in range(n):
                    io = self.io[ix[i]]
                    chroms = io.rd_chromosomes()
                    for c in chroms:
                        if (c in self.chrom) or len(self.chrom) == 0:
                            flag = (FLAG_USEMASK if self.rd_use_mask else 0) | \
                                   (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
                            if io.signal_exists(c, bin_size, "calls", flag):
                                calls = io.read_calls(c, bin_size, "calls", flag)
                                for call in calls:
                                    if in_interval(call["size"], self.size_range) \
                                            and in_interval(call["p_val"], self.p_range) \
                                            and in_interval(call["pN"], self.pN_range) \
                                            and in_interval(call["Q0"], self.Q0_range) \
                                            and in_interval(call["dG"], self.dG_range):
                                        type = "duplication" if call["type"] == 1 else "deletion"

                                        row = [self.file_title(ix[i]), caller, type, c, call["start"], call["end"],
                                               call["size"], call["cnv"], call["p_val"], call["p_val_2"],
                                               call["p_val_3"], call["p_val_4"], call["Q0"], call["pN"], call["dG"]]
                                        if self.annotate:
                                            row.append(annotator.get_info("%s:%d-%d" % (c, call["start"], call["end"])))
                                        ret.append(row)
            elif caller == "combined_mosaic":
                for i in range(n):
                    io = self.io[ix[i]]
                    chroms = io.rd_chromosomes()
                    for c in chroms:
                        if (c in self.chrom) or len(self.chrom) == 0:
                            flag = (FLAG_USEMASK if self.rd_use_mask else 0) | \
                                   (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | \
                                   (FLAG_USEMASK if self.snp_use_mask else 0) | \
                                   (FLAG_USEID if self.snp_use_id else 0) | \
                                   (FLAG_USEHAP if self.snp_use_phase else 0)
                            if io.signal_exists(c, bin_size, "calls combined", flag):
                                calls = io.read_calls(c, bin_size, "calls combined", flag)
                                for call in calls:
                                    if in_interval(call["size"], self.size_range) \
                                            and in_interval(call["p_val"], self.p_range) \
                                            and in_interval(call["pN"], self.pN_range) \
                                            and in_interval(call["Q0"], self.Q0_range) \
                                            and in_interval(call["bins"], self.bins_range) \
                                            and in_interval(abs(call["baf"]), self.baf_range):

                                        if n > 1:
                                            print("%s\t" % self.file_title(ix[i]), end="")
                                        if len(self.callers) > 1:
                                            print("%s\t" % caller, end="")
                                        keys = ["start", "end", "size", "cnv", "p_val", "lh_del", "lh_loh",
                                                "lh_dup", "Q0", "pN", "pNS", "pP", "bins", "bins", "baf",
                                                "rd_p_val", "baf_p_val", "hets", "homs"]
                                        type = {-1: "deletion", 0: "cnnloh", 1: "duplication"}[call["type"]]
                                        row = [self.file_title(i), caller, type, c] + [call[k] for k in keys]
                                        row[16] = bin_size
                                        for m in range(2):
                                            row += call["models"][m]

                                        if self.annotate:
                                            row.append(annotator.get_info("%s:%d-%d" % (data[3], data[4], data[5])))
                                        ret.append(row)
        return ret

    def print_calls_file(self):
        format = self.print_filename.split(".")[-1]
        calls = self.get_calls()
        if self.print_filename == "":
            for call in calls:
                print(*call, sep="\t")
        elif format == "tsv":
            with open(self.print_filename, 'w') as f:
                for call in calls:
                    print(*call, sep="\t", file=f)
        elif format == "xlsx":
            import xlsxwriter
            workbook = xlsxwriter.Workbook(self.print_filename)
            files_callers = []
            sheets = {}
            rix = {}
            for call in calls:
                caller = call[1]
                fc = call[0] + " (" + caller + ")"
                sfc = call[0][:25] + " " + ({"rd_mean_shift": "ms", "combined_mosaic": "2d"}[caller])
                if fc not in files_callers:
                    sheets[fc] = workbook.add_worksheet(sfc)
                    rix[fc] = 0
                    files_callers.append(fc)
            for call in calls:
                caller = call[1]
                fc = call[0] + " (" + caller + ")"
                cix = 0
                for f in call[2:]:
                    sheets[fc].write(rix[fc], cix, f)
                    cix += 1
                rix[fc] += 1
            workbook.close()
        elif format == "vcf":
            samples = []
            for call in calls:
                sample = call[0]
                if sample not in samples:
                    samples.append(sample)
            header = """##fileformat=VCFv4.1
##fileDate={date}
##reference={rg}
##source=CNVpytor
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=pytorRD,Number=1,Type=Float,Description="Normalized RD">
##INFO=<ID=pytorP1,Number=1,Type=Float,Description="e-val by t-test">
##INFO=<ID=pytorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">
##INFO=<ID=pytorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">
##INFO=<ID=pytorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">
##INFO=<ID=pytorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">
##INFO=<ID=pytorPN,Number=1,Type=Integer,Description="Fraction of N bases">
##INFO=<ID=pytorDG,Number=1,Type=Integer,Description="Distance to nearest gap in reference genome">
##INFO=<ID=pytorCL,Number=1,Type=Integer,Description="Caller method">
##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=LOH,Description="Copy number neutral loss of heterozygosity">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">;
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}"""
            if self.reference_genome:
                rg = self.reference_genome["name"]
            else:
                rg = "unknown"
            header = header.format(date=datetime.date.today().strftime("%Y-%m-%d"), rg=rg, samples="\t".join(samples))
            ii = 0
            with open(self.print_filename, 'w') as f:
                print(header, file=f)
                for call in calls:
                    ii += 1
                    id = "CNVpytor_" + {"deletion": "del", "duplication": "dup", "cnnloh": "loh"}[call[2]] + str(ii)
                    alt = {"deletion": "<DEL>", "duplication": "<DUP>", "cnnloh": "<LOH>"}[call[2]]
                    info = "END=" + str(int(call[5])) + ";IMPRECISE;SVLEN=" + str(int(call[6])) + ";SVTYPE=" + alt[1:4]
                    info += ";pytorRD=" + str(call[7])
                    info += ";pytorP1=" + str(call[8])
                    info += ";pytorP2=" + str(call[9])
                    info += ";pytorP3=" + str(call[10])
                    info += ";pytorP4=" + str(call[11])
                    info += ";pytorQ0=" + str(call[12])
                    info += ";pytorPN=" + str(call[13])
                    info += ";pytorDG=" + str(call[14])
                    info += ";pytorCL=" + call[1]
                    format = "GT:CN"
                    row = [call[3], int(call[4]), id, ".", alt, ".", "PASS", info, format]
                    for sample in samples:
                        if sample == call[0]:
                            if call[2] == "deletion" and call[7] < 0.25:
                                row.append("1/1:0")
                            elif call[2] == "deletion" and call[7] > 0.25:
                                row.append("0/1:0")
                            elif call[2] == "duplication" and call[7] <= 1.75:
                                row.append("0/1:2")
                            elif call[2] == "duplication" and call[7] > 1.75 and call[7] <= 2.25:
                                row.append("1/1:2")
                            elif call[2] == "duplication" and call[7] > 2.25:
                                row.append("./1:%.2f" % call[7])
                            else:
                                row.append("./.:.")
                        else:
                            row.append("./.:.")
                    print(*row, sep="\t", file=f)
        if self.plot:
            for call in calls:
                plot_start = call[4] - call[6]
                if plot_start < 1:
                    plot_start = 1
                plot_end = call[5] + call[6]
                self.multiple_regions(["%s:%d-%d" % (c, plot_start, plot_end)])

    def print_calls(self):
        bin_size = self.bin_size
        n = len(self.plot_files)
        ix = self.plot_files
        if self.annotate:
            annotator = Annotator(self.reference_genome)
        for caller in self.callers:
            if caller == "rd_mean_shift":
                for i in range(n):
                    io = self.io[ix[i]]
                    chroms = io.rd_chromosomes()
                    for c in chroms:
                        if (c in self.chrom) or len(self.chrom) == 0:
                            flag = (FLAG_USEMASK if self.rd_use_mask else 0) | \
                                   (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
                            if io.signal_exists(c, bin_size, "calls", flag):
                                calls = io.read_calls(c, bin_size, "calls", flag)
                                for call in calls:
                                    if in_interval(call["size"], self.size_range) \
                                            and in_interval(call["p_val"], self.p_range) \
                                            and in_interval(call["pN"], self.pN_range) \
                                            and in_interval(call["Q0"], self.Q0_range) \
                                            and in_interval(call["dG"], self.dG_range):
                                        type = "duplication" if call["type"] == 1 else "deletion"
                                        if n > 1:
                                            print("%s\t" % self.file_title(i), end="")
                                        if len(self.callers) > 1:
                                            print("%s\t" % caller, end="")
                                        print("%s\t%s:%d-%d\t%d\t%.4f\t%e\t%e\t%e\t%e\t%.4f\t%.4f\t%d\t" % (
                                            type, c, call["start"], call["end"], call["size"], call["cnv"],
                                            call["p_val"],
                                            call["p_val_2"], call["p_val_3"], call["p_val_4"], call["Q0"], call["pN"],
                                            call["dG"]), end="")
                                        if self.annotate:
                                            print("\t%s" % annotator.get_info(
                                                "%s:%d-%d" % (c, call["start"], call["end"])))
                                        else:
                                            print()
                                        if self.plot:
                                            plot_start = call["start"] - call["size"]
                                            if plot_start < 1:
                                                plot_start = 1
                                            plot_end = call["end"] + call["size"]
                                            self.multiple_regions(["%s:%d-%d" % (c, plot_start, plot_end)])
            elif caller == "combined_mosaic":
                for i in range(n):
                    io = self.io[ix[i]]
                    chroms = io.rd_chromosomes()
                    for c in chroms:
                        if (c in self.chrom) or len(self.chrom) == 0:
                            flag = (FLAG_USEMASK if self.rd_use_mask else 0) | \
                                   (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | \
                                   (FLAG_USEMASK if self.snp_use_mask else 0) | \
                                   (FLAG_USEID if self.snp_use_id else 0) | \
                                   (FLAG_USEHAP if self.snp_use_phase else 0)
                            if io.signal_exists(c, bin_size, "calls combined", flag):
                                calls = io.read_calls(c, bin_size, "calls combined", flag)
                                for call in calls:
                                    if in_interval(call["size"], self.size_range) \
                                            and in_interval(call["p_val"], self.p_range) \
                                            and in_interval(call["pN"], self.pN_range) \
                                            and in_interval(call["Q0"], self.Q0_range) \
                                            and in_interval(call["bins"], self.bins_range) \
                                            and in_interval(abs(call["baf"]), self.baf_range):
                                        if n > 1:
                                            print("%s\t" % self.file_title(i), end="")
                                        if len(self.callers) > 1:
                                            print("%s\t" % caller, end="")
                                        keys = ["start", "end", "size", "cnv", "p_val", "lh_del", "lh_loh",
                                                "lh_dup", "Q0", "pN", "pNS", "pP", "bins", "bins", "baf",
                                                "rd_p_val", "baf_p_val", "hets", "homs"]
                                        type = {-1: "deletion", 0: "cnnloh", 1: "duplication"}[call["type"]]
                                        data = [type, c] + [call[k] for k in keys]
                                        data[14] = bin_size
                                        for m in range(2):
                                            data += call["models"][m]

                                        print(("%s\t%s:%d-%d\t%d\t%.4f\t%e\t%e\t%e\t%e" + \
                                               "\t%.4f\t%.4f\t%.4f\t%.4f\t" + "%d\t%d\t%.4f\t%e\t%e\t%d\t%d\t%d\t" + \
                                               "CN%d/CN%d\t%e\t%.4f\t%d\tCN%d/CN%d\t%e\t%.4f") % tuple(data), end="")
                                        if self.annotate:
                                            print("\t%s" % annotator.get_info("%s:%d-%d" % (data[1], data[2], data[3])))
                                        else:
                                            print()
                                        if self.plot:
                                            plot_start = call["start"] - call["size"]
                                            if plot_start < 1:
                                                plot_start = 1
                                            plot_end = call["end"] + call["size"]
                                            _logger.debug("Plotting region %s:%d-%d" % (c, plot_start, plot_end))
                                            self.multiple_regions(["%s:%d-%d,%s:%d-%d,%s:%d-%d" \
                                                                   % (
                                                                       c, plot_start, call["start"] - 1, c,
                                                                       call["start"],
                                                                       call["end"], c, call["end"] + 1, plot_end)])

    def print_simple_merged_calls(self):

        bin_size = self.bin_size
        n = len(self.plot_files)
        if n == 0:
            return
        ix = self.plot_files
        format = self.print_filename.split(".")[-1]
        if format == "tsv":
            f = open(self.print_filename, 'w')
        elif format == "xlsx":
            import xlsxwriter
            if os.path.exists(self.print_filename):
                os.remove(self.print_filename)
            workbook = xlsxwriter.Workbook(self.print_filename)
            sheet = workbook.add_worksheet("merged_calls")
            header = ["TYPE", "REGION", "SIZE"]
            for i in range(n):
                header.append(self.file_title(ix[i]))
            if self.annotate:
                header.append("GENES")
            styleh = workbook.add_format({'bold': True, 'font_color': 'white'})
            styleh.set_pattern(1)  # This is optional when using a solid fill.
            styleh.set_bg_color('#555555')
            styleh2 = workbook.add_format({'bold': True, 'font_color': 'white'})
            styleh2.set_pattern(1)  # This is optional when using a solid fill.
            styleh2.set_bg_color('#555555')
            styleh2.set_rotation(75)
            style_r = workbook.add_format()
            style_r.set_pattern(1)  # This is optional when using a solid fill.
            style_r.set_bg_color('red')
            style_g = workbook.add_format()
            style_g.set_pattern(1)  # This is optional when using a solid fill.
            style_g.set_bg_color('green')
            style_size = workbook.add_format({'num_format': '#,##0'})
            style_cn = workbook.add_format({'num_format': '0'})
            style_cn_b = workbook.add_format({'num_format': '0', 'bold': True})
            sheet.set_column(0, 0, 10)
            sheet.set_column(1, 1, 22)
            sheet.set_column(2, 2, 10)
            if self.annotate:
                sheet.set_column(len(header) - 1, len(header) - 1, 100)

            for col, val in enumerate(header):
                if col > 2 and col < len(header) - int(self.annotate):
                    sheet.write(0, col, val, styleh2)
                else:
                    sheet.write(0, col, val, styleh)
            ri = 0
        if self.annotate:
            annotator = Annotator(self.reference_genome)
        chroms = self.io[ix[0]].rd_chromosomes()
        for c in chroms:
            if (c in self.chrom) or len(self.chrom) == 0:
                flag = (FLAG_USEMASK if self.rd_use_mask else 0) | (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
                calls = [list(filter(lambda call: in_interval(call["size"], self.size_range) \
                                                  and in_interval(call["p_val"], self.p_range) \
                                                  and in_interval(call["pN"], self.pN_range) \
                                                  and in_interval(call["Q0"], self.Q0_range) \
                                                  and in_interval(call["dG"], self.dG_range),
                                     self.io[ix[i]].read_calls(c, bin_size, "calls", flag))) for i in range(n)]
                pointers = [0] * n
                while any([pointers[i] < len(calls[i]) for i in range(n)]):
                    starts = [calls[i][pointers[i]]["start"] if pointers[i] < len(calls[i]) else np.inf for i in
                              range(n)]
                    mini = starts.index(min(starts))
                    maxend = 0
                    toupdate = []
                    minend = calls[mini][pointers[mini]]["end"]
                    maxstart = 0
                    files = []
                    types = []
                    cns = []
                    for i in range(n):
                        if (pointers[i] < len(calls[i])) and ((min(calls[i][pointers[i]]["end"],
                                                                   calls[mini][pointers[mini]]["end"]) -
                                                               calls[i][pointers[i]]["start"]) > (
                                                                      0.5 * calls[mini][pointers[mini]]["size"])) \
                                and ((min(calls[i][pointers[i]]["end"],
                                          calls[mini][pointers[mini]]["end"]) -
                                      calls[i][pointers[i]]["start"]) > (
                                             0.5 * (calls[i][pointers[i]]["end"] - calls[i][pointers[i]]["start"]))):
                            toupdate.append(i)
                            call = calls[i][pointers[i]]
                            if call["end"] > maxend:
                                maxend = call["end"]
                            if call["end"] < minend:
                                minend = call["end"]
                            if call["start"] > maxstart:
                                maxstart = call["start"]
                            type = "duplication" if call["type"] == 1 else "deletion"
                            types.append(type)
                            files.append(i)
                            cns.append(int(call["cnv"] * 2))
                    type = max(set(types), key=types.count)
                    data = [type, c, maxstart, minend, minend - maxstart + 1]
                    genotypes = [
                        self.genotype([bin_size], "%s:%d-%d" % (c, maxstart, minend), file_index=ix[i], p_val=True)[0]
                        for i
                        in range(n)]
                    copynumbers = [c[3] for c in genotypes]
                    if np.all([np.abs(c - np.round(c)) < 0.25 for c in copynumbers]) or True:
                        if self.print_filename == "":
                            print(("%s\t%s:%d-%d\t%d" + n * "\t%.2f") % tuple(data + copynumbers), end="")
                            print("\t%s" % str(files), end="")
                            if self.annotate:
                                print("\t%s" % annotator.get_info("%s:%d-%d" % (c, maxstart, minend)))
                            else:
                                print()
                        elif format == "tsv":
                            print(("%s\t%s:%d-%d\t%d" + n * "\t%.2f") % tuple(data + copynumbers), end="", file=f)
                            print("\t%s" % str(files), end="", file=f)
                            if self.annotate:
                                print("\t%s" % annotator.get_info("%s:%d-%d" % (c, maxstart, minend)), file=f)
                            else:
                                print(file=f)
                        elif format == "xlsx":
                            ri += 1
                            if type == "deletion":
                                sheet.write(ri, 0, data[0], style_r)
                            else:
                                sheet.write(ri, 0, data[0], style_g)
                            sheet.write(ri, 1, "%s:%d-%d" % (c, maxstart, minend))
                            sheet.write(ri, 2, data[4], style_size)
                            for col, val in enumerate(copynumbers):
                                if col in files:
                                    sheet.write(ri, 3 + col, val, style_cn_b)
                                else:
                                    sheet.write(ri, 3 + col, val, style_cn)
                            if self.annotate:
                                sheet.write(ri, 3 + len(copynumbers),
                                            annotator.get_info("%s:%d-%d" % (c, maxstart, minend)))

                        if self.plot:
                            plot_start = maxstart - (minend - maxstart)
                            if plot_start < 1:
                                plot_start = 1
                            plot_end = minend + (minend - maxstart)
                            self.multiple_regions(["%s:%d-%d" % (c, plot_start, plot_end)])
                    for i in toupdate:
                        pointers[i] += 1
        if format == "tsv":
            f.close()
        elif format == "xlsx":
            sheet.conditional_format(1, 3, ri, len(header) - int(self.annotate), {'type': '3_color_scale',
                                                                                  'min_color': "#FF0000",
                                                                                  'mid_color': "#FFFFFF",
                                                                                  'max_color': "#00FF00",
                                                                                  'min_type': 'num',
                                                                                  'min_value': 0,
                                                                                  'mid_type': 'num',
                                                                                  'mid_value': 2,
                                                                                  'max_type': 'num',
                                                                                  'max_value': 4
                                                                                  })
            workbook.close()

    def manhattan(self, plot_type="rd"):
        bin_size = self.bin_size
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for manhattan.")
            return
        n = len(self.plot_files)
        ix = self.plot_files

        self.new_figure(panel_count=n, grid=(1, n), panel_size=(24, 2))
        for i in range(n):
            ax = self.next_panel()
            io = self.io[ix[i]]
            ax.set_title(self.file_title(ix[i]), position=(0.01, 1.01),
                         fontdict={'verticalalignment': 'bottom', 'horizontalalignment': 'left'})

            if plot_type == "rd":
                chroms = []
                for c, (l, t) in self.reference_genome["chromosomes"].items():
                    rd_chr = io.rd_chromosome_name(c)
                    if len(self.chrom) == 0 or (rd_chr in self.chrom) or (c in self.chrom):
                        if io.signal_exists(rd_chr, bin_size, "RD", 0) and \
                                io.signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                                (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                            chroms.append((rd_chr, l))

                apos = 0
                xticks = [0]

                max_m, stdev = io.rd_normal_level(bin_size, FLAG_GC_CORR)
                for c, l in chroms:
                    flag_rd = (FLAG_USEMASK if self.rd_use_mask else 0)
                    his_p = io.get_signal(c, bin_size, "RD", flag_rd)
                    his_p_corr = io.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
                    if self.rd_manhattan_call:
                        his_p_call = io.get_signal(c, bin_size, "RD call", flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_seg = io.get_signal(c, bin_size, "RD mosaic segments",
                                                         flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_seg = segments_decode(his_p_mosaic_seg)
                        his_p_mosaic_call = io.get_signal(c, bin_size, "RD mosaic call",
                                                          flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_seg_2d = io.get_signal(c, bin_size, "RD mosaic segments 2d",
                                                            flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_seg_2d = segments_decode(his_p_mosaic_seg_2d)
                        his_p_mosaic_call_2d = io.get_signal(c, bin_size, "RD mosaic call 2d",
                                                             flag_rd | FLAG_GC_CORR)
                        his_p_mosaic = np.zeros_like(his_p) * np.nan
                        if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and (
                                "rd_mosaic" in self.callers):
                            for seg, lev in zip(list(his_p_mosaic_seg), list(his_p_mosaic_call[0])):
                                for segi in seg:
                                    his_p_mosaic[segi] = lev
                        his_p_mosaic_2d = np.zeros_like(his_p) * np.nan
                        if his_p_mosaic_call_2d is not None and len(
                                his_p_mosaic_call_2d) > 0 and ("combined_mosaic" in self.callers):
                            for seg, lev in zip(list(his_p_mosaic_seg_2d), list(his_p_mosaic_call_2d[0])):
                                for segi in seg:
                                    his_p_mosaic_2d[segi] = lev
                    pos = range(apos, apos + len(his_p))
                    ax.text(apos + len(his_p) // 2, max_m // 10, Genome.canonical_chrom_name(c),
                            fontsize=8, verticalalignment='bottom', horizontalalignment='center', )
                    if self.markersize == "auto":
                        plt.plot(pos, his_p_corr, ls='', marker='.')
                    else:
                        plt.plot(pos, his_p_corr, ls='', marker='.', markersize=self.markersize)
                    if self.rd_manhattan_call:
                        if his_p_call is not None and len(his_p_call) > 0 and ("rd_mean_shift" in self.callers):
                            plt.step(pos, his_p_call, "r")
                        if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and (
                                "rd_mosaic" in self.callers):
                            plt.plot(pos, his_p_mosaic, "k")
                        if his_p_mosaic_call_2d is not None and len(
                                his_p_mosaic_call_2d) > 0 and ("combined_mosaic" in self.callers):
                            plt.plot(pos, his_p_mosaic_2d, "k")
                    apos += len(his_p)
                    xticks.append(apos)
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 15, 0.5) * max_m, minor=[])
                ax.xaxis.set_ticks(xticks, minor=[])
                ax.set_ylim([self.rd_manhattan_range[0] * max_m, self.rd_manhattan_range[1] * max_m])
                n_bins = apos
                ax.set_xlim([0, n_bins])
                ax.grid()

            elif plot_type == "baf_mosaic":
                chroms = []
                snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                    FLAG_USEHAP if self.snp_use_phase else 0)
                for c, (l, t) in self.reference_genome["chromosomes"].items():
                    snp_chr = io.snp_chromosome_name(c)
                    if len(self.chrom) == 0 or (snp_chr in self.chrom) or (c in self.chrom):
                        if io.signal_exists(snp_chr, bin_size, "SNP likelihood call", snp_flag) and \
                                io.signal_exists(snp_chr, bin_size, "SNP likelihood segments", snp_flag) and \
                                (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                            chroms.append((snp_chr, l))

                apos = 0
                xticks = [0]

                cix = 0
                cmap = list(map(colors.to_rgba, plt.rcParams['axes.prop_cycle'].by_key()['color']))
                for c, l in chroms:
                    likelihood = io.get_signal(c, bin_size, "SNP likelihood call", snp_flag)
                    segments = segments_decode(io.get_signal(c, bin_size, "SNP likelihood segments", snp_flag))
                    call_pos = []
                    call_baf = []
                    call_c = []
                    for s, lh in zip(segments, likelihood):
                        b, p = likelihood_baf_pval(lh)
                        if b > 0 and len(s) > self.min_segment_size:
                            alpha = -np.log(p + 1e-40) / self.contrast
                            if alpha > 1:
                                alpha = 1
                            for pos in s:
                                call_pos.append(apos + pos)
                                call_baf.append(b)
                                color = cmap[cix % len(cmap)]
                                color = (color[0], color[1], color[2], alpha)
                                call_c.append(color)

                    ax.text(apos + l // bin_size // 2, 0.4, Genome.canonical_chrom_name(c),
                            fontsize=8, verticalalignment='bottom', horizontalalignment='center', )
                    plt.scatter(call_pos, call_baf, s=20, color=np.array(call_c), edgecolors='face', marker='|')
                    apos += l // bin_size
                    xticks.append(apos)
                    cix += 1

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 0.5, 0.1), minor=[])
                ax.xaxis.set_ticks(xticks, minor=[])
                ax.set_ylim([0, 0.5])
                n_bins = apos
                ax.set_xlim([0, n_bins])
                ax.grid()

            elif plot_type == "rd_mean_shift":
                chroms = []
                flag = (FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR

                for c, (l, t) in self.reference_genome["chromosomes"].items():
                    rd_chr = io.rd_chromosome_name(c)
                    if rd_chr is not None and len(self.chrom) == 0 or (rd_chr in self.chrom) or (c in self.chrom):
                        if (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                            chroms.append((rd_chr, l))

                apos = 0
                xticks = [0]

                cix = 0
                cmap = list(map(colors.to_rgba, plt.rcParams['axes.prop_cycle'].by_key()['color']))
                for c, l in chroms:
                    call_pos = []
                    call_conc = []
                    call_c = []
                    if io.signal_exists(c, bin_size, "calls", flag):
                        calls = io.read_calls(c, bin_size, "calls", flag)

                        for call in calls:
                            if in_interval(call["size"], self.size_range) and in_interval(call["p_val"], self.p_range) \
                                    and in_interval(call["pN"], self.pN_range) \
                                    and in_interval(call["Q0"], self.Q0_range):
                                alpha = - np.log(call["p_val"] + 1e-40) / self.contrast
                                if alpha > 1:
                                    alpha = 1
                                if alpha < 0:
                                    alpha = 0
                                for pos in range(int(call["start"]) // bin_size, int(call["end"]) // bin_size + 1):
                                    call_pos.append(apos + pos)
                                    level = call["cnv"] * 2
                                    if level > 4:
                                        level = 4
                                    call_conc.append(level)
                                    if call["type"] == 1:
                                        call_c.append((0, 1, 0, alpha))
                                    elif call["type"] == -1:
                                        call_c.append((1, 0, 0, alpha))
                                    else:
                                        call_c.append((0, 0, 1, alpha))
                        ax.text(apos + l // bin_size // 2, 0.4, Genome.canonical_chrom_name(c),
                                fontsize=8, verticalalignment='bottom', horizontalalignment='center', )
                        plt.scatter(call_pos, call_conc, s=20, color=np.array(call_c), edgecolors='face', marker='|')
                        apos += l // bin_size
                        xticks.append(apos)
                        cix += 1

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 4.0, 1.0), minor=[])
                ax.xaxis.set_ticks(xticks, minor=[])
                ax.set_ylim([0, 4.0])
                n_bins = apos
                ax.set_xlim([0, n_bins])
                ax.grid()

            elif plot_type == "combined_mosaic":
                chroms = []
                flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                    FLAG_USEHAP if self.snp_use_phase else 0) | (FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR

                for c, (l, t) in self.reference_genome["chromosomes"].items():
                    snp_chr = io.snp_chromosome_name(c)
                    if snp_chr is not None and len(self.chrom) == 0 or (snp_chr in self.chrom) or (c in self.chrom):
                        if (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                            chroms.append((snp_chr, l))

                apos = 0
                xticks = [0]

                cix = 0
                cmap = list(map(colors.to_rgba, plt.rcParams['axes.prop_cycle'].by_key()['color']))
                for c, l in chroms:
                    call_pos = []
                    call_conc = []
                    call_c = []
                    if io.signal_exists(c, bin_size, "calls combined", flag):
                        calls = io.read_calls(c, bin_size, "calls combined", flag)

                        for call in calls:
                            if call["bins"] > self.min_segment_size:
                                alpha = -np.log(call["p_val"] + 1e-40) / self.contrast
                                if alpha > 1:
                                    alpha = 1
                                for pos in range(int(call["start"]) // bin_size, int(call["end"]) // bin_size + 1):
                                    call_pos.append(apos + pos)
                                    call_conc.append(call["models"][0][4])
                                    if call["type"] == 1:
                                        call_c.append((0, 1, 0, alpha))
                                    elif call["type"] == -1:
                                        call_c.append((1, 0, 0, alpha))
                                    else:
                                        call_c.append((0, 0, 1, alpha))

                        ax.text(apos + l // bin_size // 2, 0.4, Genome.canonical_chrom_name(c),
                                fontsize=8, verticalalignment='bottom', horizontalalignment='center', )
                        plt.scatter(call_pos, call_conc, s=20, color=np.array(call_c), edgecolors='face', marker='|')
                        apos += l // bin_size
                        xticks.append(apos)
                        cix += 1

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 1.0, 0.1), minor=[])
                ax.xaxis.set_ticks(xticks, minor=[])
                ax.set_ylim([0, 1.0])

                n_bins = apos
                ax.set_xlim([0, n_bins])
                ax.grid()

        self.fig_show(suffix="manhattan" if plot_type == "rd" else "snp_calls")

    def callmap(self, color="frequency", background="white", pixel_size=1700000, max_p_val=1e-20, min_freq=0.01,
                plot="cmap"):
        bin_size = self.bin_size
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for callmap.")
            return
        n = len(self.plot_files)
        ix = self.plot_files

        if plot:
            self.new_figure(panel_count=n, grid=(1, 1), panel_size=(24, 0.24 * n))

        chroms = []
        starts = []
        ends = []
        pixels = 0
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            if l > 10 * bin_size:
                if len(self.chrom) == 0 or (c in self.chrom) or (self.io[0].snp_chromosome_name(c) in self.chrom):
                    chroms.append(c)
                    starts.append(pixels)
                    pixels += l // pixel_size + 1
                    ends.append(pixels - 1)

        cmap = np.zeros((n, pixels, 3))
        cmap[:, ends, :] = 1

        for i in range(n):
            io = self.io[ix[i]]
            print(io.filename)
            flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                FLAG_USEHAP if self.snp_use_phase else 0) | (FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR
            flag_rd = FLAG_GC_CORR | (FLAG_USEMASK if self.rd_use_mask else 0)
            for c, start in zip(chroms, starts):
                snp_chr = io.snp_chromosome_name(c)
                if io.signal_exists(snp_chr, bin_size, "calls combined", flag):

                    calls = io.read_calls(snp_chr, bin_size, "calls combined", flag)
                    segments = io.get_signal(snp_chr, bin_size, "RD mosaic segments 2d", flag_rd)
                    segments = segments_decode(segments)

                    for call in calls:
                        if call["bins"] > self.min_segment_size and call["p_val"] < max_p_val and "segment" in call and \
                                call["models"][0][4] > min_freq:
                            cix = int(call["type"]) + 1
                            for b in segments[int(call["segment"])]:
                                if color == "frequency":
                                    cmap[i, start + b * bin_size // pixel_size, cix] = max(
                                        cmap[i, start + b * bin_size // pixel_size, cix], call["models"][0][4])
                                elif color == "coverage":
                                    cmap[i, start + b * bin_size // pixel_size, cix] += bin_size / pixel_size
                                else:  # model copy number
                                    if call["models"][0][0] == 0:
                                        cmap[i, start + b * bin_size // pixel_size, 0] = 1
                                    elif call["models"][0][0] == 1:
                                        cmap[i, start + b * bin_size // pixel_size, 0] = 1
                                        cmap[i, start + b * bin_size // pixel_size, 1] = 1
                                    elif call["models"][0][0] == 2:
                                        cmap[i, start + b * bin_size // pixel_size, 2] = 1
                                    else:
                                        cn = call["models"][0][0]
                                        if cn > 6:
                                            cn = 6
                                        cmap[i, start + b * bin_size // pixel_size, 1] = (2 + cn) / 8

        def b2w(pixel):
            if np.all(pixel == 1):
                pixel[:] = 0
            elif pixel[0] > pixel[1] and pixel[0] > pixel[2]:
                pixel[1] = pixel[2] = 1 - pixel[0]
                pixel[0] = 1
            elif pixel[1] > pixel[2]:
                pixel[0] = pixel[2] = 1 - pixel[1]
                pixel[1] = 1
            else:
                pixel[0] = pixel[1] = 1 - pixel[2]
                pixel[2] = 1
            return pixel

        if background == "white":
            cmap = cmap.reshape(n * pixels, 3)
            np.apply_along_axis(b2w, 1, cmap)
            cmap = cmap.reshape(n, pixels, 3)

        cmap = (255 * cmap).astype("int")
        if plot == "cmap":
            self.new_figure(panel_count=1, grid=(1, 1), panel_size=(24, 0.24 * n))
            ax = self.next_panel()
            plt.imshow(cmap, aspect='auto')
            for i in ends[:-1]:
                plt.axvline(x=i - 0.5, color='red', linewidth=0.5)
            ax.set_yticks([])
            ax.set_yticklabels([])
            ax.set_xticks((np.array(starts) + np.array(ends)) / 2)
            chroms = list(map(Genome.canonical_chrom_name, chroms))
            ax.set_xticklabels(chroms)
            self.fig_show(suffix="callmap")
        elif plot == "regions":
            self.new_figure(panel_count=1, grid=(1, 1), panel_size=(24, 24))
            ax = self.next_panel()
            corr = np.corrcoef(
                np.concatenate((cmap[:, :, 0].transpose(), cmap[:, :, 1].transpose(), cmap[:, :, 2].transpose()),
                               axis=0))
            plt.imshow(corr, aspect='auto', vmin=-1, vmax=1)
            plt.colorbar()
            starts3 = np.concatenate((np.array(starts), np.array(starts) + ends[-1], np.array(starts) + 2 * ends[-1]))
            ends3 = np.concatenate((np.array(ends), np.array(ends) + ends[-1], np.array(ends) + 2 * ends[-1]))
            for i in ends3[:-1]:
                plt.axvline(x=i - 0.5, color='red', linewidth=0.5)
                plt.axhline(y=i - 0.5, color='red', linewidth=0.5)

            ax.set_xticks((starts3 + ends3) / 2)
            ax.set_yticks((starts3 + ends3) / 2)
            chroms = list(map(Genome.canonical_chrom_name, chroms))
            ax.set_xticklabels(chroms + chroms + chroms)
            ax.set_yticklabels(chroms + chroms + chroms)
            self.fig_show(suffix="callmap")
        else:
            self.new_figure(panel_count=2, panel_size=(12, 12))
            ax = self.next_panel()
            x = np.concatenate((cmap[:, :, 0], cmap[:, :, 1], cmap[:, :, 2]),
                               axis=1)
            corr = np.corrcoef(x)
            plt.imshow(corr, aspect='auto', vmin=-1, vmax=1)
            plt.colorbar()
            ax = plt.gca()

            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax = self.next_panel()
            Z = hierarchy.linkage(x, 'average', 'correlation')
            dn = hierarchy.dendrogram(Z)

            self.fig_show(suffix="callmap")
        return cmap

    def callmap_regions(self, regions, color="frequency", background="white", pixel_size=1700000, max_p_val=1e-20,
                        min_freq=0.01, plot="cmap"):
        bin_size = self.bin_size
        n = len(self.plot_files)
        ix = self.plot_files

        if plot:
            self.new_figure(panel_count=n, grid=(1, 1), panel_size=(24, 0.24 * n))

        chroms = []
        starts = []
        ends = []
        pixels = 0
        regs = decode_region(regions, max_size=1000000000)
        for c, (pos1, pos2) in regs:
            chr_len = self.io[ix[0]].get_chromosome_length(c)
            if chr_len is not None and pos2 == 1000000000:
                pos2 = chr_len
            starts.append(pixels)
            pixels += (pos2 - pos1) // pixel_size + 1
            ends.append(pixels - 1)

        cmap = np.zeros((n, pixels, 3))
        cmap[:, ends, :] = 1

        for i in range(n):
            io = self.io[ix[i]]
            print(io.filename)
            flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                FLAG_USEHAP if self.snp_use_phase else 0) | (FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR
            flag_rd = FLAG_GC_CORR | (FLAG_USEMASK if self.rd_use_mask else 0)
            for (c, (pos1, pos2)), start in zip(regs, starts):
                chr_len = self.io[ix[0]].get_chromosome_length(c)
                if chr_len is not None and pos2 == 1000000000:
                    pos2 = chr_len
                snp_chr = io.snp_chromosome_name(c)
                if io.signal_exists(snp_chr, bin_size, "calls combined", flag):

                    calls = io.read_calls(snp_chr, bin_size, "calls combined", flag)
                    segments = io.get_signal(snp_chr, bin_size, "RD mosaic segments 2d", flag_rd)
                    segments = segments_decode(segments)

                    for call in calls:
                        if call["bins"] > self.min_segment_size and call["p_val"] < max_p_val and "segment" in call and \
                                call["models"][0][4] > min_freq:
                            cix = int(call["type"]) + 1
                            for b in segments[int(call["segment"])]:
                                if (b * bin_size > pos1) and (b * bin_size < pos2):
                                    if color == "frequency":
                                        cmap[i, start + (b * bin_size - pos1) // pixel_size, cix] = max(
                                            cmap[i, start + (b * bin_size - pos1) // pixel_size, cix],
                                            call["models"][0][4])
                                    elif color == "coverage":
                                        cmap[i, start + (
                                                b * bin_size - pos1) // pixel_size, cix] += bin_size / pixel_size
                                    else:  # model copy number
                                        if call["models"][0][0] == 0:
                                            cmap[i, start + (b * bin_size - pos1) // pixel_size, 0] = 1
                                        elif call["models"][0][0] == 1:
                                            cmap[i, start + (b * bin_size - pos1) // pixel_size, 0] = 1
                                            cmap[i, start + (b * bin_size - pos1) // pixel_size, 1] = 1
                                        elif call["models"][0][0] == 2:
                                            cmap[i, start + (b * bin_size - pos1) // pixel_size, 2] = 1
                                        else:
                                            cn = call["models"][0][0]
                                            if cn > 6:
                                                cn = 6
                                            cmap[i, start + (b * bin_size - pos1) // pixel_size, 1] = (2 + cn) / 8

        def b2w(pixel):
            if np.all(pixel == 1):
                pixel[:] = 0
            elif pixel[0] > pixel[1] and pixel[0] > pixel[2]:
                pixel[1] = pixel[2] = 1 - pixel[0]
                pixel[0] = 1
            elif pixel[1] > pixel[2]:
                pixel[0] = pixel[2] = 1 - pixel[1]
                pixel[1] = 1
            else:
                pixel[0] = pixel[1] = 1 - pixel[2]
                pixel[2] = 1
            return pixel

        if background == "white":
            cmap = cmap.reshape(n * pixels, 3)
            np.apply_along_axis(b2w, 1, cmap)
            cmap = cmap.reshape(n, pixels, 3)

        cmap = (255 * cmap).astype("int")
        if plot == "cmap":
            self.new_figure(panel_count=1, grid=(1, 1), panel_size=(24, 0.24 * n))
            ax = self.next_panel()
            plt.imshow(cmap, aspect='auto')
            for i in ends[:-1]:
                plt.axvline(x=i - 0.5, color='red', linewidth=0.5)
            ax.set_yticks([])
            ax.set_yticklabels([])
            ax.set_xticks((np.array(starts) + np.array(ends)) / 2)
            ax.set_xticklabels(regions.split(","))
            self.fig_show(suffix="callmap")
        elif plot == "regions":
            self.new_figure(panel_count=1, grid=(1, 1), panel_size=(24, 24))
            ax = self.next_panel()
            corr = np.corrcoef(
                np.concatenate((cmap[:, :, 0].transpose(), cmap[:, :, 1].transpose(), cmap[:, :, 2].transpose()),
                               axis=0))
            plt.imshow(corr, aspect='auto', vmin=-1, vmax=1)
            plt.colorbar()
            starts3 = np.concatenate((np.array(starts), np.array(starts) + ends[-1], np.array(starts) + 2 * ends[-1]))
            ends3 = np.concatenate((np.array(ends), np.array(ends) + ends[-1], np.array(ends) + 2 * ends[-1]))
            for i in ends3[:-1]:
                plt.axvline(x=i - 0.5, color='red', linewidth=0.5)
                plt.axhline(y=i - 0.5, color='red', linewidth=0.5)

            ax.set_xticks((starts3 + ends3) / 2)
            ax.set_yticks((starts3 + ends3) / 2)
            chroms = list(map(Genome.canonical_chrom_name, chroms))
            ax.set_xticklabels(regions.split(",") + regions.split(",") + regions.split(","))
            ax.set_yticklabels(regions.split(",") + regions.split(",") + regions.split(","))
            self.fig_show(suffix="callmap")
        else:
            self.new_figure(panel_count=2, panel_size=(12, 12))
            ax = self.next_panel()
            x = np.concatenate((cmap[:, :, 0], cmap[:, :, 1], cmap[:, :, 2]),
                               axis=1)
            corr = np.corrcoef(x)
            plt.imshow(corr, aspect='auto', vmin=-1, vmax=1)
            plt.colorbar()
            ax = plt.gca()

            ax.set_xticks(range(n))
            ax.set_yticks(range(n))
            ax = self.next_panel()
            Z = hierarchy.linkage(x, 'average', 'correlation')
            dn = hierarchy.dendrogram(Z)

            self.fig_show(suffix="callmap")
        return cmap

    def multiple_regions(self, regions):
        n = len(self.plot_files) * len(regions)
        self.new_figure(panel_count=n)
        j = 0
        for i in range(len(self.plot_files)):
            for r in regions:
                self.regions(self.plot_files[i], r)
                j += 1
        self.fig_show(suffix="regions")

    def regions(self, ix, region):
        panels = self.panels
        bin_size = self.bin_size
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
            FLAG_USEHAP if self.snp_use_phase else 0)
        self.new_subgrid(len(panels), hspace=0.05, wspace=0.1)
        r = decode_region(region, max_size=1000000000)
        io = self.io[ix]
        for i in range(len(panels)):
            ax = self.next_subpanel(sharex=True)
            if i == 0 and self.title:
                ax.set_title(self.file_title(ix) + ": " + region.replace(",", ", "), position=(0.01, 0.93),
                             fontdict={'verticalalignment': 'bottom', 'horizontalalignment': 'left'},
                             color='C0')

            if panels[i] == "rd":
                g_p = [0]
                g_p_corr = [0]
                g_p_seg = [0]
                g_p_call = [0]
                g_p_call_mosaic = [0]
                g_p_call_mosaic_2d = [0]
                mean, stdev = 0, 0
                borders = []
                pos_x = []
                for c, (pos1, pos2) in r:
                    if pos2 == 1000000000:
                        pos2 = io.get_chromosome_length(c)
                        if pos2 is None:
                            pos2 = 1000000000
                    flag_rd = 0
                    if self.rd_use_mask:
                        flag_rd = FLAG_USEMASK
                    mean, stdev = io.rd_normal_level(bin_size, flag_rd | FLAG_GC_CORR)
                    his_p = io.get_signal(c, bin_size, "RD", flag_rd)
                    his_p_corr = io.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
                    his_p_seg = io.get_signal(c, bin_size, "RD partition", flag_rd | FLAG_GC_CORR)
                    his_p_call = io.get_signal(c, bin_size, "RD call", flag_rd | FLAG_GC_CORR)
                    his_p_mosaic_seg = io.get_signal(c, bin_size, "RD mosaic segments",
                                                     flag_rd | FLAG_GC_CORR)
                    his_p_mosaic_seg = segments_decode(his_p_mosaic_seg)
                    his_p_mosaic_call = io.get_signal(c, bin_size, "RD mosaic call",
                                                      flag_rd | FLAG_GC_CORR)
                    if self.snp_use_phase:
                        his_p_mosaic_seg_2d = io.get_signal(c, bin_size, "RD mosaic segments 2d phased",
                                                            flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_call_2d = io.get_signal(c, bin_size, "RD mosaic call 2d phased",
                                                             flag_rd | FLAG_GC_CORR)
                    else:
                        his_p_mosaic_seg_2d = io.get_signal(c, bin_size, "RD mosaic segments 2d",
                                                            flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_call_2d = io.get_signal(c, bin_size, "RD mosaic call 2d",
                                                             flag_rd | FLAG_GC_CORR)
                    his_p_mosaic_seg_2d = segments_decode(his_p_mosaic_seg_2d)
                    his_p_mosaic = np.zeros_like(his_p) * np.nan
                    if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and ("rd_mosaic" in self.callers):
                        for seg, lev in zip(list(his_p_mosaic_seg), list(his_p_mosaic_call[0])):
                            for segi in seg:
                                his_p_mosaic[segi] = lev
                    his_p_mosaic_2d = np.zeros_like(his_p) * np.nan
                    if his_p_mosaic_call_2d is not None and len(his_p_mosaic_call_2d) > 0 and (
                            "combined_mosaic" in self.callers):
                        for seg, lev in zip(list(his_p_mosaic_seg_2d), list(his_p_mosaic_call_2d[0])):
                            for segi in seg:
                                his_p_mosaic_2d[segi] = lev
                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    bins = len(list(his_p[start_bin:end_bin]))
                    pos_x.extend(range(pos1, pos2 + bin_size, bin_size)[0:bins])

                    g_p.extend(list(his_p[start_bin:end_bin]))
                    g_p_corr.extend(list(his_p_corr[start_bin:end_bin]))
                    if his_p_seg is not None and len(his_p_seg) > 0 and self.rd_partition:
                        g_p_seg.extend(list(his_p_seg[start_bin:end_bin]))
                    if his_p_call is not None and len(his_p_call) > 0 and self.rd_call and (
                            "rd_mean_shift" in self.callers):
                        g_p_call.extend(list(his_p_call[start_bin:end_bin]))
                    if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call and (
                            "rd_mosaic" in self.callers):
                        g_p_call_mosaic.extend(list(his_p_mosaic[start_bin:end_bin]))
                    if his_p_mosaic_call_2d is not None and len(his_p_mosaic_call_2d) > 0 and self.rd_call and (
                            "combined_mosaic" in self.callers):
                        g_p_call_mosaic_2d.extend(list(his_p_mosaic_2d[start_bin:end_bin]))
                    borders.append(len(g_p) - 1)

                def format_func(value, tick_number):
                    ix = int(value)
                    if ix + 1 < len(pos_x):
                        p = pos_x[ix] + (pos_x[ix + 1] - pos_x[ix]) * (value - ix)
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    elif ix < len(pos_x):
                        p = pos_x[ix]
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    else:
                        return ""

                l = len(g_p)
                if i == len(panels) - 1:
                    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                    ax.set_xlim([-l * 0.0, (l - 1) * 1.0])
                    ax.xaxis.grid()
                else:
                    plt.setp(ax.get_xticklabels(), visible=False)

                if (self.rd_range[1] - self.rd_range[0]) < 30:
                    ax.yaxis.set_ticks(np.arange(int(self.rd_range[0]), int(self.rd_range[1] + 1), 1) * mean / 2,
                                       minor=[])
                    ax.yaxis.set_ticklabels([str(i) for i in range(int(self.rd_range[0]), int(self.rd_range[1] + 1))])
                ax.set_ylim([self.rd_range[0] * mean / 2, self.rd_range[1] * mean / 2])
                ax.set_ylabel("Read depth")
                ax.yaxis.grid()

                if self.rd_raw:
                    ax.step(g_p, self.rd_colors[0], label="raw")
                if self.rd_corrected:
                    ax.step(g_p_corr, self.rd_colors[1], label="GC corrected")
                if len(g_p_seg) > 1:
                    plt.step(g_p_seg, self.rd_colors[2], label="partitioning")
                if len(g_p_call) > 1:
                    plt.step(g_p_call, self.rd_colors[3], label="cnv calls")
                if len(g_p_call_mosaic) > 1:
                    plt.step(g_p_call_mosaic, self.rd_colors[4], label="mosaic cnv calls")
                if len(g_p_call_mosaic_2d) > 1:
                    plt.step(g_p_call_mosaic_2d, self.rd_colors[5], label="combined cnv calls")
                for i in borders[:-1]:
                    ax.axvline(i, color="g", lw=1)
                if self.legend:
                    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=2)
                self.fig.add_subplot(ax)

            elif panels[i] == "snp":
                borders = []
                hpos = []
                baf = []
                color = []
                alpha = 0.7
                start_pos = 0
                pos_x = []
                for c, (pos1, pos2) in r:
                    if pos2 == 1000000000:
                        pos2 = io.get_chromosome_length(c)
                        if pos2 is None:
                            pos2 = 1000000000
                    pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c)
                    ix = 0
                    mdp = 0
                    while ix < len(pos) and pos[ix] <= pos2:
                        if pos[ix] >= pos1 and (nref[ix] + nalt[ix]) != 0 and ((not self.snp_use_id) or (flag[ix] & 1)):
                            hpos.append((start_pos + pos[ix] - pos1) / bin_size)
                            if pos[ix] - pos1 > mdp:
                                mdp = pos[ix] - pos1
                            if gt[ix] != 5:
                                baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            else:
                                baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                            if self.snp_alpha_P:
                                alpha = None
                                color.append(
                                    colors.to_rgba(self.snp_colors[(gt[ix] % 4) * 2 + 1], (flag[ix] >> 1) * 0.4))
                            else:
                                color.append(self.snp_colors[(gt[ix] % 4) * 2 + (flag[ix] >> 1)])
                        ix += 1
                    start_pos += pos2 - pos1
                    pos_x.extend(range(pos1, pos2 + bin_size, bin_size))
                    borders.append(start_pos / bin_size)

                def format_func(value, tick_number):
                    ix = int(value)
                    if ix + 1 < len(pos_x):
                        p = pos_x[ix] + (pos_x[ix + 1] - pos_x[ix]) * (value - ix)
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    elif ix < len(pos_x):
                        p = pos_x[ix]
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    else:
                        return ""

                l = len(pos_x)
                if i == len(panels) - 1:
                    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                    ax.set_xlim([-l * 0.0, (l - 1) * 1.0])
                    ax.xaxis.grid()
                else:
                    plt.setp(ax.get_xticklabels(), visible=False)

                # ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], minor=[])
                ax.yaxis.set_ticklabels(["0", "1/4", "1/2", "3/4", "1"])
                ax.set_ylabel("Allele frequency")
                # l = max(hpos)
                ax.set_ylim([-0.05, 1.05])
                # ax.set_xlim([0, borders[-1]])
                ax.yaxis.grid()
                if self.markersize == "auto":
                    ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=alpha)
                else:
                    ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=alpha)

                for i in borders[:-1]:
                    ax.axvline(i, color="g", lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "snv" or panels[i][:4] == "snv:":
                callset = "default"
                if panels[i][:4] == "snv:":
                    callset = panels[i].split(":")[1]
                borders = []
                hpos = []
                baf = []
                color = []
                alpha = 0.7
                start_pos = 0
                pos_x = []
                for c, (pos1, pos2) in r:
                    if pos2 == 1000000000:
                        pos2 = io.get_chromosome_length(c)
                        if pos2 is None:
                            pos2 = 1000000000
                    pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c, callset=callset)
                    ix = 0
                    mdp = 0
                    while ix < len(pos) and pos[ix] <= pos2:
                        if pos[ix] >= pos1 and (nref[ix] + nalt[ix]) != 0:
                            hpos.append((start_pos + pos[ix] - pos1) / bin_size)
                            if pos[ix] - pos1 > mdp:
                                mdp = pos[ix] - pos1
                            if gt[ix] % 4 != 2:
                                baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            else:
                                baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                            if self.snp_alpha_P:
                                alpha = None
                                color.append(
                                    colors.to_rgba(self.snp_colors[(gt[ix] % 4) * 2 + 1], (flag[ix] >> 1) * 0.4))
                            else:
                                color.append(self.snp_colors[(gt[ix] % 4) * 2 + (flag[ix] >> 1)])
                        ix += 1
                    start_pos += pos2 - pos1
                    pos_x.extend(range(pos1, pos2 + bin_size, bin_size))
                    borders.append(start_pos / bin_size)

                def format_func(value, tick_number):
                    ix = int(value)
                    if ix + 1 < len(pos_x):
                        p = pos_x[ix] + (pos_x[ix + 1] - pos_x[ix]) * (value - ix)
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    elif ix < len(pos_x):
                        p = pos_x[ix]
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    else:
                        return ""

                l = len(pos_x)
                if i == len(panels) - 1:
                    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                    ax.set_xlim([-l * 0.0, (l - 1) * 1.0])
                else:
                    plt.setp(ax.get_xticklabels(), visible=False)
                ax.xaxis.grid()
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], minor=[])
                ax.yaxis.set_ticklabels(["0", "1/4", "1/2", "3/4", "1"])
                ax.set_ylabel("Allele frequency")
                ax.set_ylim([0., 1.])
                ax.yaxis.grid()
                if self.markersize == "auto":
                    ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=alpha)
                else:
                    ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=alpha)

                for i in borders[:-1]:
                    ax.axvline(i, color="g", lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "baf":
                g_baf, g_maf, g_i1, g_i2 = [0], [0], [0], [0]
                borders = []
                pos_x = []

                for c, (pos1, pos2) in r:
                    if pos2 == 1000000000:
                        pos2 = io.get_chromosome_length(c)
                        if pos2 is None:
                            pos2 = 1000000000

                    flag_snp = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                        FLAG_USEHAP if self.snp_use_phase else 0)
                    baf = io.get_signal(c, bin_size, "SNP baf", flag_snp)
                    maf = io.get_signal(c, bin_size, "SNP maf", flag_snp)
                    i1 = io.get_signal(c, bin_size, "SNP i1", flag_snp)
                    i2 = io.get_signal(c, bin_size, "SNP i2", flag_snp)

                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    bins = len(list(baf[start_bin:end_bin]))
                    pos_x.extend(range(pos1, pos2 + bin_size, bin_size)[0:bins])

                    g_baf.extend(list(baf[start_bin:end_bin]))
                    g_maf.extend(list(maf[start_bin:end_bin]))
                    g_i1.extend(list(i1[start_bin:end_bin]))
                    g_i2.extend(list(i2[start_bin:end_bin]))
                    borders.append(len(g_baf) - 1)

                def format_func(value, tick_number):
                    ix = int(value)
                    if ix + 1 < len(pos_x):
                        p = pos_x[ix] + (pos_x[ix + 1] - pos_x[ix]) * (value - ix)
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    elif ix < len(pos_x):
                        p = pos_x[ix]
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    else:
                        return ""

                l = len(g_baf)
                if i == len(panels) - 1:
                    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                    ax.set_xlim([-l * 0.0, (l - 1) * 1.0])
                    ax.xaxis.grid()

                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], minor=[])
                ax.yaxis.set_ticklabels(["0", "1/4", "1/2", "3/4", "1"])
                ax.set_ylabel("Allele frequency")

                ax.set_ylim(self.baf_range)
                # ax.set_xlim([-l * 0.0, l * 1.0])

                ax.yaxis.grid()
                # ax.xaxis.grid()
                if self.plot_baf:
                    ax.step(g_baf, self.baf_colors[0], label="BAF")
                if self.plot_maf:
                    ax.step(g_maf, self.baf_colors[1], label="MAF")
                if self.plot_dbaf:
                    ax.step(g_i1, self.baf_colors[2], label="I1")
                if self.legend:
                    ax.legend()
                for i in borders[:-1]:
                    ax.axvline(i, color="g", lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "likelihood":
                borders = []
                gl = []
                call_pos = []
                call_i1 = []
                call_i2 = []
                call_c = []
                call_pos_2d = []
                call_i1_2d = []
                call_i2_2d = []
                call_c_2d = []
                tlen = 0
                tlen_2d = 0
                pos_x = []
                for c, (pos1, pos2) in r:
                    if pos2 == 1000000000:
                        pos2 = io.get_chromosome_length(c)
                        if pos2 is None:
                            pos2 = 1000000000
                    likelihood = io.get_signal(c, bin_size, "SNP likelihood", snp_flag)
                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    bins = len(list(likelihood[start_bin:end_bin]))
                    pos_x.extend(range(pos1, pos2 + bin_size, bin_size)[0:bins])
                    gl.extend(list(likelihood[start_bin:end_bin]))
                    borders.append(len(gl) - 1)
                    if self.snp_call and ("baf_mosaic" in self.callers):
                        likelihood_call = io.get_signal(c, bin_size, "SNP likelihood call", snp_flag)
                        segments = segments_decode(io.get_signal(c, bin_size, "SNP likelihood segments", snp_flag))

                        for s, lh in zip(segments, likelihood_call):
                            i1, i2, p = likelihood_pixels_pval(lh)
                            if i1 != i2 and len(s) > self.min_segment_size:
                                alpha = -np.log(p + 1e-40) / self.contrast
                                if alpha > 1:
                                    alpha = 1
                                for pos in s:
                                    if pos >= start_bin and pos < end_bin:
                                        call_pos.append(pos - start_bin + tlen)
                                        call_i1.append(min(i1, i2))
                                        call_i2.append(max(i1, i2))
                                        color = colors.to_rgb(self.lh_colors[0]) + (alpha,)
                                        call_c.append(color)
                        tlen += end_bin - start_bin
                    if self.snp_call and ("combined_mosaic" in self.callers):
                        likelihood_call = io.get_signal(c, bin_size, "SNP likelihood call 2d", snp_flag)
                        segments = segments_decode(io.get_signal(c, bin_size, "SNP likelihood segments 2d", snp_flag))

                        for s, lh in zip(segments, likelihood_call):
                            i1, i2, p = likelihood_pixels_pval(lh)
                            if i1 != i2 and len(s) > self.min_segment_size:
                                alpha = -np.log(p + 1e-40) / self.contrast
                                if alpha > 1:
                                    alpha = 1
                                for pos in s:
                                    if pos >= start_bin and pos < end_bin:
                                        call_pos_2d.append(pos - start_bin + tlen_2d)
                                        call_i1_2d.append(min(i1, i2))
                                        call_i2_2d.append(max(i1, i2))
                                        color = colors.to_rgb(self.lh_colors[1]) + (alpha,)
                                        call_c_2d.append(color)
                        tlen_2d += end_bin - start_bin

                def format_func(value, tick_number):
                    ix = int(value)
                    if ix + 1 < len(pos_x):
                        p = pos_x[ix] + (pos_x[ix + 1] - pos_x[ix]) * (value - ix)
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    elif ix < len(pos_x):
                        p = pos_x[ix]
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    else:
                        return ""

                img = np.array(gl).transpose()
                l = img.shape[1]
                if i == len(panels) - 1:
                    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                    ax.set_xlim([-l * 0.0, (l - 1) * 1.0])
                    # ax.xaxis.grid()
                else:
                    plt.setp(ax.get_xticklabels(), visible=False)

                ax.imshow(img, aspect='auto', extent=[0, l, 0, img.shape[0] - 1])
                # ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticks([0, img.shape[0] / 4, img.shape[0] / 2, 3 * img.shape[0] / 4, img.shape[0] - 1],
                                   minor=[])
                ax.yaxis.set_ticklabels(["1", "3/4", "1/2", "1/4", "0"])
                ax.set_ylabel("Allele frequency")
                # ax.xaxis.set_ticks(np.arange(0, len(gl), 50), minor=[])
                # ax.set_xlim([-0.5, img.shape[1] - 0.5])
                ax.set_ylim([self.baf_range[0] * img.shape[0], self.baf_range[1] * img.shape[0]])
                if self.snp_call and ("baf_mosaic" in self.callers):
                    plt.scatter(call_pos, call_i1, s=self.lh_markersize, color=np.array(call_c), edgecolors='face',
                                marker=self.lh_marker)
                    plt.scatter(call_pos, call_i2, s=self.lh_markersize, color=np.array(call_c), edgecolors='face',
                                marker=self.lh_marker)
                if self.snp_call and ("combined_mosaic" in self.callers):
                    plt.scatter(call_pos_2d, call_i1_2d, s=self.lh_markersize, color=np.array(call_c_2d),
                                edgecolors='face', marker=self.lh_marker)
                    plt.scatter(call_pos_2d, call_i2_2d, s=self.lh_markersize, color=np.array(call_c_2d),
                                edgecolors='face', marker=self.lh_marker)

                for i in borders[:-1]:
                    ax.axvline(i + 1, color="g", lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "CN":
                borders = []
                gh1 = []
                gh2 = []
                tlen = 0
                tlen_2d = 0
                for c, (pos1, pos2) in r:
                    if pos2 == 1000000000:
                        pos2 = io.get_chromosome_length(c)
                        if pos2 is None:
                            pos2 = 1000000000

                    his_p = io.get_signal(c, bin_size, "RD p", flag_rd)
                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    if end_bin > len(his_p):
                        end_bin = len(his_p)
                    h1 = np.array([0] * (end_bin - start_bin))
                    h2 = np.array([0] * (end_bin - start_bin))
                    h1[his_p[start_bin:end_bin] != 0] = 1
                    h2[his_p[start_bin:end_bin] != 0] = 1

                    flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                        FLAG_USEHAP if self.snp_use_phase else 0) | (
                               FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR
                    flag_rd = FLAG_GC_CORR | (FLAG_USEMASK if self.rd_use_mask else 0)
                    if io.signal_exists(c, bin_size, "calls combined", flag):
                        calls = io.read_calls(c, bin_size, "calls combined", flag)
                        segments = io.get_signal(c, bin_size, "RD mosaic segments 2d phased", FLAG_GC_CORR)
                        print(segments)
                        segments = segments_decode(segments)

                        for call in calls:
                            for b in segments[int(call["segment"])]:
                                if b < end_bin and b >= start_bin:
                                    h1[b - start_bin] = call["models"][0][1]
                                    h2[b - start_bin] = call["models"][0][2]
                    gh1.extend(list(h1))
                    gh2.extend(list(h2))
                    borders.append(len(gh1) - 1)
                x = range(len(gh1))
                plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
                plt.stackplot(x, gh1, gh2, baseline='sym')

                def format_func(value, tick_number):
                    ix = int(value)
                    if ix + 1 < len(pos_x):
                        p = pos_x[ix] + (pos_x[ix + 1] - pos_x[ix]) * (value - ix)
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    elif ix < len(pos_x):
                        p = pos_x[ix]
                        return "{0} Mbp".format(int(p / 100) / 10000)
                    else:
                        return ""

                l = len(gh1)
                if i == len(panels) - 1:
                    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
                    ax.set_xlim([-l * 0.0, (l - 1) * 1.0])
                    ax.xaxis.grid()

                for i in borders[:-1]:
                    ax.axvline(i + 0.5, color="g", lw=1)
                self.fig.add_subplot(ax)

    def global_plot(self):
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_files[0]].rd_chromosome_name(c)
            if (len(self.chrom) == 0 or (rd_chr in self.chrom) or (c in self.chrom)) and rd_chr is not None:
                if (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                    chroms.append((rd_chr, l))
        panels = self.panels
        bin_size = self.bin_size
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
            FLAG_USEHAP if self.snp_use_phase else 0)
        rd_flag = (FLAG_USEMASK if self.rd_use_mask else 0) | (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
        n = len(self.plot_files)
        self.new_figure(panel_count=n)
        for ii in range(len(self.plot_files)):
            ix = self.plot_files[ii]
            self.new_subgrid(len(panels), hspace=0.05, wspace=0.05)
            io = self.io[ix]
            for i in range(len(panels)):
                ax = self.next_subpanel(sharex=True)
                if i == 0 and self.title:
                    ax.set_title(self.file_title(ix), position=(0.01, 0.9),
                                 fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                                 color='C0')

                if panels[i] == "rd":
                    start = 0
                    xticks = [0]
                    xticks_minor = []
                    xticks_labels = []
                    for c, l in chroms:
                        mean, stdev = io.rd_normal_level(bin_size, rd_flag | FLAG_GC_CORR)
                        his_p = io.get_signal(c, bin_size, "RD", rd_flag)
                        pos = range(start, start + len(his_p))
                        if self.markersize == "auto":
                            plt.plot(pos, his_p, ls='', marker='.', markersize=1)
                        else:
                            plt.plot(pos, his_p, ls='', marker='.', markersize=self.markersize)
                        xticks_minor.append(start + len(his_p) // 2)
                        xticks_labels.append(Genome.canonical_chrom_name(c))
                        start += l // bin_size + 1
                        xticks.append(start)

                    ax.set_xlim([0, start])
                    ax.xaxis.set_ticks(xticks)
                    ax.xaxis.set_ticklabels([""] * len(xticks))
                    if i == (len(panels) - 1):
                        ax.xaxis.set_ticks(xticks_minor, minor=True)
                        ax.xaxis.set_ticklabels(xticks_labels, minor=True)
                    else:
                        plt.setp(ax.get_xticklabels(which="both"), visible=False)
                    yticks = np.arange(self.rd_manhattan_range[0], self.rd_manhattan_range[1], 0.5)
                    ax.yaxis.set_ticklabels([str(int(2 * t)) for t in yticks])
                    ax.yaxis.set_ticks(yticks * mean)
                    ax.set_ylabel("RD [CN]")
                    ax.set_ylim([self.rd_manhattan_range[0] * mean, self.rd_manhattan_range[1] * mean])
                    ax.grid()
                    self.fig.add_subplot(ax)

                elif panels[i] == "snp":
                    start = 0
                    xticks = []
                    xticks_minor = []
                    xticks_labels = []
                    pos_x = []
                    for c, l in chroms:
                        pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c)
                        ix = 0
                        hpos = []
                        color = []
                        alpha = 0.7
                        baf = []
                        while ix < len(pos):
                            if (nref[ix] + nalt[ix]) != 0 and ((not self.snp_use_id) or (flag[ix] & 1)):
                                hpos.append(start + (pos[ix] / bin_size))
                                if gt[ix] % 4 != 2:
                                    baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                                else:
                                    baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                                if self.snp_alpha_P:
                                    alpha = None
                                    color.append(
                                        colors.to_rgba(self.snp_colors[(gt[ix] % 4) * 2 + 1], (flag[ix] >> 1) * 0.4))
                                else:
                                    color.append(self.snp_colors[(gt[ix] % 4) * 2 + (flag[ix] >> 1)])
                            ix += 1
                        if self.markersize == "auto":
                            ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=0.1, alpha=alpha)
                        else:
                            ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=alpha)
                        xticks_minor.append(start + l // bin_size // 2)
                        xticks_labels.append(Genome.canonical_chrom_name(c))
                        start += l // bin_size + 1
                        xticks.append(start)
                    ax.set_xlim([0, start])
                    ax.xaxis.set_ticks(xticks)
                    ax.xaxis.set_ticklabels([""] * len(xticks))
                    if i == (len(panels) - 1):
                        ax.xaxis.set_ticks(xticks_minor, minor=True)
                        ax.xaxis.set_ticklabels(xticks_labels, minor=True)
                    else:
                        plt.setp(ax.get_xticklabels(minor=True), visible=False)
                    ax.grid()
                    ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
                    ax.yaxis.set_ticklabels(["0", "1/4", "1/2", "3/4", "1"])
                    ax.set_ylabel("BAF")
                    ax.set_ylim([-0.05, 1.05])
                    ax.yaxis.grid()
                    self.fig.add_subplot(ax)

                elif panels[i] == "snv" or panels[i][:4] == "snv:":
                    callset = "default"
                    if panels[i][:4] == "snv:":
                        callset = panels[i].split(":")[1]
                    start = 0
                    xticks = []
                    xticks_minor = []
                    xticks_labels = []
                    pos_x = []
                    for c, l in chroms:
                        pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c, callset=callset)
                        ix = 0
                        hpos = []
                        color = []
                        alpha = 0.7
                        baf = []
                        while ix < len(pos):
                            if (nref[ix] + nalt[ix]) != 0 and ((not self.snp_use_id) or (flag[ix] & 1)):
                                hpos.append(start + (pos[ix] / bin_size))
                                if gt[ix] % 4 != 2:
                                    baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                                else:
                                    baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                                if self.snp_alpha_P:
                                    alpha = None
                                    color.append(
                                        colors.to_rgba(self.snp_colors[(gt[ix] % 4) * 2 + 1], (flag[ix] >> 1) * 0.4))
                                else:
                                    color.append(self.snp_colors[(gt[ix] % 4) * 2 + (flag[ix] >> 1)])
                            ix += 1
                        if self.markersize == "auto":
                            ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=0.1, alpha=alpha)
                        else:
                            ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=alpha)
                        xticks_minor.append(start + l // bin_size // 2)
                        xticks_labels.append(Genome.canonical_chrom_name(c))
                        start += l // bin_size + 1
                        xticks.append(start)
                    ax.set_xlim([0, start])
                    ax.xaxis.set_ticks(xticks)
                    ax.xaxis.set_ticklabels([""] * len(xticks))
                    if i == (len(panels) - 1):
                        ax.xaxis.set_ticks(xticks_minor, minor=True)
                        ax.xaxis.set_ticklabels(xticks_labels, minor=True)
                    else:
                        plt.setp(ax.get_xticklabels(minor=True), visible=False)
                    ax.grid()
                    ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
                    ax.yaxis.set_ticklabels(["0", "1/4", "1/2", "3/4", "1"])
                    ax.set_ylabel("BAF")
                    ax.set_ylim([-0.05, 1.05])
                    ax.yaxis.grid()
                    self.fig.add_subplot(ax)


                elif panels[i] == "likelihood":
                    start = 0
                    xticks = [0]
                    xticks_minor = []
                    xticks_labels = []
                    gl = []
                    for c, l in chroms:
                        likelihood = io.get_signal(c, bin_size, "SNP likelihood", snp_flag)
                        lh = list(likelihood)
                        size = l // bin_size + 1
                        if len(lh) < size:
                            if len(lh) > 0:
                                lh.extend([lh[-1] for jj in range(size - len(lh))])
                            elif len(gl) > 0:
                                lh.extend([gl[-1] for jj in range(size - len(lh))])

                        gl.extend(lh)
                        xticks_minor.append(start + l // bin_size // 2)
                        xticks_labels.append(Genome.canonical_chrom_name(c))
                        start += l // bin_size + 1
                        xticks.append(start)

                    img = np.array(gl).transpose()
                    img[0, :] = 0
                    img[-1, :] = 0
                    ax.imshow(img, aspect='auto')
                    ax.yaxis.set_ticks([0, img.shape[0] / 4, img.shape[0] / 2, 3 * img.shape[0] / 4, img.shape[0] - 1],
                                       minor=[])
                    ax.yaxis.set_ticklabels(["1", "3/4", "1/2", "1/4", "0"])
                    ax.set_ylabel("BAF")
                    ax.set_xlim([0, start])
                    ax.xaxis.set_ticks(xticks)
                    ax.xaxis.set_ticklabels([""] * len(xticks))
                    if i == (len(panels) - 1):
                        ax.xaxis.set_ticks(xticks_minor, minor=True)
                        ax.xaxis.set_ticklabels(xticks_labels, minor=True)
                    else:
                        plt.setp(ax.get_xticklabels(minor=True), visible=False)
                    ax.xaxis.grid()
                    self.fig.add_subplot(ax)

        self.fig_show(suffix="global")

    def circular(self):
        chroms = self.chrom
        bin_size = self.bin_size
        n = len(self.plot_files)
        ix = self.plot_files
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
            FLAG_USEHAP if self.snp_use_phase else 0)
        rd_flag = FLAG_GC_CORR | (FLAG_USEMASK if self.rd_use_mask else 0)
        self.new_figure(panel_count=n)
        for i in range(n):
            ax = self.next_polar_panel()
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            rainbow = ax._get_lines
            io = self.io[ix[i]]
            plot_len = 0
            plot_chroms = []
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                rd_chr = io.rd_chromosome_name(c)
                if rd_chr is not None and (len(chroms) == 0 or (rd_chr in chroms) or (c in chroms)) and (
                        Genome.is_autosome(c) or Genome.is_sex_chrom(c)
                ) and io.signal_exists(rd_chr, bin_size, "SNP maf", snp_flag) and io.signal_exists(
                    rd_chr, bin_size, "RD", rd_flag):
                    plot_chroms.append((rd_chr, l))
                    plot_len += l // bin_size + 1
            rd_mean, stdev = io.rd_normal_level(bin_size, rd_flag)
            tl = 0
            dt = 2.0 * np.pi / plot_len
            theta = np.arange(0, 2.0 * np.pi, dt)
            angles = []
            labels = []
            for j in range(len(plot_chroms)):
                c, l = plot_chroms[j]
                next_color = rainbow.get_next_color()
                rd_color = self.rd_circular_colors[j % len(self.rd_circular_colors)]
                snp_color = self.snp_circular_colors[j % len(self.snp_circular_colors)]
                rd = io.get_signal(c, bin_size, "RD", rd_flag)
                maf = io.get_signal(c, bin_size, "SNP maf", snp_flag)
                c01 = io.get_signal(c, bin_size, "SNP bin count 0|1", snp_flag)
                c10 = io.get_signal(c, bin_size, "SNP bin count 1|0", snp_flag)
                hets = c01 + c10
                np.warnings.filterwarnings('ignore')
                maf[hets < (bin_size / 10000)] = 0
                # plt.polar(theta[tl:tl + maf.size], 1 - maf / 2, color=snp_color, linewidth=0.3)
                # plt.fill_between(theta[tl:tl + maf.size], 1 - maf / 2, np.ones_like(maf), color=snp_color, alpha=0.8)
                plt.polar(theta[tl:tl + maf.size], 1 - maf / 2, linewidth=0.3, color=next_color)
                plt.fill_between(theta[tl:tl + maf.size], 1 - maf / 2, np.ones_like(maf), alpha=1, color=next_color)
                markersize = 5
                if self.markersize != "auto":
                    markersize = self.markersize
                ax.scatter(theta[tl:tl + rd.size], np.ones_like(rd) / 10. + 0.7 * rd / (self.rd_range[1] * rd_mean / 2),
                           s=markersize, alpha=0.7, color=next_color)

                # plt.polar(theta[tl:tl + rd.size], np.ones_like(rd) / 10. + 0.7 * rd / (self.rd_range[1] * rd_mean / 2),
                #          color=rd_color, linewidth=0.3)
                # plt.fill_between(theta[tl:tl + rd.size], np.ones_like(rd) / 10.,
                #                 np.ones_like(rd) / 10. + 0.7 * rd / (self.rd_range[1] * rd_mean / 2),
                #                 color=rd_color,
                #                 alpha=0.8)

                # ax.text(theta[tl + maf.size // 3], 0.8, c, fontsize=8)
                labels.append(Genome.canonical_chrom_name(c))
                angles.append(180 * theta[tl + rd.size // 2] / np.pi)
                tl += l // bin_size + 1
            for cn in range(int(self.rd_range[1])):
                plt.polar(theta, np.ones_like(theta) * (0.1 + 0.7 * (cn / self.rd_range[1])), color="k", linewidth=0.1)
            ax.set_rmax(1.0)
            ax.set_rticks([])
            ax.set_thetagrids(angles, labels=labels, fontsize=10, weight="bold", color="black")
            if self.title:
                ax.set_title(self.file_title(ix[i]), loc="left", fontsize=10, weight="bold", color="black")
            ax.grid(False)
        self.fig_show(suffix="circular")

    def rd_baf(self, hist=True):
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(12, 8), facecolor='w', edgecolor='k')
        n = len(self.plot_files)
        ix = self.plot_files
        if self.grid == "auto":
            sx, sy = self._panels_shape(n)
        else:
            sx, sy = tuple(self.grid)
        grid = gridspec.GridSpec(sy, sx, wspace=0.2, hspace=0.2)
        bin_size = self.bin_size
        for i in range(n):
            ax = self.fig.add_subplot(grid[i])
            io = self.io[ix[i]]

            chroms = []
            snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                FLAG_USEHAP if self.snp_use_phase else 0)
            rd_flag = FLAG_GC_CORR | (FLAG_USEMASK if self.rd_use_mask else 0)
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = io.snp_chromosome_name(c)
                if len(self.chrom) == 0 or (snp_chr in self.chrom) or (c in self.chrom):
                    if io.signal_exists(snp_chr, bin_size, "SNP likelihood call", snp_flag) and \
                            io.signal_exists(snp_chr, bin_size, "SNP likelihood segments", snp_flag) and \
                            io.signal_exists(snp_chr, bin_size, "RD mosaic call", rd_flag) and \
                            io.signal_exists(snp_chr, bin_size, "RD mosaic segments", rd_flag) and \
                            Genome.is_autosome(c):
                        chroms.append((snp_chr, l))
            x = []
            y = []
            for c, l in chroms:
                flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO

                likelihood = io.get_signal(c, bin_size, "SNP likelihood call", snp_flag)
                segments_baf = segments_decode(io.get_signal(c, bin_size, "SNP likelihood segments", snp_flag))
                rd = io.get_signal(c, bin_size, "RD mosaic call", rd_flag)
                segments_rd = segments_decode(io.get_signal(c, bin_size, "RD mosaic segments", rd_flag))

                mbaf = {}
                mrd = {}
                for s, lh in zip(segments_baf, likelihood):
                    b, p = likelihood_baf_pval(lh)
                    for pos in s:
                        mbaf[pos] = 0.5 - b
                for s, r in zip(segments_rd, rd[0]):
                    for pos in s:
                        mrd[pos] = r
                for p in mbaf:
                    if p in mrd:
                        x.append(mbaf[p])
                        y.append(mrd[p])

            if hist:
                from matplotlib.colors import LogNorm
                ax.hist2d(x, y, bins=[np.arange(0, 0.51, 0.01), np.arange(0, max(y), max(y) / 100.)], norm=LogNorm())
            else:
                ax.scatter(x, y, marker=".", alpha=0.5)

        if self.output_filename != "":
            plt.savefig(self._image_filename("rd_baf"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def dispersion(self, legend=True):
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(8)
            self.fig.set_figwidth(12)
        grid = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0.2)

        ax = self.fig.add_subplot(grid[0])
        for i in self.io:
            bin_sizes = sorted(set([int(x[1]) for x in i.chromosomes_bin_sizes_with_signal("RD")]))
            rd = []
            drd = []
            for bs in bin_sizes:
                if i.signal_exists(None, bs, "RD stat", flags=FLAG_AUTO):
                    stat = i.get_signal(None, bs, "RD stat", flags=FLAG_AUTO)
                    rd.append(stat[4])
                    drd.append(stat[5])
            ax.set_yscale("log")
            ax.set_xscale("log")
            ax.grid(True)
            ax.set_xlabel("mean RD")
            ax.set_ylabel("stdev RD")
            if legend:
                ax.legend(loc="upper left")
            ax.plot(rd, drd, "*-", label=i.filename)

        ax = self.fig.add_subplot(grid[1])
        for i in self.io:
            bin_sizes = sorted(set([int(x[1]) for x in i.chromosomes_bin_sizes_with_signal("RD")]))
            rd = []
            drd = []
            for bs in bin_sizes:
                if i.signal_exists(None, bs, "RD stat", flags=FLAG_AUTO | FLAG_GC_CORR):
                    stat = i.get_signal(None, bs, "RD stat", flags=FLAG_AUTO | FLAG_GC_CORR)
                    rd.append(stat[4])
                    drd.append(stat[5])
            ax.set_yscale("log")
            ax.set_xscale("log")
            ax.grid(True)
            ax.set_xlabel("mean RD (GC corr)")
            ax.set_ylabel("stdev RD (GC corr)")
            if legend:
                ax.legend(loc="upper left")
            ax.plot(rd, drd, "*-", label=i.filename)

        if self.output_filename != "":
            plt.savefig(self._image_filename("dispersion"), dpi=200)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def region_rd_stat(self, region, n_bins=21, plot=False, legend=True):
        n = len(self.plot_files)
        ix = self.plot_files
        if plot:
            plt.clf()
            plt.rcParams["font.size"] = 8
            if self.grid == "auto":
                sx, sy = self._panels_shape(n)
            else:
                sx, sy = tuple(self.grid)
            self.fig = plt.figure(1, dpi=200, facecolor='w', edgecolor='k')
            if self.output_filename != "":
                self.fig.set_figheight(3 * sy)
                self.fig.set_figwidth(4 * sx)
            grid = gridspec.GridSpec(sy, sx, wspace=0.2, hspace=0.2)
        for i in range(n):
            io = self.io[ix[i]]
            if plot:
                ax = self.fig.add_subplot(grid[i])
                ax.set_title(self.file_title(ix[i]), position=(0.01, 1.07),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
            regs = decode_region(region)
            data = []
            for c, (pos1, pos2) in regs:
                flag_rd = (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | (FLAG_USEMASK if self.rd_use_mask else 0)
                his_p = io.get_signal(c, self.bin_size, "RD", flag_rd)
                bin1 = (pos1 - 1) // self.bin_size
                bin2 = (pos2 - 1) // self.bin_size
                data += list(his_p[bin1:bin2 + 1][np.isfinite(his_p[bin1:bin2 + 1])])

            data = np.array(data)
            dmin = np.min(data)
            dmax = np.max(data)
            p1 = np.percentile(data, 1)
            p99 = np.percentile(data, 99)
            data = data[data > p1]
            data = data[data < p99]
            mean = np.mean(data)
            std = np.std(data)

            rd_min = mean - 5 * std
            rd_max = mean + 5 * std
            bins = np.linspace(rd_min, rd_max, n_bins)

            hist, binsr = np.histogram(data, bins=bins)

            fitn, fitm, fits = fit_normal(bins[:-1], hist)[0]

            print("%s\t%s\t%.4f\t%.4f\t%e\t%e\t%.4f\t%.4f\t%.4f\t%.4f" % (
                io.filename, region, fitm, fits, dmin, dmax, p1, p99, mean, std))

            if plot:
                x = np.linspace(bins[0], bins[-1], 1001)
                plt.plot(x, normal(x, fitn, fitm, fits), "g-", label=region)
                plt.plot(bins[:-1], hist, "b*")
                if legend:
                    plt.legend()

        if plot:
            if self.output_filename != "":
                plt.savefig(self._image_filename("comp"), dpi=200)
                plt.close(self.fig)
            elif self.interactive:
                plt.show(block=False)
                plt.draw()
            else:
                plt.show()

    def compare(self, region1, region2, n_bins=21, plot=False, stdout=True, legend=True):
        n = len(self.plot_files)
        ix = self.plot_files
        ret = []

        if plot:
            plt.clf()
            plt.rcParams["font.size"] = 8
            if self.grid == "auto":
                sx, sy = self._panels_shape(n)
            else:
                sx, sy = tuple(self.grid)
            self.fig = plt.figure(1, dpi=200, facecolor='w', edgecolor='k')
            if self.output_filename != "":
                self.fig.set_figheight(3 * sy)
                self.fig.set_figwidth(4 * sx)
            grid = gridspec.GridSpec(sy, sx, wspace=0.2, hspace=0.2)
        for i in range(n):
            io = self.io[ix[i]]
            if plot:
                ax = self.fig.add_subplot(grid[i])
                ax.set_title(self.file_title(ix[i]), position=(0.01, 1.07),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
            regs1 = decode_region(region1)
            regs2 = decode_region(region2)
            data1 = []
            data2 = []
            for c, (pos1, pos2) in regs1:
                flag_rd = (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | (FLAG_USEMASK if self.rd_use_mask else 0)
                his_p = io.get_signal(c, self.bin_size, "RD", flag_rd)
                bin1 = (pos1 - 1) // self.bin_size
                bin2 = (pos2 - 1) // self.bin_size
                data1 += list(his_p[bin1:bin2 + 1][np.isfinite(his_p[bin1:bin2 + 1])])
            for c, (pos1, pos2) in regs2:
                flag_rd = (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | (FLAG_USEMASK if self.rd_use_mask else 0)
                his_p = io.get_signal(c, self.bin_size, "RD", flag_rd)
                bin1 = (pos1 - 1) // self.bin_size
                bin2 = (pos2 - 1) // self.bin_size
                data2 += list(his_p[bin1:bin2 + 1][np.isfinite(his_p[bin1:bin2 + 1])])

            data1 = np.array(data1)
            p1_1 = np.percentile(data1, 1)
            p99_1 = np.percentile(data1, 99)
            data1 = data1[data1 > p1_1]
            data1 = data1[data1 < p99_1]
            mean1 = np.mean(data1)
            std1 = np.std(data1)

            data2 = np.array(data2)
            p1_2 = np.percentile(data2, 1)
            p99_2 = np.percentile(data2, 99)
            data2 = data2[data2 > p1_2]
            data2 = data2[data2 < p99_2]
            mean2 = np.mean(data2)
            std2 = np.std(data2)

            rd_min = min(mean1 - 5 * std1, mean2 - 5 * std2)
            rd_max = max(mean1 + 5 * std1, mean2 + 5 * std2)
            bins = np.linspace(rd_min, rd_max, n_bins)

            hist1, binsr = np.histogram(data1, bins=bins)
            hist2, binsr = np.histogram(data2, bins=bins)

            fitn1, fitm1, fits1 = fit_normal(bins[:-1], hist1)[0]
            fitn2, fitm2, fits2 = fit_normal(bins[:-1], hist2)[0]

            pval = t_test_2_samples(fitm1, fits1, sum(hist1), fitm2, fits2, sum(hist2))

            if stdout:
                print("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%e\t%.4f\t%.4f" % (
                    io.filename, region1, region2, fitm1, fits1, fitm2, fits2, pval, fitm1 / fitm2,
                    fitm1 / fitm2 * (fits1 / fitm1 / np.sqrt(sum(hist1)) + fits2 / fitm2 / np.sqrt(sum(hist2)))))
            ret.append([io.filename, region1, region2, fitm1, fits1, fitm2, fits2, pval, fitm1 / fitm2,
                        fitm1 / fitm2 * (fits1 / fitm1 / np.sqrt(sum(hist1)) + fits2 / fitm2 / np.sqrt(sum(hist2)))])

            if plot:
                x = np.linspace(bins[0], bins[-1], 1001)
                plt.plot(x, normal(x, fitn1, fitm1, fits1), "g-", label=region1)
                plt.plot(x, normal(x, fitn2, fitm2, fits2), "b-", label=region2)
                plt.plot(bins[:-1], hist1, "g*")
                plt.plot(bins[:-1], hist2, "b*")
                if legend:
                    plt.legend()

        if plot:
            if self.output_filename != "":
                plt.savefig(self._image_filename("comp"), dpi=200)
                plt.close(self.fig)
            elif self.interactive:
                plt.show(block=False)
                plt.draw()
            else:
                plt.show()

        return ret

    def compare_baf(self, region1, region2, plot=False, stdout=True, legend=True):
        n = len(self.plot_files)
        ix = self.plot_files
        ret = []

        if plot:
            plt.clf()
            plt.rcParams["font.size"] = 8
            if self.grid == "auto":
                sx, sy = self._panels_shape(n)
            else:
                sx, sy = tuple(self.grid)
            self.fig = plt.figure(1, dpi=200, facecolor='w', edgecolor='k')
            if self.output_filename != "":
                self.fig.set_figheight(3 * sy)
                self.fig.set_figwidth(4 * sx)
            grid = gridspec.GridSpec(sy, sx, wspace=0.2, hspace=0.2)
        for i in range(n):
            io = self.io[ix[i]]
            if plot:
                ax = self.fig.add_subplot(grid[i])
                ax.set_title(self.file_title(ix[i]), position=(0.01, 1.07),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
            regs1 = decode_region(region1)
            regs2 = decode_region(region2)
            data1 = []
            data2 = []
            for c, (pos1, pos2) in regs1:
                flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
                his_p = io.get_signal(c, self.bin_size, "SNP likelihood", flag)
                bin1 = (pos1 - 1) // self.bin_size
                bin2 = (pos2 - 1) // self.bin_size
                data1 += list(his_p[bin1:bin2 + 1])
            for c, (pos1, pos2) in regs2:
                flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
                his_p = io.get_signal(c, self.bin_size, "SNP likelihood", flag)
                bin1 = (pos1 - 1) // self.bin_size
                bin2 = (pos2 - 1) // self.bin_size

                data2 += list(his_p[bin1:bin2 + 1])

            d1 = np.array(data1)
            d2 = np.array(data2)
            h1 = np.ones_like(d1[0])
            h2 = np.ones_like(d2[0])
            for i in range(len(d1)):
                if sum(d1[i]) != 0:
                    h1 *= d1[i]
                h1 /= sum(h1)
            for i in range(len(d2)):
                if sum(d2[i]) != 0:
                    h2 *= d2[i]
                h2 /= sum(h2)

            b1, p1 = likelihood_baf_pval(h1)
            b2, p2 = likelihood_baf_pval(h2)

            if stdout:
                print("%s\t%s\t%s\t%.4f\t%e\t%.4f\t%e" % (
                    io.filename, region1, region2, b1, p1, b2, p2))
            ret.append([io.filename, region1, region2, b1, p1, b2, p2])

            if plot:
                plt.plot(h1, "g")
                plt.plot(h2, "b")

        if plot:
            if self.output_filename != "":
                plt.savefig(self._image_filename("comp_baf"), dpi=200)
                plt.close(self.fig)
            elif self.interactive:
                plt.show(block=False)
                plt.draw()
            else:
                plt.show()

        return ret

    def single_cell_allelic_dropout(self, callset=None, res=1000, n_bins=100, threshold=0.1, snp_threshold=0.01,
                                    neigh=False, plot=False, stdout=True, title=None):
        """
        Function used to identify regions without allelic dropout in the case of single cell amplification.
        It requires baf data for bin size. It will filter out all bins with at least one SNP bellow snp_threshold and
        all bins with collective maximum baf likelihood bellow threshold parameter.

        Parameters
        ----------
        callset : str or None
            Name of callset if not default.
        res : int
            Resolution in bins used to calculate percentage of dropouts in region.
        n_bins : int
            Number of bins in histograms.
        threshold : float
            Collective threshold of AF for allelic dropout
        snp_threshold : float
            Single SNP threshold of AF for allelic dropout
        neigh : bool
            Remove neighbouring bins also.
        plot : bool
            Make plots.
        stdout : bool
            Print out good regions

        """

        if plot:
            self.new_figure(panel_count=2, panel_size=(16, 6), title=title)
            ax = self.next_panel()
            bafG = []
            baf = []
            cpos = 0
            sizeG = []
            sizeB = []
        for c in self.io[self.plot_file].snp_chromosomes():
            if len(self.chrom) == 0 or (c in self.chrom):
                snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                    FLAG_USEHAP if self.snp_use_phase else 0)
                i1 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP i1", snp_flag)
                pos, ref, alt, nref, nalt, gt, flag, qual = self.io[self.plot_file].read_snp(c, callset=callset)
                c00 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 0|0", snp_flag)
                c11 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 1|1", snp_flag)
                homs = c00 + c11
                c01 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 0|1", snp_flag)
                c10 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 1|0", snp_flag)
                hets = c01 + c10
                count = c01 + c10 + c00 + c11
                mask = np.zeros_like(i1)
                density = np.zeros(len(mask) // res)
                # mask[hets == 0] = 1
                mask[hets == 0] = 2
                mask[i1 > (0.5 - threshold)] = 1
                for ix in range(len(pos)):
                    if (nref[ix] + nalt[ix]) != 0 and ((gt[ix] % 4) in [1, 2]):
                        b = 1.0 * nalt[ix] / (nref[ix] + nalt[ix])
                        if (b < snp_threshold) or (b > (1 - snp_threshold)):
                            mask[(pos[ix] - 1) // self.bin_size] = 1

                if neigh:
                    ada = mask == 1
                    ada1 = np.roll(ada, 1)
                    ada2 = np.roll(ada, -1)
                    ada1[0] = False
                    ada2[-1] = False
                    mask[ada1] = 1
                    mask[ada2] = 1
                ix = 0
                while ix < len(mask):
                    if mask[ix] == 2:
                        adan = 0
                        if ix > 0 and mask[ix - 1] == 1:
                            adan = 1
                        jx = ix
                        while jx < len(mask) and mask[jx] == 2:
                            jx += 1
                        if jx < len(mask) and mask[jx] == 1:
                            adan = 1
                        mask[ix:jx] = adan
                        ix = jx
                    else:
                        ix += 1
                ix = 0
                ojx = 0
                while ix < len(mask):
                    if mask[ix] == 0:
                        jx = ix
                        while jx < len(mask) and mask[jx] == 0:
                            jx += 1
                        if stdout:
                            print("%s\t%d\t%d" % (c, ix * self.bin_size + 1, jx * self.bin_size))
                        sizeG.append((jx - ix) * self.bin_size)
                        if ix > ojx:
                            sizeB.append((ix - ojx) * self.bin_size)
                        ojx = jx
                        ix = jx
                    else:
                        ix += 1
                if plot:
                    for ix in range(len(density)):
                        density[ix] = np.mean(mask[res * ix:res * (ix + 1)])
                    ax.plot(np.arange(cpos, cpos + len(density)) * res, density)
                    cpos += len(density)
                    for ix in range(len(pos)):
                        if (nref[ix] + nalt[ix]) != 0 and ((gt[ix] % 4) in [1, 2]):
                            baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            if mask[(pos[ix] - 1) // self.bin_size] == 0:
                                bafG.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
        ax.set_xlabel("Position in genome [bins]")
        ax.set_ylabel("Percentage of allelic dropout")
        ax.grid(True)
        if plot:
            self.new_subgrid(2, grid="horizontal", hspace=0.05, wspace=0.2)
            ax = self.next_subpanel()
            ms = 5 * max(np.mean(sizeG), np.mean(sizeB))
            ax.hist(sizeB, bins=np.arange(1, ms, self.bin_size), histtype="step", log=True,
                    label="Allelic dropout regions", linewidth=3)
            ax.hist(sizeG, bins=np.arange(1, ms, self.bin_size), histtype="step", log=True,
                    label="Region with both alleles", linewidth=3)
            ax.legend()
            ax.grid(True)
            ax.set_xlabel("Size [bp]")
            ax.set_ylabel("Number of regions")
            self.fig.add_subplot(ax)

            ax = self.next_subpanel()
            ax.hist(baf, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)),
                    label="All heterozygous variants")
            ax.hist(bafG, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)),
                    label="Region with both alleles")
            ax.legend()
            ax.grid(True)
            ax.set_xlabel("VAF")
            ax.set_ylabel("Distribution")
            self.fig.add_subplot(ax)

            self.fig_show(suffix="allelic_dropout")

    def single_cell_allelic_dropout_2(self, callset=None, res=1000, n_bins=100, threshold=0.1, snp_threshold=0.01,
                                      neigh=False, plot=True, stdout=False, title=None):
        """
        Function used to identify regions without allelic dropout in the case of single cell amplification.
        It requires baf data for bin size. It will filter out all bins with at least one SNP bellow snp_threshold and
        all bins with collective maximum baf likelihood bellow threshold parameter.

        Parameters
        ----------
        callset : str or None
            Name of callset if not default.
        res : int
            Resolution in bins used to calculate percentage of dropouts in region.
        n_bins : int
            Number of bins in histograms.
        threshold : float
            Collective threshold of AF for allelic dropout
        snp_threshold : float
            Single SNP threshold of AF for allelic dropout
        neigh : bool
            Remove neighbouring bins also.
        plot : bool
            Make plots.
        stdout : bool
            Print out good regions

        """

        if plot:
            self.new_figure(panel_count=2, panel_size=(16, 6), title=title)
            ax = self.next_panel()
            bafG = []
            baf = []
            cpos = 0
            sizeG = []
            sizeB = []
            start = 0
            xticks = [0]
            xticks_minor = []
            xticks_labels = []
        for c in self.io[self.plot_file].snp_chromosomes():
            if len(self.chrom) == 0 or (c in self.chrom):
                snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                    FLAG_USEHAP if self.snp_use_phase else 0)
                i1 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP i1", snp_flag)
                pos, ref, alt, nref, nalt, gt, flag, qual = self.io[self.plot_file].read_snp(c, callset=callset)
                c00 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 0|0", snp_flag)
                c11 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 1|1", snp_flag)
                homs = c00 + c11
                c01 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 0|1", snp_flag)
                c10 = self.io[self.plot_file].get_signal(c, self.bin_size, "SNP bin count 1|0", snp_flag)
                hets = c01 + c10
                count = c01 + c10 + c00 + c11
                mask = np.zeros_like(i1)
                density = np.zeros(len(mask) // res)
                # mask[hets == 0] = 1
                mask[hets == 0] = 2
                mask[i1 > (0.5 - threshold)] = 1
                for ix in range(len(pos)):
                    if (nref[ix] + nalt[ix]) != 0 and ((gt[ix] % 4) in [1, 2]):
                        b = 1.0 * nalt[ix] / (nref[ix] + nalt[ix])
                        if (b < snp_threshold) or (b > (1 - snp_threshold)):
                            mask[(pos[ix] - 1) // self.bin_size] = 1

                if neigh:
                    ada = mask == 1
                    ada1 = np.roll(ada, 1)
                    ada2 = np.roll(ada, -1)
                    ada1[0] = False
                    ada2[-1] = False
                    mask[ada1] = 1
                    mask[ada2] = 1
                ix = 0
                while ix < len(mask):
                    if mask[ix] == 2:
                        adan = 0
                        if ix > 0 and mask[ix - 1] == 1:
                            adan = 1
                        jx = ix
                        while jx < len(mask) and mask[jx] == 2:
                            jx += 1
                        if jx < len(mask) and mask[jx] == 1:
                            adan = 1
                        mask[ix:jx] = adan
                        ix = jx
                    else:
                        ix += 1
                ix = 0
                ojx = 0
                while ix < len(mask):
                    if mask[ix] == 0:
                        jx = ix
                        while jx < len(mask) and mask[jx] == 0:
                            jx += 1
                        if stdout:
                            print("%s\t%d\t%d" % (c, ix * self.bin_size + 1, jx * self.bin_size))
                        sizeG.append((jx - ix) * self.bin_size)
                        if ix > ojx:
                            sizeB.append((ix - ojx) * self.bin_size)
                        ojx = jx
                        ix = jx
                    else:
                        ix += 1

                if plot:
                    for ix in range(len(density)):
                        density[ix] = np.mean(mask[res * ix:res * (ix + 1)]) * 100

                    pos1 = np.arange(cpos, cpos + len(density))  # * res
                    if self.markersize == "auto":
                        plt.plot(pos1, density, ls='', marker='o', markersize=1)
                    else:
                        plt.plot(pos1, density, ls='', marker='o', markersize=self.markersize)
                    cpos += len(density)

                    xticks_minor.append(start + len(density) // 2)
                    xticks_labels.append(Genome.canonical_chrom_name(c))
                    start += len(density)
                    xticks.append(start)

                    for ix in range(len(pos)):
                        if (nref[ix] + nalt[ix]) != 0 and ((gt[ix] % 4) in [1, 2]):
                            baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            if mask[(pos[ix] - 1) // self.bin_size] == 0:
                                bafG.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))

        ax.set_xlim([0, start])
        ax.xaxis.set_ticks(xticks)
        ax.xaxis.set_ticklabels([""] * len(xticks))
        ax.minorticks_on()
        ax.xaxis.set_ticks(xticks_minor, minor=True)
        ax.xaxis.set_ticklabels(xticks_labels, minor=True)
        print(xticks_minor, xticks_labels)
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Percentage of allelic dropout")
        ax.grid(True)
        if plot:
            ax = self.next_panel()
            ax.hist(baf, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)),
                    label="All heterozygous variants")
            ax.hist(bafG, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)),
                    label="Region with both alleles")
            ax.legend()
            ax.grid(True)
            ax.set_xlabel("VAF")
            ax.set_ylabel("Distribution")

            self.fig_show(suffix="allelic_dropout")

    def compare_rd_dist(self, regions):
        self.new_figure(panel_count=1)
        ax = self.next_panel()
        ax.set_ylabel("Normalised distribution")
        ax.set_xlabel("Difference in copy number")
        regs = decode_region(regions)
        io1 = self.io[self.plot_files[0]]
        io2 = self.io[self.plot_files[1]]
        bin_size = self.bin_size
        drd = []
        for c, (pos1, pos2) in regs:
            flag_rd = 0
            if self.rd_use_mask:
                flag_rd = FLAG_USEMASK
            mean1, stdev = io1.rd_normal_level(bin_size, flag_rd | FLAG_GC_CORR)
            mean2, stdev = io2.rd_normal_level(bin_size, flag_rd | FLAG_GC_CORR)
            his_p_corr1 = io1.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            his_p_corr2 = io2.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            for i in range(len(his_p_corr1)):
                drd.append(his_p_corr1[i] * 2 / mean1 - his_p_corr2[i] * 2 / mean2)

        # for i in range(n):
        #     io = self.io[ix[i]]
        #     stat = self.io[self.plot_file].get_signal(None, self.bin_size, "RD stat", FLAG_AUTO)
        #     his_p = io.get_signal(None, self.bin_size, "RD p dist", FLAG_AUTO)
        #     bin_size = int(stat[1])
        #     max_rd = int(stat[0])
        #     lim_rd = int(max(2 * stat[4], stat[4] + 3 * stat[5]))
        #     ax.set_xlim([0, lim_rd])
        #     bins = range(0, 2*max_rd + 5*bin_size, bin_size)
        #     x = np.arange(0, max_rd // bin_size * bin_size, 0.1 * bin_size)
        #     #plt.plot(x, normal(x, 1, stat[4], stat[5]), "g-")
        #     x = np.array(bins)
        #     plt.plot(x[1:len(his_p)], his_p[1:] / stat[3],label = io.filename)
        ax.hist(drd, bins=np.linspace(-0.5, 0.5, 100))
        # ax.legend()
        ax.set_yticklabels([])
        ax.grid()
        self.fig_show(suffix="compare_rd")

    def qc(self, snp_qc=True):
        n = len(self.plot_files)
        if snp_qc:
            labels = ["filename", "meanRL", "dRL", "meanFL", "dFL", "meanRD", "dRD", "NdRD", "meanC", "dC", "NdC",
                      "NdBAF",
                      "ALT/REF"]
        else:
            labels = ["filename", "meanRL", "dRL", "meanFL", "dFL", "meanRD", "dRD", "NdRD"]
        print(("{:}\t" * len(labels)).format(*tuple(labels)))
        for i in range(n):
            io = self.io[self.plot_files[i]]
            flag_rd = (FLAG_USEMASK if self.rd_use_mask else 0) | (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
            mean, stdev = io.rd_normal_level(self.bin_size, flag_rd | FLAG_GC_CORR)
            nstdev = stdev / np.sqrt(mean) if mean > 0 else 0
            rfd = io.get_signal(None, None, "read frg dist")
            if rfd is None or len(rfd) == 0:
                mrl, stdrl, mfl, stdfl = 0, 0, 0, 0
            else:
                read_size = np.sum(rfd, axis=1)
                frag_size = np.sum(rfd, axis=0)
                mrl = np.sum(read_size * np.arange(read_size.size)) / np.sum(read_size)
                mfl = np.sum(frag_size * np.arange(frag_size.size)) / np.sum(frag_size)
                mrl2 = np.sum(read_size * np.arange(read_size.size) * np.arange(read_size.size)) / np.sum(read_size)
                mfl2 = np.sum(frag_size * np.arange(frag_size.size) * np.arange(frag_size.size)) / np.sum(frag_size)
                stdrl = np.sqrt(mrl2 - mrl * mrl)
                stdfl = np.sqrt(mfl2 - mfl * mfl)

            if snp_qc:
                snp_chrs = io.snp_chromosomes()
                cbafs = {}
                counts = []
                bafs = []
                talt, tref = 0, 0

                for c in snp_chrs:
                    if c in self.chrom or len(self.chrom) == 0:
                        pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c)
                        maxc = max(nref) + max(nalt)
                        for ix in range(len(pos)):
                            if gt[ix] % 4 in [1, 2] and (not self.snp_use_id or (flag[ix] & 1)) and (
                                    not self.snp_use_mask or (flag[ix] & 2)):
                                tc = nref[ix] + nalt[ix]
                                if tc > 0:
                                    baf = nalt[ix] / tc
                                    if tc not in cbafs:
                                        cbafs[tc] = []
                                    cbafs[tc].append(baf)
                                    counts.append(tc)
                                    bafs.append(baf)
                                    talt += nalt[ix]
                                    tref += nref[ix]
                alt_ref_ratio = talt / tref if tref > 0 else 0
                if len(counts) > 0:
                    meanc = np.mean(counts)
                    stdc = np.std(counts)
                    nstdc = stdc / np.sqrt(meanc) if meanc > 0 else 0
                    maxc = 2 * np.mean(counts)
                    bins = 20
                    x = [[] for j in range(bins)]
                    dbaf = []

                    for count, baf in zip(counts, bafs):
                        if count > (meanc / 2) and count < (3 * meanc / 2):
                            dbaf.append(2 * (baf - 0.5) * np.sqrt(count))
                    std_dbaf = np.std(dbaf)
                else:
                    meanc, stdc, nstdc, std_dbaf = 0, 0, 0, 0
                print("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (io.filename.replace(".pytor", ""), mrl,
                                                                              stdrl, mfl, stdfl, mean, stdev, nstdev,
                                                                              meanc, stdc, nstdc,
                                                                              std_dbaf, alt_ref_ratio))
            else:
                print("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (
                    io.filename.replace(".pytor", ""), mrl, stdrl, mfl, stdfl, mean, stdev, nstdev))

    def snp_dist(self, regions, callset=None, n_bins=100, gt_plot=[0, 1, 2, 3], titles=None, beta_distribution=False,
                 log_scale=False):
        nf = len(self.plot_files)
        regions = regions.split(" ")
        nr = len(regions)
        n = nf * nr
        self.new_figure(panel_count=n)
        for ii in range(nf):
            for i in range(nr):
                ax = self.next_panel()
                if titles is None:
                    ax.set_title(self.file_title(self.plot_files[ii]) + ": " + regions[i], position=(0.01, 1.10),
                                 fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
                else:
                    ax.set_title(titles[i], position=(0.01, 1.10),
                                 fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
                regs = decode_region(regions[i])
                baf = []
                bafP = []
                bafNP = []
                mean_rd = 0
                for c, (pos1, pos2) in regs:
                    pos, ref, alt, nref, nalt, gt, flag, qual = self.io[self.plot_files[ii]].read_snp(c,
                                                                                                      callset=callset)
                    ix = 0
                    while ix < len(pos) and pos[ix] <= pos2:
                        if pos[ix] >= pos1 and (nref[ix] + nalt[ix]) != 0 and ((gt[ix] % 4) in gt_plot):
                            if gt[ix] % 4 != 2:
                                baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                                if flag[ix] & 2:
                                    bafP.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                                    mean_rd += nref[ix] + nalt[ix]
                                else:
                                    bafNP.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            else:
                                baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                                if flag[ix] & 2:
                                    bafP.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                                    mean_rd += nref[ix] + nalt[ix]
                                else:
                                    bafNP.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                        ix += 1
                mean_rd /= len(bafP)
                x_bins = np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1))
                ax.hist(baf, bins=x_bins, label="All heterozygous variants")
                ax.hist(bafP, bins=x_bins, label="P bases only")
                # ax.hist(bafNP, bins=x_bins, label="non-P bases only", histtype=u'step')
                if log_scale:
                    plt.yscale('log', nonposy='clip')

                if beta_distribution:
                    xx = np.linspace(0.2, 0.8, 200)
                    ax.plot(xx, beta_fun.pdf(xx, mean_rd / 2, mean_rd / 2) * len(bafP) / n_bins, c="black",
                            label="Beta distribution")
                ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=3)
                ax.set_xlabel("VAF")
                ax.set_ylabel("Distribution")

        self.fig_show(suffix="snp_dist")

    def baf_dist(self, regions, n_bins=100, titles=None, log_scale=False, minbaf=0.4, maxbaf=0.6):
        bin_size = self.bin_size
        nf = len(self.plot_files)
        regions = regions.split(" ")
        nr = len(regions)
        n = nf * nr
        self.new_figure(panel_count=nf)
        for ii in range(nf):
            ax = self.next_panel()
            for i in range(nr):
                if titles is None:
                    ax.set_title(self.file_title(self.plot_files[ii]), position=(0.01, 1.10),
                                 fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
                else:
                    ax.set_title(titles[i], position=(0.01, 1.10),
                                 fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
                regs = decode_region(regions[i])
                baf = []
                mean_rd = 0
                for c, (pos1, pos2) in regs:
                    flag_snp = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                        FLAG_USEHAP if self.snp_use_phase else 0)
                    cbaf = self.io[self.plot_files[ii]].get_signal(c, bin_size, "SNP baf", flag_snp)
                    baf = baf + list(cbaf[pos1 // bin_size:pos2 // bin_size])
                x_bins = np.arange(minbaf, maxbaf + 1. / (n_bins + 1), 1. / (n_bins + 1))
                ax.hist(baf, bins=x_bins, label=regions[i], lw=3, alpha=0.3)
            if log_scale:
                plt.yscale('log', nonposy='clip')

            ax.set_xlabel("AF")
            ax.set_ylabel("Distribution")
            ax.legend()

        self.fig_show(suffix="baf_dist")

    def phased_baf(self, regions, callset=None, stdout=False):
        regions = regions.split(" ")
        n = len(regions)
        ret = []
        for i in range(n):
            regs = decode_region(regions[i])
            talt = 0
            tref = 0
            taltP = 0
            trefP = 0
            for c, (pos1, pos2) in regs:
                pos, ref, alt, nref, nalt, gt, flag, qual = self.io[self.plot_file].read_snp(c, callset=callset)
                ix = 0
                while ix < len(pos) and pos[ix] <= pos2:
                    if pos[ix] >= pos1 and (nref[ix] + nalt[ix]) != 0:
                        if gt[ix] == 5:
                            talt += nalt[ix]
                            tref += nref[ix]
                            if flag[ix] & 2:
                                taltP += nalt[ix]
                                trefP += nref[ix]
                        elif gt[ix] == 6:
                            tref += nalt[ix]
                            talt += nref[ix]
                            if flag[ix] & 2:
                                trefP += nalt[ix]
                                taltP += nref[ix]
                    ix += 1
            baf = talt / (tref + talt)
            bafP = taltP / (trefP + taltP)
            ret.append([baf, bafP])
            if stdout:
                print("%s\t%f\t%f" % (regions[i], baf, bafP))
        return ret

    def sample2sample(self, s1, s2, callset=None):
        io = self.io[s1]
        io2 = self.io[s2]
        bin_size = self.bin_size
        chroms = io.rd_chromosomes()
        ret = []
        for c in chroms:
            _logger.info("Processing sample2sample %d vs %d, chrom %s" % (s1, s2, c))
            if (c in self.chrom) or len(self.chrom) == 0:
                flag = (FLAG_USEMASK if self.rd_use_mask else 0) | \
                       (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | \
                       (FLAG_USEMASK if self.snp_use_mask else 0) | \
                       (FLAG_USEID if self.snp_use_id else 0) | \
                       (FLAG_USEHAP if self.snp_use_phase else 0)
                snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                    FLAG_USEHAP if self.snp_use_phase else 0)

                if io.signal_exists(c, bin_size, "calls combined", flag):
                    calls = io.read_calls(c, bin_size, "calls combined", flag)
                    if len(calls) > 0:
                        segments = segments_decode(
                            io.get_signal(c, bin_size, "SNP likelihood segments 2d", snp_flag))
                        pos1, ref1, alt1, nref1, nalt1, gt1, flag1, qual1 = io.read_snp(c, callset=callset)
                        pos2, ref2, alt2, nref2, nalt2, gt2, flag2, qual2 = io2.read_snp(c, callset=callset)
                        phase = {}
                        for i in range(len(pos1)):
                            if (nref1[i] + nalt1[i]) > 0 and ((gt1[i] % 4) in [1, 2]):
                                phase[str(pos1[i]) + ref1[i] + alt1[i]] = 0
                                if nalt1[i] > nref1[i]:
                                    phase[str(pos1[i]) + ref1[i] + alt1[i]] = 1
                        binned_snps = {}
                        binned_snps_P = {}
                        for i in range(len(pos2)):
                            snpkey = str(pos2[i]) + ref2[i] + alt2[i]
                            if ((gt2[i] % 4) in [1, 2]) and (snpkey in phase):
                                bin = pos2[i] // bin_size

                                if bin not in binned_snps:
                                    binned_snps[bin] = []
                                if phase[snpkey] == 0:
                                    binned_snps[bin].append((nref2[i], nalt2[i]))
                                else:
                                    binned_snps[bin].append((nalt2[i], nref2[i]))
                                if flag2[i] & 2:
                                    if bin not in binned_snps_P:
                                        binned_snps_P[bin] = []
                                    if phase[snpkey] == 0:
                                        binned_snps_P[bin].append((nref2[i], nalt2[i]))
                                    else:
                                        binned_snps_P[bin].append((nalt2[i], nref2[i]))

                        for call in calls:
                            if in_interval(call["size"], self.size_range) \
                                    and in_interval(call["p_val"], self.p_range) \
                                    and in_interval(call["pN"], self.pN_range) \
                                    and in_interval(call["Q0"], self.Q0_range) \
                                    and in_interval(call["bins"], self.bins_range) \
                                    and in_interval(abs(call["baf"]), self.baf_range):
                                type = {-1: "deletion", 0: "cnnloh", 1: "duplication"}[call["type"]]

                                th1, th2 = 0, 0
                                th1P, th2P = 0, 0
                                for bin in segments[int(call["segment"])]:
                                    if bin in binned_snps:
                                        for h1, h2 in binned_snps[bin]:
                                            th1 += h1
                                            th2 += h2
                                    if bin in binned_snps_P:
                                        for h1, h2 in binned_snps_P[bin]:
                                            th1P += h1
                                            th2P += h2
                                baf = th1 / (th1 + th2) if (th1 + th2) > 0 else 0
                                bafP = th1P / (th1P + th2P) if (th1P + th2P) > 0 else 0
                                ret.append(
                                    [type, c, int(call["start"]), int(call["end"]), th1, th2, baf, th1P, th2P, bafP])
        return ret

    def sample2all(self, s):
        l = []
        for i in self.plot_files:
            l.append(self.sample2sample(s, i))
        for j in range(len(l[0])):
            x = l[0][j]
            print("%s\t%s:%d-%d" % (x[0], x[1], x[2], x[3]), end="")
            for i in range(len(self.plot_files)):
                print("\t%s\t%d\t%d\t%.5f\t%d\t%d\t%.5f" % tuple(
                    [self.io[self.plot_files[i]].filename.replace(".pytor", "")] + l[i][j][4:]), end="")
            print()

    def snp_compare(self, regions, ix1, ix2, callset=None, n_bins=100, titles=None, test_loh=False):
        regions = regions.split(" ")
        n = len(regions)
        self.new_figure(panel_count=n)
        for i in range(n):
            ax = self.next_panel()
            if titles is None:
                ax.set_title(regions[i], position=(0.01, 1.07),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
            else:
                ax.set_title(titles[i], position=(0.01, 1.07),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
            regs = decode_region(regions[i])
            oval = []
            for c, (pos_start, pos_end) in regs:
                pos1, ref1, alt1, nref1, nalt1, gt1, flag1, qual1 = self.io[ix1].read_snp(c, callset=callset)
                pos2, ref2, alt2, nref2, nalt2, gt2, flag2, qual2 = self.io[ix2].read_snp(c, callset=callset)

                counts1 = {}
                counts2 = {}
                ix = 0
                while ix < len(pos1) and pos1[ix] <= pos_end:
                    if pos1[ix] >= pos_start and (nref1[ix] + nalt1[ix]) != 0:
                        counts1[pos1[ix]] = (nref1[ix] / np.sqrt(nref1[ix] ** 2 + nalt1[ix] ** 2),
                                             nalt1[ix] / np.sqrt(nref1[ix] ** 2 + nalt1[ix] ** 2))
                    ix += 1
                ix = 0
                xx = []
                yy = []
                cc = []
                hist1 = []
                hist2 = []
                while ix < len(pos2) and pos2[ix] <= pos_end:
                    if pos2[ix] >= pos_start and (nref2[ix] + nalt2[ix]) != 0:
                        counts2[pos2[ix]] = (nref2[ix], nalt2[ix])
                    ix += 1
                for p in counts1:
                    if p in counts2:
                        xx.append(p)
                        yy.append(counts1[p][1] / (counts1[p][0] + counts1[p][1]))
                        cc.append("green")
                        xx.append(p)
                        yy.append(counts2[p][1] / (counts2[p][0] + counts2[p][1]))
                        cc.append("blue")
                        if counts2[p][1] / (counts2[p][0] + counts2[p][1]) > 0.8:
                            t = counts1[p][1] / (counts1[p][0] + counts1[p][1])
                            if t > 0.2 and t < 0.8:
                                hist1.append(t)
                        else:
                            t = counts1[p][1] / (counts1[p][0] + counts1[p][1])
                            if t > 0.2 and t < 0.8:
                                hist2.append(t)
                    else:
                        xx.append(p)
                        yy.append(counts1[p][1] / (counts1[p][0] + counts1[p][1]))
                        cc.append("red")
                        t = counts1[p][1] / (counts1[p][0] + counts1[p][1])
                        if t > 0.2 and t < 0.8:
                            hist2.append(t)
                for p in counts2:
                    if not (p in counts1):
                        xx.append(p)
                        yy.append(counts2[p][1] / (counts2[p][0] + counts2[p][1]))
                        cc.append("orange")

            if test_loh:
                ax.hist(hist1, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)), histtype='step')
                ax.hist(hist2, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)), histtype='step')
                print("H1:", np.mean(hist1), np.std(hist1), len(hist1))
                print("H2:", np.mean(hist2), np.std(hist2), len(hist2))
                ax.set_xlabel("baf")
                ax.set_ylabel("distribnution")
            else:
                ax.scatter(xx, yy, marker=".", s=0.1, c=cc)
                # ax.hist(oval, bins=np.arange(0, 1.0 + 1. / (n_bins + 1), 1. / (n_bins + 1)))
                ax.set_xlabel("position")
                ax.set_ylabel("baf")

        self.fig_show(suffix="snp_dist")

    def denovo_calls(self, sample, reference, call_type="mosaic"):
        bin_size = self.bin_size
        io = self.io[sample]
        if call_type == "mosaic":
            chroms = io.rd_chromosomes()
            for c in chroms:
                if (c in self.chrom) or len(self.chrom) == 0:
                    flag = (FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR
                    if io.signal_exists(c, bin_size, "calls", flag):
                        calls = io.read_calls(c, bin_size, "calls", flag)
                        for call in calls:
                            if in_interval(call["size"], self.size_range) \
                                    and in_interval(call["p_val"], self.p_range) \
                                    and in_interval(call["pN"], self.pN_range) \
                                    and in_interval(call["Q0"], self.Q0_range):
                                type = "duplication" if call["type"] == 1 else "deletion"
                                region = "%s:%d-%d" % (c, call["start"], call["end"])

                                cn0 = self.genotype([bin_size], region, file_index=sample)[0][3]
                                cref = list(
                                    map(lambda x: self.genotype([bin_size], region, file_index=x)[0][3], reference))
                                if (((sum(map(lambda x: 0 if (cn0 - x) > 0.5 else 1, cref)) == 0) and cn0 > 2.5) \
                                    or ((sum(map(lambda x: 0 if (x - cn0) > 0.5 else 1, cref)) == 0) and cn0 < 1.5)) \
                                        and (sum(map(lambda x: 0 if np.abs(x - 2.) < 0.5 else 1, cref)) == 0):
                                    print(type, region, call["cnv"], cn0, cref)

                                # if n > 1:
                                #     print("%s\t" % self.file_title(i), end="")
                                # print("%s\t%s:%d-%d\t%d\t%.4f\t%e\t%e\t%e\t%e\t%.4f\t%.4f\t" % (
                                #     type, c, call["start"], call["end"], call["size"], call["cnv"], call["p_val"],
                                #     call["p_val_2"], call["p_val_3"], call["p_val_4"], call["Q0"], call["pN"]))

    def genotype(self, bin_sizes, region, p_val=False, interactive=False, file_index=None):
        if file_index is None:
            file_index = self.plot_file
        ret = []
        regs = decode_region(region, max_size=1000000000)
        for c, (pos1, pos2) in regs:
            chr_len = self.io[file_index].get_chromosome_length(c)
            if chr_len is not None and pos2 == 1000000000:
                pos2 = chr_len
            if interactive:
                print(c + ":" + str(pos1) + "-" + str(pos2), end="")
            ret.append([c, pos1, pos2])
            for bs in bin_sizes:
                flag_rd = (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | (FLAG_USEMASK if self.rd_use_mask else 0)
                stat = self.io[file_index].get_signal(c, bs, "RD stat", flag_rd | FLAG_AUTO)
                if stat is None or len(stat) == 0:
                    stat = self.io[file_index].get_signal(c, bs, "RD stat", flag_rd | FLAG_SEX)
                his_p = self.io[file_index].get_signal(c, bs, "RD", flag_rd)
                bin1 = (pos1 - 1) // bs
                bin2 = (pos2 - 1) // bs
                rc = 0
                rc2 = 0
                if bin1 == bin2:
                    try:
                        rc = (pos2 - pos1 + 1) * his_p[bin1] / bs
                        rc2 = (pos2 - pos1 + 1) * his_p[bin1] * his_p[bin1] / bs
                    except IndexError:
                        pass
                else:
                    try:
                        rc += (bin1 * bs - pos1 + 1 + bs) * his_p[bin1] / bs
                        rc += (pos2 - bin2 * bs) * his_p[bin2] / bs
                        rc2 += (bin1 * bs - pos1 + 1 + bs) * his_p[bin1] * his_p[bin1] / bs
                        rc2 += (pos2 - bin2 * bs) * his_p[bin2] * his_p[bin2] / bs
                    except IndexError:
                        pass
                    for ix in range(bin1 + 1, bin2):
                        try:
                            rc += his_p[ix]
                            rc2 += his_p[ix] * his_p[ix]
                        except IndexError:
                            pass
                e2 = 0
                if p_val:
                    e1 = getEValue(stat[4], stat[5], his_p, bin1, bin2 + 1) * 2.9e9 / bs
                    e2 = gaussianEValue(stat[4], stat[5], his_p, bin1, bin2 + 1) * 2.9e9
                if interactive:
                    print("\t%f" % (2. * rc / (stat[4] * (pos2 - pos1 + 1) / bs)), end="")
                    if p_val:
                        print("\t%e\t%e" % (e1, e2), end="")

                ret[-1].append(2. * rc / (stat[4] * (pos2 - pos1 + 1) / bs))
                if p_val:
                    ret[-1].append(e2)
            if interactive:
                print()

        return ret

    def genotype_all(self, bin_sizes, regions, interactive=False, file_index=None):
        if file_index is None:
            file_index = self.plot_file
        rd_gc_chromosomes = {}
        for c in self.io_gc.gc_chromosomes():
            rd_name = self.io[file_index].rd_chromosome_name(c)
            if not rd_name is None:
                rd_gc_chromosomes[rd_name] = c
        ret = {}
        for bs in bin_sizes:
            oc = ""
            ret[bs] = []
            for region in regions:
                regs = decode_region(region, max_size=1000000000)
                c, (pos1, pos2) = regs[0]
                if oc != c:
                    chr_len = self.io[file_index].get_chromosome_length(c)
                    if chr_len is not None and pos2 == 1000000000:
                        pos2 = chr_len
                    flag_rd = (FLAG_GC_CORR if self.rd_use_gc_corr else 0) | (FLAG_USEMASK if self.rd_use_mask else 0)
                    stat = self.io[file_index].get_signal(c, bs, "RD stat", flag_rd | FLAG_AUTO)
                    if stat is None or len(stat) == 0:
                        stat = self.io[file_index].get_signal(c, bs, "RD stat", flag_rd | FLAG_SEX)
                    his_p = self.io[file_index].get_signal(c, bs, "RD", flag_rd)
                    qrd_p = self.io[file_index].get_signal(c, bs, "RD")
                    qrd_u = self.io[file_index].get_signal(c, bs, "RD unique")
                    gc, at, distN = False, False, False
                    if c in rd_gc_chromosomes and self.io_gc.signal_exists(rd_gc_chromosomes[c], None, "GC/AT"):
                        gcat = self.io_gc.get_signal(rd_gc_chromosomes[c], None, "GC/AT")
                        gc, at = gc_at_decompress(gcat)
                        NN = 100 - np.array(gc) - np.array(at)
                        distN = np.zeros_like(NN, dtype="long") - 1
                        distN[NN == 100] = 0
                        prev = 0
                        for Ni in range(0, distN.size):
                            if distN[Ni] == -1:
                                prev += 100
                                distN[Ni] = prev
                            else:
                                prev = 0
                        prev = 0
                        for Ni in range(distN.size - 1, -1, -1):
                            if distN[Ni] > 0:
                                prev += 100
                                if prev < distN[Ni]:
                                    distN[Ni] = prev
                            else:
                                prev = 0
                    snp = c in self.io[file_index].snp_chromosomes()
                    snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                        FLAG_USEHAP if self.snp_use_phase else 0)
                    if snp:
                        snp_likelihood = list(
                            self.io[file_index].get_signal(c, bs, "SNP likelihood", snp_flag).astype("float64"))
                        snp_hets = self.io[file_index].get_signal(c, bs, "SNP bin count 0|1", snp_flag)
                        snp_hets += self.io[file_index].get_signal(c, bs, "SNP bin count 1|0", snp_flag)
                        snp_homs = self.io[file_index].get_signal(c, bs, "SNP bin count 1|1", snp_flag)
                else:
                    if chr_len is not None and pos2 == 1000000000:
                        pos2 = chr_len
                oc = c
                ret[bs].append([c, pos1, pos2])

                bin1 = (pos1 - 1) // bs
                bin2 = (pos2 - 1) // bs
                rc = 0
                rc2 = 0
                sp = 0
                su = 0
                nansize = 0
                if bin1 == bin2:
                    try:
                        if not np.isnan(his_p[bin1]):
                            rc = (pos2 - pos1 + 1) * his_p[bin1] / bs
                            rc2 = (pos2 - pos1 + 1) * his_p[bin1] * his_p[bin1] / bs
                            sp = (pos2 - pos1 + 1) * qrd_p[bin1] / bs
                            su = (pos2 - pos1 + 1) * qrd_u[bin1] / bs
                            nansize = (pos2 - pos1 + 1)
                    except IndexError:
                        pass
                else:
                    try:
                        if not np.isnan(his_p[bin1]):
                            rc += (bin1 * bs - pos1 + 1 + bs) * his_p[bin1] / bs
                            rc2 += (bin1 * bs - pos1 + 1 + bs) * his_p[bin1] * his_p[bin1] / bs
                            sp += (bin1 * bs - pos1 + 1 + bs) * qrd_p[bin1] / bs
                            su += (bin1 * bs - pos1 + 1 + bs) * qrd_u[bin1] / bs
                            nansize += (bin1 * bs - pos1 + 1 + bs)
                        if not np.isnan(his_p[bin2]):
                            rc += (pos2 - bin2 * bs) * his_p[bin2] / bs
                            rc2 += (pos2 - bin2 * bs) * his_p[bin2] * his_p[bin2] / bs
                            sp += (pos2 - bin2 * bs) * qrd_p[bin2] / bs
                            su += (pos2 - bin2 * bs) * qrd_u[bin2] / bs
                            nansize += (pos2 - bin2 * bs)

                    except IndexError:
                        pass
                    for ix in range(bin1 + 1, bin2):
                        try:
                            if not np.isnan(his_p[ix]):
                                rc += his_p[ix]
                                rc2 += his_p[ix] * his_p[ix]
                                sp += qrd_p[ix]
                                su += qrd_u[ix]
                                nansize += bs
                        except IndexError:
                            pass
                if gc:
                    sbin1 = (pos1 - 1) // 100
                    sbin2 = (pos2 - 1) // 100
                    pN = 0
                    if bin1 == bin2:
                        try:
                            pN = (pos2 - pos1 + 1) * (gc[sbin1] + at[sbin1]) / 100
                        except IndexError:
                            pass
                    else:
                        try:
                            pN += (sbin1 * 100 - pos1 + 101) * (gc[sbin1] + at[sbin1]) / 100
                            pN += (pos2 - sbin2 * 100) * (gc[sbin2] + at[sbin2]) / 100

                        except IndexError:
                            pass
                        for ix in range(sbin1 + 1, sbin2):
                            try:
                                pN += gc[ix] + at[ix]
                            except IndexError:
                                pass

                e1 = getEValue(stat[4], stat[5], his_p, bin1, bin2 + 1) * 2.9e9 / bs
                e2 = gaussianEValue(stat[4], stat[5], his_p, bin1, bin2 + 1) * 2.9e9
                dG = -1
                if gc:
                    pN = 1 - pN / (pos2 - pos1 + 1)
                    dG = np.min(distN[sbin1:sbin2])
                else:
                    pN = -1
                    dG = -1
                if nansize == 0:
                    rc = np.nan
                else:
                    rc = 2 * rc / (stat[4] * nansize / bs)
                ret[bs][-1].append(rc)
                ret[bs][-1].append(e1)
                ret[bs][-1].append(e2)
                q0 = 0
                if sp != 0:
                    q0 = (sp - su) / sp
                ret[bs][-1].append(q0)
                ret[bs][-1].append(pN)
                ret[bs][-1].append(dG)
                ret[bs][-1].append(nansize / (pos2 - pos1 + 1))
                if snp:
                    homs = np.sum(snp_homs[bin1:bin2 + 1])
                    hets = np.sum(snp_hets[bin1:bin2 + 1])
                    log_lh = np.zeros_like(snp_likelihood[0], dtype="g")
                    for ix in range(bin1, min(bin2 + 1, len(snp_likelihood))):
                        if np.isfinite(np.sum(snp_likelihood[ix])):
                            tmp = np.log(snp_likelihood[ix] + 1e-100)
                            tmp[tmp == -np.inf] = -100
                            log_lh += tmp

                    log_lh -= np.max(log_lh)
                    lh = np.exp(log_lh)
                    lh = lh / np.sum(lh)
                    baf, baf_p = likelihood_baf_pval(lh)
                    ret[bs][-1] += [homs, hets, baf, baf_p]
                else:
                    ret[bs][-1] += [0, 0, 0, 1]
        if interactive:
            plist = []
            for bs in bin_sizes:
                if len(plist) == 0:
                    plist = ret[bs]
                else:
                    for ix in range(len(ret[bs])):
                        plist[ix] += ret[bs][ix][3:]
            for r in plist:
                print(
                    ("%s:%d-%d" + (len(bin_sizes) * "\t%.4f\t%e\t%e\t%.4f\t%.4f\t%d\t%.4f\t%d\t%d\t%.4f\t%e")) % tuple(
                        r))
        return ret

    def genotype_prompt(self, bin_sizes=[], all=False):
        done = False
        while not done:
            try:
                try:
                    line = raw_input("")
                except NameError:
                    line = input("")
            except EOFError:
                return
            if line is None or line == "":
                done = True
            else:
                if all:
                    self.genotype_all(bin_sizes, [line], interactive=True)
                else:
                    self.genotype(bin_sizes, line, interactive=True)

    def rd_baf_call_models(self, maxcn=10):
        bin_size = self.bin_size
        n = len(self.plot_files)
        ix = self.plot_files
        self.new_figure(panel_count=n)

        for i in range(n):
            ax = self.next_panel()
            io = self.io[ix[i]]
            ax.set_title(self.file_title(ix[i]), position=(0.1, 0.1),
                         fontdict={'verticalalignment': 'bottom', 'horizontalalignment': 'left'})

            chroms = []
            flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                FLAG_USEHAP if self.snp_use_phase else 0) | (FLAG_USEMASK if self.rd_use_mask else 0) | FLAG_GC_CORR

            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = io.snp_chromosome_name(c)
                if len(self.chrom) == 0 or (snp_chr in self.chrom) or (c in self.chrom):
                    if (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                        chroms.append((snp_chr, l))

            x = np.linspace(0, 1, 1000)
            master_lh = {}
            for cn in range(maxcn, -1, -1):
                for h1 in range(cn // 2 + 1):
                    h2 = cn - h1
                    mrd = 2 - 2 * x + x * cn
                    np.seterr(divide='ignore')
                    mbaf = 0.5 - (1 - x + x * h1) / (2 - 2 * x + (h1 + h2) * x)
                    plt.plot(mbaf, mrd, "-", label="%d: %d/%d" % (cn, h1, h2), zorder=6 - cn)

            cix = 0
            cmap = list(map(colors.to_rgba, plt.rcParams['axes.prop_cycle'].by_key()['color']))
            for c, l in chroms:
                call_rd = []
                call_baf = []
                call_label = []
                if io.signal_exists(c, bin_size, "calls combined", flag):
                    calls = io.read_calls(c, bin_size, "calls combined", flag)

                    for call in calls:
                        if call["bins"] > self.min_segment_size:
                            call_rd.append(call["cnv"] * 2)
                            call_baf.append(call["baf"])
                            call_label.append(c + ":" + str(call["start"]) + "-" + str(call["end"]))

                plt.scatter(call_baf, call_rd, s=20, edgecolors='face', marker='.')
                cix += 1

            ax.set_xlabel("|ΔBAF|")
            ax.set_ylabel("Relative RD level")

            ax.legend()

            ax.set_ylim([0, maxcn])
            ax.set_xlim([-0.02, 0.5])
            ax.grid()

        self.fig_show(suffix="models")

    def rd_stat_violin(self):
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_files[0]].rd_chromosome_name(c)
            if (len(self.chrom) == 0 or (rd_chr in self.chrom) or (c in self.chrom)) and rd_chr is not None:
                chroms.append(rd_chr)
        n = len(self.plot_files)
        ix = self.plot_files
        self.new_figure(panel_count=n)
        for i in range(n):
            ax = self.next_panel()
            io = self.io[ix[i]]
            allrd = []
            flag_rd = (FLAG_USEMASK if self.rd_use_mask else 0) | (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
            mean, stdev = io.rd_normal_level(self.bin_size, flag_rd)
            pchr = []
            for c in chroms:
                rd = io.get_signal(c, self.bin_size, "RD", flag_rd) * 2 / mean
                rd = rd[rd < self.rd_range[1]]
                if rd is not None and len(rd) > 5:
                    allrd.append(rd)
                    pchr.append(c)
            pos = list(range(len(allrd), 0, -1))
            ax.violinplot(allrd, pos, showmeans=False, showmedians=True, showextrema=False, vert=False)
            ax.set_yticks(pos)
            ax.set_yticklabels(pchr)
            if self.rd_range[1] < 30:
                ax.set_xticks(range(int(self.rd_range[0]), int(self.rd_range[1]), 1))
            ax.set_xlim([self.rd_range[0], self.rd_range[1]])
            ax.grid(True)
        self.fig_show(suffix="rd_stat_violin")

    def rd_stat_violin(self, regions, factor=1):
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_files[0]].rd_chromosome_name(c)
            if (len(self.chrom) == 0 or (rd_chr in self.chrom) or (c in self.chrom)) and rd_chr is not None:
                chroms.append(rd_chr)
        n = len(self.plot_files)
        ix = self.plot_files
        self.new_figure(panel_count=n)
        for i in range(n):
            ax = self.next_panel()
            io = self.io[ix[i]]
            allrd = []
            flag_rd = (FLAG_USEMASK if self.rd_use_mask else 0) | (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
            mean, stdev = io.rd_normal_level(self.bin_size, flag_rd)
            rchr = []
            for region in regions:
                regs = decode_region(region, max_size=1000000000)
                rchr.append(region)
                crd = np.array([])
                for reg in regs:
                    rd = io.get_signal(reg[0], self.bin_size, "RD", flag_rd) * 2 * factor / mean
                    crd = np.append(crd, rd[reg[1][0] // self.bin_size:reg[1][1] // self.bin_size])
                allrd.append(crd)
            pos = list(range(len(allrd), 0, -1))
            ax.violinplot(allrd, pos, showmeans=False, showmedians=True, showextrema=False, vert=False)
            ax.set_yticks(pos)
            ax.set_yticklabels(rchr)
            if self.rd_range[1] < 30:
                ax.set_xticks(range(int(self.rd_range[0]), int(self.rd_range[1]), 1))
            ax.set_xlim([self.rd_range[0], self.rd_range[1]])
            ax.grid(True)
        self.fig_show(suffix="rd_stat_violin")

    def rd_stat_fit_factor(self, regions, factor_range=[0.5, 10], df=0.01):
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_files[0]].rd_chromosome_name(c)
            if (len(self.chrom) == 0 or (rd_chr in self.chrom) or (c in self.chrom)) and rd_chr is not None:
                chroms.append(rd_chr)
        n = len(self.plot_files)
        ix = self.plot_files
        ret = []
        for i in range(n):
            io = self.io[ix[i]]
            allrd = []
            flag_rd = (FLAG_USEMASK if self.rd_use_mask else 0) | (FLAG_GC_CORR if self.rd_use_gc_corr else 0)
            mean, stdev = io.rd_normal_level(self.bin_size, flag_rd)
            rchr = []
            means = []
            for region in regions:
                regs = decode_region(region, max_size=1000000000)
                rchr.append(region)
                crd = np.array([])
                for reg in regs:
                    rd = io.get_signal(reg[0], self.bin_size, "RD", flag_rd) * 2 / mean
                    crd = np.append(crd, rd[reg[1][0] // self.bin_size:reg[1][1] // self.bin_size])
                allrd.append(crd)
                means.append(np.mean(crd))
            def min_int(f):
                t=np.tensordot(means, f, axes=0)
                return np.sum((np.round(t)-t)**2,axis=0)
            fac = min(np.linspace(factor_range[0], factor_range[1], int((factor_range[1]-factor_range[0])/df)), key=lambda x: min_int(x))
            ret.append(fac)
        return ret









    def read_fragment_dist(self):
        n = len(self.plot_files)
        ix = self.plot_files
        self.new_figure(panel_count=n)
        for i in range(n):
            self.new_subgrid(2, hspace=0.2, wspace=0.2)
            io = self.io[ix[i]]
            rfd = io.get_signal(None, None, "read frg dist")
            rd = np.sum(rfd, axis=1)
            fd = np.sum(rfd, axis=0)
            ax = self.next_subpanel()
            ax.plot(rd)
            ax = self.next_subpanel()
            ax.plot(fd)

        self.fig_show(suffix="read_fragment_dist")


def anim_plot_likelihood(likelihood, segments, n, res, iter, prefix, maxp, minp):
    mm = [[0] * res] * n
    for i in range(len(segments)):
        for b in segments[i]:
            mm[b] = list(likelihood[i])
    fig = plt.figure(1, figsize=(16, 9), dpi=120, facecolor='w', edgecolor='k')
    fig.suptitle(
        "Iter: " + str(iter) + "   /   Segments: " + str(len(segments)) + "   /   Overlap interval: (" + (
                '%.4f' % minp) + "," + (
                '%.4f' % maxp) + ")", fontsize='large')
    plt.subplot(211)
    plt.ylabel("BAF")
    plt.imshow(np.transpose(np.array(mm)), aspect='auto')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.yticks([0, 50.5, 101, 151.5, 201], ("1.00", "0.75", "0.50", "0.25", "0.00"))
    # plt.grid(True,color="w")
    plt.subplot(212)
    plt.xlabel("BAF")
    plt.ylabel("Likelihood")
    plt.xticks([0, 0.25, 0.50, 0.75, 1.0])
    plt.grid(True, color="b")
    for i in range(len(likelihood)):
        plt.plot(np.linspace(1. / (res + 1), 1. - 1. / (res + 1), res), likelihood[i])
    plt.savefig(prefix + "_" + str(iter).zfill(4), dpi=150)
    plt.close(fig)


def anim_plot_rd(level, error, segments, n, iter, prefix, maxp, minp, mean):
    rd = [np.nan] * n
    for i in range(len(segments)):
        for b in segments[i]:
            rd[b] = level[i]

    fig = plt.figure(1, figsize=(16, 9), dpi=120, facecolor='w', edgecolor='k')
    fig.suptitle(
        "Iter: " + str(iter) + "   /   Segments: " + str(len(segments)) + "   /   Overlap interval: (" + (
                '%.4f' % minp) + "," + (
                '%.4f' % maxp) + ")", fontsize='large')
    plt.subplot(211)
    plt.ylabel("RD")
    plt.step(range(n), rd, "k")
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.yticks(np.arange(0, 3, 0.5) * mean, [])
    plt.ylim([0, 3 * mean])
    plt.grid(True, color="grey")

    plt.subplot(212)
    plt.xlabel("RD")
    plt.ylabel("Likelihood")
    plt.xticks(np.arange(0, 3, 0.5) * mean, [])
    plt.xlim([0, 3 * mean])
    plt.grid(True, color="grey")
    for i in range(len(level)):
        xx = np.linspace(0, 3 * mean, 300)
        yy = normal(xx, 1, level[i], error[i])
        plt.plot(xx, yy)
    plt.savefig(prefix + "_" + str(iter).zfill(4), dpi=150)
    plt.close(fig)


def anim_plot_rd_likelihood(level, error, likelihood, segments, n, res, iter, prefix, maxp, mean):
    rd = [np.nan] * n
    for i in range(len(segments)):
        for b in segments[i]:
            rd[b] = level[i]

    mm = [[0] * res] * n
    for i in range(len(segments)):
        for b in segments[i]:
            mm[b] = list(likelihood[i])

    fig, ax = plt.subplots(2, 2, figsize=(16, 9), dpi=120, facecolor='w', edgecolor='k',
                           gridspec_kw={
                               'width_ratios': [2, 1],
                               'height_ratios': [1, 1]})
    fig.suptitle(
        "Iter: " + str(iter) + "   /   Segments: " + str(len(segments)) + "   /   Maximal overlap: " + (
                '%.4f' % maxp), fontsize='large')
    ax[0][0].set_ylabel("RD")
    ax[0][0].step(range(n), rd, "k")
    ax[0][0].set_xlim([0, n])
    ax[0][0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax[0][0].set_yticks(np.arange(0, 3, 0.5) * mean)
    ax[0][0].set_yticklabels([])

    ax[0][0].set_ylim([0, 3 * mean])
    ax[0][0].grid(True, color="grey")

    ax[0][1].set_ylabel("")
    ax[0][1].set_xlabel("")
    ax[0][1].set_yticks(np.arange(0, 3, 0.5) * mean)
    ax[0][1].set_yticklabels([])
    ax[0][1].set_xticklabels([])

    ax[0][1].set_ylim([0, 3 * mean])
    ax[0][1].grid(True, color="grey")
    for i in range(len(level)):
        xx = np.linspace(0, 3 * mean, 300)
        yy = normal(xx, 1, level[i], error[i])
        ax[0][1].plot(yy, xx)

    ax[1][0].set_ylabel("BAF")
    ax[1][0].imshow(np.transpose(np.array(mm)), aspect='auto')
    ax[1][0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax[1][0].set_yticks([0, 50.5, 101, 151.5, 201])
    ax[1][0].set_yticklabels(["1.00", "0.75", "0.50", "0.25", "0.00"])

    # plt.grid(True,color="w")

    ax[1][1].set_ylabel("")
    ax[1][1].set_xlabel("Likelihood")
    ax[1][1].set_yticks([0, 0.25, 0.50, 0.75, 1.0])
    ax[1][1].set_yticklabels([])
    ax[1][1].set_xticklabels([])

    ax[1][1].grid(True, color="b")
    for i in range(len(likelihood)):
        ax[1][1].plot(likelihood[i], np.linspace(1. / (res + 1), 1. - 1. / (res + 1), res))

    plt.subplots_adjust(bottom=0.1, left=0.1, wspace=0., hspace=0.)

    plt.savefig(prefix + "_" + str(iter).zfill(4), dpi=150)
    plt.close(fig)
