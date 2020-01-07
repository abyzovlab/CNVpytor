""" cnvpytor.viewer

Class Viewer: ploting CNVpytor data
"""
from __future__ import absolute_import, print_function, division

from .io import *
from .utils import *
from .genome import *
from .viewparams import ViewParams, HelpDescription
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors

import numpy as np
import logging
import readline
import traceback

_logger = logging.getLogger("cnvpytor.viewer")


class Reader:
    def __init__(self, files):
        self.io = [IO(f, ro=True) for f in files]


class Export(Reader):
    def jbrowse(self):
        pass

    def cnvnator(self):
        pass


class Show(Reader):
    def ls(self):
        for i in self.io:
            i.ls()

    def meta(self):
        for i in self.io:
            i.read_meta_attribute()

    def info(self, bin_sizes):
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
        print(("{:25}{:>20}{:>20}{:>20}{:>20}" + "{:>20}" * (len(labels) - 5)).format(*tuple(labels)))
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
            print("{:25}{:20.2f}{:20.2f}{:20.2f}{:20.2f}".format(i.filename, mrl, sdr, mfl, sdf), end="")
            for bs in bin_sizes:
                for flag in [FLAG_AUTO, FLAG_SEX, FLAG_MT]:
                    if bs <= 500 or not flag == FLAG_MT:
                        if i.signal_exists(None, bs, "RD stat", flags=flag):
                            stat = i.get_signal(None, bs, "RD stat", flags=flag)
                            if stat[4] > 0:
                                stat[5] /= stat[4] / 100.
                            print("{:20.2f}{:20.2f}".format(stat[4], stat[5]), end="")
                        else:
                            print("{:>20}{:>20}".format("-", "-"), end="")
                        if i.signal_exists(None, bs, "RD stat", flags=(flag | FLAG_GC_CORR)):
                            stat = i.get_signal(None, bs, "RD stat", flags=(flag | FLAG_GC_CORR))
                            if stat[4] > 0:
                                stat[5] /= stat[4] / 100.
                            print("{:20.2f}{:20.2f}".format(stat[4], stat[5]), end="")
                        else:
                            print("{:>20}{:>20}".format("-", "-"), end="")
            print()


class Viewer(Show, ViewParams, HelpDescription):

    def __init__(self, files, params):
        _logger.debug("Viewer class init: files [%s], params %s." % (", ".join(files), str(params)))
        ViewParams.__init__(self, params)
        Show.__init__(self, files)
        self.io_gc = self.io[0]
        self.io_mask = self.io[0]
        self.reference_genome = None
        self.interactive = False
        self.plot_files = list(range(len(files)))
        self.fig = None
        if self.io[0].signal_exists(None, None, "reference genome"):
            rg_name = np.array(self.io[0].get_signal(None, None, "reference genome")).astype("str")[0]
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
                if current == "rd":
                    self.rd(int(p), self.rd_use_mask)
                if current == "baf":
                    self.baf(int(p), self.rd_use_mask)
                if current == "likelihood":
                    self.likelihood(int(p))
                elif current == "manhattan":
                    self.manhattan(int(p), use_mask=self.rd_use_mask)
                elif current == "calls":
                    self.manhattan(int(p), use_mask=self.rd_use_mask, plot_type="calls")
                elif current == "stat":
                    self.stat(int(p))
                elif current == "circular":
                    self.circular(int(p), self.chrom, self.rd_use_mask)
                elif current == "regions":
                    self.multiple_regions(int(p), regions, panels=self.panels)
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

    def plot(self, command):
        self.interactive = False
        self.parse(command)

    def prompt(self):
        self.interactive = True

        for p in self.params:
            self.command_tree["set"][p] = None
            if type(self.params[p]) == type(True):
                self.command_tree["unset"][p] = None
        for c in self.param_help:
            self.command_tree["help"][c] = None
        chromosomes = set({})
        for f in self.io:
            chromosomes = chromosomes.union(set(f.rd_chromosomes()))
            chromosomes = chromosomes.union(set(f.snp_chromosomes()))
        for c in chromosomes:
            self.command_tree[c] = None
        self.command_tree["set"]["style"] = dict(zip(plt.style.available, [None] * len(plt.style.available)))

        readline.parse_and_bind("tab: complete")
        completer = PromptCompleter(self.command_tree)
        readline.set_completer(completer.complete)
        quit = False
        try:
            while not quit:
                try:
                    line = raw_input("cnvpytor> ")
                except NameError:
                    line = input("cnvpytor> ")

                pre = line.split(">")
                f = pre[0].strip().split(" ")
                n = len(f)
                if len(line) == 0:
                    continue
                elif f[0] == "quit":
                    quit = True
                elif line[0] == "|":
                    try:
                        eval(compile(line[1:], '<string>', 'single'))
                    except Exception as e:
                        print(traceback.format_exc())
                elif f[0] == "save" and n > 1:
                    plt.savefig(f[1])
                elif f[0] in ["draw", "repaint", "update"] and n == 1:
                    self.fig.canvas.draw()
                elif f[0] == "ls":
                    self.ls()
                elif f[0] == "show":
                    if n == 1:
                        self.show()
                elif f[0] == "set" and n > 1:
                    self.set(f[1], f[2:])
                elif f[0] == "help" and n > 1:
                    self.help(f[1])
                elif f[0] == "help" and n == 1:
                    self.help("help")
                elif f[0] == "unset" and n > 1:
                    self.unset(f[1])
                elif f[0] == "genotype" and n > 1:
                    for ni in range(1, n):
                        self.genotype([self.bin_size], f[ni])
                elif f[0] == "compare" and n == 3:
                    self.compare(f[1], f[2], plot=True)
                elif f[0] == "compare" and n == 4:
                    self.compare(f[1], f[2], n_bins=int(f[3]), plot=True)
                elif f[0] == "info" and n > 1:
                    self.info(list(map(binsize_type,f[1:])))
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

    def image_filename(self, sufix):
        parts = self.output_filename.split(".")
        if parts[-1] != "png" and parts[-1] != "pdf" and parts[-1] != "jpg" and parts[-1] != "eps" and parts[
            -1] != "svg":
            _logger.warning("File extension should be: .jpg, .png, .svg, .eps or .pdf")
            exit(0)
        parts[-1] = sufix + "." + parts[-1]
        return ".".join(parts)

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
                    _logger.info("RD stat for %s - number of miochondria per cell: %.2f +- %.2f" % (
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
            plt.savefig(self.image_filename("stat"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def rd(self, bin_size, use_mask):
        plt.clf()
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        plt.rcParams["font.size"] = 8
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_file].rd_chromosome_name(c)
            if self.io[self.plot_file].signal_exists(rd_chr, bin_size, "RD", 0) and \
                    self.io[self.plot_file].signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                chroms.append((rd_chr, l))
        sx, sy = self.panels_shape(len(chroms))
        self.fig = plt.figure(1, figsize=(sx, sy), dpi=200, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(8)
            self.fig.set_figwidth(12)
        ix = 1
        for c, l in chroms:
            flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
            stat = self.io[self.plot_file].get_signal(None, bin_size, "RD stat", flag)
            if stat is None:
                _logger.error(
                    "Data for bin size %d is missing in file '%s'!" % (bin_size, self.io[self.plot_file].filename))
                exit(0)
            flag_rd = 0
            if use_mask:
                flag_rd = FLAG_USEMASK
            his_p = self.io[self.plot_file].get_signal(c, bin_size, "RD", flag_rd)
            his_p_corr = self.io[self.plot_file].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            his_p_seg = self.io[self.plot_file].get_signal(c, bin_size, "RD partition", flag_rd | FLAG_GC_CORR)
            his_p_call = self.io[self.plot_file].get_signal(c, bin_size, "RD call", flag_rd | FLAG_GC_CORR)
            his_p_mosaic_seg = self.io[self.plot_file].get_signal(c, bin_size, "RD mosaic segments",
                                                                  flag_rd | FLAG_GC_CORR)
            his_p_mosaic_seg = segments_decode(his_p_mosaic_seg)
            his_p_mosaic_call = self.io[self.plot_file].get_signal(c, bin_size, "RD mosaic call",
                                                                   flag_rd | FLAG_GC_CORR)
            his_p_mosaic = np.zeros_like(his_p) * np.nan
            if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call_mosaic:
                for seg, lev in zip(list(his_p_mosaic_seg), list(his_p_mosaic_call[0])):
                    for segi in seg:
                        his_p_mosaic[segi] = lev
            ax = plt.subplot(sx, sy, ix)
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * stat[4], [])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // bin_size, 10e6 // bin_size), [])
            ax.set_ylim([0, max(3. * stat[4], stat[4] + 5. * stat[5])])
            n_bins = l // bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()

            plt.step(his_p, "grey")
            plt.step(his_p_corr, "k")
            if his_p_seg is not None and len(his_p_seg) > 0 and self.rd_partition:
                plt.step(his_p_seg, "r")
            if his_p_call is not None and len(his_p_call) > 0 and self.rd_call:
                plt.step(his_p_call, "g")
            if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call_mosaic:
                plt.step(his_p_mosaic, "b")
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.output_filename != "":
            plt.savefig(self.image_filename("gview"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def rd_diff(self, bin_size, use_mask, file1, file2):
        plt.clf()
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        plt.rcParams["font.size"] = 8
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[self.plot_file].rd_chromosome_name(c)
            if self.io[self.plot_file].signal_exists(rd_chr, bin_size, "RD", 0) and \
                    self.io[self.plot_file].signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                chroms.append((rd_chr, l))
        sx, sy = self.panels_shape(len(chroms))
        self.fig = plt.figure(1, figsize=(sx, sy), dpi=200, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(8)
            self.fig.set_figwidth(12)
        ix = 1
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
            flag_rd = 0
            if use_mask:
                flag_rd = FLAG_USEMASK
            his_p_corr1 = self.io[file1].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            his_p_corr2 = self.io[file2].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            ax = plt.subplot(sx, sy, ix)
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks(np.arange(0, 2, 0.25), [])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // bin_size, 10e6 // bin_size), [])
            ax.set_ylim([0, 1])
            n_bins = l // bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()

            plt.step(np.abs(his_p_corr1 / stat1[4] - his_p_corr2 / stat2[4]), "k")
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.output_filename != "":
            plt.savefig(self.image_filename("gview"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def likelihood(self, bin_size):
        plt.clf()
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        plt.rcParams["font.size"] = 8
        chroms = []
        if self.reference_genome is None:
            chroms = self.io[self.plot_file].snp_chromosomes()
        else:
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = self.io[self.plot_file].snp_chromosome_name(c)
                if self.io[self.plot_file].signal_exists(snp_chr, bin_size, "SNP likelihood", snp_flag) and (
                        Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                    chroms.append(snp_chr)
        sx, sy = self.panels_shape(len(chroms))
        self.fig = plt.figure(1, figsize=(sx, sy), dpi=200, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(8)
            self.fig.set_figwidth(12)
        ix = 1
        for c in chroms:
            likelihood = self.io[self.plot_file].get_signal(c, bin_size, "SNP likelihood", snp_flag)
            img = np.array(likelihood).transpose()
            ax = plt.subplot(sx, sy, ix)
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.imshow(img, aspect='auto')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticks(np.arange(0, likelihood.shape[0], 50), [])
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.output_filename != "":
            plt.savefig(self.image_filename("gview"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def snp(self, size=10, plot_gt=None, plot_pmask=None):
        if plot_pmask is None:
            plot_pmask = [0, 1]
        if plot_gt is None:
            plot_gt = [0, 1, 2, 3]
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(12, 8), dpi=90, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(8)
            self.fig.set_figwidth(12)
        chroms = []
        if self.reference_genome is None:
            chroms = self.io[self.plot_file].snp_chromosomes()
        else:
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = self.io[self.plot_file].snp_chromosome_name(c)
                if self.io[self.plot_file].signal_exists(snp_chr, None, "SNP pos", 0) and \
                        self.io[self.plot_file].signal_exists(snp_chr, None, "SNP desc", 0) and \
                        self.io[self.plot_file].signal_exists(snp_chr, None, "SNP counts", 0) and \
                        self.io[self.plot_file].signal_exists(snp_chr, None, "SNP qual", 0) and \
                        (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                    chroms.append(snp_chr)
        sx, sy = self.panels_shape(len(chroms))
        ix = 1
        for c in chroms:
            pos, ref, alt, nref, nalt, gt, flag, qual = self.io[self.plot_file].read_snp(c)
            hpos = []
            baf = []
            color = []
            for i in range(len(pos)):
                if (nref[i] + nalt[i]) != 0:
                    if (gt[i] % 4 in plot_gt) and ((flag[i] >> 1) in plot_pmask):
                        hpos.append(pos[i])
                        if gt[i] % 4 != 2:
                            baf.append(1.0 * nalt[i] / (nref[i] + nalt[i]))
                        else:
                            baf.append(1.0 * nref[i] / (nref[i] + nalt[i]))
                        color.append(self.snp_colors[(gt[i] % 4) * 2 + (flag[i] >> 1)])

            ax = plt.subplot(sx, sy, ix)
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], [])
            l = max(pos)
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6), 10e6), [])
            ax.set_ylim([0., 1.])
            ax.set_xlim([-0.05 * l, 1.05 * l])
            ax.grid()
            if self.markersize == "auto":
                ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=0.7)
            else:
                ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=0.7)
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.output_filename != "":
            plt.savefig(self.image_filename("gview"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()


    @staticmethod
    def panels_shape(n):
        sx, sy = 1, 1
        if n == 2:
            sx = 2
        elif n in [3, 4]:
            sx, sy = 2, 2
        elif n in [5, 6]:
            sx, sy = 3, 2
        elif n in [7, 8, 9]:
            sx, sy = 3, 3
        elif n in [10, 11, 12]:
            sx, sy = 4, 3
        elif n in [13, 14, 15, 16]:
            sx, sy = 4, 4
        elif n in [17, 18, 19, 20]:
            sx, sy = 5, 4
        elif n in [21, 22, 23, 24]:
            sx, sy = 6, 4
        else:
            while sx * sy < n:
                sx += 1
                sy = int(2. * sx / 3 + 1.)
        return sx, sy

    def manhattan(self, bin_size, use_mask=False, plot_type="rd"):
        plt.clf()
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        n = len(self.plot_files)
        ix = self.plot_files

        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(1.5 * n)
            self.fig.set_figwidth(12)
        grid = gridspec.GridSpec(n, 1, wspace=0.2, hspace=0.2)
        for i in range(n):
            ax = self.fig.add_subplot(grid[i])
            io = self.io[ix[i]]
            ax.set_title(io.filename, position=(0.01, 1.07),
                         fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})

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

                max_m = 0
                for c, l in chroms:
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    stat = io.get_signal(None, bin_size, "RD stat", flag)
                    if stat[4] > max_m:
                        max_m = stat[4]
                    flag_rd = 0
                    if use_mask:
                        flag_rd = FLAG_USEMASK
                    his_p = io.get_signal(c, bin_size, "RD", flag_rd)
                    his_p_corr = io.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
                    if self.rd_manhattan_call:
                        his_p_call = io.get_signal(c, bin_size, "RD call", flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_seg = io.get_signal(c, bin_size, "RD mosaic segments",
                                                         flag_rd | FLAG_GC_CORR)
                        his_p_mosaic_seg = segments_decode(his_p_mosaic_seg)
                        his_p_mosaic_call = io.get_signal(c, bin_size, "RD mosaic call",
                                                          flag_rd | FLAG_GC_CORR)
                        his_p_mosaic = np.zeros_like(his_p) * np.nan
                        if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call_mosaic:
                            for seg, lev in zip(list(his_p_mosaic_seg), list(his_p_mosaic_call[0])):
                                for segi in seg:
                                    his_p_mosaic[segi] = lev
                    pos = range(apos, apos + len(his_p))
                    ax.text(apos + len(his_p) // 2, stat[4] // 10, Genome.canonical_chrom_name(c),
                            fontsize=8, verticalalignment='bottom', horizontalalignment='center', )
                    if self.markersize == "auto":
                        plt.plot(pos, his_p_corr, ls='', marker='.')
                    else:
                        plt.plot(pos, his_p_corr, ls='', marker='.', markersize=self.markersize)
                    if self.rd_manhattan_call:
                        if his_p_call is not None and len(his_p_call) > 0 and self.rd_call:
                            plt.step(pos, his_p_call, "r")
                        if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call_mosaic:
                            plt.plot(pos, his_p_mosaic, "k")
                    apos += len(his_p)
                    xticks.append(apos)
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 15, 0.5) * max_m, [])
                ax.xaxis.set_ticks(xticks, [])
                ax.set_ylim([self.rd_manhattan_range[0] * max_m, self.rd_manhattan_range[1] * max_m])
            else:
                chroms = []
                snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
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
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO

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
                    # plt.plot(call_pos, call_baf, color=np.array(call_c), ls='', marker='.')
                    apos += l // bin_size
                    xticks.append(apos)
                    cix += 1

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 0.5, 0.1), [])
                ax.xaxis.set_ticks(xticks, [])
                ax.set_ylim([0, 0.5])

            n_bins = apos
            ax.set_xlim([0, n_bins])
            ax.grid()
        plt.subplots_adjust(bottom=0.05, top=0.95, wspace=0, hspace=0, left=0.05, right=0.95)

        if self.output_filename != "":
            plt.savefig(self.image_filename("manhattan" if plot_type == "rd" else "snp_calls"), dpi=200)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def multiple_regions(self, bin_size, regions, panels=["rd"], sep_color="g"):
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(12, 8), facecolor='w', edgecolor='k')
        grid = gridspec.GridSpec(len(self.io), len(regions), wspace=0.2, hspace=0.2)
        ix = 0
        for i in self.io:
            for r in regions:
                self.regions(i, grid[ix], bin_size, r, panels=panels, sep_color=sep_color)
                ix += 1
        plt.subplots_adjust(bottom=0.05, top=0.95, wspace=0, hspace=0, left=0.05, right=0.95)

        if self.output_filename != "":
            plt.savefig(self.image_filename("regions"), dpi=150)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def regions(self, io, element, bin_size, region, panels=["rd"], sep_color="g"):
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
        grid = gridspec.GridSpecFromSubplotSpec(len(panels), 1, subplot_spec=element, wspace=0, hspace=0.1)
        r = decode_region(region)
        for i in range(len(panels)):
            ax = self.fig.add_subplot(grid[i])
            if i == 0:
                ax.set_title(io.filename + ": " + region, position=(0.01, 0.9),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                             color='C0')

            if panels[i] == "rd":
                g_p = []
                g_p_corr = []
                g_p_seg = []
                g_p_call = []
                g_p_call_mosaic = []
                mean, stdev = 0, 0
                borders = []
                for c, (pos1, pos2) in r:
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    flag_rd = 0
                    if self.rd_use_mask:
                        flag_rd = FLAG_USEMASK
                    stat = io.get_signal(None, bin_size, "RD stat", flag)
                    mean = stat[4]
                    stdev = stat[5]
                    his_p = io.get_signal(c, bin_size, "RD", flag_rd)
                    his_p_corr = io.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
                    his_p_seg = io.get_signal(c, bin_size, "RD partition", flag_rd | FLAG_GC_CORR)
                    his_p_call = io.get_signal(c, bin_size, "RD call", flag_rd | FLAG_GC_CORR)
                    his_p_mosaic_seg = io.get_signal(c, bin_size, "RD mosaic segments",
                                                     flag_rd | FLAG_GC_CORR)
                    his_p_mosaic_seg = segments_decode(his_p_mosaic_seg)
                    his_p_mosaic_call = io.get_signal(c, bin_size, "RD mosaic call",
                                                      flag_rd | FLAG_GC_CORR)
                    his_p_mosaic = np.zeros_like(his_p) * np.nan
                    if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call_mosaic:
                        for seg, lev in zip(list(his_p_mosaic_seg), list(his_p_mosaic_call[0])):
                            for segi in seg:
                                his_p_mosaic[segi] = lev

                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    g_p.extend(list(his_p[start_bin:end_bin]))
                    g_p_corr.extend(list(his_p_corr[start_bin:end_bin]))
                    if his_p_seg is not None and len(his_p_seg) > 0 and self.rd_partition:
                        g_p_seg.extend(list(his_p_seg[start_bin:end_bin]))
                    if his_p_call is not None and len(his_p_call) > 0 and self.rd_call:
                        g_p_call.extend(list(his_p_call[start_bin:end_bin]))
                    if his_p_mosaic_call is not None and len(his_p_mosaic_call) > 0 and self.rd_call_mosaic:
                        g_p_call_mosaic.extend(list(his_p_mosaic[start_bin:end_bin]))
                    borders.append(len(g_p) - 1)

                # ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                l = len(g_p)
                # ax.xaxis.set_ticks(np.arange(0, l, 10), [])
                ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * mean, [])
                ax.set_ylim([0, max(3. * mean, mean + 5. * stdev)])
                ax.set_xlim([-l * 0.0, l * 1.0])

                ax.yaxis.grid()
                ax.step(g_p, "grey")
                ax.step(g_p_corr, "k")
                if len(g_p_seg) > 0:
                    plt.step(g_p_seg, "r")
                if len(g_p_call) > 0:
                    plt.step(g_p_call, "g")
                if len(g_p_call_mosaic) > 0:
                    plt.step(g_p_call_mosaic, "b")
                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "snp":
                borders = []
                hpos = []
                baf = []
                color = []
                start_pos = 0
                for c, (pos1, pos2) in r:
                    pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c)
                    ix = 0
                    mdp = 0
                    while ix < len(pos) and pos[ix] <= pos2:
                        if pos[ix] >= pos1 and (nref[ix] + nalt[ix]) != 0:
                            hpos.append(start_pos + pos[ix] - pos1)
                            if pos[ix] - pos1 > mdp:
                                mdp = pos[ix] - pos1
                            if gt[ix] % 4 != 2:
                                baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            else:
                                baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                            color.append(self.snp_colors[(gt[ix] % 4) * 2 + (flag[ix] >> 1)])
                        ix += 1
                    start_pos += mdp
                    borders.append(start_pos)

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], [])
                l = max(hpos)
                # ax.xaxis.set_ticks(np.arange(0, (l + 10e6), 10e6), [])
                ax.set_ylim([0., 1.])
                ax.set_xlim([0, borders[-1]])
                ax.yaxis.grid()
                if self.markersize == "auto":
                    ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=0.7)
                else:
                    ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=self.markersize, alpha=0.7)

                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "baf":
                g_baf, g_maf, g_i1, g_i2 = [], [], [], []
                borders = []
                for c, (pos1, pos2) in r:
                    flag_snp = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0) | (
                        FLAG_USEHAP if self.snp_use_phase else 0)
                    baf = io.get_signal(c, bin_size, "SNP baf", flag_snp)
                    maf = io.get_signal(c, bin_size, "SNP maf", flag_snp)
                    i1 = io.get_signal(c, bin_size, "SNP i1", flag_snp)
                    i2 = io.get_signal(c, bin_size, "SNP i2", flag_snp)

                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    g_baf.extend(list(baf[start_bin:end_bin]))
                    g_maf.extend(list(maf[start_bin:end_bin]))
                    g_i1.extend(list(i1[start_bin:end_bin]))
                    g_i2.extend(list(i2[start_bin:end_bin]))
                    borders.append(len(g_baf) - 1)

                # ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                l = len(g_baf)
                # ax.xaxis.set_ticks(np.arange(0, l, 10), [])
                ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], [])
                ax.set_ylim([0, 1])
                ax.set_xlim([-l * 0.0, l * 1.0])

                ax.yaxis.grid()
                ax.step(g_baf, "grey")
                ax.step(g_maf, "k")
                ax.step(g_i1, "r")
                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                self.fig.add_subplot(ax)

            elif panels[i] == "likelihood":
                borders = []
                gl = []
                for c, (pos1, pos2) in r:
                    likelihood = io.get_signal(c, bin_size, "SNP likelihood", snp_flag)
                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    gl.extend(list(likelihood[start_bin:end_bin]))
                    borders.append(len(gl) - 1)
                img = np.array(gl).transpose()
                ax.imshow(img, aspect='auto')
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.xaxis.set_ticks(np.arange(0, len(gl), 50), [])

                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                self.fig.add_subplot(ax)

    def circular(self, bin_size, chroms=[], use_mask_rd=True):
        n = len(self.plot_files)
        ix = self.plot_files
        snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
        rd_flag = FLAG_GC_CORR | (FLAG_USEMASK if use_mask_rd else 0)
        if self.grid == "auto":
            sx, sy = self.panels_shape(n)
        else:
            sx, sy = self.grid
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(sy * 8)
            self.fig.set_figwidth(sx * 8)
        grid = gridspec.GridSpec(sy, sx, wspace=0.2, hspace=0.2)
        for i in range(n):
            ax = self.fig.add_subplot(grid[i], projection='polar')
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
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
            rd_mean = io.get_signal(None, bin_size, "RD stat", FLAG_AUTO)[4]
            tl = 0
            dt = 2.0 * np.pi / plot_len
            theta = np.arange(0, 2.0 * np.pi, dt)
            angles = []
            labels = []
            for j in range(len(plot_chroms)):
                c, l = plot_chroms[j]
                rd_color = self.palette1[j % len(self.palette1)]
                snp_color = self.palette2[j % len(self.palette2)]
                rd = io.get_signal(c, bin_size, "RD", rd_flag)
                maf = io.get_signal(c, bin_size, "SNP maf", snp_flag)
                plt.polar(theta[tl:tl + maf.size], 1 - maf, color=snp_color, linewidth=0.3)
                plt.fill_between(theta[tl:tl + maf.size], 1 - maf, np.ones_like(maf), color=snp_color, alpha=0.8)
                plt.polar(theta[tl:tl + rd.size], rd / (3. * rd_mean), color=rd_color, linewidth=0.3)
                plt.fill_between(theta[tl:tl + rd.size], np.ones_like(rd) / 10., rd / (3. * rd_mean), color=rd_color,
                                 alpha=0.8)
                # ax.text(theta[tl + maf.size // 3], 0.8, c, fontsize=8)
                labels.append(Genome.canonical_chrom_name(c))
                angles.append(180 * theta[tl + rd.size // 2] / np.pi)
                tl += l // bin_size + 1
            ax.set_rmax(0.9)
            ax.set_rticks([])
            ax.set_thetagrids(angles, labels=labels, fontsize=10, weight="bold", color="black")
            ax.set_title(io.filename.split("/")[-1], loc="left", fontsize=10, weight="bold", color="black")
            ax.grid(False)
        plt.subplots_adjust(bottom=0.05, top=0.95, wspace=0.2, hspace=0.2, left=0.05, right=0.95)
        if self.output_filename != "":
            plt.savefig(self.image_filename("circular"), dpi=200)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def rd_baf(self, hist=True):
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(12, 8), facecolor='w', edgecolor='k')
        n = len(self.plot_files)
        ix = self.plot_files
        if self.grid == "auto":
            sx, sy = self.panels_shape(n)
        else:
            sx, sy = tuple(self.grid)
        grid = gridspec.GridSpec(sy, sx, wspace=0.2, hspace=0.2)
        bin_size = self.bin_size
        for i in range(n):
            ax = self.fig.add_subplot(grid[i])
            io = self.io[ix[i]]


            chroms = []
            snp_flag = (FLAG_USEMASK if self.snp_use_mask else 0) | (FLAG_USEID if self.snp_use_id else 0)
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
            x=[]
            y=[]
            for c, l in chroms:
                flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO

                likelihood = io.get_signal(c, bin_size, "SNP likelihood call", snp_flag)
                segments_baf = segments_decode(io.get_signal(c, bin_size, "SNP likelihood segments", snp_flag))
                rd = io.get_signal(c, bin_size, "RD mosaic call", rd_flag)
                segments_rd = segments_decode(io.get_signal(c, bin_size, "RD mosaic segments", rd_flag))

                mbaf={}
                mrd={}
                for s, lh in zip(segments_baf, likelihood):
                    b, p = likelihood_baf_pval(lh)
                    for pos in s:
                        mbaf[pos]=0.5-b
                for s, r in zip(segments_rd, rd[0]):
                    for pos in s:
                        mrd[pos]=r
                for p in mbaf:
                    if p in mrd:
                        x.append(mbaf[p])
                        y.append(mrd[p])

            if hist:
                from matplotlib.colors import LogNorm
                ax.hist2d(x,y,bins=[np.arange(0,0.51,0.01),np.arange(0,max(y),max(y)/100.)],norm=LogNorm())
            else:
                ax.scatter(x,y,marker=".",alpha=0.5)

        if self.output_filename != "":
            plt.savefig(self.image_filename("regions"), dpi=150)
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
            plt.savefig(self.image_filename("dispersion"), dpi=200)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def compare(self, region1, region2, n_bins=21, plot=False, legend=True):
        n = len(self.plot_files)
        ix = self.plot_files
        if plot:
            plt.clf()
            plt.rcParams["font.size"] = 8
            if self.grid == "auto":
                sx, sy = self.panels_shape(n)
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
                ax.set_title(io.filename, position=(0.01, 1.07),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'})
            stat = io.get_signal(None, self.bin_size, "RD stat", FLAG_AUTO)
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

            print("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%e\t%.4f\t%.4f" % (
                io.filename, region1, region2, fitm1, fits1, fitm2, fits2, pval, fitm1 / fitm2,
                fitm1 / fitm2 * (fits1 / fitm1 / np.sqrt(sum(hist1)) + fits2 / fitm2 / np.sqrt(sum(hist2)))))

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
                plt.savefig(self.image_filename("comp"), dpi=200)
                plt.close(self.fig)
            elif self.interactive:
                plt.show(block=False)
                plt.draw()
            else:
                plt.show()

    def genotype(self, bin_sizes, region):
        regs = decode_region(region)
        for c, (pos1, pos2) in regs:
            print(c + ":" + str(pos1) + ":" + str(pos2), end="")
            for bs in bin_sizes:
                flag_rd = FLAG_GC_CORR if self.rd_use_gc_corr else o
                stat = self.io[self.plot_file].get_signal(c, bs, "RD stat", flag_rd)
                his_p = self.io[self.plot_file].get_signal(c, bs, "RD", flag_rd)
                bin1 = (pos1 - 1) // bs
                bin2 = (pos2 - 1) // bs
                rc = 0
                if bin1 == bin2:
                    rc = (pos2 - pos1 + 1) * his_p[bin1] / bs
                else:
                    rc += (bin1 * bs - pos1 + 1 + bs) * his_p[bin1] / bs
                    rc += (pos2 - bin2 * bs) * his_p[bin1] / bs
                    for ix in range(bin1 + 1, bin2):
                        rc += his_p[ix]
                print("\t%f" % rc * stat[4] * (pos2 - pos1 + 1) / bs, end="")
        print()

    def genotype_prompt(self, bin_sizes=[]):
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
                self.genotype(bin_sizes, line)


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
