""" cnvpytor.viewer

Class Viewer: ploting CNVpytor data
"""

from .io import *
from .utils import *
from .genome import *
import matplotlib.pyplot as plt
import numpy as np
import logging

_logger = logging.getLogger("cnvpytor.viewer")


class Viewer:
    def __init__(self, files, png_prefix):
        _logger.debug("Viewer class init: files [%s], png_prefix '%s'." % (", ".join(files), png_prefix))
        self.io = [IO(f) for f in files]
        self.png_prefix = png_prefix
        self.set_style('seaborn-deep')
        self.io_gc = self.io[0]
        self.io_mask = self.io[0]
        self.reference_genome = None
        if self.io[0].signal_exists(None, None, "reference genome"):
            rg_name = np.array(self.io[0].get_signal(None, None, "reference genome")).astype("str")[0]
            self.reference_genome = Genome.reference_genomes[rg_name]
            if "mask_file" in Genome.reference_genomes[rg_name]:
                self.io_mask = IO(Genome.reference_genomes[rg_name]["mask_file"])
            if "gc_file" in Genome.reference_genomes[rg_name]:
                self.io_gc = IO(Genome.reference_genomes[rg_name]["gc_file"])

    @staticmethod
    def set_style(style):
        if style in plt.style.available:
            plt.style.use(style)

    def stat(self):
        auto = self.io[0].signal_exists(None, None, "RD stat", FLAG_AUTO)
        sex = self.io[0].signal_exists(None, None, "RD stat", FLAG_SEX)
        mt = self.io[0].signal_exists(None, None, "RD stat", FLAG_MT)
        if not (auto or sex or mt):
            return
        cond = [auto, sex, mt]
        stat_list = []
        n_cols = sum(map(int, cond))
        ix = 1
        plt.rcParams["font.size"] = 8
        fig = plt.figure(1, figsize=(4 * n_cols, 8), dpi=90, facecolor='w', edgecolor='k')
        for t, c, flag in zip(["Autosomes", "X/Y", "Mitochondria"], cond, [FLAG_AUTO, FLAG_SEX, FLAG_MT]):
            if c:
                stat = self.io[0].get_signal(None, None, "RD stat", flag)
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
                his_p = self.io[0].get_signal(None, None, "RD p dist", flag)
                his_u = self.io[0].get_signal(None, None, "RD u dist", flag)
                his_rd_gc = self.io[0].get_signal(None, None, "RD GC dist", flag)
                gc_corr = self.io[0].get_signal(None, None, "GC corr", flag)
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
                ax.set_ylim([0, 1.1])
                bins = range(0, max_rd, bin_size)
                x = np.arange(0, max_rd // bin_size * bin_size, 0.1 * bin_size)
                plt.plot(x, normal(x, 1, stat[4], stat[5]), "g-")
                x = np.array(bins)
                plt.plot(x[:len(his_u)], his_u / stat[3], "y*")
                plt.plot(x[:len(his_p)], his_p / stat[3], "b*")
                ix += 1
        plt.subplots_adjust(bottom=0.08, top=0.95, wspace=0.25, hspace=0, left=0.05 * 3 / n_cols, right=0.95)
        if self.png_prefix != "":
            plt.savefig(self.png_prefix + "_stat.png", dpi=150)
            plt.close(fig)
        else:
            plt.show()

    def rd(self, bin_size, chroms=[]):
        plt.rcParams["font.size"] = 8
        fig = plt.figure(1, figsize=(12, 8), dpi=90, facecolor='w', edgecolor='k')
        count = 0
        for c in self.io[0].chromosomes_with_signal(None, "RD p"):
            if len(chroms) == 0 or (c in chroms):
                count += 1
        sx = 1
        sy = 1
        while sx * sy < count:
            sx += 1
            sy = int(2. * sx / 3 + 1.)
        ix = 1
        for c in self.io[0].chromosomes_with_signal(None, "RD p"):
            if len(chroms) == 0 or (c in chroms):
                flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                stat = self.io[0].get_signal(None, None, "RD stat", flag)
                his_p = self.io[0].get_signal(c, bin_size, "RD", 0)
                his_p_corr = self.io[0].get_signal(c, bin_size, "RD", FLAG_GC_CORR)
                ax = plt.subplot(sx, sy, ix)
                ax.set_ylim([0, max(3. * stat[4], stat[4] + 5. * stat[5]) * bin_size / 100])
                plt.step(his_p, "grey")
                plt.step(his_p_corr, "k")
                ix += 1
        if self.png_prefix != "":
            plt.savefig(self.png_prefix + "_stat.png", dpi=150)
            plt.close(fig)
        else:
            plt.show()

    def gview(self, bin_size, use_mask):
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        plt.rcParams["font.size"] = 8
        fig = plt.figure(1, figsize=(12, 8), dpi=90, facecolor='w', edgecolor='k')
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[0].rd_chromosome_name(c)
            if self.io[0].signal_exists(rd_chr, bin_size, "RD", 0) and \
                    self.io[0].signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                chroms.append((rd_chr, l))
        sx, sy = self.panels_shape(len(chroms))
        ix = 1
        for c, l in chroms:
            flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
            stat = self.io[0].get_signal(None, None, "RD stat", flag)
            flag_rd = 0
            if use_mask:
                flag_rd = FLAG_USEMASK
            his_p = self.io[0].get_signal(c, bin_size, "RD", flag_rd)
            his_p_corr = self.io[0].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            ax = plt.subplot(sx, sy, ix)
            ax.set_title(c, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
            ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * stat[4] * bin_size / 100, [])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // bin_size, 10e6 // bin_size), [])
            ax.set_ylim([0, max(3. * stat[4], stat[4] + 5. * stat[5]) * bin_size / 100])
            n_bins = l // bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()

            plt.step(his_p, "grey")
            plt.step(his_p_corr, "k")
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.png_prefix != "":
            plt.savefig(self.png_prefix + "_stat.png", dpi=150)
            plt.close(fig)
        else:
            plt.show()

    def panels_shape(self, n):
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

    def manhattan(self, bin_size, use_mask):
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        chroms = []
        for c, (l, t) in self.reference_genome["chromosomes"].items():
            rd_chr = self.io[0].rd_chromosome_name(c)
            if self.io[0].signal_exists(rd_chr, bin_size, "RD", 0) and \
                    self.io[0].signal_exists(rd_chr, bin_size, "RD", FLAG_GC_CORR) and \
                    (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                chroms.append((rd_chr, l))

        ix = 1
        apos = 0
        xticks = [0]

        plt.rcParams["font.size"] = 8
        fig = plt.figure(1, figsize=(12, 4), dpi=90, facecolor='w', edgecolor='k')
        ax = plt.gca()
        max_m = 0
        for c, l in chroms:
            flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
            stat = self.io[0].get_signal(None, None, "RD stat", flag)
            if stat[4] > max_m:
                max_m = stat[4]
            flag_rd = 0
            if use_mask:
                flag_rd = FLAG_USEMASK
            his_p = self.io[0].get_signal(c, bin_size, "RD", flag_rd)
            his_p_corr = self.io[0].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            pos = range(apos, apos + len(his_p))
            ax.text(apos + len(his_p) // 2, int(stat[4]) * bin_size // 100 * 2 // 3, Genome.canonical_chrom_name(c),
                    fontsize=12, verticalalignment='top', horizontalalignment='center',)
            plt.plot(pos, his_p_corr, ls='', marker='.')
            apos += len(his_p)
            xticks.append(apos)

        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * max_m * bin_size / 100, [])
        ax.xaxis.set_ticks(xticks, [])
        ax.set_ylim([0, 2 * max_m * bin_size / 100])
        n_bins = apos
        ax.set_xlim([0, n_bins])
        ax.grid()
        plt.subplots_adjust(bottom=0.05, top=0.95, wspace=0, hspace=0, left=0.05, right=0.95)

        if self.png_prefix != "":
            plt.savefig(self.png_prefix + "_manhattan.png", dpi=150)
            plt.close(fig)
        else:
            plt.show()
