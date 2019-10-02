""" cnvpytor.viewer

Class Viewer: ploting CNVpytor data
"""

from .io import *
from .utils import *
from .genome import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import logging

_logger = logging.getLogger("cnvpytor.viewer")


class Viewer:
    def __init__(self, files, plot_file):
        _logger.debug("Viewer class init: files [%s], png_prefix '%s'." % (", ".join(files), plot_file))
        self.io = [IO(f) for f in files]
        self.plot_file = plot_file
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

    def image_filename(self, sufix):
        parts = self.plot_file.split(".")
        if parts[-1] != "png" and parts[-1] != "pdf" and parts[-1] != "jpg" and parts[-1] != "eps" and parts[
            -1] != "svg":
            _logger.warning("File extension should be: .jpg, .png, .svg, .eps or .pdf")
            exit(0)
        parts[-1] = sufix + "." + parts[-1]
        return ".".join(parts)

    def parse(self, args):
        current = "gview"
        regions = []
        for p in args.plot:
            if p.isdigit() and (int(p) % 100) == 0:
                if current == "gview":
                    self.gview(int(p), args.use_mask_with_rd)
                elif current == "manhattan":
                    self.manhattan(int(p), args.use_mask_with_rd)
                elif current == "stat":
                    self.stat(int(p))
                elif current == "regions":
                    self.regions(int(p), regions, panels=args.panels)
                    regions = []
            elif p == "rdstat":
                self.stat()
            elif p == "baf":
                self.baf()
            elif current == "regions":
                regions.append(p)
            else:
                current = p

    def stat(self, his_bin_size=100):
        auto = self.io[0].signal_exists(None, his_bin_size, "RD stat", FLAG_AUTO)
        sex = self.io[0].signal_exists(None, his_bin_size, "RD stat", FLAG_SEX)
        mt = self.io[0].signal_exists(None, his_bin_size, "RD stat", FLAG_MT) and (his_bin_size < 1001)
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
                stat = self.io[0].get_signal(None, his_bin_size, "RD stat", flag)
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
                his_p = self.io[0].get_signal(None, his_bin_size, "RD p dist", flag)
                his_u = self.io[0].get_signal(None, his_bin_size, "RD u dist", flag)
                his_rd_gc = self.io[0].get_signal(None, his_bin_size, "RD GC dist", flag)
                gc_corr = self.io[0].get_signal(None, his_bin_size, "GC corr", flag)
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
        if self.plot_file != "":

            plt.savefig(self.image_filename("stat"), dpi=150)
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
            stat = self.io[0].get_signal(None, bin_size, "RD stat", flag)
            if stat is None:
                _logger.error("Data for bin size %d is missing in file '%s'!" % (bin_size, self.io[0].filename))
                exit(0)
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
            ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * stat[4], [])
            ax.xaxis.set_ticks(np.arange(0, (l + 10e6) // bin_size, 10e6 // bin_size), [])
            ax.set_ylim([0, max(3. * stat[4], stat[4] + 5. * stat[5])])
            n_bins = l // bin_size
            ax.set_xlim([-n_bins * 0.05, n_bins * 1.05])
            ax.grid()

            plt.step(his_p, "grey")
            plt.step(his_p_corr, "k")
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.plot_file != "":
            plt.savefig(self.image_filename("gview"), dpi=150)
            plt.close(fig)
        else:
            plt.show()

    def baf(self):
        plt.rcParams["font.size"] = 8
        fig = plt.figure(1, figsize=(12, 8), dpi=90, facecolor='w', edgecolor='k')
        chroms = []
        if self.reference_genome is None:
            chroms = self.io[0].snp_chromosomes()
        else:
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                snp_chr = self.io[0].snp_chromosome_name(c)
                if self.io[0].signal_exists(snp_chr, None, "SNP pos", 0) and \
                        self.io[0].signal_exists(snp_chr, None, "SNP desc", 0) and \
                        self.io[0].signal_exists(snp_chr, None, "SNP counts", 0) and \
                        self.io[0].signal_exists(snp_chr, None, "SNP qual", 0) and \
                        (Genome.is_autosome(c) or Genome.is_sex_chrom(c)):
                    chroms.append(snp_chr)
        sx, sy = self.panels_shape(len(chroms))
        ix = 1
        for c in chroms:
            pos, ref, alt, nref, nalt, gt, flag, qual = self.io[0].read_snp(c)
            hpos = []
            baf = []
            for i in range(len(pos)):
                if (gt[i] == 1 or gt[i] == 5 or gt[i] == 6) and (nref[i] + nalt[i]) != 0:
                    hpos.append(pos[i])
                    if gt[i] % 4 == 1:
                        baf.append(1.0 * nalt[i] / (nref[i] + nalt[i]))
                    else:
                        baf.append(1.0 * nref[i] / (nref[i] + nalt[i]))

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
            plt.scatter(hpos, baf, marker='.', c=['k' if f // 2 else 'y' for f in flag], s=1.5, alpha=1)
            ix += 1
        plt.subplots_adjust(bottom=0., top=1., wspace=0, hspace=0, left=0., right=1.)
        if self.plot_file != "":
            plt.savefig(self.image_filename("gview"), dpi=150)
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
            stat = self.io[0].get_signal(None, bin_size, "RD stat", flag)
            if stat[4] > max_m:
                max_m = stat[4]
            flag_rd = 0
            if use_mask:
                flag_rd = FLAG_USEMASK
            his_p = self.io[0].get_signal(c, bin_size, "RD", flag_rd)
            his_p_corr = self.io[0].get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
            pos = range(apos, apos + len(his_p))
            ax.text(apos + len(his_p) // 2, stat[4] // 10, Genome.canonical_chrom_name(c),
                    fontsize=12, verticalalignment='bottom', horizontalalignment='center', )
            plt.plot(pos, his_p_corr, ls='', marker='.')
            apos += len(his_p)
            xticks.append(apos)

        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * max_m, [])
        ax.xaxis.set_ticks(xticks, [])
        ax.set_ylim([0, 2 * max_m])
        n_bins = apos
        ax.set_xlim([0, n_bins])
        ax.grid()
        plt.subplots_adjust(bottom=0.05, top=0.95, wspace=0, hspace=0, left=0.05, right=0.95)

        if self.plot_file != "":
            plt.savefig(self.image_filename("manhattan"), dpi=150)
            plt.close(fig)
        else:
            plt.show()

    def regions(self, bin_size, regions, panels=["rd"], use_mask_rd=False, sep_color="g"):
        plt.rcParams["font.size"] = 8
        fig = plt.figure(1, figsize=(12, 8), dpi=90, facecolor='w', edgecolor='k')
        grid = gridspec.GridSpec(len(self.io), len(regions), wspace=0.2, hspace=0.2)
        ix = 0
        for i in self.io:
            for r in regions:
                self.region(i, fig, grid[ix], bin_size, r, panels=panels, use_mask_rd=use_mask_rd, sep_color=sep_color)
                ix += 1
        plt.subplots_adjust(bottom=0.05, top=0.95, wspace=0, hspace=0, left=0.05, right=0.95)

        if self.plot_file != "":
            plt.savefig(self.image_filename("regions"), dpi=150)
            plt.close(fig)
        else:
            plt.show()

    def region(self, io, fig, element, bin_size, region, panels=["rd"], use_mask_rd=False, sep_color="g"):
        grid = gridspec.GridSpecFromSubplotSpec(len(panels), 1, subplot_spec=element, wspace=0, hspace=0.1)
        r = decode_region(region)
        for i in range(len(panels)):
            ax = fig.add_subplot(grid[i])
            if i==0:
                ax.set_title(region, position=(0.01, 0.9), fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
                         color='C0')
            g_p = []
            g_p_corr = []
            if panels[i] == "rd":
                mean, stdev = 0, 0
                borders = []
                for c, (pos1, pos2) in r:
                    flag = FLAG_MT if Genome.is_mt_chrom(c) else FLAG_SEX if Genome.is_sex_chrom(c) else FLAG_AUTO
                    flag_rd = 0
                    if use_mask_rd:
                        flag_rd = FLAG_USEMASK
                    stat = io.get_signal(None, bin_size, "RD stat", flag)
                    mean = stat[4]
                    stdev = stat[5]
                    his_p = io.get_signal(c, bin_size, "RD", flag_rd)
                    his_p_corr = io.get_signal(c, bin_size, "RD", flag_rd | FLAG_GC_CORR)
                    start_bin = (pos1 - 1) / bin_size
                    end_bin = pos2 / bin_size
                    g_p.extend(list(his_p[start_bin:end_bin]))
                    g_p_corr.extend(list(his_p_corr[start_bin:end_bin]))
                    borders.append(len(g_p))

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks(np.arange(0, 3, 0.5) * mean, [])
                l = len(g_p)
                ax.xaxis.set_ticks(np.arange(0, l, 10), [])
                ax.set_ylim([0, max(3. * mean, mean + 5. * stdev)])
                ax.set_xlim([-l * 0.0, l * 1.0])
                ax.grid()
                ax.step(g_p, "grey")
                ax.step(g_p_corr, "k")
                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                fig.add_subplot(ax)
            elif panels[i] == "baf":
                borders = []
                hpos = []
                baf = []
                color = []
                start_pos = 0
                for c, (pos1, pos2) in r:
                    pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c)
                    ix = 0
                    last_hpos = -1
                    while ix < len(pos) and pos[ix] <= pos2:
                        if pos[ix] >= pos1 and (gt[ix] == 1 or gt[ix] == 5 or gt[ix] == 6) and (nref[ix] + nalt[ix]) != 0:
                            hpos.append(start_pos + pos[ix] - pos1)
                            last_hpos = start_pos + pos[ix] - pos1
                            if gt[i] % 4 == 1:
                                baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            else:
                                baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                            color.append('k' if flag[i] // 2 else 'y')
                        ix += 1
                    start_pos += pos2 - pos1
                    borders.append(start_pos)

                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.yaxis.set_ticks([0, 0.25, 0.5, 0.75, 1.0], [])
                l = max(hpos)
                # ax.xaxis.set_ticks(np.arange(0, (l + 10e6), 10e6), [])
                ax.set_ylim([0., 1.])
                ax.set_xlim([0, borders[-1]])
                ax.grid()
                ax.scatter(hpos, baf, marker='.', c=['k' if f // 2 else 'y' for f in flag], s=1.5, alpha=1)

                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                fig.add_subplot(ax)

            elif panels[i] == "likelihood":
                borders = []
                gl = []
                for c, (pos1, pos2) in r:
                    likelihood = io.get_signal(c, bin_size, "SNP likelihood")
                    start_bin = (pos1 - 1) / bin_size
                    end_bin = pos2 / bin_size
                    gl.extend(list(likelihood[start_bin:end_bin]))
                    borders.append(len(gl))
                img = np.array(gl).transpose()
                ax.imshow(img, aspect='auto')
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.xaxis.set_ticks(np.arange(0, l, 50), [])
                ax.grid()

                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                fig.add_subplot(ax)
