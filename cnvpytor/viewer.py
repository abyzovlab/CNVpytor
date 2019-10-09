""" cnvpytor.viewer

Class Viewer: ploting CNVpytor data
"""
from __future__ import absolute_import, print_function, division

from .io import *
from .utils import *
from .genome import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import logging
import readline
import traceback

_logger = logging.getLogger("cnvpytor.viewer")


class Viewer:
    def __init__(self, files, output_filename):
        _logger.debug("Viewer class init: files [%s], png_prefix '%s'." % (", ".join(files), output_filename))
        self.io = [IO(f, ro=True) for f in files]
        self.output_filename = output_filename
        self.set_style('classic')
        self.io_gc = self.io[0]
        self.io_mask = self.io[0]
        self.reference_genome = None
        self.interactive = False
        self.plot_files = [True for i in files]
        self.fig = None
        self.xkcd = False
        self.grid = "auto"
        self.plot_file = 0
        self.use_mask = True
        self.use_id = False
        self.palette1 = ["#555555", "#aaaaaa"]
        self.palette2 = ["#00ff00", "#0000ff"]
        self.baf_colors = {(0, 0): "yellow", (0, 1): "orange", (1, 0): "cyan", (1, 1): "blue", (2, 0): "lime",
                           (2, 1): "green", (3, 0): "yellow", (3, 1): "orange"}
        if self.io[0].signal_exists(None, None, "reference genome"):
            rg_name = np.array(self.io[0].get_signal(None, None, "reference genome")).astype("str")[0]
            self.reference_genome = Genome.reference_genomes[rg_name]
            if "mask_file" in Genome.reference_genomes[rg_name]:
                self.io_mask = IO(Genome.reference_genomes[rg_name]["mask_file"], ro=True, buffer=True)
            if "gc_file" in Genome.reference_genomes[rg_name]:
                self.io_gc = IO(Genome.reference_genomes[rg_name]["gc_file"], ro=True, buffer=True)

    def parse(self, command, args):
        current = "regions"
        regions = []

        for p in command:
            if p.isdigit() and (int(p) % 100) == 0:
                if current == "rd":
                    self.rd(int(p), args.use_mask_with_rd)
                if current == "likelihood":
                    self.likelihood(int(p))
                elif current == "manhattan":
                    self.manhattan(int(p), args.use_mask_with_rd)
                elif current == "stat":
                    self.stat(int(p))
                elif current == "circular":
                    self.circular(int(p), args.chrom, args.use_mask_with_rd)
                elif current == "regions":
                    self.multiple_regions(int(p), regions, panels=args.panels)
                    regions = []
            elif p == "rdstat":
                self.stat()
            elif p == "baf":
                self.baf()
            elif p in ["rd", "manhattan", "stat", "regions", "likelihood", "circular"]:
                current = p
            elif current == "regions":
                regions.append(p)
            else:
                current = p

    def plot(self, args):
        self.interactive = False
        self.use_mask = not args.no_mask
        self.use_id = args.use_id
        self.parse(args.plot, args)

    def prompt(self, bin_size, args):
        self.interactive = True
        self.use_mask = not args.no_mask
        self.use_id = args.use_id
        command_tree = {"set": {"bin_size": None,
                                "panels": None,
                                "use_mask_rd": None,
                                "plot_files": None,
                                "plot_file": None,
                                "grid": None,
                                "use_mask": None,
                                "use_id": None,
                                "xkcd": None,
                                "style": {}
                                },
                        "unset": {"use_mask_rd": None,
                                  "xkcd": None
                                  },
                        "save": None, "show": None, "quit": None, "rd": None, "likelihood": None, "baf": None,
                        "stat": None, "rdstat": None, "circular": None, "manhattan": None
                        }
        chromosomes = set({})
        for f in self.io:
            chromosomes = chromosomes.union(set(f.rd_chromosomes()))
            chromosomes = chromosomes.union(set(f.snp_chromosomes()))
        for c in chromosomes:
            command_tree[c] = None
        command_tree["set"]["style"] = dict(zip(plt.style.available, [None] * len(plt.style.available)))
        readline.parse_and_bind("tab: complete")
        completer = PromptCompleter(command_tree)
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
                elif f[0] == "show":
                    if n == 1:
                        print("    Parameters")
                        print("        * bin_size:", bin_size)
                        print("        * panels:", args.panels)
                        print("        * use_mask_rd:", args.use_mask_with_rd)
                        print("        * use_mask:", self.use_mask)
                        print("        * use_id:", self.use_id)
                        print("        * plot_files:")
                        for fi, fn, fs in list(
                                zip(range(len(self.plot_files)), map(lambda x: x.filename, self.io), self.plot_files)):
                            print("            %2d" % fi, fn, fs)
                        print("        * plot_file:", self.plot_file)
                        print("        * grid:", self.grid)
                        # print("    Available plot styles:", ", ".join(plt.style.available))
                elif f[0] == "set":
                    if n > 2 and f[1] == "bin_size":
                        try:
                            bin_size = binsize_type(f[2])
                        except (ArgumentTypeError, ValueError):
                            _logger.warning("bin_size should be intiger divisible by 100")
                    elif n > 3 and f[1] == "grid" and f[2].isdigit() and f[3].isdigit():
                        self.grid = (int(f[2]), int(f[3]))
                    elif n > 2 and f[1] == "plot_file" and f[2].isdigit() and int(f[2]) < len(self.io):
                        self.plot_file = int(f[2])
                    elif n > 2 and f[1] == "plot_files":
                        select = [int(i) for i in f[2:] if i.isdigit()]
                        self.plot_files = list([ix in select for ix in range(len(self.io))])
                    elif n > 2 and f[1] == "panels":
                        args.panels = f[2:]
                    elif n > 2 and f[1] == "style":
                        self.set_style(f[2])
                        if self.fig:
                            self.fig.canvas.draw()
                    elif n > 1 and f[1] == "use_mask_rd":
                        args.use_mask_with_rd = True
                    elif n > 1 and f[1] == "use_mask":
                        self.use_mask = True
                    elif n > 1 and f[1] == "use_id":
                        self.use_id = True
                    elif n > 1 and f[1] == "xkcd":
                        self.set_style("classic")
                        from matplotlib import patheffects
                        plt.xkcd()
                        plt.rcParams["path.effects"] = [patheffects.withStroke(linewidth=0.5, foreground="w")]
                        self.xkcd = True
                    else:
                        print("Unrecognized set argument!")
                elif f[0] == "unset":
                    if n > 1 and f[1] == "use_mask_rd":
                        args.use_mask_with_rd = False
                    elif n > 1 and f[1] == "use_mask":
                        self.use_mask = False
                    elif n > 1 and f[1] == "use_id":
                        self.use_id = False
                    elif n > 1 and f[1] == "xkcd" and self.xkcd:
                        plt.rcdefaults()
                        self.set_style('classic')
                        self.xkcd = False

                    else:
                        print("Unrecognized unset argument!")
                else:

                    try:
                        if f[0] not in ["rdstat", "baf"]:
                            self.parse(f + [str(bin_size)], args)
                        else:
                            self.parse(f, args)
                        if len(pre) > 1:
                            fns = pre[1].strip().split(" ")
                            if fns[0] != "":
                                plt.savefig(fns[0], dpi=200)
                    except Exception as e:
                        print(traceback.format_exc())
        except (EOFError, KeyboardInterrupt):
            print()
            return

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

    def stat(self, his_bin_size=100):
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
                ax.set_ylim([0, 1.1])
                bins = range(0, max_rd, bin_size)
                x = np.arange(0, max_rd // bin_size * bin_size, 0.1 * bin_size)
                plt.plot(x, normal(x, 1, stat[4], stat[5]), "g-")
                x = np.array(bins)
                plt.plot(x[:len(his_u)], his_u / stat[3], "y*")
                plt.plot(x[:len(his_p)], his_p / stat[3], "b*")
                ix += 1
        plt.subplots_adjust(bottom=0.08, top=0.95, wspace=0.25, hspace=0, left=0.05 * 3 / n_cols, right=0.95)
        if self.output_filename != "":
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
        snp_flag = (FLAG_USEMASK if self.use_mask else 0) | (FLAG_USEID if self.use_id else 0)
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

    def baf(self):
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(12, 8), dpi=90, facecolor='w', edgecolor='k')
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
                    hpos.append(pos[i])
                    if gt[i] % 4 != 2:
                        baf.append(1.0 * nalt[i] / (nref[i] + nalt[i]))
                    else:
                        baf.append(1.0 * nref[i] / (nref[i] + nalt[i]))
                    color.append(self.baf_colors[(gt[i] % 4, flag[i] >> 1)])

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
            plt.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=0.7)
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

    def manhattan(self, bin_size, use_mask=False):
        plt.clf()
        if self.reference_genome is None:
            _logger.warning("Missing reference genome required for gview.")
            return
        n = self.plot_files.count(True)
        ix = [x for x in range(len(self.plot_files)) if self.plot_files[x]]

        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, facecolor='w', edgecolor='k')
        if self.output_filename != "":
            self.fig.set_figheight(2 * n)
            self.fig.set_figwidth(12)
        grid = gridspec.GridSpec(n, 1, wspace=0.2, hspace=0.2)
        for i in range(n):
            ax = self.fig.add_subplot(grid[i])
            io = self.io[ix[i]]

            chroms = []
            for c, (l, t) in self.reference_genome["chromosomes"].items():
                rd_chr = io.rd_chromosome_name(c)
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

        if self.output_filename != "":
            plt.savefig(self.image_filename("manhattan"), dpi=200)
            plt.close(self.fig)
        elif self.interactive:
            plt.show(block=False)
            plt.draw()
        else:
            plt.show()

    def multiple_regions(self, bin_size, regions, panels=["rd"], use_mask_rd=False, sep_color="g"):
        plt.clf()
        plt.rcParams["font.size"] = 8
        self.fig = plt.figure(1, figsize=(12, 8), facecolor='w', edgecolor='k')
        grid = gridspec.GridSpec(len(self.io), len(regions), wspace=0.2, hspace=0.2)
        ix = 0
        for i in self.io:
            for r in regions:
                self.regions(i, grid[ix], bin_size, r, panels=panels, use_mask_rd=use_mask_rd, sep_color=sep_color)
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

    def regions(self, io, element, bin_size, region, panels=["rd"], use_mask_rd=False, sep_color="g"):
        snp_flag = (FLAG_USEMASK if self.use_mask else 0) | (FLAG_USEID if self.use_id else 0)
        grid = gridspec.GridSpecFromSubplotSpec(len(panels), 1, subplot_spec=element, wspace=0, hspace=0.1)
        r = decode_region(region)
        for i in range(len(panels)):
            ax = self.fig.add_subplot(grid[i])
            if i == 0:
                ax.set_title(io.filename + ": " + region, position=(0.01, 0.9),
                             fontdict={'verticalalignment': 'top', 'horizontalalignment': 'left'},
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
                    start_bin = (pos1 - 1) // bin_size
                    end_bin = pos2 // bin_size
                    g_p.extend(list(his_p[start_bin:end_bin]))
                    g_p_corr.extend(list(his_p_corr[start_bin:end_bin]))
                    borders.append(len(g_p))

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
                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                self.fig.add_subplot(ax)
            elif panels[i] == "baf":
                borders = []
                hpos = []
                baf = []
                color = []
                start_pos = 0
                for c, (pos1, pos2) in r:
                    pos, ref, alt, nref, nalt, gt, flag, qual = io.read_snp(c)
                    ix = 0
                    while ix < len(pos) and pos[ix] <= pos2:
                        if pos[ix] >= pos1 and (nref[ix] + nalt[ix]) != 0:
                            hpos.append(start_pos + pos[ix] - pos1)
                            if gt[ix] % 4 != 2:
                                baf.append(1.0 * nalt[ix] / (nref[ix] + nalt[ix]))
                            else:
                                baf.append(1.0 * nref[ix] / (nref[ix] + nalt[ix]))
                            color.append(self.baf_colors[(gt[ix] % 4, flag[ix] >> 1)])
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
                ax.yaxis.grid()
                ax.scatter(hpos, baf, marker='.', edgecolor=color, c=color, s=10, alpha=0.7)

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
                    borders.append(len(gl))
                img = np.array(gl).transpose()
                ax.imshow(img, aspect='auto')
                ax.xaxis.set_ticklabels([])
                ax.yaxis.set_ticklabels([])
                ax.xaxis.set_ticks(np.arange(0, len(gl), 50), [])

                for i in borders[:-1]:
                    ax.axvline(i, color=sep_color, lw=1)
                self.fig.add_subplot(ax)

    def circular(self, bin_size, chroms=[], use_mask_rd=True):
        n = self.plot_files.count(True)
        ix = [x for x in range(len(self.plot_files)) if self.plot_files[x]]
        snp_flag = (FLAG_USEMASK if self.use_mask else 0) | (FLAG_USEID if self.use_id else 0)
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
