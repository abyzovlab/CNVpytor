"""
__main__

"""
from __future__ import print_function
from .root import *
from .viewer import *
from .version import __version__
from .fasta import *
from .export import *
from .trio import *
import sys
import os
import logging
import argparse
import matplotlib.pyplot as plt


def main():
    """ main()
    CNVpytor main commandline program.

    """
    parser = argparse.ArgumentParser(
        description="Lite version of the CNVnator written in Python.\nA tool for CNV discovery from depth of read mapping.")
    parser.add_argument('-version', '--version', action='store_true', help='show version number and exit')
    parser.add_argument('-root', '--root', type=str, nargs="+",
                        help="CNVnator hd5 file: data storage for all calculations", default=None)

    parser.add_argument('-download', '--download_resources', action='store_true', help='download resource files')

    parser.add_argument('-chrom', '--chrom', type=str, nargs="+", help="list of chromosomes to apply calculation",
                        default=[])
    parser.add_argument('-v', '--verbose', type=str,
                        choices=["none", "debug", "info", "warning", "error", "d", "e", "i", "w"],
                        help="verbose level: debug, info (default), warning, error", default="info")
    parser.add_argument('-log', '--log_file', type=str, help='log file')
    parser.add_argument('-j', '--max_cores', type=int,
                        help="maximal number of cores to use in calculation", default=8)
    parser.add_argument('-rd', '--rd', nargs="+", type=str, help="read bam/sam/cram and store read depth information")
    parser.add_argument('-T', '--reference_filename', type=str, help="reference fasta for CRAM")

    parser.add_argument('-gc', '--gc', type=str, help="read fasta file and store GC/AT content")
    parser.add_argument('-cgc', '--copy_gc', type=str, help="copy GC/AT content from another cnvnator file")
    parser.add_argument('-his', '--his', type=binsize_type, nargs="+",
                        help="create histograms for specified bin size (multiple bin sizes separate by space)")
    parser.add_argument('-snp2his', '--his_from_snp', type=binsize_type, nargs="+",
                        help="create histograms for specified bin size (multiple bin sizes separate by space)")
    parser.add_argument('-stat', '--stat', type=binsize_type, nargs="+",
                        help="calculate statistics for specified bin size (multiple bin sizes separate by space)")
    parser.add_argument('-partition', '--partition', type=binsize_type, nargs="+",
                        help="calculate segmentation for specified bin size (multiple bin sizes separate by space)")
    parser.add_argument('-call', '--call', type=str, nargs="+",
                        help="CNV caller: [baf] bin_size [bin_size2 ...] (multiple bin sizes separate by space)")
    parser.add_argument('-vcf', '-snp', '--vcf', nargs="+", type=str, help="read SNP data from vcf files")
    parser.add_argument('-somatic_snv', '--somatic_snv', nargs="+", type=str, help="read SNP data from vcf files")

    parser.add_argument('-minc', '--min_count', type=int,
                        help="minimal count of haterozygous SNPs", default=None)
    parser.add_argument('-vcf2rd', '--rd_from_vcf', type=str, help="read SNP data from vcf files")
    parser.add_argument('-noAD', '--no_snp_counts', action='store_true',
                        help="read positions of variants, not counts (AD tag)")
    parser.add_argument('-nofilter', '--no_filter', action='store_true',
                        help="read all variants (not only PASS)")
    parser.add_argument('-ad', '--ad_tag', type=str, help="counts tag (default: AD)", default="AD")
    parser.add_argument('-gt', '--gt_tag', type=str, help="genotype tag (default: GT)", default="GT")
    parser.add_argument('-dp', '--dp_tag', type=str, help="read depth tag (default: DP)", default="DP")
    parser.add_argument('-callset', '--callset', type=str, help="name for somatic VCF signal", default=None)
    parser.add_argument('-maxcn', '--max_copy_number', type=int, help="maximal copy number", default=10)
    parser.add_argument('-mindbaf', '--baf_threshold', type=float, help="threshold for change in BAF level",
                        default=0.0)
    parser.add_argument('-mincf', '--min_cell_fraction', type=float, help="minimal cell fraction", default=0.0)

    parser.add_argument('-pileup', '--pileup_bam', nargs="+", type=str, help="calculate SNP counts from bam files")
    parser.add_argument('-snp2rd', '--rd_from_snp', action='store_true', help="calculate RD from SNP counts")
    parser.add_argument('-sbin', '--s_bin_size', type=binsize_type, help="Super bin size (use with -snp2rd)",
                        default=10000)

    parser.add_argument('-mask', '--mask', type=str, help="read fasta mask file and flag SNPs in P region")
    parser.add_argument('-mask_snps', '--mask_snps', action='store_true', help="flag SNPs in P region")
    parser.add_argument('-trio_phase', '--trio_phase', action='store_true', help="Phase trio")
    parser.add_argument('-parents', '--phase_parents', action='store_true', help="Phase parents")
    parser.add_argument('-mask_snvs', '--mask_snvs', type=str, help="flag SNVs in P region")
    parser.add_argument('-idvar', '--idvar', type=str, help="read vcf file and flag SNPs that exist in database file")
    parser.add_argument('-random_phase', '--random_phase', action='store_true', help="randomly phase SNPs")
    parser.add_argument('-baf', '--baf', type=binsize_type, nargs="+",
                        help="create BAF histograms for specified bin size (multiple bin sizes separate by space)")
    parser.add_argument('-nomask', '--no_mask', action='store_true', help="do not use P mask in BAF histograms")
    parser.add_argument('-useid', '--use_id', action='store_true', help="use id flag filtering in SNP histograms")
    parser.add_argument('-usehom', '--use_hom', action='store_true', help="use hom")
    parser.add_argument('-usephase', '--use_phase', action='store_true',
                        help="use information about phase while processing SNP data")
    parser.add_argument('-reducenoise', '--reduce_noise', action='store_true',
                        help="reduce noise in processing SNP data")
    parser.add_argument('-blw', '--baf_likelihood_width', type=float,
                        help="likelihood width used in processing SNP data (default=0.8)", default=0.8)
    parser.add_argument('-altc', '--alt_corr', action='store_true',
                        help="Remove alt/ref bias")

    parser.add_argument('-plot', '--plot', type=str, nargs="+", help="plotting")
    parser.add_argument('-view', '--view', type=binsize_type,
                        help="Enters interactive ploting mode")
    parser.add_argument('-agg', '--force_agg', action='store_true', help="Force Agg matplotlib backend")

    parser.add_argument('-panels', '--panels', type=str, nargs="+", default=["rd"], choices=["rd", "baf", "likelihood"],
                        help="plot panels (with -plot regions)")

    parser.add_argument('-style', '--plot_style', type=str,
                        help="available plot styles: " + ", ".join(plt.style.available), choices=plt.style.available)
    parser.add_argument('-o', '--plot_output_file', type=str, help="output filename prefix and extension", default="")
    parser.add_argument('-anim', '--animation', type=str, help="animation folder/prefix", default="")

    parser.add_argument('-make_gc_file', '--make_gc_genome_file', action='store_true',
                        help="used with -gc will create genome gc file")
    parser.add_argument('-make_mask_file', '--make_mask_genome_file', action='store_true',
                        help="used with -mask will create genome mask file")
    parser.add_argument('-rd_use_mask', '--use_mask_with_rd', action='store_true', help="used P mask in RD histograms")
    parser.add_argument('-nogc', '--no_gc_corr', action='store_true', help="do not use GC correction in RD histograms")
    parser.add_argument('-rg', '--reference_genome', type=str, help="Manually set reference genome", default=None)
    parser.add_argument('-sample', '--vcf_sample', type=str, help="Sample name in vcf file", default="")
    parser.add_argument('-conf', '--reference_genomes_conf', type=str, help="Configuration with reference genomes",
                        default=None)

    parser.add_argument('-ls', '--ls', action='store_true', help='list pytor file(s) content')
    parser.add_argument('-gc_info', '--gc_info', action='store_true', help='list pytor file(s) gc content stat')
    parser.add_argument('-info', '--info', type=binsize_type, nargs="*", help='print statistics for pythor file(s)')
    parser.add_argument('-qc', '--qc', type=binsize_type, nargs="*", help='print quality control statistics')
    parser.add_argument('-rdqc', '--rd_qc', type=binsize_type, nargs="*",
                        help='print quality control statistics without SNP data')
    parser.add_argument('-comp', '--compare', type=str, nargs="*", help='compere two regions: -comp reg1 reg2 [n_bins]')
    parser.add_argument('-genotype', '--genotype', type=str, nargs="*")
    parser.add_argument('-a', '--all', action='store_true', help='Genotype with all columns')
    parser.add_argument('-meta', '--metadata', action='store_true', help='list Metadata')
    parser.add_argument('-fasta2rg', '--reference_genome_template', type=str,
                        help="create template for reference genome using chromosome lengths from fasta file")
    parser.add_argument('-export', '--export', type=str, nargs="*", help='Export to jbrowse and cnvnator')
    args = parser.parse_args(sys.argv[1:])

    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    if args.verbose in {"debug", "d"}:
        level = logging.DEBUG
    elif args.verbose in {"info", "i"}:
        level = logging.INFO
    elif args.verbose in {"warning", "w"}:
        level = logging.WARNING
    elif args.verbose in {"error", "e"}:
        level = logging.ERROR
    else:
        level = logging.CRITICAL

    if args.log_file:
        logging.basicConfig(filename=args.log_file, level=logging.DEBUG, format=log_format)
        logger = logging.getLogger('cnvpytor')
        ch = logging.StreamHandler()
        formatter = logging.Formatter(log_format)
        ch.setFormatter(formatter)
        ch.setLevel(level)
        logger.addHandler(ch)
    else:
        logging.basicConfig(level=level, format=log_format)
        logger = logging.getLogger('cnvpytor')
    logger.debug("Start logging...")

    if args.reference_genome_template is not None:
        Fasta(args.reference_genome_template).print_reference_genome_template()

    if args.download_resources:
        Genome.download_resources()
        return 0

    if not Genome.check_resources():
        logger.error("Some reference genome resource files are missing. "
                     "Run 'cnvpytor -download' as same user who has installed cnvpytor.")
        return 0

    if args.version:
        print('CNVpytor {}'.format(__version__))
        return 0

    if args.reference_genomes_conf:
        Genome.load_reference_genomes(args.reference_genomes_conf)
    elif os.path.exists(os.path.expanduser('~/.cnvpytor/reference_genomes_conf.py')):
        Genome.load_reference_genomes(os.path.expanduser('~/.cnvpytor/reference_genomes_conf.py'))

    if args.root is not None:

        if args.ls:
            show = Show(args.root)
            show.ls()

        if args.gc_info:
            show = Show(args.root)
            show.gc_info()

        if args.export:
            if len(args.export) > 0:
                dir_name_list = args.export[1:]
                dir_name = ''
                if len(dir_name_list) > 0:
                    dir_name = dir_name_list[0]
                export_program = args.export[0].lower()
                if export_program in ['jbrowse', 'cnvnator']:
                    if export_program == 'jbrowse':
                        export_j = ExportJBrowse(args.root, dir_name)
                        export_j.create_reference_json()
                        export_j.rd_signal()
                        export_j.snp_signal()
                        export_j.create_tracklist_json()
                    elif export_program == 'cnvnator':
                        logger.info("Under Development")
                else:
                    logger.error("Incorrect export program name")

        if args.metadata:
            show = Show(args.root)
            show.meta()

        if args.info is not None:
            show = Show(args.root)
            show.info(args.info)


        if args.genotype is not None:
            params = {"output_filename": args.plot_output_file,
                      "chrom": args.chrom,
                      "panels": args.panels,
                      "snp_use_mask": not args.no_mask,
                      "snp_use_id": args.use_id,
                      "rd_use_mask": args.use_mask_with_rd
                      }
            view = Viewer(args.root, params, force_agg=args.force_agg)
            view.genotype_prompt(list(map(binsize_type, args.genotype)), all=args.all)

        if args.qc is not None:
            params = {"bin_size": binsize_type(args.qc[-1]),
                      "chrom": args.chrom,
                      "snp_use_mask": not args.no_mask,
                      "snp_use_id": args.use_id,
                      "rd_use_mask": args.use_mask_with_rd,
                      "rd_use_gc_corr": not args.no_gc_corr
                      }
            view = Viewer(args.root, params, force_agg=args.force_agg)
            view.qc()

        if args.rd_qc is not None:
            params = {"bin_size": binsize_type(args.rd_qc[-1]),
                      "chrom": args.chrom,
                      "snp_use_mask": not args.no_mask,
                      "snp_use_id": args.use_id,
                      "rd_use_mask": args.use_mask_with_rd,
                      "rd_use_gc_corr": not args.no_gc_corr
                      }
            view = Viewer(args.root, params, force_agg=args.force_agg)
            view.qc(snp_qc=False)


        if args.compare is not None:
            params = {"bin_size": binsize_type(args.compare[-1]),
                      "rd_use_gc_corr": not args.no_gc_corr,
                      "rd_use_mask": args.use_mask_with_rd
                      }
            view = Viewer(args.root, params, force_agg=args.force_agg)
            if len(args.compare) == 3:
                view.compare(args.compare[0], args.compare[1])
            elif len(args.compare) == 4:
                view.compare(args.compare[0], args.compare[1], int(args.compare[2]))

        if args.rd:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.rd(args.rd, chroms=args.chrom, reference_filename=args.reference_filename)

        if args.reference_genome:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.set_reference_genome(args.reference_genome)

        if args.plot:
            params = {"output_filename": args.plot_output_file,
                      "chrom": args.chrom,
                      "panels": args.panels,
                      "snp_use_mask": not args.no_mask,
                      "snp_use_id": args.use_id,
                      "rd_use_mask": args.use_mask_with_rd,
                      "rd_use_gc_corr": not args.no_gc_corr
                      }
            if args.plot_style:
                params["style"] = args.plot_style
            view = Viewer(args.root, params)
            view.plot_command(args.plot)

        if args.view:
            params = {"bin_size": args.view,
                      "output_filename": args.plot_output_file,
                      "chrom": args.chrom,
                      "panels": args.panels,
                      "snp_use_mask": not args.no_mask,
                      "snp_use_id": args.use_id,
                      "rd_use_mask": args.use_mask_with_rd,
                      "rd_use_gc_corr": not args.no_gc_corr
                      }
            if args.plot_style:
                params["style"] = args.plot_style
            view = Viewer(args.root, params, force_agg=args.force_agg)
            view.prompt()

        if args.gc:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.gc(args.gc, chroms=args.chrom, make_gc_genome_file=args.make_gc_genome_file)

        if args.copy_gc:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.copy_gc(args.copy_gc, chroms=args.chrom)

        if args.vcf:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.vcf(args.vcf, chroms=args.chrom, sample=args.vcf_sample, no_counts=args.no_snp_counts,
                    ad_tag=args.ad_tag, gt_tag=args.gt_tag, filter=not args.no_filter)

        if args.idvar:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.variant_id(args.idvar, chroms=args.chrom)

        if args.somatic_snv:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            callset = "default" if args.callset is None else args.callset
            app.vcf(args.somatic_snv, chroms=args.chrom, sample=args.vcf_sample, no_counts=args.no_snp_counts,
                    ad_tag=args.ad_tag, gt_tag=args.gt_tag, filter=not args.no_filter, callset=callset)

        if args.rd_from_vcf:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.rd_from_vcf(args.rd_from_vcf, chroms=args.chrom, sample=args.vcf_sample, ad_tag=args.ad_tag,
                            dp_tag=args.dp_tag)

        if args.pileup_bam:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.pileup(args.pileup_bam, chroms=args.chrom, reference_filename=args.reference_filename)

        if args.rd_from_snp:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.rd_from_snp(chroms=args.chrom, use_mask=not args.no_mask, use_id=args.use_id,
                            s_bin_size=args.s_bin_size)

        if args.mask:
            app = Root(args.root[0], create=True, max_cores=args.max_cores)
            app.mask(args.mask, chroms=args.chrom, make_mask_genome_file=args.make_mask_genome_file)

        if args.mask_snps:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.mask_snps()

        if args.mask_snvs:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.mask_snps(callset=args.mask_snvs)

        if args.random_phase:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.random_phase()

        if args.trio_phase:
            app = Trio(args.root)
            app.trio_phase(parents=args.phase_parents)

        if args.stat:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.rd_stat(chroms=args.chrom)

        if args.his:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.calculate_histograms(args.his, chroms=args.chrom)

        if args.his_from_snp:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.calculate_histograms_from_snp_counts(args.his_from_snp, chroms=args.chrom, use_mask=not args.no_mask,
                                                     use_id=args.use_id, callset=args.callset,
                                                     min_count=args.min_count)
        if args.baf:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.calculate_baf(args.baf, chroms=args.chrom, use_mask=not args.no_mask, use_id=args.use_id,
                              use_phase=args.use_phase, reduce_noise=args.reduce_noise, blw=args.baf_likelihood_width,
                              use_hom=args.use_hom, alt_ref_correct=args.alt_corr)
        if args.partition:
            app = Root(args.root[0], max_cores=args.max_cores)
            app.partition(args.partition, chroms=args.chrom, use_gc_corr=not args.no_gc_corr,
                          use_mask=args.use_mask_with_rd)

        if args.call:
            app = Root(args.root[0], max_cores=args.max_cores)
            if args.call[0] == "baf":
                if args.call[1] in ["mosaic", "germline"]:
                    event_type = args.call[1]
                    bins = list(map(binsize_type, args.call[2:]))
                else:
                    event_type = "both"
                    bins = list(map(binsize_type, args.call[1:]))
                if args.use_phase:
                    app.call_baf_phased(bins, chroms=args.chrom, event_type=event_type, print_calls=True,
                            use_gc_corr=not args.no_gc_corr,
                            rd_use_mask=args.use_mask_with_rd, snp_use_mask=not args.no_mask, snp_use_id=args.use_id,
                            mcount=args.min_count, max_copy_number=args.max_copy_number,
                            min_cell_fraction=args.min_cell_fraction, baf_threshold=args.baf_threshold,
                            use_hom=args.use_hom, anim=args.animation)
                else:
                    app.call_baf(bins, chroms=args.chrom, event_type=event_type, print_calls=True,
                            use_gc_corr=not args.no_gc_corr,
                            rd_use_mask=args.use_mask_with_rd, snp_use_mask=not args.no_mask, snp_use_id=args.use_id,
                            mcount=args.min_count, max_copy_number=args.max_copy_number,
                            min_cell_fraction=args.min_cell_fraction, baf_threshold=args.baf_threshold,
                            use_hom=args.use_hom, anim=args.animation)
                #app.call_baf_old([binsize_type(x) for x in args.call[1:]], chroms=args.chrom, use_id=args.use_id,
                #             use_mask=not args.no_mask, mcount=args.min_count, anim=args.animation)
            elif args.call[0] == "mosaic":
                app.call_mosaic(list(map(binsize_type, args.call[1:])), chroms=args.chrom,
                                use_gc_corr=not args.no_gc_corr,
                                use_mask=args.use_mask_with_rd, anim=args.animation)
            elif args.call[0] == "combined":
                if args.call[1] in ["mosaic", "germline"]:
                    event_type = args.call[1]
                    bins = list(map(binsize_type, args.call[2:]))
                else:
                    event_type = "both"
                    bins = list(map(binsize_type, args.call[1:]))
                if args.use_phase:
                    app.call_2d_phased(bins, chroms=args.chrom, event_type=event_type, print_calls=True,
                            use_gc_corr=not args.no_gc_corr,
                            rd_use_mask=args.use_mask_with_rd, snp_use_mask=not args.no_mask, snp_use_id=args.use_id,
                            mcount=args.min_count, max_copy_number=args.max_copy_number,
                            min_cell_fraction=args.min_cell_fraction, baf_threshold=args.baf_threshold,
                            use_hom=args.use_hom, anim=args.animation)
                else:
                    app.call_2d(bins, chroms=args.chrom, event_type=event_type, print_calls=True,
                            use_gc_corr=not args.no_gc_corr,
                            rd_use_mask=args.use_mask_with_rd, snp_use_mask=not args.no_mask, snp_use_id=args.use_id,
                            mcount=args.min_count, max_copy_number=args.max_copy_number,
                            min_cell_fraction=args.min_cell_fraction, baf_threshold=args.baf_threshold,
                            use_hom=args.use_hom, anim=args.animation)
            else:
                app.call(list(map(binsize_type, args.call)), chroms=args.chrom, print_calls=True,
                         use_gc_corr=not args.no_gc_corr, use_mask=args.use_mask_with_rd)


if __name__ == '__main__':
    main()
