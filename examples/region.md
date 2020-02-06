# Circular plot example

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/region.png">

Import RD signal for each sample from bam/sam/cram file and calculate histograms with specified bin size. 
It can be any positive integer divisible by 100. We use 100000 in this example

```
> cnvpytor -root sample.root -rd sample.bam
> cnvpytor -root sample.root -his 100000
```

Import BAF signal for each sample from vcf file and calculate histograms with specified bin size. 

```
> cnvpytor -root sample.root -snp sample.vcf.gz
> cnvpytor -root sample.root -mask_snps
> cnvpytor -root sample.root -baf 100000
```
Second line is used to filter out all SNP-s that are not in P-region of the strict 1kG mask.

Enter interactive plotting mode with all sample you want to plot listed:
```
> cnvpytor -root sample1.root sample2.root sample3.root sample4.root -view 100000

cnvpytor> set style classic
cnvpytor> set rd_use_mask
cnvpytor> set file_titles normal tubular villous cancer 
cnvpytor> 15:20M-100M
cnvpytor> save image.png
cnvpytor> quit
```

By pressing *tab* two times after typing 'set style ' you will got list of all available styles.

If you skip third line (set file_titles) file names will be used for subplot title.

Instead 15:20M-100M you can specify any genomic region or regions. 
If there are more then one region they can be comma or space separated 
(comma - to be plotted in the same subplot, space - to be separately plotted).

Instead 'quit' CTRL+D can be used.