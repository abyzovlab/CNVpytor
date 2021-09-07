# Compare plot example

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/compare.png">

Import RD signal from bam/sam/cram file and calculate histograms with specified bin size. 
It can be any positive integer divisible by 100. We use 100000 in this example

```
> cnvpytor -root sample.pytor -rd sample.bam
> cnvpytor -root sample.pytor -his 100000
```

Enter interactive plotting mode with all sample you want to plot listed:

```
> cnvpytor -root sample1.pytor sample2.pytor sample3.pytor sample4.pytor -view 100000

cnvpytor> set style classic
cnvpytor> set rd_use_mask
cnvpytor> set plot
cnvpytor> compare 16:5400k-6200k 16:17500k-18200k

sample.pytor 16:5400k-6200k 16:17500k-18200k 1751.9391 72.2750 2865.2145 96.6045 0.000000e+00 0.6115 0.0053

cnvpytor> save image.png
cnvpytor> quit
```

Output after 'compare' command contains: 
* filename (sample.pytor)
* region1 (16:5400k-6200k)
* region2 (16:17500k-18200k)
* mean read depth signal in region1 (1751.9391)
* stdev of the read depth signal in region1 (72.2750)
* mean read depth signal in region2 (2865.2145)
* stdev of the read depth signal in region2 (96.6045)
* significance: T test p-value (0.0)
* ratio between two read depth levels (0.6115)
* ratio error (0.0053)

By pressing *tab* two times after typing 'set style ' you will got list of all available styles.

If you skip third line (set file_titles) file names will be used for subplot title.

Instead 'quit' CTRL+D can be used.

## Reference

In this example we used data from following article:

[1] Oncotarget 2017 Dec 26;9(6):6780-6792. doi: [10.18632/oncotarget.23687](https://www.doi.org/10.18632/oncotarget.23687)<br>
**Inferring modes of evolution from colorectal cancer with residual polyp of origin.**<br>
Kim M, Druliner BR, Vasmatzis N, Bae T, Chia N, Abyzov A, Boardman LA.
