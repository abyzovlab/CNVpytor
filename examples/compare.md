# Manhattan plot example

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/compare.png">

Import RD signal from bam/sam/cram file and calculate histograms with specified bin size. 
It can be any positive integer divisible by 100. We use 100000 in this example

```
> cnvpytor -root sample.root -rd sample.bam
> cnvpytor -root sample.root -his 100000
```

Enter interactive plotting mode with all sample you want to plot listed:

```
> cnvpytor -root sample1.root sample2.root sample3.root sample4.root -view 100000

cnvpytor> set style classic
cnvpytor> set rd_use_mask
cnvpytor> compare 16:5400k-6200k 16:17500k-18200k
cnvpytor> save image.png
cnvpytor> quit
```

By pressing *tab* two times after typing 'set style ' you will got list of all available styles.

If you skip third line (set file_titles) file names will be used for subplot title.

Instead 'quit' CTRL+D can be used.
