# Manhattan plot example

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/manhattan.png">

Import RD signal for each sample from bam/sam/cram file and calculate histograms with specified bin size. 
It can be any positive intiger divisible by 100. We use 100000 in this example

```
> cnvpytor -root sample.root -rd sample.bam
> cnvpytor -root sample.root -his 100000
```

Enter interactive ploting mode:
```
> cnvpytor -root sample1.root sample2.root sample3.root sample4.root -view 100000

cnvpytor> set style bmh
cnvpytor> set rd_use_mask
cnvpytor> set file_titles normal tubular villous cancer 
cnvpytor> manhattan
cnvpytor> manhattan > image.png
cnvpytor> quit
```

By pressing *tab* two times after typing 'set style ' you will got list of all available styles.

If you do not specify file_titles, file names will be used for subplot title.

Instead 'quit' CTRL+D can be used.




