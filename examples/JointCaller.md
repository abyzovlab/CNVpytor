# Using joint caller (prototype)

Import RD signal from bam/sam/cram file, calculate histograms, import SNPs form VCF file and calculate BAF histograms
with specified bin size. Bin size should be larger then 5kbp to have enough SNPs per bin. We use 100kbp in this example:

```
> cnvpytor -root sample.pytor -rd sample.bam
> cnvpytor -root sample.pytor -his 100000
> cnvpytor -root sample.pytor -snp sample.vcf.gz
> cnvpytor -root sample.pytor -mask_snps
> cnvpytor -root sample.pytor -baf 100000
```

Run joint caller:

```
> cnvpytor -root sample.pytor -call combined > calls.tsv
```

Description of columns you can find in [Geting started](GettingStarted.md)

## Filtering and exporting calls

Enter interactive mode with one or multiple samples
```
> cnvpytor -root sample.pytor -view 100000

cnvpytor> set callers combined_mosaic      # set caller to be used for ploting/exporting
cnvpytor> set Q0_range -1 0.5              # filter calls with more than half not uniqly maped reads
cnvpytor> set p_range 0 0.0001             # filter non-confident calls 
cnvpytor> set size_range 500000 inf        # filter calls smaller than 500kbp
cnvpytor> set dG_range 100000 inf          # filter calls close to gaps in reference genome (<100kbp)
cnvpytor> print calls                      # printing calls on screen
...
...
cnvpytor> set print_filename calls.xlsx    # Output filename (tsv, vcf or xlsx)
cnvpytor> print calls                      # Generate Excel output
cnvpytor> quit
```
File calls.xlsx contains list of filtered calls.



## Annotating and plotting all calls

Enter interactive mode with one or multiple samples
```
> cnvpytor -root sample.pytor -view 100000

cnvpytor> set callers combined_mosaic      # set caller to be used for ploting/exporting
cnvpytor> set Q0_range -1 0.5              # filter calls with more than half not uniqly maped reads
cnvpytor> set p_range 0 0.0001             # filter non-confident calls 
cnvpytor> set size_range 500000 inf        # filter calls smaller than 500kbp
cnvpytor> set dG_range 100000 inf          # filter calls close to gaps in reference genome (<100kbp)
cnvpytor> set print_filename calls.vcf     # Output filename
cnvpytor> set output_filename calls.png    # Prefix for graphical output files
cnvpytor> set annotate                     # Turn on annotation (optional - takes a lot of time)
cnvpytor> set plot                         # Turn on ploting for each calls (optional - takes a lot of time)
cnvpytor> set panels rd likelihood         # Set plotting panels
cnvpytor> print calls                      # Generate Excel output and png files with RD/likelihood plots
cnvpytor> quit
```
File calls.vcf contains filtered calls. Files calls.regions.0000.png, calls.regions.0001.png, contain RD/likelihood 
region plots for all calls.

## Plotting with calls

Enter interactive mode with one or multiple samples:
```
> cnvpytor -root sample.pytor -view 100000
cnvpytor> set callers combined_mosaic      # set joint caller
cnvpytor> set rd_range 0 4                 # set rd range
cnvpytor> unset title                      # remove title (sample:region) from plot
cnvpytor> chr16:4m-10m                     # generate plot
```

This is an example of generated plot:

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/joint_caller_plot_example.png">




