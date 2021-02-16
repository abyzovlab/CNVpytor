# Annotating calls and merging over multiple samples

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/merging.png">

Import RD signal for each sample from bam/sam/cram file and calculate histograms with specified bin size. 
It can be any positive integer divisible by 100. We use 10000 in this example

```
> cnvpytor -root sample.pytor -rd sample.bam
> cnvpytor -root sample.pytor -his 10000
```

Enter interactive plotting mode with all sample you want to plot listed:
```
> cnvpytor -root sample1.pytor sample2.pytor sample3.pytor sample4.pytor -view 10000

cnvpytor> set Q0_range -1 0.5              # filter calls with more than half not uniqly maped reads
cnvpytor> set p_range 0 0.0001             # filter non-confident calls 
cnvpytor> set size_range 50000 inf         # filter calls smaller than 50kbp
cnvpytor> set dG_range 100000 inf          # filter calls close to gaps in reference genome (<100kbp)
cnvpytor> print joint_calls                # printing calls on screen
...
...
cnvpytor> set print_filename joint.xlsx    # Output filename
cnvpytor> set output_filename joint.png    # Prefix for graphical output files
cnvpytor> set annotate                     # Turn on annotation (optional - takes a lot of time)
cnvpytor> set plot                         # Turn on ploting for each calls (optional - takes a lot of time)
cnvpytor> print joint_calls                # Generate Excel output and png files with RD plots
cnvpytor> quit
```
File joint.xlsx contains list Excel file with list of all calls merged over samples.

Files joint.regions.0000.png to joint.regions.0047.png contain RD region plots for all 48 calls.
This could be used for manual filtering of false positive calls. This is an example of generated plot:

<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/joint.regions.0017.png">


## Reference

In this example we used data from [Illumina Platinum Genomes](https://github.com/Illumina/PlatinumGenomes).

