# Getting Started with CNVpytor

## Installing from GitHub

```
> git clone https://github.com/abyzovlab/CNVpytor.git
> cd CNVpytor
> python setup.py install
```
For single user (without admin privileges) use:
```
> python setup.py install --user
```

## Import read depth signal

Make sure that you have indexed SAM, BAM or CRAM file.

Initialize your CNVpytor project by running:

```
> cnvpytor -root file.pytor -rd file.bam
```

File file.pytor will be created and read depth data binned to 100 base pair bins will be stored 
in _pytor_ file.

CNVpytor will detect reference genome and use internal database for GC content and 1000 genome strict mask.

This works for hg19 and hg38 genomes. For other species or reference genomes you have to 
[specify reference genome](examples/AddReferenceGenome.md).

To check is reference genome detected use:

```
> cnvpytor -root file.pytor -ls

```
CNVpytor will print out details about file.pytor including line that specify which reference genome is
used and are there available GC and mask data.
```
Using reference genome: hg19 [ GC: yes, mask: yes ]

```

## Calling CNV/CNA events

## Import SNP data

### From single variant file

### Using SNP positions from variant file and counts from alignment file

## Calculating RD histograms

## Calculating BAF histograms

## Ploting

## Visualize CNVpytor data inside [JBrowse](https://github.com/GMOD/jbrowse)

## Other
