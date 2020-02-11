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

### From alignment file

Make sure that you have indexed SAM, BAM or CRAM file.

Initialize your CNVpytor project by running:

```
> cnvpytor -root file.pytor -rd file.bam
```

File file.pytor will be created and read depth data binned to 100 base pair bins will be stored in _pytor_ file.

### From variant file

## Import SNP data

### From single variant file

### Using SNP positions from variant file and counts from alignment file

## Calculating RD histograms

## Calculating BAF histograms

## Calling CNV/CNA events

## Ploting

## Other
