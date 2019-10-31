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

## Create and store read depth signal

Make sure that you have indexed SAM, BAM or CRAM file.

Initialize your CNVpytor project by running:

```
> cnvpytor -root file.pytor -rd file.bam
```

File file.pytor will be created and read depth data binned to 100 base pair bins will be stored in _pytor_ file.

