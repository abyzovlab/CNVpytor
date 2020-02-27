![CNVpytor Logo](https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/cnvpytor/imgs/cnvpytor_w_640.png)

# CNVpytor - a python extension of CNVnator

CNVpytor is a Python package and command line tool for CNV/CNA analysis from depth-of-coverage by mapped reads developed in Abyzov Lab, Mayo Clinic.

**NEW: CNVpytor release ver 1.0 is available to install!**

## Learn how to use CNVpytor in 10 minutes

* [Geting started](https://github.com/abyzovlab/CNVpytor/blob/master/GettingStarted.md) with command line interface
* [Jupyter notebook](https://github.com/abyzovlab/CNVpytor/blob/master/examples/CNVpytor.ipynb): How to use CNVpytor from Python 

## Gallery

| | |
|---|---|
| Manhattan plot ([see example](https://github.com/abyzovlab/CNVpytor/blob/master/examples/manhattan.md))| Circular plot ([see example](https://github.com/abyzovlab/CNVpytor/blob/master/examples/circular.md))|
|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/manhattan.png" width="512px">|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/circular.png" width="512px">|
| Region plot ([see example](https://github.com/abyzovlab/CNVpytor/blob/master/examples/region.md))| Compare regions ([see example](https://github.com/abyzovlab/CNVpytor/blob/master/examples/compare.md))|
|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/region.png" width="512px">|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/compare.png" width="512px">|

## Dependencies

* requests>=2.0
* gnureadline
* pathlib>=1.0
* pysam>=0.15
* numpy>=1.16
* scipy>=1.1
* matplotlib>=2.2
* h5py>=2.9

Optional:

* pyBigWig - for JBrowse export functionality
* ROOT - for CNVnator root import/export functionality
* seaborn - for additional plotting styles 

## Install

```
> pip install cnvpytor
> cnvpytor -download
```

## Using

Check [CNVpytor github page](https://github.com/abyzovlab/CNVpytor).

## License

Released under MIT licence.
