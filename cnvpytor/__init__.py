"""
CNVpytor
========
CNVpytor is a Python package and command line tool for CNV analysis
from depth-of-coverage by mapped reads developed in Abyzov Lab, Mayo Clinic.

Source::
    https://github.com/abyzovlab/CNVpytor

Bug reports::
    https://github.com/abyzovlab/CNVpytor/issues

API Documentation::
    https://abyzovlab.github.io/CNVpytor/

## Dependencies

* gnureadline
* requests
* pysam
* numpy
* scipy
* matplotlib
* h5py >= 2.9

Optional:

* ROOT - for CNVnator root import/export functionality
* seaborn - for additional plotting styles

## Install

### Install by cloning from GitHub

```
> git clone https://github.com/abyzovlab/CNVpytor.git
> cd CNVpytor
> pip install .
```
For single user (without admin privileges) use:
```
> pip install --user .
```

### Install using pip

```
> pip install cnvpytor
> cnvpytor -download
```

## Simple example

Call CNV using read depth:
```
> cnvpytor -root file.pytor -rd file.bam
> cnvpytor -root file.pytor -his 1000 10000 100000
> cnvpytor -root file.pytor -partition 1000 10000 100000
> cnvpytor -root file.pytor -call 1000 10000 100000
```

Call CNV using single nucleotide polymorphism::
```
> cnvpytor -root file.pytor -snp file.vcf
> cnvpytor -root file.pytor -pileup file.bam
> cnvpytor -root file.pytor -baf 10000 100000
> cnvpytor -root file.pytor -call baf 10000 100000
```

Plot
```
> cnvpytor -root file.pytor -plot stat
> cnvpytor -root file.pytor -plot 10000 100000
> cnvpytor -root file.pytor -plot stat manhattan 100000 -o prefix.pdf
> cnvpytor -root file.pytor -plot baf -chrom 1 2 3 4
> cnvpytor -root file.pytor -plot regions 1:10M-20M,2:20M-43M 3:10M-20M 10000
> cnvpytor -root file.pytor -plot circular 100000 -use_mask_rd -o prefix.png
```

Plot - interactive mode
```
> cnvpytor -root file.pytor -view 10000
cnvpytor> chr1:1M-50M
cnvpytor> rd
cnvpytor> set panels rd likelihood
cnvpytor> show
    Parameters
        * bin_size: 100000
        * panels: ['rd','likelihood']
        * use_mask_rd: False
        * use_mask: True
        * use_id: False
        * plot_files:
             0 file1.pytor True
             1 file2.pytor True
             2 file3.pytor True
        * plot_file: 0
        * grid: auto

cnvpytor> set bin_size 100000
cnvpytor> chr1:1M-50M chr2:60M-65M > filename.png
```


## Bugs

Please report any bugs that you find on GitHub:
https://github.com/abyzovlab/CNVpytor/issues

Or, even better, fork the repository on GitHub and create a pull request.

## License

Released under MIT licence.

"""
from __future__ import absolute_import, print_function, division
from .io import *
from .bam import *
from .utils import *
from .fasta import *
from .vcf import *
from .root import *
from .viewer import *
from .genome import Genome
from .version import __version__
