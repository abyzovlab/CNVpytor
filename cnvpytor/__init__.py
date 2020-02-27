"""
CNVpytor
========
CNVpytor is a Python package and command line tool for CNV analysis
from depth-of-coverage by mapped reads developed in Abyzov Lab, Mayo Clinic.

Source and tutorial::
    https://github.com/abyzovlab/CNVpytor

Bug reports::
    https://github.com/abyzovlab/CNVpytor/issues

API Documentation::
    https://abyzovlab.github.io/CNVpytor/

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
