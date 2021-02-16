![CNVpytor Logo](https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/cnvpytor/imgs/cnvpytor_w_640.png)

# CNVpytor - a python extension of CNVnator

CNVpytor is a Python package and command line tool for CNV/CNA analysis from depth-of-coverage by mapped reads developed in Abyzov Lab, Mayo Clinic.

**NEW: CNVpytor release ver 1.0 is available to install!**

## Learn how to use CNVpytor in 10 minutes

* [Geting started](GettingStarted.md) with command line interface
* [Jupyter notebook](examples/CNVpytor.ipynb): How to use CNVpytor from Python
* [Google Colab](examples/Colab.ipynb): With CEPH trio example dataset
* [Video Tutorial](https://www.youtube.com/watch?v=RJMQtrD0SuE)

## Gallery

| | |
|---|---|
| Manhattan plot ([see example](examples/manhattan.md))| Circular plot ([see example](examples/circular.md))|
|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/manhattan.png" width="512px">|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/circular.png" width="512px">|
| Region plot ([see example](examples/region.md))| Compare regions ([see example](examples/compare.md))|
|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/region.png" width="512px">|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/compare.png" width="512px">|
| Merging calls over multiple samples ([see example](examples/merging_cals.md))|
|<img src="https://raw.githubusercontent.com/abyzovlab/CNVpytor/master/imgs/merging.png" width="512px">|

## Install

### Dependencies

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


## Use as a command line tool

![scheme](cnvpytor/imgs/CNVpytor_v1.1.svg)

_Diagram made using [Draw.io](https://github.com/jgraph/drawio)._ 

### Call CNVs using read depth:
```
> cnvpytor -root file.pytor -rd file.bam
> cnvpytor -root file.pytor -his 1000 10000 100000
> cnvpytor -root file.pytor -partition 1000 10000 100000
> cnvpytor -root file.pytor -call 1000 10000 100000
```

### Importing and using single nucleotide polymorphism data:
```
> cnvpytor -root file.pytor -snp file.vcf -sample sample_name
> cnvpytor -root file.pytor -pileup file.bam # OPTIONAL
> cnvpytor -root file.pytor -baf 10000 100000
```

### Filtering calls form view mode
```
> cnvpytor -root file.pytor -view 100000 
print calls
set Q0_range 0 0.5
set size_range 100000 inf
print calls
set p_range 0 0.00001
set print_filename output.xls
print calls
set print_filename output.vcf
print calls
```
Annotating calls:
```
> cnvpytor -root file.pytor -view 100000 
set Q0_range 0 0.5
set size_range 100000 inf
set print_filename output.tsv
set annotate
print calls
```

### Merging calls form multiple samples
```
> cnvpytor -root file1.pytor file2.pytor ... -view 100000 
print joint_calls
set Q0_range 0 0.5
set size_range 100000 inf
set print_filename output.xls
print joint_calls
```
Plotting all merged calls:
```
> cnvpytor -root file1.pytor file2.pytor ... -view 100000 
set Q0_range 0 0.5
set size_range 100000 inf
set print_filename output.xls
set print
set output_filename prefix.png
print joint_calls
```
Annotate calls:
```
> cnvpytor -root file1.pytor file2.pytor ... -view 100000 
set Q0_range 0 0.5
set size_range 100000 inf
set print_filename output.xls
set annotate
print joint_calls
```


### Plot from command line
```
> cnvpytor -root file.pytor -plot rdstat
> cnvpytor -root file.pytor -plot rd 10000 100000
> cnvpytor -root file.pytor -plot rdstat manhattan 100000 -o prefix.pdf
> cnvpytor -root file.pytor -plot baf 100000 
> cnvpytor -root file.pytor -plot regions 1:10M-20M,2:20M-43M 3:10M-20M 10000
> cnvpytor -root file.pytor -plot circular 100000 -rd_use_mask -o prefix.png
```

### Plot from script
```
> echo "rdstat" | cnvpytor -root file.pytor -view 100000 -o prefix.png

> cnvpytor -root file.pytor -view 100000 <<ENDL
set rd_use_mask
set markersize 1
set grid vertical
set output_filename prefix.png
manhattan
circular
ENDL

> cnvpytor -root file.pytor -view 100000 < script.spytor

```

### Plot from interactive mode

CNVpytor view interactive mode is implemented with **<tab> completion** and internal documentation (**help** command).

To enter interactive mode use '-view bin_size' option:

```
> cnvpytor -root file.pytor -view 10000
cnvpytor> chr1:1M-50M
cnvpytor> rd
cnvpytor> set panels rd likelihood
cnvpytor> show
Parameters
    * baf_colors: ['gray', 'black', 'red', 'green', 'blue']
    * bin_size: 100000
    * chrom: []
    * contrast: 20
    * dpi: 200
    * file_titles: []
    * grid: auto
    * lh_colors: ['yellow']
    * markersize: auto
    * min_segment_size: 0
    * output_filename: 
    * panels: ['rd']
    * plot_file: 0
    * plot_files: [0]
            0: file.pytor
    * rd_call: True
    * rd_call_mosaic: False
    * rd_circular_colors: ['#555555', '#aaaaaa']
    * rd_colors: ['grey', 'black', 'red', 'green', 'blue']
    * rd_manhattan_call: False
    * rd_manhattan_range: [0, 2]
    * rd_partition: True
    * rd_range: [0, 3]
    * rd_raw: True
    * rd_use_gc_corr: True
    * rd_use_mask: False
    * snp_call: False
    * snp_circular_colors: ['#00ff00', '#0000ff']
    * snp_colors: ['yellow', 'orange', 'cyan', 'blue', 'lime', 'green', 'yellow', 'orange']
    * snp_use_id: False
    * snp_use_mask: True
    * snp_use_phase: False
    * style: None
    * xkcd: False

cnvpytor> help markersize

markersize
    Size of markers used in scatter like plots (e.g. manhattan, snp).

TYPE
    float or str

DEFAULT
    auto

PLOTS AFFECTS
    manhattan, snp, region plot with snp panel

EXAMPLE(s)
    set markersize 10
    set markersize auto

SEE ALSO
    rd_colors, snp_colors, baf_colors, lh_colors

cnvpytor> set bin_size 100000
cnvpytor> chr1:1M-50M chr2:60M-65M > filename.png
```

## Use as a Python package

CNVpytor is not just command line tool but also Python package. 

For more details check [API Documentation](https://abyzovlab.github.io/CNVpytor/) or 
see examples in [Jupyter notebook](examples/CNVpytor.ipynb).

## JBrowse plugin for CNVpytor
A JBrowse plugin is developed that does on-fly analysis of read depth from VCF file.
https://github.com/abyzovlab/CNVpytorVCF

## Export
### 1. CNVpytor data visualization using JBrowse
#### JBrowse version and plugins
JBrowse version: https://github.com/GMOD/jbrowse/archive/1.16.6-release.tar.gz

 Plugins: 
 - multibigwig (https://github.com/elsiklab/multibigwig )
 - multiscalebigwig (https://github.com/cmdcolin/multiscalebigwig)

**Note:** The JBrowse development version is required as integration of different jbrowse plugins are needed.

#### Usage
To generate cnvpytor file for JBrowse visualization:
```
cnvpytor -root [pytor files] -export jbrowse [optional argument: output path]

Default export directory name: 
 - For single pytor file input:  jbrowse_[pytor file name]
 - For multiple pytor file input: cnvpytor_jbrowse_export
```
The above command creates all the necessary files that are required to visualize the cnvpytor data. 

To view cnvpytor file using JBrowse, users need to install JBrowse and required plugins (See JBrowse version and plugins section).
`
http://localhost/jbrowse/?data=[export directory] 
`

``` 
# Example usage
cnvpytor -root test.pytor -export jbrowse
http://localhost/jbrowse/?data=jbrowse_test
```
[See More](GettingStarted.md#visualize-cnvpytor-data-inside-jbrowse)
## Bugs

Please report any bugs that you find on GitHub:
https://github.com/abyzovlab/CNVpytor/issues

Or, even better, fork the repository on GitHub and create a pull request.

## License

Released under MIT licence.



