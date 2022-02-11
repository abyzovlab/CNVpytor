# Configuring reference genome

For GC correction and 1000 genome strict mask filtering CNVpytor uses information 
related to the reference genome. With installation two reference genomes are
available: hg19 (GRCh37) and hg38 (GRCh38).

If we want to use CNVpytor for other species or different reference genome for human first we have 
to create GC and mask file (optional).

In this example we will configure mouse reference genome MGSCv37.

To create GC file we need sequence of the reference genome in fasta.gz file:

```
> cnvpytor -root MGSCv37_gc_file.pytor -gc ~/hg19/mouse.fasta.gz -make_gc_file
```

This command will produce _MGSCv37_gc_file.pytor_ file that contains information about 
GC content in 100-base-pair bins.

For reference genomes where we have strict mask in the same format as 100 Genomes Project 
[strict mask](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/),
we can create mask file using command:

```
> cnvpytor -root MGSCv37_mask_file.pytor -mask ~/hg19/mouse.strict_mask.whole_genome.fasta.gz -make_mask_file
```

If we do not have mask file, we can skip this step. Mask file contains information about 
regions of the genome that are more accessible to next generation sequencing methods 
using short reads. CNVpytor uses P marked positions to filter SNP-s and read depth signal.
If reference genome configuration does not contain mask file, CNVpytor will still be fully functional, 
apart from the filtering step.
You may also generate your own mask file by creating fasta file that contains character "P" if corresponding 
base pair passes the filter and any character different than "P" if not.

Now, we will create _example_ref_genome_conf.py_ file containing list of chromosomes and chromosome lengths:

```
import_reference_genomes = {
    "mm9": {
        "name": "MGSCv37",
        "species": "Mus musculus",
        "chromosomes": OrderedDict(
            [("chr1", (197195432, "A")), ("chr2", (181748087, "A")), ("chr3", (159599783, "A")),
            ("chr4", (155630120, "A")), ("chr5", (152537259, "A")), ("chr6", (149517037, "A")),
            ("chr7", (152524553, "A")), ("chr8", (131738871, "A")), ("chr9", (124076172, "A")),
            ("chr10", (129993255, "A")), ("chr11", (121843856, "A")), ("chr12", (121257530, "A")),
            ("chr13", (120284312, "A")), ("chr14", (125194864, "A")), ("chr15", (103494974, "A")),
            ("chr16", (98319150, "A")), ("chr17", (95272651, "A")), ("chr18", (90772031, "A")),
            ("chr19", (61342430, "A")), ("chrX", (166650296, "S")), ("chrY", (15902555, "S")),
            ("chrM", (16299, "M"))]),
        "gc_file":"/..PATH../MGSCv37_gc_file.pytor",
        "mask_file": "/..PATH../MGSCv37_mask_file.pytor"
    }
}
```

Last line can be skipped, if there is no mask file. Letter next to chromosome length denote type of a chromosome:
A - autosome, S - sex chromosome, M - mitochondria.


To use CNVpytor with new reference genome us -conf option in each cnvpytor command, e.g.
```
cnvpytor -conf REL_PATH/example_ref_genome_conf.py -root file.pytor -rd file.bam
```

CNVpytor will use chromosome lengths from alignment file to detect reference genome. 
However, if you configured reference genome after you had already run -rd step you 
could assign reference genome using -rg:
```
cnvpytor -conf REL_PATH/example_ref_genome_conf.py -root file.pytor -rg mm9
```

To avoid typing ```-conf REL_PATH/example_ref_genome_conf.py"``` each time you run cnvpytor, 
you can create an bash alias or make configuration permanent by copying ```example_ref_genome_conf.py``` 
to ```~/.cnvpytor/reference_genomes_conf.py```. 

We would like to encourage you to send us configuration, gc and mask file and we would be glad to include it 
into the CNVpytor code. Or, even better, fork the repository on GitHub, add configuration in cnvpytor/genome.py, 
data files in cnvpytor/data and create a pull request.


