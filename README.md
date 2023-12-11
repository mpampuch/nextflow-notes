# Nextflow Notes

Some tips / things of notes for myself while I'm learning Nextflow

## Process outputs

When you run a program, theres a very high likelihood that many output or intermediate files will be created. what the `output:` syntax specifies is the only file or files (or stdout) that your want to include in your output *channel* for the next process or processes.

## Value Channels vs Queue Channels

Value channels allow processes to consume elements infinite times.

Elements in the queue channel are consumed. You cannot consume the same element for the same process twice.

## Single Value Coersion to Value Channels

Processes can only take channels as inputs. That being said, if you pass in a regular variable that has a single value (e.g. string, number, etc), nextflow will implicitly create a value channel containing that value that you are providing. Therefore you can pass a regular value (i.e. string, number) into a process, but just be aware that this is what is going on behind the scenes

Example

```groovy
/*
 * Define the default parameters 
 */

params.genome     = "${baseDir}/data/genome.fa" 
params.variants   = "${baseDir}/data/known_variants.vcf.gz"
params.blacklist  = "${baseDir}/data/blacklist.bed"
params.reads      = "${baseDir}/data/reads/ENCSR000COQ1_{1,2}.fastq.gz" 
params.results    = "results"

/*
 * Process 1A: Create a FASTA genome index with samtools
 */

process prepare_genome_samtools {
    container 'quay.io/biocontainers/samtools:1.3.1--h0cf4675_11'

    input:
    path genome 

    output:
    path "${genome}.fai" 

    script:
    """
    samtools faidx ${genome}
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)

    prepare_genome_samtools(params.genome) 
}
```

The `prepare_genome_samtools(params.genome)` is a  valid call to a process because `params.genome` will be converted from a string into a value channel.

## Biocontainers

Each program should have its own designated container. Don't create container images with too many things or things your don't need.

Run only one process per container: In almost all cases, you should only run a single process in a single container. Decoupling applications into multiple containers makes it much easier to scale horizontally and reuse containers. If that service depends on another service, make use of container linking.

Biocontainers is a project to create a docker container for every recipe they have in bioconda.

Sometimes you'll need to have a container with more than one tool, in this case there is a multi-package container. 

You can request a multi-package container here: https://biocontainers.pro/multipackage

## get*() in Groovy

In Groovy, any method that looks like `get*()` can also be accessed as a field. For example, `myFile.getName()` is equivalent to `myFile.name`, `myFile.getBaseName()` is equivalent to `myFile.baseName`, and so on.

## Important functions for file attributes and paths

The `file()` method returns a Path, so any method defined for Path can also be used in a Nextflow script.

Additionally, the following methods are also defined for Paths in Nextflow:

`exists()`
Returns true if the file exists.

`getBaseName()`
Gets the file name without its extension, e.g. `/some/path/file.tar.gz` -> `file.tar`.

`getExtension()`
Gets the file extension, e.g. `/some/path/file.txt` -> `txt`.

`getName()`
Gets the file name, e.g. `/some/path/file.txt` -> `file.txt`.

`getSimpleName()`
Gets the file name without any extension, e.g. `/some/path/file.tar.gz` -> `file`.

`getParent()`
Gets the file parent path, e.g. `/some/path/file.txt` -> `/some/path`.

`isDirectory()`
Returns true if the file is a directory.

`isEmpty()`
Returns true if the file is empty or does not exist.

`isFile()`
Returns true if it is a regular file (i.e. not a directory).

`isLink()`
Returns true if the file is a symbolic link.
