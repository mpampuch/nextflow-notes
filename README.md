# Nextflow Notes

Some tips / things of notes for myself while I'm learning Nextflow

## What is Nextflow?

Nextflow is a framework used for writing parallel and distributed computational pipelines, common in bioinformatics, genomics, and other fields where complex data analysis pipelines are common. It simplifies the process of creating complex workflows that involve processing large volumes of data.

Nextflow pipelines are partitioned into channels, processes, and workflows. Channels represent data streams, processes perform computations on that data, and workflows orchestrate the execution of processes and manage the flow of data through the pipeline.

![Nextflow example](what-is-nextflow.png)

The central tenant of Nextflow is generating pipelines with 100% reproducibility. Most bioinformatic studies use multiple softwares and languages all with different versions. This makes reproducibility very difficult as you have to go back and figure out what softwares, versions and what dependancies were used. Worst off different softwares might not work on certain machines or even can work differently depending on what machine you are running the analysis on, resulting in different statistical calculations and outputs from the same input data (see https://doi.org/10.1038/nbt.3820 for more information). With Nextflow, if you use it correctly and follow all the best practices, your data analysis pipelines should be fully reproducible regardless of when or where it is run.

## Running Nextflow

Nextflow pipelines need to be written into files with the `.nf` extension.

Example, create a file called `main.nf` and provide it the following code:

```nextflow
#!/usr/bin/env nextflow

process FASTQC {
    input:
    path input

    output:
    path "*_fastqc.{zip, html}"

    script:
    """
    fastqc -q $input
    """
}

workflow {
    Channel.fromPath("*.fastq.gz") | FASTQC
}
```

**Note**: The shebang (`#!/usr/bin/env nextflow`) is a line that helps the operating system decide what program should be used to interpret this script code. If you always use Nextflow to call this script, *this line is optional*.


Nextflow can then run by using the following command

```bash
nextflow run main.nf
```

Only run it like this the first time. If you run it again it will begin running the entire workflow over again. It's better to not rerun processes that don't need to be rerun to save a lot of time.

When you run a nextflow pipeline, you will get an output like this in your terminal

```
N E X T F L O W  ~  version 23.10.0
Launching `main.nf` [focused_noether] DSL2 - revision: 197a0e289a
executor >  local (3)
[18/f6351b] process > SPLITLETTERS (1)   [100%] 1 of 1 ✔
[2f/007bc5] process > CONVERTTOUPPER (1) [100%] 2 of 2 ✔
WORLD!
HELLO 
```

This is what this all means

```
N E X T F L O W  ~  

# THE VERSION OF NEXTFLOW THAT WAS USED FOR THIS RUN SPECIFICALLY
version 23.10.0

# THE NAME OF THE SCRIPT FILE THAT WAS RUN WITH NEXTFLOW
Launching `main.nf` 

# A GENERATED MNEMONIC WHICH IS AN RANDOM OBJECTIVE AND A RANDOM SURNAME OF A FAMOUS SCIENTIST
[focused_noether] 

# THE VERSION OF THE NEXTFLOW LANGUAGE
DSL2

# THE REVISION HASH, WHICH IS LIKE AN ID OF YOUR PIPELINE.
 **If you change the workflow script code, this revision ID will change**
- revision: 197a0e289a

# THE EXECUTOR FOR THE PIPELINE
executor >  local 

# NEXTFLOW'S GUESS AT HOW MANY TASKS THAT WILL OCCUR IN YOUR PIPELINE
(3)

# YOUR PROCESSES AND OUTPUTS
[18/f6351b] process > SPLITLETTERS (1)   [100%] 1 of 1 ✔
[2f/007bc5] process > CONVERTTOUPPER (1) [100%] 2 of 2 ✔
WORLD!
HELLO 
```

### Rerunning Nextflow

When using the `-resume` flag, successfully completed tasks are skipped and the previously cached results are used in downstream tasks. 

```bash
nextflow run main.nf -resume
```

In practice, every execution starts from the beginning. However, when using resume, before launching a task, Nextflow uses the unique ID to check if:

- the working directory exists
- it contains a valid command exit status
- it contains the expected output files.

The mechanism works by assigning a unique ID to each task. This unique ID is used to create a separate execution directory, called the working directory, where the tasks are executed and the results stored. A task’s unique ID is generated as a 128-bit hash number obtained from a composition of the task’s:

- Inputs values
- Input files
- Command line string
- Container ID
- Conda environment
- Environment modules
- Any executed scripts in the bin directory

## Configuring Nextflow

The first thing that nextflow looks for when a workflow is run is configuration files in multiple locations. Since each configuration file can contain conflicting settings, the sources are ranked to determine which settings are applied. Possible configuration sources, in order of priority:

1. Parameters specified on the command line (`--something value`)

2. Parameters provided using the `-params-file` option

3. Config file specified using the `-c my_config` option

4. The config file named `nextflow.config` in the current directory

5. The config file named `nextflow.config` in the workflow project directory

6. The config file `$HOME/.nextflow/config`

7. Values defined within the pipeline script itself (e.g. `main.nf`)

When more than one of these options for specifying configurations are used, they are merged, so that the settings in the first override the same settings appearing in the second, and so on.

This is an example of a nextflow configuration in a file called `nextflow.config`

```nextflow
propertyOne = 'world'
propertyTwo = "Hello $propertyOne"
customPath = "$PATH:/my/app/folder"
```

*Note*: The quotes act like in bash. Variables inside single quotes remain literal. Variables inside double quotes get expanded (including environment variables)

### Config scopes
Configuration settings can be organized in different scopes by dot prefixing the property names with a scope identifier, or grouping the properties in the same scope using the curly brackets notation. For example:

```nextflow
alpha.x = 1
alpha.y = 'string value..'

beta {
     p = 2
     q = 'another string ..'
}
```

### Config params
The scope params allows the definition of workflow parameters that override the values defined in the main workflow script.

This is useful to consolidate one or more execution parameters in a separate file.

`nextflow.config`

```nextflow
params.foo = 'Bonjour'
params.bar = 'le monde!'
```

`snippet.nf`

```nextflow
params.foo = 'Hello'
params.bar = 'world!'

// print both params
println "$params.foo $params.bar"
```

Any variables that are set using `params.` can be modified through the command line when executing Nextflow using two dashes (`--`)

E.g.

```bash
nextflow run snippet.nf --foo "Bye"
```

**Note**: In nextflow, single dashes (`-`) in command line arguments refer to Nextflow commands (e.g. `-resume`), while double dashes (`--`) refer to workflow parameters

### Config env
The env scope allows the definition of one or more variables that will be exported into the environment where the workflow tasks will be executed.

### Config process scopes and process directives

Process directives allow the specification of settings for the task execution such as cpus, memory, container, and other resources in the workflow script.

This is useful when prototyping a small workflow script.

However, it’s always a good practice to decouple the workflow execution logic from the process configuration settings, i.e. ***it’s strongly suggested to define the process settings in the workflow configuration file instead of the workflow script***.

The process configuration scope allows the setting of any process directives in the Nextflow configuration file:

`nextflow.config`


```nextflow
process {
    cpus = 10
    memory = 8.GB
    container = 'biocontainers/bamtools:v2.4.0_cv3'
}
```
The above config snippet defines the cpus, memory and container directives for all processes in your workflow script. Depending on your executor these things may behave differently.

### Process selectors
The `withLabel` selectors allow the configuration of all processes annotated with a label directive as shown below:

```nextflow
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
        queue = 'long'
    }
}
```
The above configuration example assigns 16 cpus, 64 Gb of memory and the long queue to all processes annotated with the `big_mem` label.

Labels can be added using the `label` directive

e.g.

```nextflow
process bigTask {
    label 'big_mem'

    ```
    <task script>
    ```
}
```

In the same manner, the `withName` selector allows the configuration of a specific process in your pipeline by its name. For example:

```nextflow
process {
    withName: hello {
        cpus = 4
        memory = 8.GB
        queue = 'short'
    }
}
```

A process selector can also be negated prefixing it with the special character `!`. For example:

```nextflow
process {
    withLabel: 'foo' { cpus = 2 }
    withLabel: '!foo' { cpus = 4 }
    withName: '!align.*' { queue = 'long' }
}
```

The above configuration snippet sets 2 cpus for the processes annotated with the `foo` label and 4 cpus to all processes not annotated with that label. Finally it sets the use of `long` queue to all process whose name does not start with `align`.


### Config Conda execution

If you already have a conda environment in your machine that you want to use for your processes, you can use

`nextflow.config`

```
process.conda = "/home/ubuntu/miniconda2/envs/nf-tutorial"
```

You can specify the path of an existing Conda environment as either directory or the path of Conda environment YAML file.

## Processes

In nextflow, it is a best practice to always name processes in UPPERCASE. This way you can easily see what are process blocks and what are regular functions.

## The `fair` process directive

Because the processes in nextflow can finish in different orders than they were inputted as, the `fair` directive is a way that you can ensure the tasks finish in the order that they were started. This is because the `fair` process directive distributes computing resources in a "fair" way (comes from fair-threading) to ensure the first one finishes first and so on.

## Task directories

Take a look at the last part of the output when you run a nextflow pipeline

```
# YOUR PROCESSES AND OUTPUTS
[18/f6351b] process > SPLITLETTERS (1)   [100%] 1 of 1 ✔
[2f/007bc5] process > CONVERTTOUPPER (1) [100%] 2 of 2 ✔
WORLD!
HELLO 
```

Breaking down the processes section further this is what each part means

```
# THE HEXADECIMAL HASH FOR THE TASK (ALSO THE TASK DIRECTORY NAME)
[18/f6351b] 

# THE NAME OF THE PROCESS USED FOR THE TASK
process > SPLITLETTERS 

# I'M NOT SURE WHAT THIS NUMBER MEANS YET. MAYBE THE CURRENT PROCESS?
(1)   

# THE PROGRESS OF ALL THE TASKS (OR MAYBE THE CURRENT TASK?)
[100%] 

# ENUMERATED PROGRESS OF ALL THE TASKS
1 of 1 ✔

# THE HEXADECIMAL HASH FOR THE TASK (LAST TASK ONLY IF THERE ARE MORE THAN ONE)
[2f/007bc5] 

# THE NAME OF THE PROCESS USED FOR THE TASK
process > CONVERTTOUPPER 

# I'M NOT SURE WHAT THIS NUMBER MEANS YET. MAYBE THE CURRENT PROCESS?
(1) 

# THE PROGRESS OF ALL THE TASKS (OR MAYBE THE CURRENT TASK?)
[100%] 

# ENUMERATED PROGRESS OF ALL THE TASKS
2 of 2 ✔
```

The hexadecimal numbers, like `18/f6351b`, identify the unique process execution, that we call a task. These numbers are also the prefix of the directories where each task is executed. You can inspect the files produced by changing to the directory `$PWD/work` and using these numbers to find the task-specific execution path (e.g. Go to  `$PWD/work/18/f6351b46bb9f65521ea61baaaa9eff` to find all the information on the task performed using the `SPLITLETTERS` process).

### The ANSI log

If you look at the second process in the above examples, you notice that it runs twice (once for each task), executing in two different work directories for each input file. The ANSI log output from Nextflow dynamically refreshes as the workflow runs; in the previous example the work directory `2f/007bc5` is the second of the two directories that were processed (overwriting the log with the first). To print all the relevant paths to the screen, disable the ANSI log output using the `-ansi-log` flag.

Example

```bash
nextflow run hello.nf -ansi-log false
```

Will output

```
N E X T F L O W  ~  version 23.10.0
Launching `hello.nf` [boring_bhabha] DSL2 - revision: 197a0e289a
[18/f6351b] Submitted process > SPLITLETTERS (1)
[dc/e177f3] Submitted process > CONVERTTOUPPER (1) # NOTICE THE HASH FOR TASK 1 IS VISIBLE NOW
[2f/007bc5] Submitted process > CONVERTTOUPPER (2)
HELLO 
WORLD!
```

Now you can find out in which directory everything related to every task performed is stored straight from the console.

**Note:** Even if you don't use the `-ansi-log false` flag, you can still see the hashes/directories all the tasks are stored in using the `.nextflow.log` file. The task directories can be found in the `[Task monitor]` logs.

### The `.command.sh` file

Inside every task directory that was successfully run, there will be a `.command.sh` file (e.g. `$PWD/work/18/f6351b46bb9f65521ea61baaaa9eff/.command.sh`). This file contains the **final** script that was run for that task.

Example: If this is in the workflow

```nextflow
params.greeting = 'Hello world!'

process SPLITLETTERS {
    input:
    val x

    output:
    path 'chunk_*'

    """
    printf '$x' | split -b 6 - chunk_
    """
}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
}
```

The `.command.sh` file for a task run on this process will look like this

```bash
printf 'Hello world!' | split -b 6 - chunk_
```

This is very useful for troubleshooting when things don't work like you'd expect.

### The `.command.run` file

The `.command.run` file is very similar to the `.command.sh` file, except that it shows you the jobscript that Nextflow created to run the script (e.g. If you are running your scripts with SLURM, it will show you the SLURM job script Nextflow created and that was subsequently called with `sbatch`).

## Creating a configuration file
When a workflow script is launched, Nextflow looks for a file named `nextflow.config` in the current directory and in the script base directory (if it is not the same as the current directory). Finally, it checks for the file: `$HOME/.nextflow/config`.

When more than one of the above files exists, they are merged, so that the settings in the first override the same settings that may appear in the second, and so on.

The default config file search mechanism can be extended by providing an extra configuration file by using the command line option: `-c <config file>`.

Information on writing these config files can be found here https://training.nextflow.io/basic_training/config/.

### Configuring a SLURM script for each process

If you add the following code to the `nextflow.config` file

```nextflow
process.executor = 'slurm'
```

Then nextflow will write the SLURM job script for every file for you. Nextflow will manage each process as a separate job that is submitted to the cluster using the `sbatch` command.

More information on how to configure this further can be found here https://www.nextflow.io/docs/latest/executor.html#slurm

### Running Nextflow workflows on SLURM / HPC's

When running Nextflow on a HPC, it's recommended to run it as a job on a compute node. This is because a lot of computing clusters have strict rules on running processes on login nodes. Therefore, it's always advisable to create jobscripts like this for all your nextflow jobs. 

`launch_nf.sh`
```bash
#!/bin/bash
#SBATCH --partition WORK
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 12:00:00

WORKFLOW=$1
CONFIG=$2

# Use a conda environment where you have installed Nextflow
# (may not be needed if you have installed it in a different way)
conda activate nextflow

nextflow -C ${CONFIG} run ${WORKFLOW}
```

and launch the workflow using

```bash
sbatch launch_nf.sh /home/my_user/path/my_workflow.nf /home/my_user/path/my_config_file.conf
```

## The `.view` Channel Operator

`.view` is a channel operator that consumes every element of a channel and prints it to the screen.

Example

```nextflow
// main.nf
params.reads = "${baseDir}/data/reads/ENCSR000COQ1_{1,2}.fastq.gz"

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)
    reads_ch.view()
}
```

```bash
nextflow run main.nf
# Outputs
# [ENCSR000COQ1, [path/to/ENCSR000COQ1_1.fastq.gz, path/to/ENCSR000COQ1_2.fastq.gz]]
```

You can see in this case it outputs a single channel element created from the `.fromFilePairs` channel operator.
- The first item is an ID (in this case the replicate ID)
- The second item is a list with all the items (in this case the path to the reads)

## Process Inputs

The `input` block allows you to define the input channels of a process, similar to function arguments. A process may have at most one input block, and it must contain at least one input.

The input block follows the syntax shown below:

```nextflow
input:
  <input qualifier> <input name>
```

An input definition consists of a *qualifier* and a *name*. The input qualifier defines the type of data to be received. This information is used by Nextflow to apply the semantic rules associated with each qualifier, and handle it properly depending on the target execution platform (grid, cloud, etc).

The following input qualifiers are available:

`val`: Access the input value by name in the process script.

`file`: (DEPRECATED) Handle the input value as a file, staging it properly in the execution context.

`path`: Handle the input value as a path, staging the file properly in the execution context.

`env`: Use the input value to set an environment variable in the process script.

`stdin`: Forward the input value to the process stdin special file.

`tuple`: Handle a group of input values having any of the above qualifiers.

`each`: Execute the process for each element in the input collection.

More information on how each of these input qualifiers work can be found here https://www.nextflow.io/docs/latest/process.html#inputs

### The double asterisk

A double asterisk (`**`) in a glob pattern works like `*` but also searches through subdirectories. For example, imagine this is your file structure

```
data
├── 2023-11-08_upn22_monarch-hmw_100mg_rep-2
│   └── no_sample
│       └── 20231108_1310_MC-115499_FAX00407_d135d0ec
│           ├── fastq_fail
│           │   ├── FAX00407_fail_d135d0ec_b0bb43ca_0.fastq.gz
│           │   └── FAX00407_fail_d135d0ec_b0bb43ca_1.fastq.gz
│           ├── fastq_pass
│           │   ├── FAX00407_pass_d135d0ec_b0bb43ca_0.fastq.gz
│           │   └── FAX00407_pass_d135d0ec_b0bb43ca_1.fastq.gz
│           ├── other_reports
│           │   ├── pore_scan_data_FAX00407_d135d0ec_b0bb43ca.csv
│           │   └── temperature_adjust_data_FAX00407_d135d0ec_b0bb43ca.csv
│           ├── pod5_fail
│           │   ├── FAX00407_fail_d135d0ec_b0bb43ca_0.pod5
│           │   └── FAX00407_fail_d135d0ec_b0bb43ca_1.pod5
│           └── pod5_pass
│               ├── FAX00407_pass_d135d0ec_b0bb43ca_0.pod5
│               └── FAX00407_pass_d135d0ec_b0bb43ca_1.pod5
├── 2023-11-12_mhmw-dna-recovery-elude-1
│   └── no_sample
│       └── 20231112_1338_MC-115499_FAX00228_b67d08a5
│           ├── fastq_fail
│           │   ├── FAX00228_fail_b67d08a5_dc19481f_0.fastq.gz
│           │   └── FAX00228_fail_b67d08a5_dc19481f_1.fastq.gz
│           ├── fastq_pass
│           │   ├── FAX00228_pass_b67d08a5_dc19481f_0.fastq.gz
│           │   └── FAX00228_pass_b67d08a5_dc19481f_1.fastq.gz
│           ├── final_summary_FAX00228_b67d08a5_dc19481f.txt
│           ├── other_reports
│           │   ├── pore_scan_data_FAX00228_b67d08a5_dc19481f.csv
│           │   └── temperature_adjust_data_FAX00228_b67d08a5_dc19481f.csv
│           ├── pod5_fail
│           │   ├── FAX00228_fail_b67d08a5_dc19481f_0.pod5
│           │   └── FAX00228_fail_b67d08a5_dc19481f_1.pod5
│           ├── pod5_pass
│           │   ├── FAX00228_pass_b67d08a5_dc19481f_0.pod5
│           │   └── FAX00228_pass_b67d08a5_dc19481f_1.pod5
│           └── sequencing_summary_FAX00228_b67d08a5_dc19481f.txt
├── 2023-11-16
│   └── no_sample
│       └── 20231116_0945_MC-115499_FAX00393_849b7392
│           ├── barcode_alignment_FAX00393_849b7392_a554d814.tsv
│           ├── fastq_fail
│           │   ├── FAX00393_fail_849b7392_a554d814_0.fastq.gz
│           │   └── FAX00393_fail_849b7392_a554d814_1.fastq.gz
│           ├── fastq_pass
│           │   ├── FAX00393_pass_849b7392_a554d814_0.fastq.gz
│           │   └── FAX00393_pass_849b7392_a554d814_1.fastq.gz
│           ├── final_summary_FAX00393_849b7392_a554d814.txt
│           ├── other_reports
│           │   ├── pore_scan_data_FAX00393_849b7392_a554d814.csv
│           │   └── temperature_adjust_data_FAX00393_849b7392_a554d814.csv
│           ├── pod5_fail
│           │   ├── FAX00393_fail_849b7392_a554d814_0.pod5
│           │   └── FAX00393_fail_849b7392_a554d814_1.pod5
│           ├── pod5_pass
│           │   ├── FAX00393_pass_849b7392_a554d814_0.pod5
│           │   └── FAX00393_pass_849b7392_a554d814_1.pod5
│           ├── pore_activity_FAX00393_849b7392_a554d814.csv
│           ├── report_FAX00393_20231116_0945_849b7392.html
│           ├── report_FAX00393_20231116_0945_849b7392.json
│           ├── report_FAX00393_20231116_0945_849b7392.md
│           ├── sample_sheet_FAX00393_20231116_0945_849b7392.csv
│           ├── sequencing_summary_FAX00393_849b7392_a554d814.txt
│           └── throughput_FAX00393_849b7392_a554d814.csv
└── Oct-26-23-upn22-pci-rep-3
    └── no_sample
        └── 20231026_1515_MC-115499_FAW96674_9d505d15
            ├── barcode_alignment__9d505d15_1f674c3a.tsv
            ├── fastq_fail
            │   ├── FAW96674_fail_9d505d15_1f674c3a_0.fastq.gz
            │   └── FAW96674_fail_9d505d15_1f674c3a_1.fastq.gz
            ├── fastq_pass
            │   ├── FAW96674_pass_9d505d15_1f674c3a_0.fastq.gz
            │   └── FAW96674_pass_9d505d15_1f674c3a_1.fastq.gz
            ├── final_summary_FAW96674_9d505d15_1f674c3a.txt
            ├── other_reports
            │   ├── pore_scan_data_FAW96674_9d505d15_1f674c3a.csv
            │   └── temperature_adjust_data_FAW96674_9d505d15_1f674c3a.csv
            ├── pod5_fail
            │   ├── FAW96674_fail_9d505d15_1f674c3a_0.pod5
            │   └── FAW96674_fail_9d505d15_1f674c3a_1.pod5
            ├── pod5_pass
            │   ├── FAW96674_pass_9d505d15_1f674c3a_0.pod5
            │   └── FAW96674_pass_9d505d15_1f674c3a_1.pod5
            ├── pore_activity__9d505d15_1f674c3a.csv
            ├── report__20231026_1515_9d505d15.html
            ├── report__20231026_1515_9d505d15.json
            ├── report__20231026_1515_9d505d15.md
            ├── sample_sheet__20231026_1515_9d505d15.csv
            ├── sequencing_summary_FAW96674_9d505d15_1f674c3a.txt
            └── throughput__9d505d15_1f674c3a.csv

33 directories, 60 files
```

All the `.fastq.gz` files can be grabbed using the `**` wildcard to search through all the subdirectories to look for the files with the `.fastq` extensions.

```
fastq_ch = Channel.fromPath("data/**/*.fastq.gz").collect()
```

## Process outputs

When you run a program, theres a very high likelihood that many output or intermediate files will be created. what the `output:` syntax specifies is the only file or files (or stdout) that your want to include in your output *channel* for the next process or processes.

### The `.out` attribute

By using `.out`, your are getting the output channel of one process, and you can pass it in as the input channel of another process

Example

```nextflow
workflow {
    reads_ch = Channel.fromFilePairs(params.reads)

    prepare_star_genome_index(params.genome)

    rnaseq_mapping_star(params.genome, 
                        prepare_star_genome_index.out,
                        reads_ch)
}
```
When a process defines multiple output channels, each output can be accessed using the array element operator (`out[0]`, `out[1]`, etc.)

Outputs can also be accessed by name if the `emit` option is specified

```nextflow
process example_process {
    output:
    path '*.bam', emit: samples_bam

    '''
    your_command --here
    '''
}

workflow {
    example_process()
    example_process.out.samples_bam.view()
}
```

### Outputs with the same filenames

If you have a process that only has a static pathname, for example

```nextflow
process rnaseq_call_variants {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${sampleId}" 

    input:
    path genome 
    path index 
    path dict  
    tuple val(sampleId), path(bam), path(bai) 

    output:
    tuple val(sampleId), path('final.vcf') 

    script:
    """
    ...
    """
```

This won't create a name conflict for every sample that gets processed

This is because nextflow creates a new directory for every task a process performs. So if you're trying to process 10 samples (run 10 tasks from a single process), you're going to have 10 isolated folders.

If you only had the `path` variable defined and not the `tuple` with the `sampleId` then it may have caused an issue but the way it's defined here file conflicts won't be an issue because every sample will get it's own folder.

## Processing things in order / the `fair` directive

While channels do emit items in the order that they are received (FIFO structure), ***processes do not necessarily process items in the order that they are received*** (because of implicit parallelization and later processes ending before earlier ones). While this isn't an issue in most cases, it is important to know. 

For example

```nextflow
process basicExample {
  input:
  val x

  "echo process job $x"
}

workflow {
  def num = Channel.of(1,2,3)
  basicExample(num)
}
```

Can output 

```
process job 3
process job 1
process job 2
```


Notice in the above example that the value `3` was processed before the others.

The `fair` directive (new in version 22.12.0-edge), when enabled, guarantees that process outputs will be emitted in the order in which they were received. For example:

```nextflow
process EXAMPLE {
    fair true

    input:
    val x
    output:
    tuple val(task.index), val(x)

    script:
    """
    sleep \$((RANDOM % 3))
    """
}

workflow {
    channel.of('A','B','C','D') | EXAMPLE | view
}
```

The above example produces:

```
[1, A]
[2, B]
[3, C]
[4, D]
```

## Process scripts

The script block within a process is a string statement that defines the command to be executed by the process' task.

```nextflow
process CONVERTTOUPPER { 
    input: 
    path y 

    output: 
    stdout 

    script: 
    """
    cat $y | tr '[a-z]' '[A-Z]'
    """
}
```

By default, nextflow expects a shell script in the script block. If you are using another language, like R or Python, you need the shebang so that Nextflow knows which software to use to interpret this code.

Example:

```nextflow
process CONVERTTOUPPER { 
    input: 
    path y 

    output: 
    stdout 

    script: 
    """
    #!/usr/bin/env python
    with open("$y") as f:
        print(f.read().upper(), end="")
    """
}
```


## Channel Operators / performing operations on channels outside of a process

Sometimes the output channel of one process doesn't quite match the input channel of the next process and so it has to be modified slightly. This can be performed using channel operators. A full list of channel operators can be found here https://www.nextflow.io/docs/latest/operator.html.

For example, in this code:

```nextflow
rnaseq_gatk_analysis
    .out
    .groupTuple()
    .join(prepare_vcf_for_ase.out.vcf_for_ASE)
    .map { meta, bams, bais, vcf -> [meta, vcf, bams, bais] }
    .set { grouped_vcf_bam_bai_ch }
```
- `.groupTuple` groups tuples that contain a common first element
- `.join` joins two channels taking the key into consideration
- `.map` applies a function to every element of a channel
- `.set` saves this channel with a new name

Step by step this looks like:

```nextflow
rnaseq_gatk_analysis
    .out
    /* Outputs
    [ENCSR000COQ, /workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam, /workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam.bai]
    [ENCSR000COQ, /workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ2.final.uniq.bam, /workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ2.final.uniq.bam.bai]
    */
    .groupTuple()
    /* Outputs
    [ENCSR000COQ, [/workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam, /workspace/gitpod/hands-on/work/92/b1ea340ce922d13bdce2985b2930f2/ENCSR000COQ2.final.uniq.bam], [/workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam.bai, /workspace/gitpod/hands-on/work/92/b1ea340ce922d13bdce2985b2930f2/ENCSR000COQ2.final.uniq.bam.bai]]
    */
    .join(prepare_vcf_for_ase.out.vcf_for_ASE)
    /* Outputs
    [ENCSR000COQ, [/workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam, /workspace/gitpod/hands-on/work/92/b1ea340ce922d13bdce2985b2930f2/ENCSR000COQ2.final.uniq.bam], [/workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam.bai, /workspace/gitpod/hands-on/work/92/b1ea340ce922d13bdce2985b2930f2/ENCSR000COQ2.final.uniq.bam.bai], /workspace/gitpod/hands-on/work/ea/4a41fbeb591ffe48cfb471890b8f5c/known_snps.vcf]
    */
    .map { meta, bams, bais, vcf -> [meta, vcf, bams, bais] }
    /* Outputs
    [ENCSR000COQ, /workspace/gitpod/hands-on/work/ea/4a41fbeb591ffe48cfb471890b8f5c/known_snps.vcf, [/workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam, /workspace/gitpod/hands-on/work/92/b1ea340ce922d13bdce2985b2930f2/ENCSR000COQ2.final.uniq.bam], [/workspace/gitpod/hands-on/work/c9/dfd66e253754b61195a166ac7726ff/ENCSR000COQ1.final.uniq.bam.bai, /workspace/gitpod/hands-on/work/92/b1ea340ce922d13bdce2985b2930f2/ENCSR000COQ2.final.uniq.bam.bai]]
    */
    .set { grouped_vcf_bam_bai_ch }
```

Just keep in mind processes and channel operators are not guaranteed to emit items in the order that they were received, as they are executed concurrently. This can lead to unintended effects based if you use a operator that takes multiple inputs. For example, the using the `merge` channel operator can lead to different results upon different runs based on the order in which the processes finish. You should always use a matching key (e.g. sample ID) to merge multiple channels, so that they are combined in a deterministic way (using an operator like `join` instead). 

### The `.flatten` channel operator

The `.flatten` operator transforms a channel in such a way that every item of type `Collection` or `Array` is flattened so that each single entry is emitted separately by the resulting channel. For example:

```nextflow
Channel
    .of( [1,[2,3]], 4, [5,[6]] )
    .flatten()
    .view()
```

Outputs

```
1
2
3
4
5
6
```

### The `.collect` channel operator

The `.collect` channel operator is basically the oposite of the `.flatten` channel operator, where it collects all the items emitted by a channel to a `List` and return the resulting object as a sole emission. For example:

```nextflow
Channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
```

Outputs

```
[1,2,3,4]
```

### The `.buffer` channel operator

The `.buffer` channel operator gathers the items emitted by the source channel into subsets and emits these subsets separately.

There are a number of ways you can regulate how buffer gathers the items from the source channel into subsets, however one of the most convienient ways of using it is with `buffer( size: n )`. transform the source channel in such a way that it emits tuples made up of `n` elements. For example:

```nextflow
Channel
    .of( 1, 2, 3, 1, 2, 3, 1 )
    .buffer( size: 2 )
    .view()
```

Outputs

```
[1, 2]
[3, 1]
[2, 3]
```

Be aware that if there is an incomplete tuple it is discarded. To emit the last items in a tuple containing less than n elements, use `buffer( size: n, remainder: true )`. For example:

```nextflow
Channel
    .of( 1, 2, 3, 1, 2, 3, 1 )
    .buffer( size: 2, remainder: true )
    .view()
```

Outputs

```
[1, 2]
[3, 1]
[2, 3]
[1]
```

## Channel Factories

Channel factories are methods that can be used to create channels explicitly.

For example, the `of` method allows you to create a channel that emits the arguments provided to it, for example:
```nextflow
ch = Channel.of( 1, 3, 5, 7 )
ch.view { "value: $it" }
```

The first line in this example creates a variable `ch` which holds a channel object. This channel emits the arguments supplied to the of method. Thus the second line prints the following:
```nextflow
value: 1
value: 3
value: 5
value: 7
```

### Range Expansion 

Ranges of values are expanded accordingly:

```nextflow
Channel
    .of(1..23, 'X', 'Y')
    .view()
```
Prints:

```nextflow
1
2
3
4
:
23
X
Y
```

## Tuples

Inputs and outputs in nextflow need to a data type assigned before a variable name. If the data type is a tuple, all the items in the tuple need a data type as well.

Example

```nextflow
process example_process {
    input:
    tuple val(replicateId), path(reads)

    output:
    tuple val(replicateId),
          path('Aligned.sortedByCoord.out.bam'),
          path('Aligned.sortedByCoord.out.bam.bai')

    script:
    """
    ...
    """
}
```

## Value Channels vs Queue Channels

Value channels allow processes to consume elements infinite times.

Elements in the queue channel are consumed. You cannot consume the same element for the same process twice.

## Single Value Coersion to Value Channels

Processes can only take channels as inputs. That being said, if you pass in a regular variable that has a single value (e.g. string, number, etc), nextflow will implicitly create a value channel containing that value that you are providing. Therefore you can pass a regular value (i.e. string, number) into a process, but just be aware that this is what is going on behind the scenes

Example

```nextflow
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

## The Working Directory and the `publishDir` directive

By default, the task work directories are created in the directory from where the pipeline is launched. This is often a scratch storage area that can be cleaned up once the computation is completed. A different location for the execution work directory can be specified using the command line option `-w` 

Example

```bash
nextflow run <script> -w /some/scratch/dir
```

Note that if you delete or move the pipeline work directory, this will prevent to use the resume feature in subsequent runs.

Also note that the pipeline work directory is intended to be used as a temporary scratch area. The final workflow outputs are expected to be stored in a different location specified using the `publishDir` directive.

## The Golden Practice for Containers

If you just try to install your dependencies using conda or mamba, there's no guarantee that the dependency graph will be solved the same way forever, therefore there is a chance your programs will work slightly different if you try to reinstall your dependencies sometime in the future.

That is why the golden rule for dependencies is to use a program like conda or mamba to intall them, but to do so inside a container image. 

The container images are supposed to work the same way forever so doing so in this way each programs should install the same dependencies and produce the exact same enviroment for better reproducibility. 

E.g.

`env.yml`

```yml
name: nf-tutorial
channels:
    - conda-forge
    - defaults
    - bioconda
dependencies:
    - bioconda::salmon=1.5.1
    - bioconda::fastqc=0.11.9
    - bioconda::multiqc=1.12
    - conda-forge::tbb=2020.2
```

`Dockerfile`
```docker
FROM mambaorg/micromamba:0.25.1

LABEL image.author.name "Your Name Here"
LABEL image.author.email "your@email.here"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba create -n nf-tutorial

RUN micromamba install -y -n nf-tutorial -f /tmp/env.yml && \
    micromamba clean --all --yes

ENV PATH /opt/conda/envs/nf-tutorial/bin:$PATH
```

## Biocontainers

Each program should have its own designated container. Don't create container images with too many things or things your don't need.

Run only one process per container: In almost all cases, you should only run a single process in a single container. Decoupling applications into multiple containers makes it much easier to scale horizontally and reuse containers. If that service depends on another service, make use of container linking.

Biocontainers is a project to create a docker container for every recipe they have in bioconda.

Sometimes you'll need to have a container with more than one tool, in this case there is a *mulled* container. 

You can request a multi-package container here: https://biocontainers.pro/multipackage

### Mulled Container Search

You can also install `galaxy-util-tools` and search for mulled containers in your CLI using the `mulled-search` tool that comes with this tool

```bash
conda activate a-conda-env-you-already-have
conda install galaxy-tool-util
mulled-search --destination quay singularity --channel bioconda --search bowtie samtools | grep mulled
```

## `get*()` in Groovy

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

## Process Directives

Process directives are optional settings that affect the execution of the process. They are written at the top of a process block.

Example

```nextflow
/*
 * Process 1C: Create the genome index file for STAR
 */

process prepare_star_genome_index {
    container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'

    input:
    path genome

    output:
    path 'genome_dir' 

    script: 
    """
    mkdir -p genome_dir

    STAR --runMode genomeGenerate \
         --genomeDir genome_dir \
         --genomeFastaFiles ${genome} \
         --runThreadN ${task.cpus}
    """
}
```

The `container ` process directive tells nextflow that if it is using docker, then to use that specific container for this specific task.

Process directives can also be defaults and invisible. For example. The default number of CPUs to run a task is 1. You can pass this into a process parameter by using `tasks.cpus`. If you want to change this however, you can write at the top of your process block

```nextflow
cpus = 4
```

### The `.tag` process directive

Whenever you have multiple samples being processed in a process, it's useful to use the `.tag` process directive.

The `.tag` process directive doesn't change anything in the analysis, but it allows you to associate each process execution with a custom label to make it easier to identify them in the log file or the trace execution report. This is useful if you want to know more information about multiple samples that are being run (which ones passed and which ones failed). 

Example

```nextflow
/*
 * Process 3: GATK Split on N
 */

process rnaseq_gatk_splitNcigar {
    container 'quay.io/broadinstitute/gotc-prod-gatk:1.0.0-4.1.8.0-1626439571'
    tag "${replicateId}" 

    input:
    path genome 
    path index 
    path genome_dict 
    tuple val(replicateId), path(bam), path(bai) 

    output:
    tuple val(replicateId), path('split.bam'), path('split.bai') 

    script:
    """
    java -jar /usr/gitc/GATK35.jar -T SplitNCigarReads \
                                   -R ${genome} -I ${bam} \
                                   -o split.bam \
                                   -rf ReassignOneMappingQuality \
                                   -RMQF 255 -RMQT 60 \
                                   -U ALLOW_N_CIGAR_READS \
                                   --fix_misencoded_quality_scores
    """
}
```

Will give you detials on every replicate ID that is being processed.

## Only Use Positional Arguments 

As of December 2023, Nextflow does not allow for keyword arguments. Therefore if you're trying to for example pass in 2 inputs, you look at the process block and the first input will be the first positional argument and the second input will be the second positional argument (and so on and so forth if there are more inputs)

Example

```nextflow
/*
 * Process 1D: Create a file containing the filtered and recoded set of variants
 */

process prepare_vcf_file {
    container 'quay.io/biocontainers/mulled-v2-b9358559e3ae3b9d7d8dbf1f401ae1fcaf757de3:ac05763cf181a5070c2fdb9bb5461f8d08f7b93b-0'

    input:
    path variantsFile
    path blacklisted

    output: 
    tuple path("${variantsFile.baseName}.filtered.recode.vcf.gz"),
          path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")

    script:
    """
    vcftools --gzvcf ${variantsFile} -c \
             --exclude-bed ${blacklisted} \
             --recode | bgzip -c \
             > ${variantsFile.baseName}.filtered.recode.vcf.gz

    tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads)

    prepare_genome_samtools(params.genome)
    prepare_genome_picard(params.genome)
    prepare_star_genome_index(params.genome)
    // CORRECT
    prepare_vcf_file(params.variants, params.blacklist) 
    // INCORRECT
    prepare_vcf_file(variantsFile = params.variants, blacklisted = params.blacklist) 
}
```

## You can run code before you execute a script in the script part of a process

Example

```nextflow
/*
 * Process 4: GATK Recalibrate
 */

process rnaseq_gatk_recalibrate {
    container 'quay.io/biocontainers/mulled-v2-aa1d7bddaee5eb6c4cbab18f8a072e3ea7ec3969:f963c36fd770e89d267eeaa27cad95c1c3dbe660-0'
    tag "${replicateId}"

    input:
    path genome
    path index
    path dict
    tuple val(replicateId), path(bam), path(bai) 
    tuple path(prepared_variants_file),
          path(prepared_variants_file_index) 

    output: 
    tuple val(sampleId),
          path("${replicateId}.final.uniq.bam"),
          path("${replicateId}.final.uniq.bam.bai")

    script:
    // Notice how the sampleID is modified with a regular expression before the script is executed
    sampleId = replicateId.replaceAll(/[12]$/,'') 
    """
    gatk3 -T BaseRecalibrator \
          --default_platform illumina \
          -cov ReadGroupCovariate \
          -cov QualityScoreCovariate \
          -cov CycleCovariate \
          -knownSites ${prepared_variants_file} \
          -cov ContextCovariate \
          -R ${genome} -I ${bam} \
          --downsampling_type NONE \
          -nct ${task.cpus} \
          -o final.rnaseq.grp

    gatk3 -T PrintReads \
          -R ${genome} -I ${bam} \
          -BQSR final.rnaseq.grp \
          -nct ${task.cpus} \
          -o final.bam

    (samtools view -H final.bam; samtools view final.bam | \
    grep -w 'NH:i:1') | samtools view -Sb -  > ${replicateId}.final.uniq.bam

    samtools index ${replicateId}.final.uniq.bam
    """
}
```

## The `return` keyword

You can place the `return` keyword in your workflow block to aburptly stop the run at a certain point. This is really useful if you're troubleshooting and you just want to get the details up to a certain part of your workflow, the rest of it will remain untouched.

Example

```nextflow
workflow {
    reads_ch =  Channel.fromFilePairs(params.reads)
    // The return keyword helps you troubleshoot. 
    // Here you want to check what the read ares
    // The return keyword ensures the rest of the workflow doesn't get run
    reads.view()
    return 


    prepare_genome_samtools(params.genome)
    prepare_genome_picard(params.genome)
    prepare_star_genome_index(params.genome)
    prepare_vcf_file(params.variants, params.blacklist)

    rnaseq_mapping_star(params.genome, prepare_star_genome_index.out, reads_ch)

    rnaseq_gatk_splitNcigar(params.genome,
                            prepare_genome_samtools.out,
                            prepare_genome_picard.out,
                            rnaseq_mapping_star.out)

   rnaseq_gatk_recalibrate(params.genome,
                           prepare_genome_samtools.out,
                           prepare_genome_picard.out,
                           rnaseq_gatk_splitNcigar.out,
                           prepare_vcf_file.out)

    // New channel to aggregate bam from different replicates into sample level.
    rnaseq_gatk_recalibrate.out
        | groupTuple
        | set { recalibrated_samples }

    rnaseq_call_variants(params.genome,
                         prepare_genome_samtools.out,
                         prepare_genome_picard.out,
                         recalibrated_samples)

    post_process_vcf(rnaseq_call_variants.out,
                        prepare_vcf_file.out)

    prepare_vcf_for_ase(post_process_vcf.out)

    recalibrated_samples
        .join(prepare_vcf_for_ase.out.vcf_for_ASE)
        .map { meta, bams, bais, vcf -> [meta, vcf, bams, bais] }
        .set { grouped_vcf_bam_bai_ch }

    ASE_knownSNPs(params.genome,
                  prepare_genome_samtools.out,
                  prepare_genome_picard.out,
                  grouped_vcf_bam_bai_ch)
}
```

## Getting Syntax Examples using nf-core

When you don't know how a specific nextflow function or element works, a really good resource is seeing how it was implemented in nf-core. The nf-core repository contains dozens of professional and expertly curated pipelines. By going to https://github.com/nf-core and typing into the search bar a specific function or operator, you can get tons and tons of examples of how it is supposed to be used.

### Nextflow Tower

Nextflow tower provides a useful tool to monitor your pipelines and work on teams. It offered through a web app and a CLI. To learn more about how to use it you can click here:

- Web app:  
    
    - https://youtu.be/ERbTqLtAkps?feature=shared&t=5852

- CLI:

    - https://youtu.be/ERbTqLtAkps?feature=shared&t=6556

