# RiboViz and Nextflow

* [Nextflow](https://www.nextflow.io/)
* [GitHub](https://github.com/nextflow-io/nextflow)
* [Documentation](https://www.nextflow.io/docs/latest/index.html)

See notes on:

* [Nextflow Documentation Notes](./NextflowDocNotes.md) from reading through the Nextflow documentation.

---

## Run the RiboViz example

### Install Nextflow using conda (recommended)

Create a new conda environment from the current RiboViz one and activate it:

```console
$ conda create --name riboviz-nextflow --clone base
$ conda activate riboviz-nextflow
```

Install Nextflow:

```console
$ conda install -c bioconda nextflow
```

Check install:

```console
$ which javac
/home/ubuntu/miniconda3/envs/riboviz-nextflow/bin/javac
$ javac -version
javac 1.8.0_152-release
$ which java
/home/ubuntu/miniconda3/envs/riboviz-nextflow/bin/java
$ java -version
openjdk version "1.8.0_152-release"
OpenJDK Runtime Environment (build 1.8.0_152-release-1056-b12)
OpenJDK 64-Bit Server VM (build 25.152-b12, mixed mode)
$ which nextflow
/home/ubuntu/miniconda3/envs/riboviz-nextflow/bin/nextflow
$ nextflow -version

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (02:14 PDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

### Install Nextflow "by hand" (alternative)

Install [OpenJDK](https://openjdk.java.net) 1.8:

* CentOS 7 users:

```console
$ sudo yum install -y openjdk-8-jdk-headless
```

* Ubuntu 18 users:

```console
$ sudo apt-get install -y openjdk-8-jdk-headless
```

Check install:

```console
$ which javac
/usr/bin/javac
$ javac -version
javac 1.8.0_242
$ which java
/usr/bin/java
$ java -version
openjdk version "1.8.0_242"
OpenJDK Runtime Environment (build 1.8.0_242-8u242-b08-0ubuntu3~18.04-b08)
OpenJDK 64-Bit Server VM (build 25.242-b08, mixed mode)
```

Install Nextflow:

```console
$ curl -s https://get.nextflow.io | bash
$ export PATH=$HOME/nextflow:$PATH
$ which nextflow
/home/ubuntu/nextflow/nextflow
$ nextflow -version

      N E X T F L O W
      version 20.01.0 build 5264
      created 12-02-2020 10:14 UTC (02:14 PDT)
      cite doi:10.1038/nbt.3820
      http://nextflow.io
```

Remember to set `PATH` before running:

```console
$ export PATH=$HOME/nextflow:$PATH
```

### Run Nextflow "hello" example

```console
$ nextflow run hello
N E X T F L O W  ~  version 20.01.0
Pulling nextflow-io/hello ...
downloaded from https://github.com/nextflow-io/hello.git
Launching `nextflow-io/hello` [spontaneous_magritte] - revision: 1d43afc0ec [master]
WARN: The use of `echo` method is deprecated
executor >  local (4)
[1d/bb459e] process > sayHello [100%] 4 of 4 /
Hola world!

Bonjour world!

Ciao world!

Hello world!
```

This runs https://github.com/nextflow-io/hello/blob/master/main.nf from https://github.com/nextflow-io/hello.git. The Nextflow script has a `sayHello` task that is run on 4 inputs - the multi-lingual greetings - and runs bash `echo` to print the greetings. Each iteration of a task on each input is run as a separate process.

### Run the RiboViz example

[riboviz.nf](./riboviz.nf) is a subset of the RiboViz workflow which implements the following steps:

1. Read configuration information from YAML configuration file.
2. Build hisat2 indices if requested.
3. Process each sample ID-sample file pair in turn:
   1. Cut out sequencing library adapters using `cutadapt`.
   2. Remove rRNA or other contaminating reads by alignment to rRNA index files using `hisat2`.
   3. Align remaining reads to ORFs index files using `hisat2`.
   4. Trim 5' mismatches from reads and remove reads with more than 2 mismatches using `riboviz.tools.trim_5p_mismatch`.

Copy over Nextflow script:

```console
$ cd riboviz
$ cp ~/workflows/nextflow/riboviz.nf .
```

Run:

```console
$ PYTHONPATH=$HOME/riboviz nextflow run riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false -with-report report.html -with-timeline timeline.html -with-dag workflow.svg
N E X T F L O W  ~  version 20.01.0
Launching `riboviz.nf` [furious_bell] - revision: 134fcb8f6f
Missing file (NotHere): example_missing_file.fastq.gz
[5c/676372] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[f7/6fcdf1] Submitted process > cutAdapters (WT3AT)
[e5/ccf3e6] Submitted process > buildIndicesrRNA (yeast_rRNA)
[5b/ff17ba] Submitted process > cutAdapters (WTnone)
[ad/1e7c54] Submitted process > hisat2rRNA (WTnone)
[e4/deff35] Submitted process > hisat2rRNA (WT3AT)
[1b/642455] Submitted process > hisat2ORF (WTnone)
[b9/aaf172] Submitted process > trim5pMismatches (WTnone)
[e8/8514a6] Submitted process > hisat2ORF (WT3AT)
[23/a980bc] Submitted process > trim5pMismatches (WT3AT)
[b4/837e33] Submitted process > summarise
Processed samples: WTnone WT3AT 
```

Note: `PYTHONPATH=$HOME/riboviz` is required so that, when Nextflow invokes `python -m riboviz.tools.trim_5p_mismatch`, Python can find the `riboviz` module.

Every invocation of a task - every process - has its own subdirectory within Nextflow's `work/` directory named after the process identifiers (e.g. `ad/1e7c54`). These subdirectories have:

* Input files. These are symbolic links to the input files for the task which, depending on the task, can be:
  - Output files in other `work/` subdirectories. For example, the directory for an `hisat2rRNA` proces will have input files which are symbolic links to the output files produced by a `cutAdapters` process,
  - Input files for the workflow. For example, the directory for a `cutAdapters` process will have an input file which is a symbolic link to a sample file in `vignettte/input`.
* Output files, from the invocation of the task.
* Log files and bash scripts (see below).

For example, for the process `ad/1e7c54`, an invocation of task `hisat2rRNA` for sample `WTnone`, the `work/` directory includes:

```console
$ find work/ad/1e7c54a889f21451cb07d29655e0be/ -printf '%P\t%l\n' | sort
.command.begin	
.command.err	
.command.log	
.command.out	
.command.run	
.command.sh	
.command.trace	
.exitcode	
nonrRNA.fq	
rRNA_map.sam	
trim.fq	/home/ubuntu/riboviz/work/5b/ff17ba3e19c5d659e54b04b08dab85/trim.fq
yeast_rRNA.1.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.1.ht2
yeast_rRNA.2.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.2.ht2
yeast_rRNA.3.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.3.ht2
yeast_rRNA.4.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.4.ht2
yeast_rRNA.5.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.5.ht2
yeast_rRNA.6.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.6.ht2
yeast_rRNA.7.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.7.ht2
yeast_rRNA.8.ht2	/home/ubuntu/riboviz/work/e5/ccf3e6388cde7038658d88a79e81d1/yeast_rRNA.8.ht2
```

The `.ht2` files are symbolic links to the outputs of process `e5/ccf3e6`, an invocation of task `buildIndicesrRNA`.

`riboviz.nf` uses Nextflow's [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directive which allows files to be published to a specific directory. By default, the files in this directory are symlinked to those in `work/`. `riboviz.nf` uses this to publishes the outputs to `vignette/index/` and `vignette/tmp/`. 

### `-with-report`, `-with-timeline`, `-with-dag`

These flag generate reports on the run and an image of the task execution workflow: See, for example:

* [report.html](./report.html)
* [timeline.html](./timeline.html)
* [workflow.svg](./workflow.svg)

### `-ansi-log`

The `-ansi-log` parameter is enabled by default. When Nextflow is run it displays and updates progress information in-place within the bash window. For example, for the above, the output, on completion, would look like:

```console
[5c/676372] process > buildIndicesORF  [100%] 1 of 1 /
[e5/ccf3e6] process > buildIndicesrRNA [100%] 1 of 1 /
[5b/ff17ba] process > cutAdapters      [100%] 2 of 2 /
[e4/deff35] process > hisat2rRNA       [100%] 2 of 2 /
[e8/8514a6] process > hisat2ORF        [100%] 2 of 2 /
[23/a980bc] process > trim5pMismatches [100%] 2 of 2 /
[b4/837e33] process > summarise        [100%] 1 of 1 /
```

To find out which `work/` subdirectory holds the outputs for each step requires inspecting the Nextflow log file. For example:

```console
$ grep INFO .nextflow.log
Mar-27 05:32:04.711 [main] INFO  nextflow.cli.CmdRun - N E X T F L O W  ~  version 20.01.0
Mar-27 05:32:04.722 [main] INFO  nextflow.cli.CmdRun - Launching `riboviz.nf` [furious_bell] - revision: 134fcb8f6f
Mar-27 05:32:06.377 [Task submitter] INFO  nextflow.Session - [5c/676372] Submitted process > buildIndicesORF (YAL_CDS_w_250)
Mar-27 05:32:06.386 [Task submitter] INFO  nextflow.Session - [f7/6fcdf1] Submitted process > cutAdapters (WT3AT)
Mar-27 05:32:06.397 [Task submitter] INFO  nextflow.Session - [e5/ccf3e6] Submitted process > buildIndicesrRNA (yeast_rRNA)
Mar-27 05:32:06.404 [Task submitter] INFO  nextflow.Session - [5b/ff17ba] Submitted process > cutAdapters (WTnone)
Mar-27 05:32:15.571 [Task submitter] INFO  nextflow.Session - [ad/1e7c54] Submitted process > hisat2rRNA (WTnone)
Mar-27 05:32:17.530 [Task submitter] INFO  nextflow.Session - [e4/deff35] Submitted process > hisat2rRNA (WT3AT)
Mar-27 05:32:24.110 [Task submitter] INFO  nextflow.Session - [1b/642455] Submitted process > hisat2ORF (WTnone)
Mar-27 05:32:27.106 [Task submitter] INFO  nextflow.Session - [b9/aaf172] Submitted process > trim5pMismatches (WTnone)
Mar-27 05:32:29.389 [Task submitter] INFO  nextflow.Session - [e8/8514a6] Submitted process > hisat2ORF (WT3AT)
Mar-27 05:32:31.368 [Task submitter] INFO  nextflow.Session - [23/a980bc] Submitted process > trim5pMismatches (WT3AT)
Mar-27 05:32:31.916 [Task submitter] INFO  nextflow.Session - [b4/837e33] Submitted process > summarise
```

Setting `-ansi-log` to `false`, as above, explicitly shows the progress as separate steps in a bash window.

---

## Assessment

### Ease of download, install, tutorials, initial use

Nextflow and its dependencies are straightforward to download and install (whether using conda or a native package manager). The "hello" example and [Get started](https://www.nextflow.io/docs/latest/getstarted.html) tutorial - writing a workflow with two tasks - are straightforward to follow.

### Ease of implementation of key RiboViz steps

The user documentation took ~1 day to read through to get a full feel of Nextflow's capabilities. An expanded tutorial, to implement a multi-task workflow, analogous to that provided by Snakemake, would have been useful. I did, later on, find Nextflow 2017 workshop [tutorial](https://github.com/nextflow-io/nf-hack17-tutorial) but my understanding of Nextflow had passed beyond its content when I found it.

Implementation of a subset of the RiboViz workflow took ~2 days. Despite having read the documentation, a lot of searching was required. The most challenging aspect was trying to have the index files reused when processing each sample file, but the solution, once the command was found, was straightforward (the `each` [input repeater](https://www.nextflow.io/docs/latest/process.html#input-repeaters)). I think progress to complete the workflow would be more rapid as more experience of using, as opposed to reading about, Nextflow's constructs is gained.

There are myriad examples to draw upon, for example:

* Examples menu on [Nextflow](https://www.nextflow.io/index.html).
* Documentation [Example](https://www.nextflow.io/docs/latest/example.html).
* Documentation [FAQ](https://www.nextflow.io/docs/latest/faq.html).
* Collection of Nextflow implentation [patterns](https://github.com/nextflow-io/patterns).

The Nextflow scripts are readable and makes clear the key steps in a workflow.

### Sample-specific sub-directories

Every invocation of a task - every process - has its own subdirectory within Nextflow's `work/` directory.

`work/` subdirectory names are cryptic as they are auto-generated. For example

```
work/1b/ff0a43332a3feaca02c998cc354878
work/20/2975a5ab45a6707685a85f5aa219f3
work/37/b11ee1d2fb315a1b72adb65c151b44
work/3e/98ec992925cc16885b3dd12967a532
```

The [tag](https://www.nextflow.io/docs/latest/process.html#tag) directive can be used to associate meaningful names with tasks and processes. These, in conjunction with the Nextflow log file, can make it easier to find out which `work/` subdirectory has the files for which process.

Nextflow supports a [publishDir](https://www.nextflow.io/docs/latest/process.html#publishdir) directive which allows files to be published to a specific directory (e.g. one with a more readable name). By default, the files in this directory are symlinked to those in `work/`, but they can also be, for example, copied or moved there.

### Parse YAML configuration files

Nextflow supports a `-params-file` command-line option which allows a YAML or JSON file to be provided. The parameters can then be accessed within a Nextflow script via a `params.` prefix e.g. `params.rrna_fasta_file`.

### Dry run option, validating configuration and input file existence

A "dry run" option was suggested in a Nextflow issue [suggestion: run -dry #31](https://github.com/nextflow-io/nextflow/issues/31) of 12/02/15. A Nextflow author commented that:

> The problem is that it's pretty [difficult?] to implement in the dataflow model on which Nextflow is based, because you need always an event i.e. a piece of data to trigger a process/operator execution.
>
> No plan to implement it due to the reason explained above. We may take in consideration if we find a way to handle a dry-run execution on top of the dataflow paradigm.

No "dry run" option is supported at present. On the issue above, the Nextflow author comments that:

> In my opinion the easiest workaround is to have always in your pipeline a tiny dataset that can be used to validate the script execution.

### Tool/step-specific log files

Standard output and standard error (as well as exit codes) for each process are automatically captured and placed in subdirectories of Nextflow's `work/` directory. For example, for a process `37/b11ee1`, an invocation of task `hisat2rRNA` for sample `WTnone`, the `work/` directory includes:

```console
$ ls -A work/37/b11ee1d2fb315a1b72adb65c151b44/
.command.begin	
.command.err	
.command.log	
.command.out	
.command.run	
.command.sh	
.command.trace	
.exitcode	
...

$ cat work/37/b11ee1d2fb315a1b72adb65c151b44/.command.out 
/home/ubuntu/hisat2-2.1.0/hisat2-align-s version 2.1.0
64-bit
Built on login-node03
Wed Jun  7 15:53:42 EDT 2017
Compiler: gcc version 4.8.2 (GCC) 
Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}

$ cat work/37/b11ee1d2fb315a1b72adb65c151b44/.command.err 
952343 reads; of these:
  952343 (100.00%) were unpaired; of these:
    467194 (49.06%) aligned 0 times
    485149 (50.94%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
50.94% overall alignment rate

$ cat work/37/b11ee1d2fb315a1b72adb65c151b44/.exitcode 
0
```

### Output bash script, that can be rerun

Nextflow does not output a bash script for a whole workflow. However each processes subdirectory within the `work/` directory has a bash script which runs that specific operation of the workflow. For example:

```console
$ cd work/37/b11ee1d2fb315a1b72adb65c11b44/
$ cat .command.sh 
#!/bin/bash -ue
hisat2 --version
hisat2 -p 1 -N 1 -k 1 --un nonrRNA.fq -x yeast_rRNA -S rRNA_map.sam -U trim.fq
```

As each subdirectory has symbolic links to any input files it requires, this means that the step can be run standalone. This can be used for debugging. For example:

```console
$ bash .command.sh 
/home/ubuntu/hisat2-2.1.0/hisat2-align-s version 2.1.0
64-bit
Built on login-node03
Wed Jun  7 15:53:42 EDT 2017
Compiler: gcc version 4.8.2 (GCC) 
Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}
952343 reads; of these:
  952343 (100.00%) were unpaired; of these:
    467194 (49.06%) aligned 0 times
    485149 (50.94%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
50.94% overall alignment rate
```

### Access configuration file name from within a Nextflow script

There doesn't seem to be a way to access the configuration file itself from within a Nextflow script i.e. the value of `-params-file`.

`riboviz.tools.count_reads` needs the RiboViz configuration file for sample ID-file name values in `fq_files`. These could be written to another file from within the Nextflow script, or via another script invoked by Nextflow, and that file provided as the input to `riboviz.tools.count_reads`.

### If processing of one sample fails will the rest be processed?

Each task can have an [errorStrategy](https://www.nextflow.io/docs/latest/process.html#errorstrategy) defined which indicates how a process, and the workflow, should behave if an error is encountered. `ignore` allows for other processes to continue if one fails. This can be used to ensure that other samples are processed if processing of one sample fails. For example, if the file for sample `WTnone` was corrupt in some way:

```console
$ PYTHONPATH=$HOME/riboviz nextflow run riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `riboviz.nf` [mad_heisenberg] - revision: 134fcb8f6f
Missing file (NotHere): example_missing_file.fastq.gz
[a4/e3fc5d] Submitted process > cutAdapters (WT3AT)
[03/f1cfc3] Submitted process > cutAdapters (WTnone)
[9d/43b3e0] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[ec/1aafe5] Submitted process > buildIndicesrRNA (yeast_rRNA)
[03/f1cfc3] NOTE: Process `cutAdapters (WTnone)` terminated with an error exit status (1) -- Error is ignored
[35/1c7ac8] Submitted process > hisat2rRNA (WT3AT)
[1a/af8ff1] Submitted process > hisat2ORF (WT3AT)
[c4/2f18a2] Submitted process > trim5pMismatches (WT3AT)
[df/28719b] Submitted process > summarise
Processed samples: WT3AT 
```

### Conditional behaviour

Nextflow supports a `when` declaration in tasks to allow tasks to only be executed if certain conditions hold.

### Aggregation of sample-specific results.

Nextflow supports a number of [combining operators](https://www.nextflow.io/docs/latest/operator.html#combining-operators) which can be used to aggregate together results. 

The [collect](https://www.nextflow.io/docs/latest/operator.html#collect) operator can be used to collect all items from a channel (e.g. the names of all successfully processed samples) and output the list as a single item. The `summarise` task in `riboviz.nf` has an example of its use.

### Containers, clusters and cloud

Nextflow has support for:

* Running jobs in [Docker](https://www.nextflow.io/docs/latest/docker.html?highlight=docker) and [Singularity](https://www.nextflow.io/docs/latest/singularity.html) containers.
* [executors](https://www.nextflow.io/docs/latest/executor.html) which manage execution of a Nextflow script on underyling execution systems. Executors are provided for:
  - Sun Grid Engine/Open Grid Engine/Univa Grid Engine - (Open) Grid Engine is used on the University of Edinburgh's [ECDF Linux Compute Cluster (Eddie)](https://www.ed.ac.uk/information-services/research-support/research-computing/ecdf/high-performance-computing).
  - PBS Pro - as used by [Cirrus](https://www.cirrus.ac.uk) and [ARCHER](https://www.archer.ac.uk).
  - SLURM resource manager
  - Kubernetes
  - AWS Batch
  - Many others.
* Cloud support, "out-of-the-box", for:
  - [Amazon Cloud](https://www.nextflow.io/docs/latest/awscloud.html)
  - [Google Cloud](https://www.nextflow.io/docs/latest/google.html)

### Other

Nextflow requires Java JDK 1.8+ to run.

Writing scripts requires knowledge of Groovy - a Python-esque language based on Java. I think anyone familiar with Python or R would not find Groovy too challenging. I hadn't used it myself before looking at Nextflow and though I have a number of years of experience with Java, I didn't really need to draw upon that knowledge when writing the Nextflow script.

Nextflow has built-in functions for FASTA and FASTQ files e.g. count records, extract records, query the NCBI SRA database for specific FASTQ files, which may be of interest or use in future.

If a Nextflow workflow fails then a `-resume` option allows it to be rerun from the point at which it failed:

```console
$ nextflow run riboviz.nf -params-file vignette/vignette_config.yaml -resume
```

This also supports incremental build. For example, given a `vignette_config.yaml` which specifies only sample `WTnone`, running Nextflow gives:

```console
$ PYTHONPATH=$HOME/riboviz nextflon riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false
N E X T F L O W  ~  version 20.01.0
Launching `riboviz.nf` [serene_noether] - revision: 26bb98ec7b
Missing file (NotHere): example_missing_file.fastq.gz
[bf/ddd7de] Submitted process > cutAdapters (WTnone)
[0d/f82601] Submitted process > buildIndicesrRNA (yeast_rRNA)
[61/5a1186] Submitted process > buildIndicesORF (YAL_CDS_w_250)
[02/dac1f8] Submitted process > hisat2rRNA (WTnone)
[15/e5f9ed] Submitted process > hisat2ORF (WTnone)
[ed/fc50b9] Submitted process > trim5pMismatches (WTnone)
[c0/29b60a] Submitted process > summarise
Processed samples: WTnone 
```

If `WT3AT` is then added to `vignette_config.yaml` and Nextflow is run with the `-resume` option, then only the processing for `WT3AT` is done, the cached outputs to date for `WTnone` being reused:

```console
$ PYTHONPATH=$HOME/riboviz nextflow run riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false -resume
N E X T F L O W  ~  version 20.01.0
Launching `riboviz.nf` [desperate_coulomb] - revision: 26bb98ec7b
Missing file (NotHere): example_missing_file.fastq.gz
[80/79a5f4] Submitted process > cutAdapters (WT3AT)
[bf/ddd7de] Cached process > cutAdapters (WTnone)
[61/5a1186] Cached process > buildIndicesORF (YAL_CDS_w_250)
[0d/f82601] Cached process > buildIndicesrRNA (yeast_rRNA)
[02/dac1f8] Cached process > hisat2rRNA (WTnone)
[15/e5f9ed] Cached process > hisat2ORF (WTnone)
[ed/fc50b9] Cached process > trim5pMismatches (WTnone)
[e3/9c39d7] Submitted process > hisat2rRNA (WT3AT)
[a6/49addd] Submitted process > hisat2ORF (WT3AT)
[46/d5db2c] Submitted process > trim5pMismatches (WT3AT)
[a5/5c6c05] Submitted process > summarise
Processed samples: WTnone WT3AT 
```
