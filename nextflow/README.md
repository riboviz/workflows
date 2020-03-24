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

Remeber to set `PATH` before running:

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
[1d/bb459e] process > sayHello [100%] 4 of 4 ?
Hola world!

Bonjour world!

Ciao world!

Hello world!
```

This runs https://github.com/nextflow-io/hello/blob/master/main.nf from https://github.com/nextflow-io/hello.git. The Nextflow script has a `sayHello` task that is run on 4 inputs - the multi-lingual greetings - and runs bash `echo` to print the greetings. Each iteration of task on each input is run as a separate process.

As an example of what goes on "under the hood", Nextflow caches the standard output into 4 directories - one per process - within a `.command.out` file:

```console
$ ls work/*/*/.command.out
work/4b/0f3c4b0c0d453c3177ca2ede1efa14/.command.out
work/8f/c8dd05e613586574e434a4506de6f9/.command.out
work/b4/2fab9eb37fe3901671fa5eee8aaaf0/.command.out
work/f7/bb1e3a1f84a2b6f03c94ec6932b760/.command.out
$ cat work/*/*/.command.out 
Ciao world!
Hello world!
Hola world!
Bonjour world!
```

### Run the RiboViz example

Copy over Nextflow script:

```console
$ cd riboviz
$ cp ~/workflows/nextflow/riboviz.nf .
```

TODO add riboviz.nf to repository.

Run:

```console
$ nextflow run riboviz.nf -params-file vignette/vignette_config.yaml
```

Create HTML report, timeline, pipeline graph:

```console
$ nextflow run riboviz.nf -params-file vignette/vignette_config.yaml -with-report report.html -with-timeline timeline.html -with-dag workflow.svg
```

See, for example:

* [report.html](./report.html)
* [timeline.html](./timeline.html)
* [workflow.svg](./workflow.svg)

TODO add files to repository.

---

## Assessment

### Ease of download, install, tutorials, initial use

Nextflow and its dependencies are straightforward to download and install (whether using conda or a native package manager). The "hello" example and [Get started](https://www.nextflow.io/docs/latest/getstarted.html) tutorial - writing a workflow with two tasks - are straightforward to follow.

### Ease of implementation of key RiboViz steps

The user documentation took ~1 day to read through to get a full feel of Nextflow's capabilities. An expanded tutorial, to implement a multi-task workflow, analogous to that provided by Snakemake, would have been useful.

TODO

### Sample-specific sub-directories

Every invocation of a task - every process - has its own subdirectory within Nextflow's `work/` directory into which the input files for that task are copied and into which the output files from that task are written. As a result, sample-specific invocations of tasks will have their own directories.

TODO is there a way to easily identify which directories hold the files for which samples?

### Parse YAML configuration files

Nextflow supports a `-params-file` command-line option which allows a YAML or JSON file to be provided. The parameters can then be accessed within a Nextflow script via a `params.` prefix e.g. `params.rrna_fasta_file`.

### Dry run option, validating configuration and input file existence

TODO

### Tool/step-specific log files

Standard output and standard error (as well as exit codes) for each task invocation are automatically captured and placed in process-specific subdirectories of Nextflow's `work/` directory.

For example:

```console
$ nextflow run riboviz.nf -params-file vignette/vignette_config.yaml
N E X T F L O W  ~  version 20.01.0
Launching `riboviz.nf` [special_goodall] - revision: 78a68e034f
executor >  local (2)
[b3/92880f] process > buildIndicesRrna [  0%] 0 of 1
[52/edae02] process > buildIndicesOrf  [  0%] 0 of 1
...
[b3/92880f] process > buildIndicesRrna [100%] 1 of 1 ?
[52/edae02] process > buildIndicesOrf  [100%] 1 of 1 ?
$ ls -1A work/b3/92880fd4c2dffe97e20b16f363182a/
.command.begin
.command.err
.command.log
.command.out
.command.run
.command.sh
.exitcode
yeast_rRNA.1.ht2
yeast_rRNA.2.ht2
yeast_rRNA.3.ht2
yeast_rRNA.4.ht2
yeast_rRNA.5.ht2
yeast_rRNA.6.ht2
yeast_rRNA.7.ht2
yeast_rRNA.8.ht2
yeast_rRNA_R64-1-1.fa
```

TODO

### Output bash script, that can be rerun

TODO

### Access configuration file name from within a Nextflow script

There doesn't seem to be a way to access the configuration file itself from within a Nextflow script i.e. the value of `-params-file`.

`riboviz.tools.count_reads` needs the RiboViz configuration file for sample ID-file name values in `fq_files`. These could be written to another file from within the Nextflow script, or via another script invoked by Nextflow, and that file provided as the input to `riboviz.tools.count_reads`.

TODO: The current Snakefile dumps `fq_files` into an output `sample_sheet.yaml` file that is then provided as input to `riboviz.tools.count_reads`.

### If processing of one sample fails will the rest be processed?

TODO

### Conditional behaviour

Nextflow supports a `when` declaration in tasks to allow tasks to only be executed if certain conditions hold.

### Other

Nextflow requires Java JDK 1.8+ to run.

Writing scripts requires knowledge of Groovy - a Python-esque language based on Java. I think anyone familiar with Python or R would not find Groovy too challenging. I hadn't used it myself before looking at Nextflow and though I have a number of years of experience with Java, I didn't really draw upon that knowledge.

Nextflow has built-in functions for FASTA and FASTQ files e.g. count records, extract records, query the NCBI SRA database for specific FASTQ files.

---

## TODO

* Ease of implementation of key RiboViz steps, for example:
  - index => [cutadapt => align (rRNA)]* => count_reads
  - index => cutadapt => demultiplex => [align (rRNA)]* => count_reads
  - Requires:
    - Iteration over samples.
    - Aggregation of sample-specific results.

* More readable directory names? For samples? for log files?
* Is there a way to easily identify which directories hold the files for which samples?

* Best practices on naming?
* Dry run option, validating configuration and input file existence.
* Output bash script, that can be rerun?

* Incremental build?

```
$ nextflow run riboviz.nf -params-file vignette/vignette_config.yaml -resume
```

* If processing of one sample fails can/will the rest be processed?

* If any criteria above are not met, then could it be added easily. It is OK to do development to extend a tool if necessary.

Read...

Julian Mazzitelli, [NGS Workflows](https://jmazz.me/blog/NGS-Workflows), @thejmazz, 8 June 2016

Brian Naughton, [Comparing bioinformatics workflow systems: nextflow, snakemake, reflow](http://blog.booleanbiotech.com/nextflow-snakemake-reflow.html), Boolean Biotech, 2 June 2019

Examples...

https://github.com/nf-core/rnaseq

RNA-Seq pipeline:

> The example below shows how put together a RNAseq pipeline with basic functionality. It maps a collection of read-pairs to a given reference genome and outputs the respective transcript model.

https://www.nextflow.io/example4.html

* https://www.nextflow.io/docs/latest/example.html

* FAQ
* How do I process multiple input files in parallel?
* How do I get a unique ID based on the file name?
* How do I use the same channel multiple times?
* How do I invoke custom scripts and tools?
* How do I iterate over a process n times?
* How do I iterate over nth files from within a process?
