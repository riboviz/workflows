# RiboViz and Nextflow

* [Nextflow](https://www.nextflow.io/)
* [GitHub](https://github.com/nextflow-io/nextflow)
* [Documentation](https://www.nextflow.io/docs/latest/index.html)

See notes on:

* [Nextflow Documentation Notes](./NextflowUserGuideNotes.md) from reading through the Nextflow documentation.

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

TODO

---

## Assessment

### Ease of download, install, tutorials, initial use

### Ease of implementation of key RiboViz steps

A lot of preparatory reading is required to go beyond the "hello world!" example. The User Guide took about 4 hours to work through, running all the examples.

### Sample-specific sub-directories

### Parse YAML configuration files

### Dry run option, validating configuration and input file existence

### Tool/step-specific log files

### Output bash script, that can be rerun

### Access configuration file name from within CWL

### If processing of one sample fails will the rest be processed?

### Conditional behaviour

### Other

---

## TODO

Julian Mazzitelli, [NGS Workflows](https://jmazz.me/blog/NGS-Workflows), @thejmazz, 8 June 2016

Brian Naughton, [Comparing bioinformatics workflow systems: nextflow, snakemake, reflow](http://blog.booleanbiotech.com/nextflow-snakemake-reflow.html, Boolean Biotech, 2 June 2019

* Ease of implementation of key RiboViz steps, for example:
  - index => [cutadapt => align (rRNA)]* => count_reads
  - index => cutadapt => demultiplex => [align (rRNA)]* => count_reads
  - Requires:
    - Ability to handle implicit naming of index files.
    - Iteration over samples.
    - Aggregation of sample-specific results.
    - Conditional behaviour e.g. indexing, UMI groups.

* Tool/step-specific log files - try Snakemake style.
* Parse YAML configuration files.
* Dry run option, validating configuration and input file existence.
* Output bash script, that can be rerun.
* If processing of one sample fails will the rest be processed?
* If any criteria above are not met, then could it be added easily. It is OK to do development to extend a tool if necessary.

Remember docs have:

* Examples
* FAQ
* How do I process multiple input files in parallel?
* How do I get a unique ID based on the file name?
* How do I use the same channel multiple times?
* How do I invoke custom scripts and tools?
* How do I iterate over a process n times?
* How do I iterate over nth files from within a process?

Notes:

* Tutorial consists of a simple two process (step) workflow. The rest of docs are more of a reference guide which took a while to read through to get a feel of the capabilities. A more complex example, as provided by Snakemake, would be nice.
* Requires Java JDK 1.8+ to run.
* Requires knowledge of Groovy (Pythonesque language underpinned by Java) for scripting.
* Built-in functions to count FASTQ, FASTA records, to query NCBI SRA database and emit FASTQ files matching specified criteria.
* Similar data flow model to Snakemake.
* Feeling that features are far richer.
* HTML reports and workflow images.
