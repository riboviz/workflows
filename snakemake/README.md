# RiboViz and Snakemake

* [GitHub](https://github.com/snakemake/snakemake)
* [Documentation](https://snakemake.readthedocs.io/en/stable/)

See notes on [Snakemake](./Snakemake.md).

---

## Run the RiboViz example

[Snakefile](./Snakefile) contains a Snakemake Snakefile with steps mimicing those in `riboviz.tools.prep_riboviz` to run the RiboViz vignette (i.e. no UMI extraction, deduplication or demultiplexing).

Create a new conda environment from the current RiboViz one and activate it:

```console
$ conda create --name riboviz-snakemake --clone base
$ conda activate riboviz-snakemake
```

Install Snakemake (and pygraphviz and networkx for HTML report generation):

```console
$ conda install -c bioconda -c conda-forge snakemake-minimal
$ conda install -c bioconda -c conda-forge pygraphviz
$ conda install -c bioconda -c conda-forge networkx
```

(opt for minimal install to keep dependencies down)

Copy over Snakefile:

```console
$ cd riboviz
$ cp ~/workflows/snakemake/Snakefile .
```

Dry run Snakemake with `vignette/vignette_config.yaml`:

```console
$ snakemake --configfile vignette/vignette_config.yaml -n
```

Run Snakemake:

```console
$ snakemake --configfile vignette/vignette_config.yaml
```

Create SVG of workflow:

```console
$ snakemake --configfile vignette/vignette_config.yaml --dag | dot -Tsvg > workflow.svg
```

See, for example, [workflow.svg](./workflow.svg).

Create an HTML report with a graphical representation of the workflow:

```console
$ snakemake --configfile vignette/vignette_config.yaml --report report.html
```

See, for example, [report.html](./report.html).

Remove RiboViz index, temporary, output and log directories:

```console
$ snakemake clean --configfile vignette/vignette_config.yaml
```

---

## Assessment

### Ease of download, install, tutorials, initial use.

* Pure Python package.
* Very easy to download, install and complete tutorial.

### Ease of implementation of key RiboViz steps

The Snakefile was implemented in a day.

### Sample-specific sub-directories

This is supported.

### Parse YAML configuration files

https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#standard-configuration

Hard-code configuration file:

```
configfile: "vignette/vignette_config.yaml"

rule cut_adapters:
    input:
        get_sample_file
    output:
        os.path.join(config['dir_tmp'], "{sample}", "trim.fq")
    log: os.path.join(config['dir_logs'], "{sample}", "cutadapt.log")
    shell:
        "cutadapt --trim-n -O 1 -m 5 -a {config[adapters]} -o {output} {input}  -j 0 >> {log}"
```

Specify/override via command-line:

```console
$ snakemake --configfile vignette/vignette_config.yaml ...
```

YAML can be validated against a custom schema, see https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation.

### Dry run option, validating configuration and input file existence

Use dry run, `-n`, flag. For example:

For example:

```console
$ mv vignette/input/yeast_rRNA_R64-1-1.fa .
$ snakemake vignette/index/yeast_rRNA.1.ht2 --configfile vignette/vignette_config.yaml -n
Building DAG of jobs...
MissingInputException in line 20 of /home/ubuntu/riboviz/Snakefile:
Missing input files for rule build_indices_rrna:
vignette/input/yeast_rRNA_R64-1-1.fa

$ snakemake vignette/index/YAL_CDS_w_250.1.ht2 --configfile vignette/vignette_config.yaml -n
Building DAG of jobs...
Job counts:
	count	jobs
	1	build_indices_orf
	1

[Tue Mar  3 02:05:03 2020]
rule build_indices_orf:
    input: vignette/input/yeast_YAL_CDS_w_250utrs.fa
    output: vignette/index/YAL_CDS_w_250.1.ht2, vignette/index/YAL_CDS_w_250.2.ht2, vignette/index/YAL_CDS_w_250.3.ht2, vignette/index/YAL_CDS_w_250.4.ht2, vignette/index/YAL_CDS_w_250.5.ht2, vignette/index/YAL_CDS_w_250.6.ht2, vignette/index/YAL_CDS_w_250.7.ht2, vignette/index/YAL_CDS_w_250.8.ht2
    log: vignette/logs/20200303-020503/hisat2_build_orf.log
    jobid: 0

Job counts:
	count	jobs
	1	build_indices_orf
	1
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

### Tool/step-specific log files

https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files

Example:

```
rule build_indices_rrna:
    input:
        config['rrna_fasta_file']
    output:
        directory(os.path.join(config['dir_index'], "rrna"))
    params:
        prefix=os.path.join(config['dir_index'], "rrna", config['rrna_index_prefix'])
    log: os.path.join(config['dir_logs'], "hisat2_build_r_rna.log")
    shell:
       "mkdir {output}; hisat2-build --version &> {log}; hisat2-build {input} {params.prefix} &>> {log}"
```

`&>` is a bash 4.0+ operator to capture both standard output and standard error.

Log files are not deleted upon error.

### Output bash script, that can be rerun

`--detailed-summary` produces a TSV file including bash commands run:

```console
$ snakemake --configfile vignette/vignette_config.yaml --detailed-summary
Building DAG of jobs...
output_file	date	rule	version	log-file(s)	input-file(s)	shellcmdstatus	plan
vignette/output/TPMs_collated.tsv	Mon Mar  2 06:02:52 2020	collate_tpms	-	vignette/logs/collate_tpms.log	vignette/output/WT3AT/tpms.tsv,vignette/output/WTnone/tpms.tsv	Rscript --vanilla rscripts/collate_tpms.R --sample-subdirs=True --output-dir=vignette/output WTnone WT3AT >> vignette/logs/collate_tpms.log	set of input files changed	no update
vignette/output/read_counts.tsv	Mon Mar  2 06:08:06 2020	count_reads	-vignette/logs/count_reads.log		python -m riboviz.tools.count_reads -c vignette/vignette_config.yaml -i vignette/input -t vignette/tmp -o vignette/output -r vignette/output/read_counts.tsv >> vignette/logs/count_reads.log	ok	no update
vignette/output/WTnone/3nt_periodicity.tsv	Mon Mar  2 06:01:52 2020	generate_stats_figs	-	vignette/logs/WTnone/generate_stats_figs.log	data/yeast_codon_pos_i200.RData,data/yeast_features.tsv,data/yeast_standard_asite_disp_length.txt,data/yeast_tRNAs.tsv,vignette/input/yeast_YAL_CDS_w_250utrs.fa,vignette/input/yeast_YAL_CDS_w_250utrs.gff3,vignette/output/WTnone/WTnone.h5	Rscript --vanilla rscripts/generate_stats_figs.R --num-processes=1 --min-read-length=10 --max-read-length=50 --buffer=250 --primary-id=Name --dataset=vignette --hd-file=vignette/output/WTnone/WTnone.h5 --orf-fasta-file=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True --output-dir=vignette/output/WTnone --do-pos-sp-nt-freq=True --t-rna-file=data/yeast_tRNAs.tsv --codon-positions-file=data/yeast_codon_pos_i200.RData --features-file=data/yeast_features.tsv --orf-gff-file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --asite-disp-length-file=data/yeast_standard_asite_disp_length.txt --count-threshold=64 > vignette/logs/WTnone/generate_stats_figs.log	ok	no update
vignette/output/WTnone/3nt_periodicity.pdf	Mon Mar  2 06:01:52 2020	generate_stats_figs	-	vignette/logs/WTnone/generate_stats_figs.log	data/yeast_codon_pos_i200.RData,data/yeast_features.tsv,data/yeast_standard_asite_disp_length.txt,data/yeast_tRNAs.tsv,vignette/input/yeast_YAL_CDS_w_250utrs.fa,vignette/input/yeast_YAL_CDS_w_250utrs.gff3,vignette/output/WTnone/WTnone.h5	Rscript --vanilla rscripts/generate_stats_figs.R --num-processes=1 --min-read-length=10 --max-read-length=50 --buffer=250 --primary-id=Name --dataset=vignette --hd-file=vignette/output/WTnone/WTnone.h5 --orf-fasta-file=vignette/input/yeast_YAL_CDS_w_250utrs.fa --rpf=True --output-dir=vignette/output/WTnone --do-pos-sp-nt-freq=True --t-rna-file=data/yeast_tRNAs.tsv --codon-positions-file=data/yeast_codon_pos_i200.RData --features-file=data/yeast_features.tsv --orf-gff-file=vignette/input/yeast_YAL_CDS_w_250utrs.gff3 --asite-disp-length-file=data/yeast_standard_asite_disp_length.txt --count-threshold=64 > vignette/logs/WTnone/generate_stats_figs.log	ok	no update
...
```

To create a bash script would require parsing column `shellcmd` column and executing the commands (in reverse order). Directory creation commands would also need to be derived.

```console
$ snakemake --configfile vignette/vignette_config.yaml --printshellcmds -n
Building DAG of jobs...
Job counts:
	count	jobs
	1	all
	2	bam_to_h5
	1	build_indices_orf
	1	build_indices_rrna
	1	collate_tpms
	1	count_reads
	2	cut_adapters
	2	generate_stats_figs
	2	index_bam
	2	map_to_orf
	2	map_to_r_rna
	2	sort_bam
	2	trim_5p_mismatches
	21

[Mon Mar  2 06:16:48 2020]
rule cut_adapters:
    input: vignette/input/SRR1042855_s1mi.fastq.gz
    output: vignette/tmp/WTnone/trim.fq
    log: vignette/logs/WTnone/cutadapt.log
    jobid: 18
    wildcards: sample=WTnone

cutadapt --trim-n -O 1 -m 5 -a	CTGTAGGCACC -o vignette/tmp/WTnone/trim.fq vignette/input/SRR1042855_s1mi.fastq.gz -j 0 >> vignette/logs/WTnone/cutadapt.log

[Mon Mar  2 06:16:48 2020]
rule build_indices_orf:
    input: vignette/input/yeast_YAL_CDS_w_250utrs.fa
    output: vignette/index/orf
    log: vignette/logs/hisat2_build_orf.log
    jobid: 16

mkdir vignette/index/orf; hisat2-build --version > vignette/logs/hisat2_build_orf.log; hisat2-build vignette/input/yeast_YAL_CDS_w_250utrs.fa vignette/index/orf/YAL_CDS_w_250 >> vignette/logs/hisat2_build_orf.log

[Mon Mar  2 06:16:48 2020]
rule cut_adapters:
    input: vignette/input/SRR1042864_s1mi.fastq.gz
    output: vignette/tmp/WT3AT/trim.fq
    log: vignette/logs/WT3AT/cutadapt.log
    jobid: 20
    wildcards: sample=WT3AT

cutadapt --trim-n -O 1 -m 5 -a	CTGTAGGCACC -o vignette/tmp/WT3AT/trim.fq vignette/input/SRR1042864_s1mi.fastq.gz -j 0 >> vignette/logs/WT3AT/cutadapt.log

[Mon Mar  2 06:16:48 2020]
rule count_reads:
    output: vignette/output/read_counts.tsv
    log: vignette/logs/count_reads.log
    jobid: 2

python -m riboviz.tools.count_reads -c vignette/vignette_config.yaml -i vignette/input -t vignette/tmp -o vignette/output -r vignette/output/read_counts.tsv >> vignette/logs/count_reads.log

[Mon Mar  2 06:16:48 2020]
rule build_indices_rrna:
    input: vignette/input/yeast_rRNA_R64-1-1.fa
    output: vignette/index/rrna
    log: vignette/logs/hisat2_build_r_rna.log
    jobid: 19

mkdir vignette/index/rrna; hisat2-build --version > vignette/logs/hisat2_build_r_rna.log; hisat2-build vignette/input/yeast_rRNA_R64-1-1.fa vignette/index/rrna/yeast_rRNA >> vignette/logs/hisat2_build_r_rna.log
...
```

To create a bash script would require parsing this (semi-structured) output and executing the commands (in reverse order). Directory creation commands would also need to be derived.

Another approach would be to explicitly have commands within Snakemake rules to export the shell commands to a specific file, but this would lead to a lot of duplicated code. If one has the Snakemake file, why bother with a bash file?

### Access configuration file from within Snakefile

There doesn't seem to be a way to access the configuration file itself from within a Snakefile i.e. the value of `--configfile`:

```console
$ snakemake --configfile vignette/vignette_config.yaml
```

For `riboviz.tools.count_reads`, which needs the configuration file, a workaround is to provide this as an additional configuration value e.g.

```console
$ snakemake --configfile vignette/vignette_config.yaml --config config_file=vignette/vignette_config.yaml
```

and to reference it as follows:

```
        "python -m riboviz.tools.count_reads -c {config[config_file]} -i {config[dir_in]} -t {config[dir_tmp]} -o {config[dir_out]} -r {output} >> {log}"
```

However, this feels hacky. Also, a comment in [Retrieve value of --configfile parameter in rule params](https://bitbucket.org/snakemake/snakemake/issues/594/retrieve-value-of-configfile-parameter-in) (link now broken) in the Snakemake issue tracker recommends not accessing the configuration file within a Snakefile:

> The reason for not documenting is that --configfile is not an obligatory parameter, and hence it does not make too much sense to access it from the Snakemake. E.g., the user could omit it, and in fact it is recommended to store workflows in a way such that no additional arguments are needed (see here), in order to maximize reproducibility. Further, the user could provide --configfile and further overwrite values via --config.

`riboviz.tools.count_reads` only uses the configuration file for sample ID-file name values in `fq_files`. These could be written to another file from within the Snakefile, or another script, and that file provided as the input to `riboviz.tools.count_reads`. The current Snakefile dumps `fq_files` into an output `sample_sheet.yaml` file that is then provided as input to `riboviz.tools.count_reads`.

### If processing of one sample fails will the rest be processed?

If Snakemake is run with the vignette configuration, it fails:

```console
$ snakemake  --configfile vignette/vignette_config.yaml 
Building DAG of jobs...
MissingInputException in line 45 of /home/ubuntu/riboviz/Snakefile:
Missing input files for rule cut_adapters:
vignette/input/example_missing_file.fastq.gz
```

The `-k` / `--keep-going` flag allows Snakemake to run independent jobs (rules to generate output files) if a specific job fails. However, in this scenario, the workflow still fails as above.

To avert this, a `SAMPLES` variable is created from `fq_files` with only those sample name-file pairs for which a sample file exists:

```
SAMPLES = {tag: filename for (tag, filename) in config['fq_files'].items()
        if os.path.exists(os.path.join(config['dir_in'], filename))
        and os.path.isfile(os.path.join(config['dir_in'], filename))}
```

If a sample file exists but causes problems (e.g. is an invalid format) then the workflow will fail. For example, if there is an `fq_files` entry: `Invalid: invalid.fastq.gz` and `invalid.fastq.gz` has content `XXX`, the workflow will fail:

```console
$ snakemake --configfile vignette/vignette_config.yaml
...
[Tue Mar  3 04:12:02 2020]
rule cut_adapters:
    input: vignette/input/invalid.fastq.gz
    output: vignette/tmp/Invalid/trim.fq
    log: vignette/logs/20200303-041202/Invalid/cutadapt.log
    jobid: 35
    wildcards: sample=Invalid

[Tue Mar  3 04:12:02 2020]
Error in rule cut_adapters:
    jobid: 35
    output: vignette/tmp/Invalid/trim.fq
    log: vignette/logs/20200303-041202/Invalid/cutadapt.log (check log file(s) for error message)
    shell:
        cutadapt --trim-n -O 1 -m 5 -a	CTGTAGGCACC -o vignette/tmp/Invalid/trim.fq vignette/input/invalid.fastq.gz -j 0 &> vignette/logs/20200303-041202/Invalid/cutadapt.log
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
...
```

Using the `-k` option allows the workflow to continue and process the other samples. However, it will still fail before invoking the rule for `collate_tpms.R` as this has a dependency on the `tpms.tsv` files for all samples, which, following dependencies back, has a dependency on the invalid input file.

I tried to implement workarounds (e.g. checking for existing `tpms.tsv` files in the rule for `collate_tpms.R`) but these felt hacky and undermined the way in which Snakemake builds the workflow based on the available files.

In such cases a user _can_ edit their configuration file to remove the `fq_files` entry for the problematic sample, and rerun Snakemake. As incremental build is supported, only the rule for `collate_tpms.R` will be run.

### Conditional behaviour

For example, for selecting to create index file, bedgraphs, UMI group summary files, to deduplicate or count reads, a response to [Add complex conditional file dependency](https://bitbucket.org/snakemake/snakemake/issues/37/add-complex-conditional-file-dependency) (link now broken) comments:

> You can put rules in conditional statements. These are evaluated before the workflow is executed though, so they cannot be data dependent.

For example:

```
if config['count_reads']:
    rule count_reads:
        input:
            os.path.join(config['dir_out'], "sample_sheet.yaml")
	output:
            os.path.join(config['dir_out'], "read_counts.tsv")
        log: os.path.join(config['dir_logs'], "count_reads.log")
        shell:
            "python -m riboviz.tools.count_reads -c {input} -i {config[dir_in]} -t {config[dir_tmp]} -o {config[dir_out]} -r {output} >> {log}"
```

The syntax for specifying optional files as inputs look a bit contrived:

```
rule all:
    input:
        [os.path.join(config['dir_out'], sample, "plus.bedgraph") for sample in SAMPLES] if config['make_bedgraph'] else [],
        [os.path.join(config['dir_out'], sample, "minus.bedgraph") for sample in SAMPLES] if config['make_bedgraph'] else [],
	os.path.join(config['dir_out'], "TPMs_collated.tsv"),
        [os.path.join(config['dir_out'], "read_counts.tsv")] if config['count_reads'] else []
```

### Common Workflow Language

https://snakemake.readthedocs.io/en/stable/executing/interoperability.html

Snakemake has a `-export-cwl` flag which exports a Common Workflow Language document:

```console
$ snakemake --configfile=vignette/vignette_config.yaml --export-cwl snakemake.cwl
```

This is very verbose. For example it includes the following entries for _every_ file in RiboViz whether or not these are used in the workflow e.g.:

```
{
    "class": "File",
    "location": "website/Data/F2_2014_Gerashchenko_complete_SD_media.tsv"
},
```

I'll try running this file when looking at CWL execution environments.

### Containers, clusters and cloud

Snakemake has support for:

* [Running jobs in containers](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#running-jobs-in-containers) - Docker and Singularity.
* [Cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) via cluster engines that support shell scripts and use a common file system. An example is provided for using Sun Grid Engine to submit a job.
* [Cloud support](https://snakemake.readthedocs.io/en/stable/executing/cloud.html). Examples are provided for:
  - Google Cloud Engine via [Kubernetes](https://kubernetes.io/) "open-source system for automating deployment, scaling, and management of containerized applications."
  - Amazon Web Services via [Tibanna](https://tibanna.readthedocs.io/en/latest/) tool to run genomic pipelines.
  - Input and output files are staged from/to remote locations and the Snakemake workflow is expected to be stored within a Git repository. There is no requirement for a shared network file system.

### Other

Snakemake Wrapper Repository includes reusable wrappers for using popular bioinformatics tools within Snakemake, https://snakemake-wrappers.readthedocs.io/en/stable/index.html.

Arbitrary Python code can be embedded anywhere within the Snakefile or it can be imported.
