# RiboViz workflow options

Relates to issue Assess workflow language options ([riboviz#48](https://github.com/riboviz/riboviz/issues/48)).

## Workflow options

[Notes on Workflow Options 28/02/2020](./workflows.md) - survey on available workflow options.

"Overchoice" is an issue, so focus is on the following that have been mentioned in meetings and seem most prevalent and popular.

Snakemake:

* [RiboViz and Snakemake](./snakemake/README.md): Discussion on Snakemake for RiboViz and how to run an example.
* [Snakefile](./snakemake/Snakefile): Example Snakefile for RiboViz.
* [workflow.svg](./snakemake/workflow.svg): Example workflow from above.
* [report.html](./snakemake/report.html): Example report from above.
* [Snakemake](./snakemake/Snakemake.md): Snakemake notes.

Common Workflow Language (CWL):

* [RiboViz and CWL](./cwl/README.md): Discussion on CWL for RiboViz and how to run an example.
* [Common Workflow Language (CWL)](./cwl/CommonWorkflowLanguage.md): CWL notes.
* [Common Workflow Language User Guide Notes](./cwl/CwlUserGuideNotes.md) from working through the examples in the CWL user guide.
* [Snakemake and CWL](./cwl/SnakemakeCwl.md): Exploring Snakemake's `--export-cwl` option.
* [cutadapt.cwl](./cwl/cutadapt.cwl): tool wrapper for `cutadapt`.
* [cutadapt-job.yml](./cwl/cutadapt-job.yml): job configuration for above.
* [riboviz-workflow.cwl](./cwl/riboviz-workflow.cwl): workflow invoking `hisat2-build` and `cutadapt`.
* [riboviz-job.yml](./cwl/riboviz-job.yml): job configuration for the above.

Nextflow:

* TODO

---

## Assessment criteria

* Ease of download, install, tutorials, initial use.
* Ease of implementation of key RiboViz steps, for example:
  - index => [cutadapt => align (rRNA)]* => count_reads
  - index => cutadapt => demultiplex => [align (rRNA)]* => count_reads
  - Requires:
    - Ability to handle implicit naming of index files.
    - Iteration over samples.
    - Aggregation of sample-specific results.
    - Conditional behaviour e.g. indexing, UMI groups.
* Tool/step-specific log files.
* Parse YAML configuration files.
* Dry run option, validating configuration and input file existence.
* Output bash script, that can be rerun.
* If processing of one sample fails will the rest be processed?

If any criteria above are not met, then could it be added easily. It is OK to do development to extend a tool if necessary.

---

## General

[Docker](./Docker.md)
