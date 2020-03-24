# RiboViz workflow options

## About

Options for reimplementing the RiboViz workflow, `riboviz.tools.prep_riboviz`, using workflow technologies popular within the bioinformatics community.

Relates to issue Assess workflow language options ([riboviz#48](https://github.com/riboviz/riboviz/issues/48)).

---

## Workflow options

General:

* [Notes on Workflow Options 28/02/2020](./workflows.md): survey on available workflow options.
* [Docker](./Docker.md): notes on Docker as some workflow technogies use Docker.

"Overchoice" is an issue, so focus is on the following that have been mentioned in meetings and seem most prevalent and popular: Snakemake, CWL and Nextflow.

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

* [RiboViz and Nextflow](./nextflow/README.md): Discussion on Nextflow for RiboViz and how to run an example.
* [Nextflow Documentation Notes](./NextflowDocNotes.md) from reading through the Nextflow documentation.

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

## Summary of observations/recommendations to date (12/03/2020)

Snakemake:

* Straightforward to learn.
* Bulk of RiboViz workflow was prototyped in Snakemake in a day.
* Prototype would take 2-3d to complete as an implementation of the RiboViz workflow. Adding support for `riboviz.tools.demultiplex_fastq` would be the most significant addition.

Snakemake and CWL:

* A Snakemake "hello world" example in a *Nature* article was exported to CWL, using the command in the article, but this could not be run using cwltool or Toil. It is unclear as to how this can be made to run.
* A CWL file exported from Snakemake is a simple workflow that invokes Snakemake itself, rather than being the CWL equivalent of the workflow steps implemented within a Snakefile.

CWL:

* More effort to learn.
* 3 steps of the RiboViz workflow were prototyped in a day.
* Prototype would take 2-3 weeks to complete.
* No support for conditional execution of steps as yet.
* "edit-compile-run" development cycle iss very slow.
* Developers would need to know JavaScript.
* Tool wrappers can be verbose for what are single command-line invocations e.g. the basic `cutadapt` wrapper in `cutadapt.cwl` is 47 lines.
* Workflows can be verbose e.g. the 3 step CWL workflow is 58 lines, whereas the 15 step Snakemake workflow is 220 lines (which are elegant in their conciseness).

General:

* Our Python workflow has custom error handling. Users using Snakemake and CWL will need to rely on their error reporting which may be more verbose, or, in the case of Toil, a possibly impenetrable stack trace. We can extend our documentation with the symptoms of common errors (e.g. missing files, missing configuration, step fails) and how to resolve these.

Recommendation to date:

* Adopt Snakemake - the work requires little effort to complete.
* Snakemake is a popular tool within the community.
* Rather than try and guess what the user community might want, ask the user community about their views on CWL.
