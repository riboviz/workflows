# RiboViz workflow options

Mike Jackson, EPCC, The University of Edinburgh, March 2020

## About

Options for reimplementing the [RiboViz](https://github.com/riboviz/riboviz) workflow, `riboviz.tools.prep_riboviz`, using workflow technologies popular within the bioinformatics community.

Relates to issue Assess workflow language options ([riboviz#48](https://github.com/riboviz/riboviz/issues/48)).

## Workflow options

* [Notes on Workflow Options 28/02/2020](./workflows.md): survey on available workflow options.
* [Docker](./Docker.md): notes on Docker as some workflow technogies use Docker.

"Overchoice" is an issue, so focus is on the following that have been mentioned in meetings and seem most prevalent and popular: Snakemake, CWL and Nextflow.

## Assessment criteria

Ease of download, install, tutorials, initial use.

Ease of implementation of key RiboViz steps, for example:

1. Read configuration information from YAML configuration file.
2. Build hisat2 indices if requested.
3. Process each sample ID-sample file pair in turn:
   1. Cut out sequencing library adapters using `cutadapt`.
   2. Remove rRNA or other contaminating reads by alignment to rRNA index files using `hisat2`.
   3. Align remaining reads to ORFs index files using `hisat2`.
   4. Trim 5' mismatches from reads and remove reads with more than 2 mismatches using `riboviz.tools.trim_5p_mismatch`.

Other necessary and useful features:

* Iteration over samples, and processing remaining samples if processing of one sample fails.
* Aggregation of sample-specific results.
* Conditional behaviour e.g. indexing, UMI groups.
* Tool/step-specific log files.
* Accept YAML configuration files.
* Dry run option, validating configuration and input file existence.
* Output a bash script, that can be rerun.
* If processing of one sample fails will the rest be processed?

## Snakemake

* [RiboViz and Snakemake](./snakemake/README.md): Discussion on Snakemake for RiboViz and how to run an example.
* [Snakefile](./snakemake/Snakefile): Example Snakefile for RiboViz.
* [Snakemake](./snakemake/Snakemake.md): Snakemake notes.

## Common Workflow Language (CWL)

* [RiboViz and CWL](./cwl/README.md): Discussion on CWL for RiboViz and how to run an example.
* [Common Workflow Language (CWL)](./cwl/CommonWorkflowLanguage.md): CWL notes.
* [Common Workflow Language User Guide Notes](./cwl/CwlUserGuideNotes.md) from working through the examples in the CWL user guide.
* [Snakemake and CWL](./cwl/SnakemakeCwl.md): Exploring Snakemake's `--export-cwl` option.
* [cutadapt.cwl](./cwl/cutadapt.cwl): tool wrapper for `cutadapt`.
* [cutadapt-job.yml](./cwl/cutadapt-job.yml): job configuration for above.
* [riboviz-workflow.cwl](./cwl/riboviz-workflow.cwl): workflow invoking `hisat2-build` and `cutadapt`.
* [riboviz-job.yml](./cwl/riboviz-job.yml): job configuration for the above.

## Nextflow

* [RiboViz and Nextflow](./nextflow/README.md): Discussion on Nextflow for RiboViz and how to run an example.
* [Nextflow Documentation Notes](./nextflow/NextflowDocNotes.md) from reading through the Nextflow documentation.
* [riboviz.nf](./nextflow/riboviz.nf): Example Nextflow script for RiboViz.

## Summary of observations and recommendations

Snakemake:

* Straightforward to learn.
* Bulk of RiboViz workflow was prototyped in Snakemake in a day.
* Prototype would take 2-3d to complete as an implementation of the RiboViz workflow. Adding support for `riboviz.tools.demultiplex_fastq` would be the most significant addition.
* Supports Docker, Singularity, cluster and cloud execution.

Snakemake and CWL:

* A Snakemake "hello world" example in a *Nature* article was exported to CWL, using the command in the article, but this could not be run using cwltool or Toil. It is unclear as to how this can be made to run.
* A CWL file exported from Snakemake is a simple workflow that invokes Snakemake itself, rather than being the CWL equivalent of the workflow steps implemented within a Snakefile.

CWL:

* More effort to learn than Snakemake.
* 3 steps of the RiboViz workflow were prototyped in a day.
* "edit-compile-run" development cycle is very slow.
* Prototype would take 2-3 weeks to implement.
* No support for conditional execution of steps as yet.
* Developers would need to know, or learn, JavaScript.
* Tool wrappers can be verbose for what are single command-line invocations e.g. the basic `cutadapt` wrapper in `cutadapt.cwl` is 47 lines.
* Workflows can be very verbose e.g. the 3 step CWL workflow is 58 lines, whereas the 15 step Snakemake workflow is 220 lines (which are elegant in their conciseness).
* Don't feel there is any point pursuing this option further at this time. Rather than try and guess whether the user community might want CWL, ask the user community about their views on CWL.

Nextflow:

* More effort to learn than Snakemake but less effort than CWL.
* 3 steps of the RiboViz workflow were prototyped in 2 days.
* "edit-compile-run" development cycle was initially slow, partly due to reference-style nature of documentation, but accelerated.
* It would take < 2 weeks to complete an implementation of RiboViz in Nextflow.
* Developers would need to know, or learn, Groovy. As Groovy is akin to Python I don't think this is a concern.
* Nextflow scripts are a similar structure and verbosity to Snakefiles.
* Every invocation of a task takes place in its own isolated directory with symbolic links to the required input files and a small bash script with the command that is to be invoked. This would be very useful for debugging and feels elegant.
* Supports Docker, Singularity, cluster and cloud execution.

General:

* Our Python workflow has custom error handling. Users using Snakemake, CWL or Nextflow will need to rely on their error reporting which may be more verbose, or, in the case of Toil, a possibly impenetrable stack trace. We can extend our documentation with the symptoms of common errors (e.g. missing files, missing configuration, step fails) and how to resolve these.

## Recommendation

I'd recommend adopting Nextflow over Snakemake for the following reasons:

* While it would require more time to complete a port of the RiboViz workflow to Nextflow than to Snakemake, the time required is not significant (i.e. < 2 weeks).
* Nextflow feels far richer in terms of its features and expressivity.
* The execution of each task in isolated directories is very useful for debugging.
* Nextflow's built-in support for, and documentation around, Docker, Singularity, cluster and cloud execution, seems more comprehensive than that of Snakemake.
* While Nextflow and Snakemake are both very popular tools within the bioinformatics community, Nextflow seems to have the slight edge (based on [Notes on Workflow Options 28/02/2020](./workflows.md)).
