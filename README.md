# RiboViz workflow options

Relates to issue Assess workflow language options ([riboviz#48](https://github.com/riboviz/riboviz/issues/48)).

## Workflow options

[Notes on Workflow Options 28/02/2020](./workflows-202002.md) - survey on available workflow options.

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



Observations/recommendations to date:

* Snakemake was straightforward to learn.
* Snakemake will take 2-3d to complete, with "dedup" being the challenge.
* Nature Snakemake example, export to CWL, fails to run with cwltool and Toil, can't find file, can't discover why.
* Snakemake CWL is a single step workflow with the invocation of Snakemake.
* CWL more ramp up effort.
* CWL YAML file per tool, YAML file for workflow.
* CWL 2-3 weeks, 4-8 weeks timescale.
* Adopt Snakemake now, "in" to community tool now.
* Ask about CWL at June meeting, push the decision to the user community (YAGNI). Also run by Edward's contacts (e.g. Wellcome) to see if they can advise based on our experiences and recommendations.
* Could have both, maintainability overhead, but potentially wider user appeal.
* Python workflow has custom error handling, Snakemake/CWL runners, rely on their error reporting which may be more verbose, or, in the case of Toil, an inpenetrable stack trace. Extend our documentation with common errors (e.g. missing files, missing configuration, step fails) and how they appear.



47 cwl/cutadapt.cwl

Mike Jackson@7390MJ MINGW64 ~/EPCC/RiboViz/workflows (master)
$ wc -l cwl/cutadapt-job.yml
9 cwl/cutadapt-job.yml

Mike Jackson@7390MJ MINGW64 ~/EPCC/RiboViz/workflows (master)
$ wc -l cwl/riboviz-workflow.cwl
58 cwl/riboviz-workflow.cwl

Mike Jackson@7390MJ MINGW64 ~/EPCC/RiboViz/workflows (master)
$ wc -l snakemake/Snakefile
220 snakemake/Snakefile


Mike Jackson@7390MJ MINGW64 ~/EPCC/RiboViz/workflows (master)
$ wc -l cwl/cutadapt-job.yml
9 cwl/cutadapt-job.yml

Mike Jackson@7390MJ MINGW64 ~/EPCC/RiboViz/workflows (master)
$ wc -l cwl/riboviz-workflow.cwl
58 cwl/riboviz-workflow.cwl

Mike Jackson@7390MJ MINGW64 ~/EPCC/RiboViz/workflows (master)
$ wc -l snakemake/Snakefile
220 snakemake/Snakefile



Add to "workflows" repository:

Strozzi F. et al. (2019) Scalable Workflows and Reproducible Data Analysis for Genomics. In: Anisimova M. (eds) Evolutionary Genomics. Methods in Molecular Biology, vol 1910. Humana, New York, NY. doi: [10.1007/978-1-4939-9074-0_24](https://doi.org/10.1007/978-1-4939-9074-0_24). Discussion of Snakemake, Nextflow, CWL and GWL complemented with examples at https://github.com/EvolutionaryGenomics/scalability-reproducibility-chapter

Marijn van Vliet, Guidelines for data analysis scripts, v2, 9 Aug 2019, [arXiv:1904.06163](https://arxiv.org/abs/1904.06163) [cs.SE].

