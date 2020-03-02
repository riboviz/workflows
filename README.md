# RiboViz workflow options

Relates to issue Assess workflow language options ([riboviz#48](https://github.com/riboviz/RiboViz/issues/48)).

---

## Content

* Notes on Workflow Options 28/02/2020 [Markdown](./workflows-202002.md)
* RiboViz and Common Workflow Language (CWL) (30/08/2018) [Markdown](./SsiRiboVizCwl-201808.md) (from Software Sustainability Institute review).

---

## Workflow shortlist

"Overchoice" is an issue, so focus on the following that have been mentioned in meetings and seem most prevalent and popular:

* SnakeMake
* NextFlow
* CWL with cwltool
* Others that support CWL:
  - Pegasus
  - Galaxy
  - Cromwell
  - Toil
  - (use both their own formats and reuse CWL from cwltool)

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
