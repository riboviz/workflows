# Options for RiboViz workflow management

Edward Wallace, Felicity Anderson: [The Wallace Lab](https://ewallace.github.io/), The University of Edinburgh.

Kostas Kavoussanakis, Mike Jackson: [EPCC](https://www.epcc.ed.ac.uk/), The University of Edinburgh.

March-July 2020.

## About

This repository contains notes related to exploring options for reimplementing the [RiboViz](https://github.com/riboviz/riboviz) workflow (implemented at that time as a custom Python script) using workflow technologies popular within the bioinformatics community.

This repository also contains Snakemake, Nextflow and CWL workflows rapidly prototyped during these explorations. These prototype workflows are unsupported and not designed for actual use, but may prove of interest to others as snapshots of how a subset of the RiboViz workflow was reimplemented using these technologies.

For our notes on options for a workflow management system for RiboViz and links to the prototype workflows within this repository, see [Exploring options for a workflow management system for RiboViz](./WorkflowOptions.md).

The `presentations/` directory has the sources for two presentations where these explorations were discussed with colleagues. The online versions of these are hosted on GitHub:


* [Choosing a workflow management system for RiboViz](https://riboviz.github.io/workflows/RiboVizWorkflowsPresentation-20200429.html), 29 April 2020. Presentation for RiboViz team plus guests. ([./presentations/RiboVizWorkflowsPresentation-20200429.html](./presentations/RiboVizWorkflowsPresentation-20200429.html))
* [Choosing a workflow management system for RiboViz](https://riboviz.github.io/workflows/RiboVizWorkflowsPresentation-20200603.html), 03 June 2020.  Presentation for EPCC seminar ([./presentations/RiboVizWorkflowsPresentation-20200603.html](./presentations/RiboVizWorkflowsPresentation-20200603.html))

## RiboViz

These notes relate to the RiboViz development ticket Assess workflow language options ([riboviz#48](https://github.com/riboviz/riboviz/issues/48)).

RiboViz release [2.0](https://github.com/riboviz/riboviz/releases/tag/2.0) of 8th July 2020 contains the first complete implementation of the RiboViz workflow in Nextflow, our chosen workflow management system, which can be compared to the Python implementation within the same release:

* Nextflow workflow ([prep_riboviz.nf](https://github.com/riboviz/riboviz/blob/2.0/prep_riboviz.nf)).
* Python workflow ([riboviz.tools.prep_riboviz](https://github.com/riboviz/riboviz/blob/2.0/riboviz/tools/prep_riboviz.py)).

## Copyright and License

The contents of this repository are Copyright (2020) The University of Edinburgh.

The contents of this repository are released under the [Apache License 2.0](./LICENSE).

## Acknowledgements

This work was supported by the following:

* [BBSRC-NSF/BIO Lead Agency](https://bbsrc.ukri.org/funding/filter/2018-nsfbio-lead-agency-scheme/) collaboration funding:
  - BBSRC grant number BB/S018506/1, "RiboViz for reliable, reproducible and rigorous quantification of protein synthesis from ribosome profiling data", 05/2019-04/2022, awarded to Edward Wallace,
  - NSF/BIO grant number 1936046, "Collaborative Research: BBSRC: RiboViz for reliable, reproducible and rigorous quantification of protein synthesis from ribosome profiling data", 07/2019-06/2022, awarded to Premal Shah and Liana Lareau.
* Edward Wallace is funded by a Sir Henry Dale Fellowship jointly funded by the Wellcome Trust and the Royal Society (Grant Number 208779/Z/17/Z).

## Citation

To cite this repository, please use:

Jackson, Michael; Wallace, Edward; Kavoussanakis, Kostas; Anderson, Felicity (2020): Options for RiboViz workflow management. figshare. Software. doi: [10.6084/m9.figshare.13147979](https://doi.org/10.6084/m9.figshare.13147979)

To cite **RiboViz**, please use both of the following references:

riboviz: analysis and visualization of ribosome profiling datasets, Carja et al., BMC Bioinformatics 2017. doi:[10.1186/s12859-017-1873-8](https://doi.org/10.1186/s12859-017-1873-8).

Wallace, Edward; Anderson, Felicity; Kavoussanakis, Kostas; Jackson, Michael; Shah, Premal; Lareau, Liana; et al. (2020): riboviz: software for analysis and visualization of ribosome profiling datasets. figshare. Software. doi: [10.6084/m9.figshare.12624200](https://doi.org/10.6084/m9.figshare.12624200)
