# Notes on Workflow Options 28/02/2020

## A review of bioinformatic pipeline frameworks

Leipzig, Jeremy (2017) A review of bioinformatic pipeline frameworks. Briefings in Bioinformatics, Volume 18, Issue 3, 1 May 2017, Pages 530-536, https://doi.org/10.1093/bib/bbw020.

Requirements:

* Accommodate production pipelines.
* Serial and parallel steps.
* Complex dependencies.
* Myriad software and data file types.
* Fixed and user-defined parameters and deliverables.

Extras:

* Visualise progress in real time.
* Instantiate containerized tools.
* Run over distributed clusters or cloud.
* Create workflows via GUIs.

Dimensions:

* Implicit or explicit syntax for specifying workflow.
* Configuration, convention or class-based design paradigm.
* Command-line or workbench/GUI.

Implicit convention frameworks:

* Make-style implicit order of execution, based on target naming.
* Leverage scripting languages to implement logic inside and outside of rules.
* Snakemake, Nextflow, BigDataScript.

Explicit convention frameworks:

* Script-style, fixed, explicit order of execution.
* Tasks can reference other tasks.
* Ruffus, bpipe.

Explicit configuration frameworks:

* Describe individual job run instances and dependencies (using job IDs).
* Pegasus:
  - XML file describes job run instances and dependencies.
  - Explicit job IDs identify task antecedents.

Explicit class-based frameworks:

* Closely bound to existing code libraries rather than various executables.
* Often contain many lines of code implementing domain logic. 
* Queue, bound to Genome Analysis Toolkit (GATK) (Scala), Toil (Python)

Implicit class-based frameworks:

* Closely bound to existing code libraries rather than various executables.
* Luigi (Python).

Explicit configuration open source erver workbenches:

* Link preconfigured modular tools together via drag-and-drop GUI. 
* Require exacting specifications of inputs and outputs.
* Galaxy:
  - Web-based interface for command-line tools.
  - New components require 20 lines of configuration.
  - Wrappers less robust than Taverna.
* Taverna standalone clients, and access tools over Internet.
  - New components require XML specification.

Explicit configuration commercial cloud workbenches:

* DNAnexus, SevenBridges (SBGenomics, Illumina's BaseSpace.

Explicit configuration open source cloud API:

* Curoverse's Arvados, iPlant Collaborative's Agave.

Common Workflow Language Specification (CWL):

* Consistent means of distributing popular tools across myriad frameworks.
* Toil, Galaxy, Taverna, SevenBridges. Arvados.

Containerization of bioinformatic tools:

* Enable frameworks to accommodate tools with complex software dependencies.
* Docker.

Choosing a pipeline framework:

* Choosing between implicit or explicit syntax is a question of personal preference.
* Implicit syntax: can seem unintuitive (-), idiomatic style (+), convenience (+).
* Convention-based: encourage high level of internal business logic (+), enable polished deliverables (e.g. HTML, PDF) (+), seem less reproducible than configuration-based pipelines (-).
* Configuration-based: cluster schedulers can anticipate load and allocate memory and compute resources (+), often require dynamic tools to create static configurations (-).
* Workbenches and class-based: ease-of-use (+), performance (+), heavyweight (-), require highly-skilled developers (-), performance improvements not guaranteed to justify development time (-), HPC favours class-based (+)
* Cloud-based platforms: scalability (+), collaborative research advantages (+).

Consider "return on investment"¬

* Conduct large-scale, highly repetitive research, require high degree of data provenance and versioning => configuration-based pipelines
* Exploratory proofs-of-concept => explicit DSL-based pipelines.

"Once a script-based pipeline is implemented in one framework, transitioning to a different one is relatively simple should priorities change."

---

## biocontainers project

https://biocontainers.pro/#/

Felipe da Veiga Leprevost, Björn A Grüning, Saulo Alves Aflitos, Hannes L Röst, Julian Uszkoreit, Harald Barsnes, Marc Vaudel, Pablo Moreno, Laurent Gatto, Jonas Weber, Mingze Bai, Rafael C Jimenez, Timo Sachsenberg, Julianus Pfeuffer, Roberto Vera Alvarez, Johannes Griss, Alexey I Nesvizhskii, Yasset Perez-Riverol, BioContainers: an open-source and community-driven framework for software standardization, Bioinformatics, Volume 33, Issue 16, 15 August 2017, Pages 2580–2582, https://doi.org/10.1093/bioinformatics/btx192

* Open source, community drive.
* Platform independent executable environments.
* Container-based technologies:
  - Docker, https://www.docker.com/
  - rkt, https://coreos.com/rkt/
* Integrate containers into piplelines and architectures.
* Infrastructure and guidelines to create, manage and distribute containers.
* Integrated with BioConda, https://bioconda.github.io/
  - Automatically generate containers for each BioConda recipe.
* 2700+ containers.
* Container provider with open source projects:
  - Galaxy, https://galaxyproject.org/
    - Uses Docker containers to solve workflow dependencies.
  - PhenoMeNal H2020, https://phenomenal-h2020.eu/home/
    - Adopted and implemented BioContainers guidelines.
    - Deploying containers into the BioContainers architecture.
* Workflow systems:
  - OpenMS, https://www.openms.de/, runs within Galaxy etc.
  - Taverna
  - Galaxy

GitHub organisation:

* https://github.com/BioContainers/
* Docker files.
* Specification.
* Tools to create/manage containers.

BioContainers registries and Registry-UI:

* https://biocontainers.pro/registry/
* Containers built and made available for download and use by Docker or rkt.
* Search, tag and find BioContainers independently of where they have been deployed. 

Request new container:

* Open issue in https://github.com/BioContainers/containers/.
* Specify software (name, URL or binary to be packaged).
* Member of the BioContainers community:
  - Picks up the issue.
  - Generate the specific container.
  - Configures/deploys automated build system.
* New container available within hours.

Create and build container:

* EITHER create BioConda recipe for the software:
  - BioConda guidelines, https://bioconda.github.io/guidelines.html.
  - Container generation tool, https://github.com/BioContainers/auto-mulled/, automatically creates a "mulled" container.
    - Uses involucro, https://github.com/involucro/involucro, to create containers without any Dockerfile definition, reusing already existing recipes from other package managers e.g. Conda, Alpine.
    - involucro installs Conda package into build-time container which has preferred package manager already installed.
    - Copies resulting new image layer on top of runtime environment defined by BioContainers (busybox).
  - Container is pushed to BioContainers quay.io registry.
* OR create Dockerfile recipe in https://github.com/BioContainers/containers.
  - Instructions to create complete container.
  - Template, https://github.com/BioContainers/specs/blob/master/container-specs.md
  - Name, version, licence, web page, maintainer.
  - Consistent with metadata needed for BioConda recipie.

Documentation:

* https://biocontainers-edu.biocontainers.pro/en/latest/
* Docker Containers: Dockerfile recipes to automatically build containers in BioContainers.
* Conda based Containers: Conda recipes to automatically build first a conda package and then a Docker container.
* Biocontainers Registry: hosted registry of all BioContainers images that are ready to be used.
* Specifications: guidelines and rules to contribute with BioContainers.

Integration with workflow engines

* https://biocontainers-edu.biocontainers.pro/en/latest/galaxy_container_dependencies.html#integration-with-workflow-engines
* Galaxy tools (also called wrappers) can use Conda packages and Docker containers as dependency resolvers.
* Recommends Conda packages as Docker is not available on every (HPC-) system.
* Conda can be installed by Galaxy and maintained entirely in user space.
* Run jobs in container runtimes:
  - Docker, requires `sudo` with no password prompt.
  - Singularity, https://www.sylabs.io/.

---

## Developing reproducible bioinformatics analysis workflows for heterogeneous computing environments to support African genomics

Choice determined by:

* Existing community workflow standards or workflow systems with languages commonly in use by the bioinformatics community.
* Existing skills within H3ABioNet and its collaborators.

Selections:

* Common workflow language (CWL).
* Nextflow.

---

## Given the experience of others writing bioinformatic pipelines, what are the pros/cons of Toil vs Snakemake vs Nextflow?

[Given the experience of others writing bioinformatic pipelines, what are the pros/cons of Toil vs Snakemake vs Nextflow?](https://www.reddit.com/r/bioinformatics/comments/a4fq4i/given_the_experience_of_others_writing/), December 2018

samiwillbe:

* Snakemake:
  - Don't like backwards make-style of thinking (results -> inputs).
* Toil:
  - Could never get it to behave how I expected.
* Nextflow:
  - Made sense and worked straight away.
  - Once dataflow concepts were grasped, wrote fairly complicated pipelines that ran locally, on clusters, and in the cloud, transparently.
  - Team ofcomputational biologists were able to come up to speed quickly and write their own pipelines.
* WDL/Cromwell
  - Nice features e.g. subworkflows.
  - CWL/WDL aims for declarative separation of workflow description and execution.
  - Works great for 90% of cases.
  - Remaining 10% of edge cases require more custom code and compromises on performance.
  - Unlike CWL and WDL, Nextflow is a proper programming language witth the tools needed to handle edge cases in much cleaner ways.
* If you work with multiple, complicated pipelines day-to-day, tweaking, optimizing, and experimenting then go with Nextflow.
* If you work on a single, relatively simple, canonical pipeline then WDL/Cromwell is probably fine.

geoffjentry (Cromwell development lead, WDL, CWL stakeholder):

* Cromwell/WDL
  - WDL and CWL help create reproducible and shareable workflows which can be easily shared across different environments, platforms, etc.
  - Immensely valuable for those who distribute, and use, commonly used tools and workflows.
* Snakemake/Nextflow:
  - Good for developing in house, individual use, workflows.
  - Abstractions required for "run anywhere" can become irritating (if they're not necessary, why bother?)

theHM:

* Nextflow:
  - Makefiles => Snakemake => Nextflow.
  - 90% of time that used to be spent optimising parallelization, file format wrangling and HPC submission now trivial with minimal boilerplate code.
  - Fail-safe hassles, retrying jobs, and containerisation. 
  - Seamless embedding of R scripts, Python snippets, Bash snippets.
  - Groovy-based, so seems daunting at first. 
* Snakemake:
  - Has added some of these features but feels more kludgey and less intuitive.

sayerskt:

* Toil/CWL:
  - CWL painful if workflow requires anything above basic logic.
  - Need to write JavaScript to achieve things that other workflow managers handle by default (e.g. rearranging output folders).
* Snakemake:
  - Nice for working with files, very concise, largely just works.
  - Documentation can be bit lacking.
  - Doesn't offer much flexibility
* Nextflow:
  - Offers incredible flexibility.
  - Easily handle anything from simple workflows to those with considerable logic.
  - Documentation and community are quite good.

davidmasp:

* Bash:
  - Simple stuff.
* Snakemake:
  - Couple of frustrating experiences.
* Nextflow:
  - Mostly write in this.
  - Groovy.
  - Error messages can be cryptic sometimes although it returns software output stream directly. Experience with Snakemake was much worse.
  - Portability and version control.
  - conda environments, Docker/Singularity images.
  - User friendly framework.
  - Reporting, logging, resuming, emails.
  - Great community.
  - HPC integration.

samuellampa:

* Snakemake:
  - Great when you want to reproduce a specific output that you know the desired filename for.
  - Really complicated to write complex workflows that nest a lot of looping and scatter/gather etc - need to mentally map the computational graph to file name patterns.
* Nextflow:
  - Much more straight-forward and more natural to reason about
  - Define steps data should go through, and their order, input data into the initial components, let data flow through the dataflow network.
* SciPipe
  - Reproduce specific outputs in a workflow system.
  - http://scipipe.org/howtos/partial_workflows
  - https://www.biorxiv.org/content/early/2018/10/06/380808)

agapow:

* Snakemake:
  - Big fan, but for a particular use - analysis, hacking at data, "tactical" pipeline.
  - Python, so works well with analysis in Python.
  - Other pipelines have other needs and focuses (reliability, operationalising, HPC control) so may not be the best choice there.

gumbos:

* Toil:
  - Most flexibility and features.
  - Steep learning curve and code overhead. 
* Snakemake:
  - Opposite of the above.

Summary:

* Nextflow seems to be viewed very highly.
* CWL/WDL seem to be viewed OK.
* Snakemake seems to be viewed OK.

---

## Poll of workflow managers, December 2018

https://twitter.com/AlbertVilella/status/1069635987427532800

* 44.4% NextFlow
* 22.2% SevenBridges (SBGenomics)
* 19.4% CWL
* Snakemake
* Galaxy

36 responses.

---

## Tools

Implicit convention frameworks:

* SnakeMake
  - https://github.com/snakemake/snakemake
  - https://snakemake.readthedocs.io/en/stable/
* NextFlow
  - https://github.com/nextflow-io/nextflow
  - https://www.nextflow.io/
* BigDataScript
  - https://github.com/pcingola/BigDataScript
  - https://pcingola.github.io/BigDataScript/

Explicit convention frameworks:

* Ruffus
  - http://www.ruffus.org.uk/
  - https://github.com/cgat-developers/ruffus
* bpipe
  - http://docs.bpipe.org/
  - https://github.com/ssadedin/bpipe

Explicit configuration frameworks:

* Pegasus
  - https://pegasus.isi.edu/
  - https://github.com/pegasus-isi/pegasus

Explicit class-based frameworks:

* Queue, bound to Genome Analysis Toolkit (GATK)
  - A response to [GATK best practices pipelinee written in scala for Queue](https://gatkforums.broadinstitute.org/wdl/discussion/8221/gatk-best-practices-pipelinee-written-in-scala-for-queue) (August 2016) comments that "As we are moving away from Queue at this point, I would like to recommend WDL to you. WDL, or Workflow Definition Language is designed to be more intuitive than other pipelining solutions (like Queue/scala)" See, for example, [Pipelining GATK with WDL and Cromwell](https://gatk.broadinstitute.org/hc/en-us/articles/360035889771-Pipelining-GATK-with-WDL-and-Cromwell).
* Toil
  - http://toil.ucsc-cgl.org/
  - https://github.com/bd2kgenomics/toil

Implicit class-based frameworks:

* Luigi
  - https://github.com/spotify/luigi

Explicit configuration open source server workbenches:

* Galaxy
  - https://galaxyproject.org/
  - https://github.com/galaxyproject/galaxy
* Taverna
  - https://taverna.incubator.apache.org/
  - https://taverna.incubator.apache.org/download/code/
  - https://taverna.incubator.apache.org/download/commandline/
    - Run workflows from a command prompt, essentially Taverna Workbench stripped of its GUI.

Explicit configuration commercial cloud workbenches:

* DNAnexus, http://dnanexus.com
* SevenBridges http://sbgenomics.com
* Illumina's BaseSpace, http://basespace.illumina.com

Explicit configuration open source cloud API:

* Curoverse's Arvados, https://curoverse.com
* iPlant Collaborative's Agave,
  - http://developer.agaveapi.co/
  - https://github.com/smirarab/iplant-agave-sdk/blob/master/README.md
  - Not sure how current this is.

Others:

* Cromwell
  - https://github.com/broadinstitute/cromwell
  - http://cromwell.readthedocs.io/

CWL's list of [Computational Data Analysis Workflow Systems](https://github.com/common-workflow-language/common-workflow-language/wiki/Existing-Workflow-systems) and [Awesome Pipeline](https://github.com/pditommaso/awesome-pipeline) refer to 100s of tools.

---

## Workflow languages

Common Workflow Language (CWL)

* https://www.commonwl.org/
* https://www.commonwl.org/#Implementations
* Production implementations:
  - cwltool - reference implementation (Python)
  - Toil
  - Cromwell
* Development implementations:
  - Galaxy
  - Taverna

Workflow Description Language (WDL)

* https://openwdl.org/
* https://github.com/openwdl/wdl
* Cromwell
* dxWDL, https://github.com/dnanexus/dxWDL
* miniwdl, https://github.com/chanzuckerberg/miniwdl

---

## Summary (28/02/20)

"Google" is hits from `"<NAME>" "bioinformatics"`.

| Product | Licence | Start | Last update | Releases | Latest | Contributors | Forks | Google | Workflow Language |
| ------- | ------- | ----- | ----------- | -------- | ------ | ------------ | ----- | ------ | ---------------- |
| SnakeMake (*) | MIT | 2013 | this week | 87 | 5.10.0 | 122 | 73 | 23,000 | - |
| NextFlow (*) | Apache 2.0 | 2013 | this week | 171 | 20.01.0 | 230 | 60 | 23,800 | - |
| BigDataScript | Apache 2.0 | 2014 | last month | 40 | v2.1a | 8 | 24 | 534 | - |
| Ruffus | MIT | 2013 | 6 months | 4 | v2.8.1 | 17 | 30 | 10,000 | - |
| bpipe | New BSD | 2013 | last month | 16 | 0.9.9.8 | 11 | 45 | 2,440 | - |
| Pegasus (*) | Apache 2.0 | 2007 | today | 40 | 4.9.3 | 19 | 55 | 750,000 | - |
| Toil | Apache 2.0 | 2011 | today | 50 | 3.24.0 | 81 | 183 | 20,700 | CWL |
| Luigi | Apache 2.0 | 2014 | last week | 55 | 2.8.12 | 419 | 2100 | 49,000 | - |
| Galaxy (*) | Academic Free 3.0 | 2005 | today | 147 | v20.01 | 225 | 594 | 655,000 | CWL |
| Taverna | Apache 2.0 | 2003 | 2019 | ? | 3.1.0 | ? | ? | 104,000 | CWL | 
| Cromwell | MIT | 2015 | this week | 57 | 49 | 89 | 184 | 80,700 | CWL, WDL |

(*) are tools that have neem mentioned in RiboViz meetings.
