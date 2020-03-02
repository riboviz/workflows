# RiboViz and Common Workflow Language (CWL)

Mike Jackson, The Software Sustainability Institute

30 August 2018

## Introduction

This report discusses the [Common Workflow Language](http://www.commonwl.org) (CWL) and its suitability for use with [RiboViz](https://riboviz.org).

---

## Common Workflow Language Specification (CWL)

[Common Workflow Language](http://www.commonwl.org) is:

> "a specification for describing analysis workflows and tools in a way that makes them portable and scalable across a variety of software and hardware environments, from workstations to cluster, cloud, and high performance computing (HPC) environments."

It is used in a number of domains including bioinformatics, medical imaging, astronomy, physics and chemistry.

The current version of the Common Workflow Language Specifications is [v1.0.2](https://www.commonwl.org/v1.0/). CWL documents are written in [YAML](http://yaml.org/), [JSON](http://json.org/) or a mix of the two.

### Tool wrappers

CWL documents called "tool wrappers" or "tool descriptions" provide a specification of a command-line tool. A simple example, wrapping the Linux `echo` command, is:

```cwl
cwlVersion: v1.0       # CWL specification version
class: CommandLineTool # States document describes command-line tool
baseCommand: echo      # Command to run
stdout: output.txt     # Capture stdout to output file
inputs:                # Input parameters
  message:             # Identifier
    type: string       # Data type
    inputBinding:      # Binding (optional), how input parameter should appear
      position: 1      # Where input parameter appears
outputs:               # Outputs
  output:              # Identifier
    type: stdout       # Data type
```

Tool wrappers can specify input flags, parameters, whether these are mandatory or optional and their data types (e.g. string, int, float, double, null, array, record, file, directory, any).

Output parameters can be the output files themselves derived from examining the name or content of the output files (e.g. checking that an output file matching a given pattern now exists)

File formats can also be specified, using existing ontologies (e.g. [EDAM](http://edamontology.org/) ontology for bioinformatics or [IANA](https://www.iana.org/assignments/media-types/media-types.xhtml) media types). For example:

```cwl
inputs:
  aligned_sequences:
    type: File
    label: Aligned sequences in BAM format
    format: edam:format_2572
    inputBinding:
      position: 1

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
```

Specifying file formats both acts as a form of documentation on how to use the tool described by the wrapper and also allows for simple type-checking of tool wrappers.

More generally, tool wrappers allow for metadata to be added. For example, documentation of inputs and outputs, citation information, keywords, code repository locations etc.

Tool wrappers support "hints" to specify software dependencies. A workflow framework that can run CWL workflows might use these hints to theck that these dependencies are available.

Tool wrappers can also list "requirements" to specify both resource requirements (e.g. minimum RAM and number of cores) and requirements on the workflow framework that itself runs the CWL workflows (not all workflow frameworks that can run CWL workflows will support all CWL features).

### Inputs file

Inputs to a specific run of a workflow are described using a YAML document. For example:

```yaml
message: Hello world!
```

### Executing a workflow

A tool wrapper and an inputs file form an executable workflow to run the tool. For example, these can be run using the Common Workflow Language reference implementation, [cwltool](https://github.com/common-workflow-language/cwltool), as follows:

```bash
cwltool echo-tool.cwl echo-job.yml 
```

This outputs information about the workflow's execution, including a description about the output file:

```
/home/ubuntu/miniconda2/bin/cwltool 1.0.20180518123035
Resolved 'echo-tool.cwl' to 'file:///home/ubuntu/cwl/echo-tool.cwl'
[job echo-tool.cwl] /tmp/tmpnyEoBe$ echo \
    'Hello world!' > /tmp/tmpnyEoBe/output.txt
[job echo-tool.cwl] completed success
{
    "output": {
        "checksum": "sha1$47a013e660d408619d894b20806b1d5086aab03b", 
        "basename": "output.txt", 
        "location": "file:///home/ubuntu/cwl/output.txt", 
        "path": "/home/ubuntu/cwl/output.txt", 
        "class": "File", 
        "size": 13
    }
}
Final process status is success
```

The contents of the output file is:

```bash
cat output.txt
```
```
Hello World
```

### Workflows 

A workflow consisting of the invocation of a number of tools is written as a CWL document that describes:

* Input parameters for the workflow (as above).
* Output parameters for the workflow (as above).
* Each step in the workflow, including the name of the tool wrapper file that describes the tool invokes that step, its inputs and its outputs.

The execution order of the steps is determined by the inputs and outputs specified by each step.

See, for example, the extract-Java-source-from-tar-file-then-compile workflow in [Writing workflows](https://www.commonwl.org/user_guide/21-1st-workflow/) in the [Common Workflow Language User Guide](https://www.commonwl.org/user_guide/).

Workflows, complemented by inputs files, are executed in the same way as tool descriptions.

### Implementations

[cwltool](https://github.com/common-workflow-language/cwltool) is a reference implementation, in Python, of CWL.

The CWL web site has a list of other CWL [Implementations](https://www.commonwl.org/#Implementations). This list includes related tools which include editors, viewers, converters and code generators.

Implementations of CWL in workflow frameworks vary in the extent to which they support execution of CWL workflows. For example, Toil and cwl-tes are noted as "build passing", Arvados as "build unstable", and Galaxy and Apache Taverna as "alpha".

[Apache Taverna](https://taverna.incubator.apache.org/) comment on [Download Taverna Command-line Tool](https://taverna.incubator.apache.org/download/commandline/) that:

> Support for executing Common Workflow Language workflows is planned. Please contact the dev@taverna mailing list for details.

Their issue tracker ([877](https://jira.apache.org/jira/browse/TAVERNA-877), [900](https://jira.apache.org/jira/browse/TAVERNA-900)) notes some work as part of Google Summer of Code 2016, but it appears to have stalled.

[Galaxy](https://galaxyproject.org/) has been [forked](https://github.com/common-workflow-language/galaxy) into the [common-workflow-language](https://github.com/common-workflow-language/) GitHub project where work is underlway to implement CWL within Galaxy.

### CWL and Docker

[Docker](https://www.docker.com/) is a platform for building and distributing applications along with all their dependencies. Unlike a virtual machine, Docker does not include a separate operating system. Rather, it exploits Linux kernel resource isolation (CPU, memory, block I/O, network) and namespaces to allow independent "containers" to run within a single Linux instance. As it does not include the operating system itself, Docker files (called "images") are far smaller than virtual machine image files.

CWL tool wrappers can specify tools that are bundled as Docker images. The framework that runs the CWL workflow takes care of starting Docker, loading the Docker image, and sharing input and output files between the running Docker image and the host upon which the CWL workflow is running.

For more information, see [Running Tools Inside Docker](https://www.commonwl.org/user_guide/07-containers/) in the [Common Workflow Language User Guide](https://www.commonwl.org/user_guide/).

### Dockstore

[Dockstore](https://dockstore.org/) is an open platform used by the [Global Alliance for Genomics and Health](https://www.ga4gh.org/) (GA4GH). It is provided to allow researchers to share Docker-based tools described in CWL. To do this Dockstore expects the following:

* Source code to be hosted on [GitHub](https://github.com), [Bitbucket](https://bitbucket.org) or [GitLab](https://gitlab.com).
* A `Dockerfile` in the Git repository that describes how to create a Docker image bundling the software.
* Docker images hosted on [Quay.io](https://quay.io/) or [Docker Hub](https://hub.docker.com/), automatically built from the source code.
* A `Dockstore.cwl` file in the Git repository that describes how to call the tools inside the Docker image.

Dockstore provide tools to allow resarchers to invoke tools registered with them.

See [Getting started with Docker and CWL](https://docs.dockstore.org/docs/prereqs/getting-started-with-docker/) and [Launching Tools and Workflows](https://docs.dockstore.org/docs/user-tutorials/launch/) for more information.

---

## Ease of use

The [Common Workflow Language User Guide](https://www.commonwl.org/user_guide/) provides a readable, step-by-step, [Software Carpentry](https://software-carpentry.org)-style, introduction to the languange.

I tried to run both the "hello world!" example using both the Common Workflow Language reference implementation, [cwltool](https://github.com/common-workflow-language/cwltool), and the implementation in the [Toil](https://toil.readthedocs.io) Python workflow engine. To download and run these took less than 30 minutes.

See [Appendix: running a "hello world!" CWL workflow](#appendix-running-a-hello-world-cwl-workflow).

---

## Using CWL for RiboViz

RiboViz's workflow is currently implemented as a shell script (`prepRiboViz.sh`). Reimplementing this in Python has been proposed as an option to improve its maintainability. However, implementing it as a CWL workflow and using the CWL reference implementation's `cwltool` is another.

The advantages of using CWL to implement RiboViz's workflow are as follows:

* The specification can describe both command-line tools (and so be used to describe RiboViz's Python and R scripts) and those bundled in Docker images (should RiboViz's Python and R scripts be bundled into Docker images to avoid researchers having to worry about their myriad dependencies).
* The specification is supported, or in the process of being supported, by a number of workflow frameworks.
* It is an open, plain-text format, based on YAML and JSON. Migrating to another specification, or even back to an implementation in a programming language, at a later date, would not be a major issue.

To refactor RiboViz to use CWL would require:

* Writing tool wrappers for RiboViz's Python (`trim_5p_mismatch.py`) and R (`bam_to_h5.R`, `generate_stats_figs.R`) scripts. These can be as simple or as rich as desired.
* Getting, or writing, tool wrappers, for each third-party tool invoked by `prepRiboViz.sh` (`samtools view|sort|index`, bedtools genomecov`, `hisat2-build`, `hisat2`):
  - A [tools](https://github.com/common-workflow-language/workflows/tree/master/tools) directory in the [workflows](https://github.com/common-workflow-language/workflows) repository of the [common-workflow-language](https://github.com/common-workflow-language/workflows/tree/master/tools) GitHub project has a number of examples, including CWL documents for `samtools view|sort|index` and `bedtools genomecov` (all of which use Docker).
  - Alternatively, simple tool wrappers for these could be initially written for command-line use.
* Migrating the `prepRiboViz.sh` into a CWL workflow.
* It would be useful to refactor the Python and R scripts so they always return non-zero codes when problematic conditions arise (i.e. they only exit with a return code of zero if they have successfully completed) and provide informative error messages. This would allow CWL workflow executor to identify when such errors have arisen. This is a desirable quality for the Python and R scripts regardless of whether CWL is adopted.

---

## Additional references

Leipzig, Jeremy (2017) "A review of bioinformatic pipeline frameworks". Briefings in Bioinformatics, Volume 18, Issue 3, 1 May 2017, Pages 530-536, doi:[10.1093/bib/bbw020](https://doi.org/10.1093/bib/bbw020).

---

## Appendix: running a "hello world!" CWL workflow

### cwltool

Common Workflow Language reference implementation, [cwltool](https://github.com/common-workflow-language/cwltool)

Get latest release:

```bash
wget https://github.com/common-workflow-language/cwltool/archive/1.0.20180809224403.tar.gz
tar -xvf 1.0.20180809224403.tar.gz
cd cwltool-1.0.20180809224403/
pip install .
which cwltool
```
```
/home/ubuntu/miniconda2/bin/cwltool
```

Check version:

```bash
cwltool --version
```
```
/home/ubuntu/miniconda2/bin/cwltool 1.0
```

Create `echo-tool.cwl`:

```cwl
cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
stdout: output.txt
inputs:
  message:
    type: string
    inputBinding:
      position: 1
outputs:
  output:
    type: stdout
```

Create `echo-job.yml`:

```yml
message: Hello world!
```

Run:

```bash
cwltool echo-tool.cwl echo-job.yml 
```
```
/home/ubuntu/miniconda2/bin/cwltool 1.0.20180518123035
Resolved 'echo-tool.cwl' to 'file:///home/ubuntu/cwl/echo-tool.cwl'
[job echo-tool.cwl] /tmp/tmpnyEoBe$ echo \
    'Hello world!' > /tmp/tmpnyEoBe/output.txt
[job echo-tool.cwl] completed success
{
    "output": {
        "checksum": "sha1$47a013e660d408619d894b20806b1d5086aab03b", 
        "basename": "output.txt", 
        "location": "file:///home/ubuntu/cwl/output.txt", 
        "path": "/home/ubuntu/cwl/output.txt", 
        "class": "File", 
        "size": 13
    }
}
Final process status is success
```

Check output:

```bash
cat output.txt
```
```
Hello world!
```

### Ubuntu cwltool package

It was discovered that Ubuntu has a `cwltool` package ([bionic (18.04LTS) > science > cwltool](https://packages.ubuntu.com/bionic/cwltool))which is installed as a dependency of `samtools`, used by RiboViz


```bash
cwltool --version
```
```
/usr/bin/cwltool 1.0.20180302231433
```

Running this gave the same results as above.

### Toil

[Toil](https://toil.readthedocs.io) is a workflow engine written in Python.

See [Running a basic CWL workflow](https://toil.readthedocs.io/en/3.15.0/gettingStarted/quickStart.html#running-a-basic-cwl-workflow).

Install Toil and `cwl` extra:

```bash
pip install toil[cwl]
toil-cwl-runner --version
```
```
DEBUG:rdflib:RDFLib Version: 4.2.2
3.17.0
```

Run:

```bash
toil-cwl-runner --defaultMemory=100M --defaultDisk=100M echo-tool.cwl  echo-job.yml 
```
```
...
Resolved 'echo-tool.cwl' to 'file:///home/ubuntu/cwl/echo-tool.cwl'
...
2018-08-29 09:15:03,034 - toil.leader - INFO - Finished toil run successfully.
{
    "output": {
        "checksum": "sha1$47a013e660d408619d894b20806b1d5086aab03b", 
        "basename": "output.txt", 
        "nameext": ".txt", 
        "nameroot": "output", 
        "http://commonwl.org/cwltool#generation": 0, 
        "location": "file:///home/ubuntu/cwl/output.txt", 
        "class": "File", 
        "size": 13
    }
...
```

Check output:

```bash
cat output.txt
```
```
Hello world!
```

**Troubleshooting: `AssertionError: The job ... is requesting 2147483648 bytes of memory...`**

If you get:

```
AssertionError: The job file:///home/ubuntu/cwl/echo-tool.cwl is requesting 2147483648 bytes of memory, more than the maximum of 2065702912 this batch system was configured with.
```

Then explictly specify the default memory required e.g. `--defaultMemory=100M`

**Troubleshooting: `toil.batchSystems.abstractBatchSystem.InsufficientSystemResources`**

If you get:

```
toil.batchSystems.abstractBatchSystem.InsufficientSystemResources: Requesting more disk than either physically available, or enforced by --maxDisk. Requested: 3221225472, Available: 2719440896
```

Then explictly specify the default disk space required e.g. `--defaultDisk=100M echo-tool.cwl`.
