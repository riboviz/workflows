# RiboViz and Common Workflow Language (CWL)

Updated from [RiboViz and Common Workflow Language (CWL)](https://github.com/riboviz/bbsrc-nsf/blob/master/docs/reports/201908-ssi-review/SsiRiboVizCwl.md), Mike Jackson, The Software Sustainability Institute, 30 August 2018.

## Common Workflow Language Specification (CWL)

[Common Workflow Language](http://www.commonwl.org) is:

> an open standard for describing analysis workflows and tools in a way that makes them portable and scalable across a variety of software and hardware environments, from workstations to cluster, cloud, and high performance computing (HPC) environments.

It is used in a number of domains including bioinformatics, medical imaging, astronomy, physics and chemistry.

The current version of the Common Workflow Language Specifications is [v1.1](https://www.commonwl.org/v1.1/).

Peter Amstutz, Michael R. Crusoe, Nebojša Tijanic (editors), Brad Chapman, John Chilton, Michael Heuer, Andrey Kartashov, Dan Leehr, Hervé Ménager, Maya Nedeljkovich, Matt Scales, Stian Soiland-Reyes, Luka Stojanovic (2016): Common Workflow Language, v1.0. Specification, Common Workflow Language working group. https://w3id.org/cwl/v1.0/ doi:10.6084/m9.figshare.3115156.v2

---

## CWL documents

CWL documents are written in [YAML](http://yaml.org/), [JSON](http://json.org/) or a mix of the two.

---

## Tool wrappers

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

---

## Inputs file

Inputs to a specific run of a workflow are described using a YAML document. For example:

```yaml
message: Hello world!
```

---

## Executing a workflow

A tool wrapper and an inputs file form an executable workflow to run the tool. For example, these can be run using the Common Workflow Language reference implementation, [cwltool](https://github.com/common-workflow-language/cwltool), as follows:

```console
$ cwltool echo-tool.cwl echo-job.yml 
```

This outputs information about the workflow's execution, including a description about the output file:

```
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool/bin/cwltool 2.0.20200224214940
INFO Resolved 'echo-tool.cwl' to 'file:///home/ubuntu/cwl-tutorial/echo-tool.cwl'
INFO [job echo-tool.cwl] /tmp/wuha0o42$ echo \
    'Hello world!' > /tmp/wuha0o42/output.txt
INFO [job echo-tool.cwl] completed success
{
    "output": {
        "location": "file:///home/ubuntu/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$47a013e660d408619d894b20806b1d5086aab03b",
        "size": 13,
        "path": "/home/ubuntu/cwl-tutorial/output.txt"
    }
}
INFO Final process status is success
```

The contents of the output file is:

```console
$ cat output.txt
Hello World
```

---

## Workflows 

A workflow consisting of the invocation of a number of tools is written as a CWL document that describes:

* Input parameters for the workflow (as above).
* Output parameters for the workflow (as above).
* Each step in the workflow, including the name of the tool wrapper file that describes the tool invokes that step, its inputs and its outputs.

The execution order of the steps is determined by the inputs and outputs specified by each step.

See, for example, the "extract Java source from tar file then compile" workflow in [Writing workflows](https://www.commonwl.org/user_guide/21-1st-workflow/) in the [Common Workflow Language User Guide](https://www.commonwl.org/user_guide/).

Workflows, complemented by inputs files, are executed in the same way as tool descriptions.

---

## Implementations

[cwltool](https://github.com/common-workflow-language/cwltool) is a reference implementation, in Python, of CWL.

The CWL web site has a list of other CWL [Implementations](https://www.commonwl.org/#Implementations). This list includes related tools which include editors, viewers, converters and code generators.

Implementations of CWL in workflow frameworks vary in the extent to which they support execution of CWL workflows. For example, Toil and Cromwell are listed as "In production" whereas Taverna and Galasy are listed as "In development".

[Apache Taverna](https://taverna.incubator.apache.org/) comment on [Download Taverna Command-line Tool](https://taverna.incubator.apache.org/download/commandline/) that:

> Support for executing Common Workflow Language workflows is planned. Please contact the dev@taverna mailing list for details.

This is unchanged from the previous version of this document, written back in August 2018.

Their issue tracker ([877](https://jira.apache.org/jira/browse/TAVERNA-877), [900](https://jira.apache.org/jira/browse/TAVERNA-900)) notes some work as part of Google Summer of Code 2016 and interest in 2018, but no active development.

[Galaxy](https://galaxyproject.org/) has been [forked](https://github.com/common-workflow-language/galaxy) into the [common-workflow-language](https://github.com/common-workflow-language/) GitHub project where work is underlway to implement CWL within Galaxy. The last activity was in November 2019.

---

## CWL and Docker

[Docker](https://www.docker.com/) is a platform for building and distributing applications along with all their dependencies. Unlike a virtual machine, Docker does not include a separate operating system. Rather, it exploits Linux kernel resource isolation (CPU, memory, block I/O, network) and namespaces to allow independent "containers" to run within a single Linux instance. As it does not include the operating system itself, Docker files (called "images") are far smaller than virtual machine image files.

CWL tool wrappers can specify tools that are bundled as Docker images. The framework that runs the CWL workflow takes care of starting Docker, loading the Docker image, and sharing input and output files between the running Docker image and the host upon which the CWL workflow is running.

For more information, see [Running Tools Inside Docker](https://www.commonwl.org/user_guide/07-containers/) in the [Common Workflow Language User Guide](https://www.commonwl.org/user_guide/).

---

## Dockstore

[Dockstore](https://dockstore.org/) is an open platform used by the [Global Alliance for Genomics and Health](https://www.ga4gh.org/) (GA4GH). It is provided to allow researchers to share Docker-based tools described in CWL. To do this Dockstore expects the following:

* Source code to be hosted on [GitHub](https://github.com), [Bitbucket](https://bitbucket.org) or [GitLab](https://gitlab.com).
* A `Dockerfile` in the Git repository that describes how to create a Docker image bundling the software.
* Docker images hosted on [Quay.io](https://quay.io/) or [Docker Hub](https://hub.docker.com/), automatically built from the source code.
* A `Dockstore.cwl` file in the Git repository that describes how to call the tools inside the Docker image.

Dockstore provide tools to allow resarchers to invoke tools registered with them.

See [Getting Started with Docker](https://docs.dockstore.org/en/develop/getting-started/getting-started-with-docker.html), [Getting started with CWL](https://docs.dockstore.org/en/develop/getting-started/getting-started-with-cwl.html) and [Launching Tools and Workflows](https://docs.dockstore.org/en/develop/launch-with/launch.html) for more information.
