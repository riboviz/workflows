# Common Workflow Language User Guide Notes

https://www.commonwl.org/user_guide/

## Running a "hello world!" CWL workflow

### cwltool

Common Workflow Language reference implementation, [cwltool](https://github.com/common-workflow-language/cwltool)

Note: though Ubuntu 18.04 (my virtual machine) has a `cwltool` package ([bionic (18.04LTS) > science > cwltool](https://packages.ubuntu.com/bionic/cwltool)), I opted to use a `pip` install for consistency with other Python dependencies used by RiboViz.

Create a new conda environment from the current RiboViz one and activate it:

```console
$ conda create --name riboviz-cwltool --clone base
$ conda activate riboviz-cwltool
```

Install CWL reference implementation and node.js:

```console
$ pip install cwlref-runner
$ cwl-runner --version
/home/ubuntu/miniconda3/envs/riboviz-cwltool/bin/cwl-runner 2.0.20200224214940
$ cwltool --version
/home/ubuntu/miniconda3/envs/riboviz-cwltool/bin/cwltool 2.0.20200224214940

$ sudo apt-get -y install nodejs
```

`cwl-runner` is an alias for the default CWL interpreter installed on a host.

Create tool wrapper, `echo-tool.cwl`:

```cwl
#!/usr/bin/env cwl-runner

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

Create input object, `echo-job.yml`:

```yml
message: Hello world!
```

Run:

```console
$ cwl-runner [tool-or-workflow-description] [input-job-settings]
```
```console
$ cwltool echo-tool.cwl echo-job.yml 
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool/bin/cwltool 2.0.20200224214940
INFO Resolved 'echo-tool.cwl' to 'file:///home/ubuntu/cwl-tutorial/echo-tool.cwl'
INFO [job echo-tool.cwl] /tmp/0wg8bch2$ echo \
    'Hello world!' > /tmp/0wg8bch2/output.txt
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
```
$ cat output.txt 
Hello world!
```

### Toil

[Toil](https://toil.readthedocs.io) is a workflow engine written in Python. It is listed by CWL as having an "in production" implementation of CWL.

See [Running a basic CWL workflow](https://toil.readthedocs.io/en/latest/gettingStarted/quickStart.html#running-a-basic-cwl-workflow):

Create a new conda environment from the current RiboViz one and activate it:

```console
$ conda create --name riboviz-toil --clone base
$ conda activate riboviz-toil
```

Install Toil and `cwl` extra:

```console
$ pip install toil[cwl]
$ cwltool --version
/home/ubuntu/miniconda3/envs/riboviz-toil/bin/cwltool 1.0.20190906054215
$ toil-cwl-runner --version
3.24.0
```

Note that `cwltool` is `1.0.x`, not `2.0.x`.

Run:

```console
$ toil-cwl-runner echo-tool.cwl echo-job.yml 
ubuntu 2020-03-03 08:51:05,402 MainThread INFO cwltool: Resolved 'echo-tool.cwl' to 'file:///home/ubuntu/cwl-tutorial/echo-tool.cwl'
ubuntu 2020-03-03 08:51:05,976 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxCores to CPU count of system (4).
ubuntu 2020-03-03 08:51:05,976 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxMemory to physically available memory (8340123648).
ubuntu 2020-03-03 08:51:05,976 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxDisk to physically available disk (19029254144).
ubuntu 2020-03-03 08:51:05,989 MainThread INFO toil: Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
ubuntu 2020-03-03 08:51:05,992 MainThread INFO toil.leader: Issued job 'file:///home/ubuntu/cwl-tutorial/echo-tool.cwl' echo kind-file_home_ubuntu_cwl-tutorial_echo-tool.cwl/instancegehkd1py with job batch system ID: 0 and cores: 1, disk: 3.0 G, and memory: 2.0 G
INFO:toil.worker:Redirecting logging to /tmp/toil-8e1d4067-cf05-45ed-94f9-f02a5548da7f-453138a4fa34406f996060e86264aa60/tmp3jb592ti/worker_log.txt
ubuntu 2020-03-03 08:51:06,588 MainThread INFO toil.leader: Job ended: 'file:///home/ubuntu/cwl-tutorial/echo-tool.cwl' echo kind-file_home_ubuntu_cwl-tutorial_echo-tool.cwl/instancegehkd1py
ubuntu 2020-03-03 08:51:08,007 MainThread INFO toil.leader: Finished toil run successfully.
{
    "output": {
        "location": "file:///home/ubuntu/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "nameroot": "output",
        "nameext": ".txt",
        "class": "File",
        "checksum": "sha1$47a013e660d408619d894b20806b1d5086aab03b",
        "size": 13
    }
}ubuntu 2020-03-03 08:51:08,017 MainThread INFO toil.common: Successfully deleted the job store: FileJobStore(/tmp/tmp5_bvvekm)
```
```console
$ cat output.txt 
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

---

## Essential Input Parameters

`inp.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  example_flag:
    type: boolean
    inputBinding:
      position: 1
      prefix: -f
  example_string:
    type: string
    inputBinding:
      position: 3
      prefix: --example-string
  example_int:
    type: int
    inputBinding:
      position: 2
      prefix: -i
      separate: false
  example_file:
    type: File?
    inputBinding:
      prefix: --file=
      separate: false
      position: 4

outputs: []
```

* `baseCommand`: command to run.
* `type: boolean`: flag is present or absent.
* `inputBinding`: if missing then the parameter does not appear on the command-line..
* `separate: false`: no space between parameter and argument.
* `type: <TYPE>?`: optional parameter.
* CWL runners create temporary directories with symbolic (soft) links to input files.
* Input files are read-only.
* Input files must be copied to the output directory if they are to be updated.
* Positions are relative, default 0.

`inp-job.yml`:

```
example_flag: true
example_string: hello
example_int: 42
example_file:
  class: File
  path: whale.txt
```

```console
$ touch whale.txt
$ cwl-runner inp.cwl inp-job.yml
INFO [job inp.cwl] /tmp/tmpzrSnfX$ echo \
    -f \
    -i42 \
    --example-string \
    hello \
    --file=/tmp/tmpRBSHIG/stg979b6d24-d50a-47e3-9e9e-90097eed2cbc/whale.txt
-f -i42 --example-string hello --file=/tmp/tmpRBSHIG/stg979b6d24-d50a-47e3-9e9e-90097eed2cbc/whale.txt
```

---

## Returning Output Files

`tar.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, --extract]
inputs:
  tarfile:
    type: File
    inputBinding:
      prefix: --file
outputs:
  example_out:
    type: File
    outputBinding:
      glob: hello.txt
```

* Working directory is the designated output directory.
* Tool must write files to output directory.
* Output parameters:
  - Output files.
  - Content of output files.
* `baseCommand`: can be an array with mandatory command-line parameters.
* `glob`: name of output file. If unknown, use a wildcard e.g. `'*.txt'`.

`tar-job.yml`:

```
tarfile:
  class: File
  path: hello.tar
```

```console
$ touch hello.txt && tar -cvf hello.tar hello.txt
$ rm hello.txt
$ cwl-runner tar.cwl tar-job.yml
INFO [job tar.cwl] /tmp/tmpqOeawQ$ tar \
    --extract --file \
    /tmp/tmpGDk8Y1/stg80bbad20-494d-47af-8075-dffc32df03a3/hello.tar
INFO [job tar.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/hello.txt",
        "basename": "hello.txt",
        "class": "File",
        "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
        "size": 0,
        "path": "/home/ubuntu/workflows/cwl-tutorial/hello.txt"
    }
}
$ ls hello.txt 
hello.txt
```

---

## Capturing Standard Output

`stdout.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  message:
    type: string
    inputBinding:
      position: 1
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* Output `type: stdout` specifies standard output.

`echo-job.yml`:

```
message: Hello world!
```

```console
$ cwl-runner stdout.cwl echo-job.yml
INFO [job stdout.cwl] /tmp/earzxjq6$ echo \
    'Hello world!' > /tmp/earzxjq6/output.txt
INFO [job stdout.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$47a013e660d408619d894b20806b1d5086aab03b",
        "size": 13,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
Hello world!
```

---

## Parameter References

`tar-param.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, --extract]
inputs:
  tarfile:
    type: File
    inputBinding:
      prefix: --file
  extractfile:                    # Added to tar.cwl
    type: string
    inputBinding:
      position: 1
outputs:
  extracted_file:
    type: File
    outputBinding:
      glob: $(inputs.extractfile) # Changed from tar.cwl
```

* `inputs`: input object provided when the CWL tool is invoked.
* Re-use parameter values in another location, via JavaScript notation:
   - `$(inputs.extractfile)`
   - `$(inputs["extractfile"])`
   - `$(inputs['extractfile'])`
* `$(inputs.<FILE>.path)`: reference a path to an input file e.g. `$(inputs.tarfile.path)`.
* Parameter references are only allowed within certain fields. See [Parameter References](https://www.commonwl.org/user_guide/06-params/index.html).

`tar-param-job.yml`:

```
tarfile:
  class: File
  path: hello.tar
extractfile: goodbye.txt          # Added to tar-job.yml
```

```console
$ rm hello.tar
$ touch hello.txt goodbye.txt
$ tar -cvf hello.tar goodbye.txt
$ rm hello.txt goodbye.txt
$ cwl-runner tar-param.cwl tar-param-job.yml
INFO [job tar-param.cwl] /tmp/oxpxkhb7$ tar \
    --extract \
    --file \
    /tmp/tmpyq8up2vw/stg3ad2c230-5309-4e59-b0b5-4a0bb19348d0/hello.tar \
    goodbye.txt
INFO [job tar-param.cwl] completed success
{
    "extracted_file": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/goodbye.txt",
        "basename": "goodbye.txt",
        "class": "File",
        "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
        "size": 0,
        "path": "/home/ubuntu/workflows/cwl-tutorial/goodbye.txt"
    }
}
```

---

## Running Tools Inside Docker

* Ensure input files are accessible within container.
* Ensure output files are recovered from the container.

node.js program:

```console
$ echo "console.log(\"Hello World\");" > hello.js
$ node hello.js 
Hello World
```

`docker.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: node
hints:
  DockerRequirement:      # Docker is required.
    dockerPull: node:slim # Argument to `docker pull` command i.e. container image.
inputs:
  src:
    type: File
    inputBinding:
      position: 1
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* `hints`: hints to run tool.
* `DockerRequirement`: Docker is required.
* `dockerPull`: Argument to `docker pull` command i.e. container image.

`docker-job.yml`:

```
src:
  class: File
  path: hello.js
```

```console
$ cwl-runner docker.cwl docker-job.yml
INFO ['docker', 'pull', 'node:slim']
slim: Pulling from library/node
...
Status: Downloaded newer image for node:slim
docker.io/library/node:slim
INFO [job docker.cwl] /tmp/3vdcbzza$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/3vdcbzza,target=/oukzmQ \
    --mount=type=bind,source=/tmp/bhs2v1zp,target=/tmp \
    --mount=type=bind,source=/home/ubuntu/workflows/cwl-tutorial/hello.js,target=/var/lib/cwl/stga4b6a901-fb3a-4ea7-8a5f-1c561713d981/hello.js,readonly \
    --workdir=/oukzmQ \
    --read-only=true \
    --log-driver=none \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/oukzmQ \
    --cidfile=/tmp/wkq8n7vn/20200305084719-143341.cid \
    node:slim \
    node \
    /var/lib/cwl/stga4b6a901-fb3a-4ea7-8a5f-1c561713d981/hello.js > /tmp/3vdcbzza/output.txt
INFO [job docker.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$648a6a6ffffdaa0badb23b8baf90b6168dd16b3a",
        "size": 12,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
Hello World
```

* `/home/ubuntu/workflows/cwl-tutorial/hello.js,target=/var/lib/cwl/stga4b6a901-fb3a-4ea7-8a5f-1c561713d981/hello.js` maps path to script outside and inside container.

Run outwith container:

```console
$ cwl-runner --no-container docker.cwl  docker-job.yml 
INFO [job docker.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$648a6a6ffffdaa0badb23b8baf90b6168dd16b3a",
        "size": 12,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
```

---

## Additional Arguments and Parameters

`arguments.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: Example trivial wrapper for Java 9 compiler
hints:
  DockerRequirement:
    dockerPull: openjdk:9.0.1-11-slim
baseCommand: javac
arguments: ["-d", $(runtime.outdir)]
inputs:
  src:
    type: File
    inputBinding:
      position: 1
outputs:
  classfile:
    type: File
    outputBinding:
      glob: "*.class"
```

* `arguments`: command line options that do not correspond exactly to input parameters.
* Can also be used for fixed, not user-configurable, arguments.
* `$(runtime.outdir)`: runtime parameter, path to designated output directory.
* `$(runtime.tmpdir)`
* `$(runtime.ram)`
* `$(runtime.cores)`
* `$(runtime.outdirSize)`
* `$(runtime.tmpdirSize)`

`arguments-job.yml`:

```
src:
  class: File
  path: Hello.java
```

```console
$ echo "public class Hello {}" > Hello.java
$ cwl-runner arguments.cwl arguments-job.yml
 javac \
 -d \
 /var/spool/cwl \
 /var/lib/cwl/stg8939ac04-7443-4990-a518-1855b2322141/Hello.java
```

---

## Array Inputs

* Multiple values for a single argument.

`array-inputs.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  filesA:
    type: string[]
    inputBinding:
      prefix: -A
      position: 1
  filesB:
    type:
      type: array
      items: string
      inputBinding:
        prefix: -B=
        separate: false
    inputBinding:
      position: 2
  filesC:
    type: string[]
    inputBinding:
      prefix: -C=
      itemSeparator: ","
      separate: false
      position: 4
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* `inputBinding` associated with `<TYPE>[]`: array parameter, single parameter whose value is all items of the array.
* `inputBinding` associated with `type:array` and `items:<TYPE>`: array element, N parameters, one for each item of the array.
* `itemSeparator`: how items should be separated.

`array-inputs-job.yml`

```
filesA: [one, two, three]
filesB: [four, five, six]
filesC: [seven, eight, nine]
```

* Either:

```
[<ITEM>, <ITEM>]
```

* Or:

```
- <ITEM>
  <ITEM>
```

```console
$ cwl-runner array-inputs.cwl array-inputs-job.yml 
INFO [job array-inputs.cwl] /tmp/inw6o2et$ echo \
    -A \
    one \
    two \
    three \
    -B=four \
    -B=five \
    -B=six \
    -C=seven,eight,nine > /tmp/inw6o2et/output.txt
$ cat output.txt 
-A one two three -B=four -B=five -B=six -C=seven,eight,nine
```

---

## Array Outputs

* Tools can produce more than one output file.

`array-outputs.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: touch
inputs:
  touchfiles:
    type:
      type: array
      items: string
    inputBinding:
      position: 1
outputs:
  output:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.txt"
```

`array-outputs-job.yml`:

```
touchfiles:
  - foo.txt
  - bar.dat
  - baz.txt
```

```console
$ cwl-runner array-outputs.cwl array-outputs-job.yml 
INFO [job array-outputs.cwl] /tmp/d7xwtbxz$ touch \
    foo.txt \
    bar.dat \
    baz.txt
INFO [job array-outputs.cwl] completed success
{
    "output": [
        {
            "location": "file:///home/ubuntu/workflows/cwl-tutorial/baz.txt",
            "basename": "baz.txt",
            "class": "File",
            "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
            "size": 0,
            "path": "/home/ubuntu/workflows/cwl-tutorial/baz.txt"
        },
        {
            "location": "file:///home/ubuntu/workflows/cwl-tutorial/foo.txt",
            "basename": "foo.txt",
            "class": "File",
            "checksum": "sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709",
            "size": 0,
            "path": "/home/ubuntu/workflows/cwl-tutorial/foo.txt"
        }
    ]
}
```

* Use `glob` to capture multiple output files.

---

## Advanced Inputs

* Constraints on parameter combinations.

`record.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  dependent_parameters:
    type:
      type: record
      name: dependent_parameters
      fields:
        itemA:
          type: string
          inputBinding:
            prefix: -A
        itemB:
          type: string
          inputBinding:
            prefix: -B
  exclusive_parameters:
    type:
      - type: record
        name: itemC
        fields:
          itemC:
            type: string
            inputBinding:
              prefix: -C
      - type: record
        name: itemD
        fields:
          itemD:
            type: string
            inputBinding:
              prefix: -D
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* `type: record`: group parameters with `fields` for dependent parameters.
* `type:` list of `type: record`: `fields` for exclusive parameters.

`record-job1.yml`:

```
dependent_parameters:
  itemA: one
exclusive_parameters:
  itemC: three
```

```console
$ cwl-runner record.cwl record-job1.yml 
ERROR Workflow error, try again with --debug for more information:
Invalid job input record:
record-job1.yml:1:1: the `dependent_parameters` field is not valid because
                       missing required field `itemB`
```

Expected!

`record-job2.yml`:

```
dependent_parameters:
  itemA: one
  itemB: two
exclusive_parameters:
  itemC: three
  itemD: four
```

```console
$ cwl-runner record.cwl record-job2.yml 
record-job2.yml:6:3: invalid field `itemD`, expected one of: 'itemC'
INFO [job record.cwl] /tmp/wm7u5thj$ echo \
    -A \
    one \
    -B \
    two \
    -C \
    three > /tmp/wm7u5thj/output.txt
INFO [job record.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$329fe3b598fed0dfd40f511522eaf386edb2d077",
        "size": 23,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
-A one -B two -C three
```

`itemC` is added and `itemD` is ignored.

`record-job3.yml`:

```
dependent_parameters:
  itemA: one
  itemB: two
exclusive_parameters:
  itemD: four
```

```console
$ cwl-runner record.cwl record-job3.yml 
INFO [job record.cwl] /tmp/yis6xv8k$ echo \
    -A \
    one \
    -B \
    two \
    -D \
    four > /tmp/yis6xv8k/output.txt
INFO [job record.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$77f572b28e441240a5e30eb14f1d300bcc13a3b4",
        "size": 22,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
-A one -B two -D four
```

---

## Environment Variables

`env.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: env
requirements:
  EnvVarRequirement:
    envDef:
      HELLO: $(inputs.message)
inputs:
  message: string
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* `requirements`: requirements to run tool.
* Tools run within restricted environments.
* Minimal environment variables are inherited from parent process.
* `EnvVarRequirement`: required environment variables.
* `envDef`: define environment variable for restricted environment.

```console
$ cwl-runner env.cwl echo-job.yml 
INFO [job env.cwl] /tmp/7es908on$ env > /tmp/7es908on/output.txt
INFO [job env.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$a0f3d81933fafe564b321aac8a4d5e7f602cccb3",
        "size": 379,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
HELLO=Hello world!
PATH=/home/ubuntu/miniconda3/bin:/home/ubuntu/miniconda3/bin:/home/ubuntu/bowtie-1.2.2-linux-x86_64:/home/ubuntu/hisat2-2.1.0:/home/ubuntu/miniconda3/envs/riboviz-cwltool-snakemake/bin:/home/ubuntu/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
HOME=/tmp/7es908on
TMPDIR=/tmp/tmp_9pr4nyx
```

---

## JavaScript Expressions

`expression.cwl`:

```console
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
requirements:
  InlineJavascriptRequirement: {}
inputs: []
arguments:
  - prefix: -A
    valueFrom: $(1+1)
  - prefix: -B
    valueFrom: $("/foo/bar/baz".split('/').slice(-1)[0])
  - prefix: -C
    valueFrom: |
      ${
        var r = [];
        for (var i = 10; i >= 1; i--) {
          r.push(i);
        }
        return r;
      }
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* Dynamic value creation.
* `InlineJavascriptRequirement`: using inline JavaScript expressions.
* Provide JavaScript anywhere where a parameter reference is legal.
* Only use when absolutely necessary, where no CWL solution exists.
* Built-in File properties e.g. `basename`, `nameroot`, `nameext` etc.

`empty.yml`:

```
{}
```

* No inputs, empty JSON object (not YAML!)

```console
$ cwl-runner expression.cwl empty.yml
INFO [job expression.cwl] /tmp/mvzfkb1k$ echo \
    -A \
    2 \
    -B \
    baz \
    -C \
    10 \
    9 \
    8 \
    7 \
    6 \
    5 \
    4 \
    3 \
    2 \
    1 > /tmp/mvzfkb1k/output.txt
INFO [job expression.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$a739a6ff72d660d32111265e508ed2fc91f01a7c",
        "size": 36,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
-A 2 -B baz -C 10 9 8 7 6 5 4 3 2 1
```

---

## Creating Files at Runtime

`createfile.cwl`:

```
class: CommandLineTool
cwlVersion: v1.0
baseCommand: ["sh", "example.sh"]
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: example.sh
        entry: |-
          PREFIX='Message is:'
          MSG="\${PREFIX} $(inputs.message)"
          echo \${MSG}
inputs:
  message: string
outputs:
  example_out:
    type: stdout
stdout: output.txt
```

* `InitialWorkDirRequirement`: listing of files to create in working directory.
* `|-`: YAML quoting syntax.

```console
$ cwl-runner createfile.cwl echo-job.yml 
INFO [job createfile.cwl] /tmp/0dunjry7$ sh \
    example.sh > /tmp/0dunjry7/output.txt
INFO [job createfile.cwl] completed success
{
    "example_out": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$9045abe4bd04dd8ccfe50c6ff61820b784b64aa7",
        "size": 25,
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
$ cat output.txt 
Message is: Hello world!
```

---

## Staging Input Files

`linkfile.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: openjdk:9.0.1-11-slim
baseCommand: javac
requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.src)
inputs:
  src:
    type: File
    inputBinding:
      position: 1
      valueFrom: $(self.basename)
outputs:
  classfile:
    type: File
    outputBinding:
      glob: "*.class"
```

* Input files are considered to be in a read-only directory separate from the output directory.
* `InitialWorkDirRequirement:` with `listing:` with `$(inputs.src)`: stage input files into the output directory.
* `valueFrom: $(self.basename)`: extract base name of input file from leading directory path (this can be ommitted and the behaviour is the same)

```console
$ cwl-runner linkfile.cwl arguments-job.yml
INFO [job linkfile.cwl] /tmp/hn4ynat9$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/hn4ynat9,target=/RuiFYY \
    --mount=type=bind,source=/tmp/_vrhdkez,target=/tmp \
    --mount=type=bind,source=/home/ubuntu/workflows/cwl-tutorial/Hello.java,target=/RuiFYY/Hello.java,readonly \
    --workdir=/RuiFYY \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/RuiFYY \
    --cidfile=/tmp/lfzax8a9/20200306032615-750056.cid \
    openjdk:9.0.1-11-slim \
    javac \
    Hello.java
INFO [job linkfile.cwl] Max memory used: 0MiB
INFO [job linkfile.cwl] completed success
{
    "classfile": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/Hello.class",
        "basename": "Hello.class",
        "class": "File",
        "checksum": "sha1$fdb876b40ad9ebc7fee873212e02d5940588642e",
        "size": 184,
        "path": "/home/ubuntu/workflows/cwl-tutorial/Hello.class"
    }
}
```

---

## File Formats

`metadata_example.cwl`:

```
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: An example tool demonstrating metadata.
baseCommand: [ wc, -l ]
inputs:
  aligned_sequences:
    type: File
    label: Aligned sequences in BAM format
    format: edam:format_2572
    inputBinding:
      position: 1
outputs:
  report:
    type: stdout
    format: edam:format_1964
    label: A text file that contains a line count
stdout: output.txt

$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl
```

* Document expected input and output file types.
* Enables type checking when creating parameter files.
* Recommended ontologies:
  - [EDAM](https://www.ebi.ac.uk/ols/ontologies/edam/terms?iri=http%3A%2F%2Fedamontology.org%2Fformat_1915)
  - [IANA](https://www.iana.org/assignments/media-types/media-types.xhtml)
  - Community or local ontologies.
* `format`: file format.
* `$namespaces`: see below.
* `$schemas`: see below.

`sample.yml`:

```
aligned_sequences:
    class: File
    format: http://edamontology.org/format_2572
    path: file-formats.bam
```

```console
$ wget https://github.com/common-workflow-language/user_guide/raw/gh-pages/_includes/cwl/16-file-formats/file-formats.bam
$ cwl-runner metadata_example.cwl sample.yml 
INFO [job metadata_example.cwl] /tmp/3cfzinmg$ wc \
    -l \
    /tmp/tmp2jh5__ev/stgaeb40084-86be-45c5-a55d-2863bd363cd9/file-formats.bam > /tmp/3cfzinmg/output.txt
INFO [job metadata_example.cwl] completed success
{
    "report": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/output.txt",
        "basename": "output.txt",
        "class": "File",
        "checksum": "sha1$8761be5bd7d26d6355f16c3dfdda557586f268b0",
        "size": 80,
        "format": "http://edamontology.org/format_1964",
        "path": "/home/ubuntu/workflows/cwl-tutorial/output.txt"
    }
}
```

---

## Metadata and Authorship

`metadata_example3.cwl`:

```
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: An example tool demonstrating metadata.
doc: Note that this is an example and the metadata is not necessarily consistent.
hints:
  ResourceRequirement:
    coresMin: 4
baseCommand: [ wc, -l ]
inputs:
  aligned_sequences:
    type: File
    label: Aligned sequences in BAM format
    format: edam:format_2572
    inputBinding:
      position: 1
outputs:
  report:
    type: stdout
    format: edam:format_1964
    label: A text file that contains a line count
stdout: output.txt

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-6130-1021
    s:email: mailto:dyuen@oicr.on.ca
    s:name: Denis Yuen
s:contributor:
  - class: s:Person
    s:identifier: http://orcid.org/0000-0002-7681-6415
    s:email: mailto:briandoconnor@gmail.com
    s:name: Brian O'Connor
s:citation: https://dx.doi.org/10.6084/m9.figshare.3115156.v2
s:codeRepository: https://github.com/common-workflow-language/common-workflow-language
s:dateCreated: "2016-12-13"
s:license: https://spdx.org/licenses/Apache-2.0 
s:keywords: edam:topic_0091 , edam:topic_0622
s:programmingLanguage: C

$namespaces:
 s: https://schema.org/
 edam: http://edamontology.org/

$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
 - http://edamontology.org/EDAM_1.18.owl
```

* Metadata about a tool or workflow.
* Ciration information.
* Extensions are not required for execution.
* Extensions fields can be prefixed by namespace prefixes.
* label`: Tool/input/output description.
* doc`: Tool documentation.
* `$namespaces`: namespace prefixes and URLs.
* `$schemas`: metadata schema URLs.

---

## Custom Types

`custom-types.cwl`:

```
#!/usr/bin/env cwl-runner 
cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMax: 1
    ramMin: 100  # just a default, could be lowered
  SchemaDefRequirement:
    types:
      - $import: biom-convert-table.yaml
hints:
  DockerRequirement:
    dockerPull: 'quay.io/biocontainers/biom-format:2.1.6--py27_0'
  SoftwareRequirement:
    packages:
      biom-format:
        specs: [ "https://doi.org/10.1186/2047-217X-1-7" ]
        version: [ "2.1.6" ]
baseCommand: [ biom, convert ]
arguments:
  - valueFrom: $(inputs.biom.nameroot).hdf5  
    prefix: --output-fp
  - --to-hdf5
inputs:
  biom:
    type: File
    format: edam:format_3746  # BIOM
    inputBinding:
      prefix: --input-fp
  table_type:
    type: biom-convert-table.yaml#table_type
    inputBinding:
      prefix: --table-type
  header_key:
    type: string?
    doc: |
      The observation metadata to include from the input BIOM table file when
      creating a tsv table file. By default no observation metadata will be
      included.
    inputBinding:
      prefix: --header-key
outputs:
  result:
    type: File
    outputBinding: { glob: "$(inputs.biom.nameroot)*" }

$namespaces:
  edam: http://edamontology.org/
  s: https://schema.org/
$schemas:
  - http://edamontology.org/EDAM_1.16.owl
  - https://schema.org/docs/schema_org_rdfa.html
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
```

* `SoftwareRequirement`: required software (metadata).
* `ResourceRequirement`: required resources (metadata).
* `SchemaDefRequirement:`, `types:`, `$import: biom-convert-table.yaml`: import custom types.
* `type: biom-convert-table.yaml#table_type`: custom type reference (name of file, name of object).

`biom-convert-table.yaml`:

```
type: enum
name: table_type
label: The type of the table to produce
symbols:
  - OTU table
  - Pathway table
  - Function table
  - Ortholog table
  - Gene table
  - Metabolite table
  - Taxon table
  - Table
```

`custom-types.yml`:

```
biom:
    class: File
    format: http://edamontology.org/format_3746
    path: rich_sparse_otu_table.biom
table_type: OTU table
```

```console
$ wget https://raw.githubusercontent.com/common-workflow-language/user_guide/gh-pages/_includes/cwl/19-custom-types/rich_sparse_otu_table.biom
```

Command run would be:

```console
$ biom convert --output-fp rich_sparse_otu_table.hdf5 --to-hdf5 --input-fp rich_sparse_otu_table.biom --table-type "OTU table" --header-key "The obs..."
```

---

## Specifying Software Requirements

```
cwlVersion: v1.0
class: CommandLineTool
label: "InterProScan: protein sequence classifier"
doc: |
      Version 5.21-60 can be downloaded here:
      https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload

      Documentation on how to run InterProScan 5 can be found here:
      https://github.com/ebi-pf-team/interproscan/wiki/HowToRun
requirements:
  ResourceRequirement:
    ramMin: 10240
    coresMin: 3
  SchemaDefRequirement:
    types:
      - $import: InterProScan-apps.yml
hints:
  SoftwareRequirement:
    packages:
      interproscan:
        specs: [ "https://identifiers.org/rrid/RRID:SCR_005829" ]
        version: [ "5.21-60" ]
inputs:
  proteinFile:
    type: File
    inputBinding:
      prefix: --input
  applications:
    type: InterProScan-apps.yml#apps[]?
    inputBinding:
      itemSeparator: ','
      prefix: --applications
baseCommand: interproscan.sh
arguments:
 - valueFrom: $(inputs.proteinFile.nameroot).i5_annotations
   prefix: --outfile
 - valueFrom: TSV
   prefix: --formats
 - --disable-precalc
 - --goterms
 - --pathways
 - valueFrom: $(runtime.tmpdir)
   prefix: --tempdir
outputs:
  i5Annotations:
    type: File
    format: iana:text/tab-separated-values
    outputBinding:
      glob: $(inputs.proteinFile.nameroot).i5_annotations
$namespaces:
 iana: https://www.iana.org/assignments/media-types/
 s: https://schema.org/
$schemas:
 - https://schema.org/docs/schema_org_rdfa.html
s:license: https://spdx.org/licenses/Apache-2.0
s:copyrightHolder: "EMBL - European Bioinformatics Institute"
```

* `hints`, `SoftwareRequirement` to indicate tool version.
* Some CWL runners can use these to check that the software is available.
* `cwl-runner` has a dependency resolvers configuration to do this.
* References:
  - Tool research resource identifiers (RRIDs). [SciCrunch](https://scicrunch.org/) registry - find, track, add ([Adding a Resource](https://scicrunch.org/page/tutorials/336)), refer to scientific resources consistently.
  - DOIs.
  - URLs.

---

## Writing Workflows

`extract-compile-workflow.cwl`:

```
#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow
inputs:
  tarball: File
  name_of_file_to_extract: string
outputs:
  compiled_class:
    type: File
    outputSource: compile/classfile        # Connects to "compile" "out".
steps:
  untar:
    run: tar-param.cwl                     # Tool description
    in:
      tarfile: tarball                     # Connects to "inputs".
      extractfile: name_of_file_to_extract # Connects to "inputs".
    out: [extracted_file]
  compile:
    run: arguments.cwl                     # Tool description
    in:
      src: untar/extracted_file            # Connects to "untar" "out".
    out: [classfile]
```

* Every step has a CWL description.
* Execution order is determined by the dependencies between steps.
* Workflow steps that do not depend on each other may run in parallel.

`extract-compile-workflow-job.yml`:

```
tarball:
  class: File
  path: hello-java.tar
name_of_file_to_extract: Hello.java
```

```console
$ tar -cf hello-java.tar Hello.java
$ cwl-runner extract-compile-workflow.cwl extract-compile-workflow-job.yml 
INFO [workflow ] start
INFO [workflow ] starting step untar
INFO [step untar] start
INFO [job untar] /tmp/kzd4wlqy$ tar \
    --extract \
    --file \
    /tmp/tmp05x8xyxy/stg81d6ee81-705c-42cc-ae3d-33495579033f/hello-java.tar \
    Hello.java
INFO [job untar] completed success
INFO [step untar] completed success
INFO [workflow ] starting step compile
INFO [step compile] start
INFO [job compile] /tmp/enpym3ef$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/enpym3ef,target=/gdJBEn \
    --mount=type=bind,source=/tmp/57oskxo6,target=/tmp \
    --mount=type=bind,source=/tmp/kzd4wlqy/Hello.java,target=/var/lib/cwl/stgc267a61c-18a6-4e66-a6c2-128611a8e067/Hello.java,readonly \
    --workdir=/gdJBEn \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/gdJBEn \
    --cidfile=/tmp/fqeegch7/20200306042038-147727.cid \
    openjdk:9.0.1-11-slim \
    javac \
    -d \
    /gdJBEn \
    /var/lib/cwl/stgc267a61c-18a6-4e66-a6c2-128611a8e067/Hello.java
INFO [job compile] Max memory used: 0MiB
INFO [job compile] completed success
INFO [step compile] completed success
INFO [workflow ] completed success
{
    "compiled_class": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/Hello.class",
        "basename": "Hello.class",
        "class": "File",
        "checksum": "sha1$fdb876b40ad9ebc7fee873212e02d5940588642e",
        "size": 184,
        "path": "/home/ubuntu/workflows/cwl-tutorial/Hello.class"
    }
}
```

---

## Nested Workflows

`nested-workflows.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}
inputs: []
outputs:
  classout:
    type: File
    outputSource: compile/compiled_class
steps:
  compile:
    run: extract-compile-workflow.cwl
    in:
      tarball: create-tar/tar_compressed_java_file
      name_of_file_to_extract:
        default: "Hello.java"
    out: [compiled_class]
  create-tar:
    in: []
    out: [tar_compressed_java_file]
    run:
      class: CommandLineTool
      requirements:
        InitialWorkDirRequirement:
          listing:
            - entryname: Hello.java
              entry: |
                public class Hello {
                  public static void main(String[] argv) {
                      System.out.println("Hello from Java");
                  }
                }
      inputs: []
      baseCommand: [tar, --create, --file=hello.tar, Hello.java]
      outputs:
        tar_compressed_java_file:
          type: File
          streamable: true
          outputBinding:
            glob: "hello.tar"
```

* `SubworkflowFeatureRequirement`: use workflow as step in another workflow.
* Conditional on CWL runner supporting this.

* CWL workflows can be used just like tool descriptions.

```console
$ cwl-runner nested-workflows.cwl empty.yml 
INFO [workflow ] start
INFO [workflow ] starting step create-tar
INFO [step create-tar] start
INFO [job create-tar] /tmp/_8m7m_6v$ tar \
    --create \
    --file=hello.tar \
    Hello.java
INFO [job create-tar] completed success
INFO [step create-tar] completed success
INFO [workflow ] starting step compile
INFO [step compile] start
INFO [workflow compile] start
INFO [workflow compile] starting step untar
INFO [step untar] start
INFO [job untar] /tmp/llg4ff98$ tar \
    --extract \
    --file \
    /tmp/tmp1ih74ioz/stg2faa9fe7-6a61-41fa-a173-554342a894a2/hello.tar \
    Hello.java
INFO [job untar] completed success
INFO [step untar] completed success
INFO [workflow compile] starting step compile_2
INFO [step compile_2] start
INFO [job compile] /tmp/05u3hn3_$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/05u3hn3_,target=/zMtEbQ \
    --mount=type=bind,source=/tmp/ibo168ho,target=/tmp \
    --mount=type=bind,source=/tmp/llg4ff98/Hello.java,target=/var/lib/cwl/stgace2d18a-4a08-4036-a364-d1dc1352e3ed/Hello.java,readonly \
    --workdir=/zMtEbQ \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/zMtEbQ \
    --cidfile=/tmp/h89i7f_z/20200306042901-149120.cid \
    openjdk:9.0.1-11-slim \
    javac \
    -d \
    /zMtEbQ \
    /var/lib/cwl/stgace2d18a-4a08-4036-a364-d1dc1352e3ed/Hello.java
INFO [job compile] Max memory used: 0MiB
INFO [job compile] completed success
INFO [step compile_2] completed success
INFO [workflow compile] completed success
INFO [step compile] completed success
INFO [workflow ] completed success
{
    "classout": {
        "location": "file:///home/ubuntu/workflows/cwl-tutorial/Hello.class",
        "basename": "Hello.class",
        "class": "File",
        "checksum": "sha1$39e3219327347c05aa3e82236f83aa6d77fe6bfd",
        "size": 419,
        "path": "/home/ubuntu/workflows/cwl-tutorial/Hello.class"
    }
}
INFO Final process status is success
```

* `>`: block, similar to `|`, but within which newlines are stripped - specify single line command across multiple lines.
* `default`: specify default value for a field, which can be overwritten by a value in the input object.

---

## Scattering Workflows

`scatter-workflow.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
inputs:
  message_array: string[]    # Input is array of string.
steps:
  echo:
    run: echo-tool.cwl
    scatter: message         # Scatter/repeat/parallelise over this input.
    in:
      message: message_array
    out: []
outputs: []
```

`echo-tool.cwl`:

```
#!/usr/bin/env cwl-runner

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

* Run tools or workflows in parallel.
* `ScatterFeatureRequirement`: run a tool or workflow multiple times over a list of inputs.
* Conditional on CWL runner supporting this.
* Specify `scatter` field for each step to scatter, referencing step inputs, not workflow inputs.
* `scatter` runs on each step independently.

`scatter-job.yml`:

```
message_array: 
  - Hello world!
  - Hola mundo!
  - Bonjour le monde!
  - Hallo welt!
```

```console
$ cwl-runner scatter-workflow.cwl scatter-job.yml
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool-snakemake/bin/cwl-runner 2.0.20200303141624
INFO Resolved 'scatter-workflow.cwl' to 'file:///home/ubuntu/workflows/cwl-tutorial/scatter-workflow.cwl'
INFO [workflow ] start
INFO [workflow ] starting step echo
INFO [step echo] start
INFO [job echo] /tmp/gw9e003g$ echo \
    'Hello world!' > /tmp/gw9e003g/output.txt
INFO [job echo] completed success
INFO [step echo] start
INFO [job echo_2] /tmp/hthkhl87$ echo \
    'Hola mundo!' > /tmp/hthkhl87/output.txt
INFO [job echo_2] completed success
INFO [step echo] start
INFO [job echo_3] /tmp/ulvj0791$ echo \
    'Bonjour le monde!' > /tmp/ulvj0791/output.txt
INFO [job echo_3] completed success
INFO [step echo] start
INFO [job echo_4] /tmp/sd76qvxf$ echo \
    'Hallo welt!' > /tmp/sd76qvxf/output.txt
INFO [job echo_4] completed success
INFO [step echo] completed success
INFO [workflow ] completed success
{}
```

`wc-tool.cwl`:

```
#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: wc
arguments: ["-c"]
inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
outputs: []
```

Edit `scatter-workflow.cwl`:

```
steps:
  echo:
    run: echo-tool.cwl
    scatter: message
    in:
      message: message_array
    out: [output]
  wc:
    run: wc-tool.cwl
    scatter: input_file
    in:
      input_file: echo/echoput
    out: []
outputs: []
```

```console
$ cwl-runner scatter-workflow.cwl scatter-job.yml 
INFO [workflow ] start
INFO [workflow ] starting step echo
...
INFO [step echo] completed success
INFO [workflow ] starting step wc
INFO [step wc] start
INFO [job wc] /tmp/u3idx1dj$ wc \
    -c \
    /tmp/tmpycq0f9a1/stgaad4b5b7-871c-4d80-81c1-ef41719456cf/output.txt
13 /tmp/tmpycq0f9a1/stgaad4b5b7-871c-4d80-81c1-ef41719456cf/output.txt
INFO [job wc] completed success
INFO [step wc] start
INFO [job wc_2] /tmp/1toxjh62$ wc \
...
INFO [job wc_4] completed success
INFO [step wc] completed success
INFO [workflow ] completed success
{}
```

Use `scatter` over a subworkflow...

`scatter-nested-workflow.cwl`:

```
!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
 ScatterFeatureRequirement: {}
 SubworkflowFeatureRequirement: {}
inputs:
  message_array: string[] 
steps:
  subworkflow:
    run: 
      class: Workflow
      inputs: 
        message: string
      outputs: []
      steps:
        echo:
          run: 1st-tool-mod.cwl
          in:
            message: message
          out: [echo_out]
        wc:
          run: wc-tool.cwl
          in:
            input_file: echo/echo_out
          out: []
    scatter: message
    in: 
      message: message_array
    out: []
outputs: []
```

```console
$ cwl-runner scatter-nested-workflow.cwl scatter-job.yml
INFO [workflow ] start
INFO [workflow ] starting step subworkflow
INFO [step subworkflow] start
INFO [workflow subworkflow] start
INFO [workflow subworkflow] starting step echo
INFO [step echo] start
INFO [job echo] /tmp/sf6_cxkg$ echo \
    'Hello world!' > /tmp/sf6_cxkg/output.txt
INFO [job echo] completed success
INFO [step echo] completed success
INFO [workflow subworkflow] starting step wc
INFO [step wc] start
INFO [job wc] /tmp/wpelyjt1$ wc \
    -c \
    /tmp/tmphn9bqeez/stg763fb294-e989-479c-9563-3f5d4b12c455/output.txt
13 /tmp/tmphn9bqeez/stg763fb294-e989-479c-9563-3f5d4b12c455/output.txt
INFO [job wc] completed success
INFO [step wc] completed success
INFO [workflow subworkflow] completed success
INFO [step subworkflow] start
INFO [workflow subworkflow_2] start
...
{}
```
