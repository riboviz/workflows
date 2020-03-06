# Snakemake and CWL

Experimenting with Snakemake's `--export-cwl` option.

## Nature's Snakemake example

Perkel, J.M. "Workflow systems turn raw data into scientific knowledge", Nature 573, 149-150 (2019) doi: [10.1038/d41586-019-02619-z](https://doi.org/10.1038/d41586-019-026).

```console
$ git clone https://github.com/jperkel/Snakemake_example snakemake-nature
$ cd snakemake-nature
$ cat name.txt 
Jane Q. Scientist
$ snakemake -s hello_world.smk
$ cat hello-world.txt
HELLO, WORLD, MY NAME IS JANE Q. SCIENTIST!
$ ls -1
chunk_aa
chunk_ab
chunk_ac
chunk_ad
chunk_ae
chunk_af
chunk_ag
chunk_ah
hello_world.smk
hello-world.tmp
hello-world.txt
name.txt
README.md
$ snakemake -s hello_world.smk clean
```
```console
$ snakemake -s hello_world.smk --export-cwl hello_world.cwl
$ cat hello_world.cwl 
```

Inspecting `hello_world.cwl` shows that it looks like the exported CWL is a wrapper for Snakemake, rather than the workflow implemented within Snakemake converted into CWL. For example, the three shell commands run by Snakemake are not present in the CWL, only a single invocation of `snakemake` itself.

This was confirmed by the Snakemake docs, in [CWL export](https://snakemake.readthedocs.io/en/stable/executing/interoperability.html#cwl-export), which comments that:

> Since, CWL is less powerful for expressing workflows than Snakemake (most importantly Snakemake offers more flexible scatter-gather patterns, since full Python can be used), export works such that every Snakemake job is encoded into a single step in the CWL workflow.

and cautions that:

> when exporting keep in mind that the resulting CWL file can become huge, depending on the number of jobs in your workflow.

In addition the CWL file has references to every file in the example Git repository. Snakemake, in [CWL export](https://snakemake.readthedocs.io/en/stable/executing/interoperability.html#cwl-export) comment that:

> ...due to limitations in CWL, it seems currently impossible to avoid that all target files (output files of target jobs), are written directly to the workdir, regardless of their relative paths in the Snakefile.

## Execute Nature CWL using `cwltool`

Create conda environment, derived from my `riboviz-snakemake` conda environment (RiboViz dependencies plus Snakemake):

```console
$ conda create --name riboviz-cwltool-snakemake --clone riboviz-snakemake
$ conda activate riboviz-cwltool-snakemake
$ pip install cwlref-runner
$ cwl-runner --version
/home/ubuntu/miniconda3/envs/riboviz-cwltool-snakemake/bin/cwl-runner 2.0.20200303141624
$ cwltool --version
/home/ubuntu/miniconda3/envs/riboviz-cwltool-snakemake/bin/cwltool 2.0.20200303141624
```

Run:

```console
$ cwltool hello_world.cwl 
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool/bin/cwltool 2.0.20200224214940
INFO Resolved 'hello_world.cwl' to 'file:///home/ubuntu/explore-workflows/snakemake-nature/hello_world.cwl'
ERROR I'm sorry, I couldn't load this CWL file, try again with --debug for more information.
The error was: cwltool requires Node.js engine to evaluate and validate Javascript expressions, but couldn't find it.  Tried nodejs, node, docker run node:slim
```
```console
$ sudo apt-get -y install nodejs
```
```console
$ cwltool hello_world.cwl 
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool/bin/cwltool 2.0.20200224214940
INFO Resolved 'hello_world.cwl' to 'file:///home/ubuntu/explore-workflows/snakemake-nature/hello_world.cwl'
INFO [workflow ] start
INFO [workflow ] starting step job-3
INFO [step job-3] start
ERROR Workflow error, try again with --debug for more information:
Docker is not available for this tool, try --no-container to disable Docker, or install a user space Docker replacement like uDocker with --user-space-docker-cmd.: docker executable is not available
```

Install Docker (see [Docker](../Docker.md)) and try again:

```console
$ cwltool hello_world.cwl
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool-snakemake/bin/cwltool 2.0.20200303141624
INFO Resolved 'hello_world.cwl' to 'file:///home/ubuntu/explore-workflows/snakemake-nature/hello_world.cwl'
INFO [workflow ] start
INFO [workflow ] starting step job-3
INFO [step job-3] start
INFO [job job-3] /tmp/u8wed2h7$ docker \
    run \
    -i \
    --mount=type=bind,source=/tmp/u8wed2h7,target=/ozwKUF \
    --mount=type=bind,source=/tmp/pq_3of28,target=/tmp \
    --mount=type=bind,source=/home/ubuntu/explore-workflows/snakemake-nature/hello_world.smk,target=/var/lib/cwl/stgec578333-bd21-41f5-b9c4-228cf1f29b74/hello_world.smk,readonly \
    --mount=type=bind,source=/home/ubuntu/explore-workflows/snakemake-nature/README.md,target=/var/lib/cwl/stg962064c4-8a7e-45ba-9889-b286bc31d89b/README.md,readonly \
    --mount=type=bind,source=/home/ubuntu/explore-workflows/snakemake-nature/name.txt,target=/var/lib/cwl/stgceae048e-c1a6-47ec-80c9-ece052f45e92/name.txt,readonly \
    --workdir=/ozwKUF \
    --read-only=true \
    --user=1000:1000 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/ozwKUF \
    --cidfile=/tmp/o0fj9jpr/20200305061551-959859.cid \
    snakemake/snakemake:v5.10.0 \
    snakemake \
    --force \
    --keep-target-files \
    --keep-remote \
    --force-use-threads \
    --wrapper-prefix \
    https://github.com/snakemake/snakemake-wrappers/raw/ \
    --notemp \
    --quiet \
    --use-conda \
    --no-hooks \
    --nolock \
    --mode \
    1 \
    --cores \
    1 \
    --allowed-rules \
    helloworld \
    --snakefile \
    /var/lib/cwl/stgec578333-bd21-41f5-b9c4-228cf1f29b74/hello_world.smk \
    hello-world.tmp
FileNotFoundError in line 17 of /var/lib/cwl/stgec578333-bd21-41f5-b9c4-228cf1f29b74/hello_world.smk:
[Errno 2] No such file or directory: 'name.txt'
  File "/var/lib/cwl/stgec578333-bd21-41f5-b9c4-228cf1f29b74/hello_world.smk", line 17, in <module>
INFO [job job-3] Max memory used: 0MiB
WARNING [job job-3] completed permanentFail
WARNING [step job-3] completed permanentFail
INFO [workflow ] completed permanentFail
{
    "job-0": null
}
WARNING Final process status is permanentFail
```

Run locally, without Docker:

```console
$ cwltool --no-container hello_world.cwl
...
FileNotFoundError in line 17 of /tmp/tmp8bn_jfse/stg35099ebb-24e9-4946-b354-e74e177997d8/hello_world.smk:
[Errno 2] No such file or directory: 'name.txt'
  File "/tmp/tmp8bn_jfse/stg35099ebb-24e9-4946-b354-e74e177997d8/hello_world.smk", line 17, in <module>
WARNING [job job-3] completed permanentFail
WARNING [step job-3] completed permanentFail
```

Set temporary directory so can see content:

```console
$ cwltool --no-container --tmpdir-prefix xxx --leave-tmpdir hello_world.cwl
...
FileNotFoundError in line 17 of /tmp/tmp_ggzm5j3/stg1a19e36c-d01b-4082-877f-218ceab6328d/hello_world.smk:
...
$ ls /tmp/tmp_ggzm5j3
ls: cannot access '/tmp/tmp_ggzm5j3': No such file or directory
```

There is an issue for this - [--tmpdir-prefix is ignored #1084](https://github.com/common-workflow-language/cwltool/issues/1084):

> --tmpdir-prefix is ignored when running a CWL tool that doesn't require a docker image and instead uses the prefix /tmp/tmp. It works as intended for tools with a DockerRequirement.

This should not be a difficult example!

I raised an issue on 06/03/2020, https://github.com/jperkel/Snakemake_example/issues/2.

## Execute Nature CWL using Toil

Create conda environment, derived from my `riboviz-snakemake` conda environment (RiboViz dependencies plus Snakemake):

```console
$ conda create --name riboviz-toil-snakemake --clone riboviz-snakemake
$ conda activate riboviz-toil-snakemake
$ pip install toil[cwl]
$ cwltool --version
/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/bin/cwltool 1.0.20190906054215
$ toil-cwl-runner --version
3.24.0
```
```console
$ toil-cwl-runner --no-container hello_world.cwl
ubuntu 2020-03-06 02:34:47,947 MainThread INFO cwltool: Resolved 'hello_world.cwl' to 'file:///home/ubuntu/workflows/snakemake-nature/hello_world.cwl'
ubuntu 2020-03-06 02:34:50,184 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxCores to CPU count of system (4).
ubuntu 2020-03-06 02:34:50,185 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxMemory to physically available memory (8340119552).
ubuntu 2020-03-06 02:34:50,185 MainThread WARNING toil.batchSystems.singleMachine: Limiting maxDisk to physically available disk (17484890112).
ubuntu 2020-03-06 02:34:50,208 MainThread INFO toil: Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
INFO:toil.worker:Redirecting logging to /tmp/toil-04eb1503-a19b-4586-9a48-77f0d909825b-453138a4fa34406f996060e86264aa60/tmpfnwdixxx/worker_log.txt
ubuntu 2020-03-06 02:34:50,837 MainThread WARNING toil.leader: The job seems to have left a log file, indicating failure: 'CWLWorkflow' kind-CWLWorkflow/instancesdwew19e
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    INFO:toil.worker:---TOIL WORKER OUTPUT LOG---
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    INFO:toil:Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    Traceback (most recent call last):
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/worker.py", line 366, in workerScript
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore, defer=defer)
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/job.py", line 1392, in _runner
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        returnValues = self._run(jobGraph, fileStore)
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/job.py", line 1329, in _run
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        return self.run(fileStore)
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 945, in run
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        if isinstance(jobobj[key][1], Promise):
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    TypeError: 'MergeInputsFlattened' object is not subscriptable
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    ERROR:toil.worker:Exiting the worker because of a failed job on host ubuntu
ubuntu 2020-03-06 02:34:50,838 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'CWLWorkflow' kind-CWLWorkflow/instancesdwew19e with ID kind-CWLWorkflow/instancesdwew19e to 1
INFO:toil.worker:Redirecting logging to /tmp/toil-04eb1503-a19b-4586-9a48-77f0d909825b-453138a4fa34406f996060e86264aa60/tmpoamqcbr2/worker_log.txt
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: The job seems to have left a log file, indicating failure: 'CWLWorkflow' kind-CWLWorkflow/instancesdwew19e
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    INFO:toil.worker:---TOIL WORKER OUTPUT LOG---
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    INFO:toil:Running Toil version 3.24.0-de586251cb579bcb80eef435825cb3cedc202f52.
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    Traceback (most recent call last):
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/worker.py", line 366, in workerScript
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        job._runner(jobGraph=jobGraph, jobStore=jobStore, fileStore=fileStore, defer=defer)
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/job.py", line 1392, in _runner
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        returnValues = self._run(jobGraph, fileStore)
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/job.py", line 1329, in _run
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        return self.run(fileStore)
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e      File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 945, in run
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e        if isinstance(jobobj[key][1], Promise):
ubuntu 2020-03-06 02:34:51,439 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    TypeError: 'MergeInputsFlattened' object is not subscriptable
ubuntu 2020-03-06 02:34:51,440 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    ERROR:toil.worker:Exiting the worker because of a failed job on host ubuntu
ubuntu 2020-03-06 02:34:51,440 MainThread WARNING toil.leader: kind-CWLWorkflow/instancesdwew19e    WARNING:toil.jobGraph:Due to failure we are reducing the remaining retry count of job 'CWLWorkflow' kind-CWLWorkflow/instancesdwew19e with ID kind-CWLWorkflow/instancesdwew19e to 0
ubuntu 2020-03-06 02:34:51,440 MainThread WARNING toil.leader: Job 'CWLWorkflow' kind-CWLWorkflow/instancesdwew19e with ID kind-CWLWorkflow/instancesdwew19e is completely failed
ubuntu 2020-03-06 02:34:55,227 MainThread INFO toil.leader: Finished toil run with 1 failed jobs.
ubuntu 2020-03-06 02:34:55,227 MainThread INFO toil.leader: Failed jobs at end of the run: 'CWLWorkflow' kind-CWLWorkflow/instancesdwew19e
Traceback (most recent call last):
  File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/bin/toil-cwl-runner", line 10, in <module>
    sys.exit(main())
  File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/cwl/cwltoil.py", line 1357, in main
    result = toil.start(wf1)
  File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/common.py", line 800, in start
    return self._runMainLoop(rootJobGraph)
  File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/common.py", line 1070, in _runMainLoop
    jobCache=self._jobCache).run()
  File "/home/ubuntu/miniconda3/envs/riboviz-toil-snakemake/lib/python3.7/site-packages/toil/leader.py", line 246, in run
    raise FailedJobsException(self.config.jobStore, self.toilState.totalFailedJobs, self.jobStore)
toil.leader.FailedJobsException
```
```console
$ toil-cwl-runner hello_world.cwl
...same error...
```

## Execute RiboViz CWL using `cwltool`

```console
$ snakemake --configfile vignette/vignette_config.yaml --export-cwl riboviz.cwl
$ cwltool --no-container riboviz.cwl
INFO /home/ubuntu/miniconda3/envs/riboviz-cwltool-snakemake/bin/cwltool 2.0.20200303141624
INFO Resolved 'riboviz.cwl' to 'file:///home/ubuntu/riboviz/riboviz.cwl'
INFO [workflow ] start
INFO [workflow ] starting step job-23
INFO [step job-23] start
INFO [job job-23] /tmp/biopa37t$ snakemake \
    --force \
    --keep-target-files \
    --keep-remote \
    --force-use-threads \
    --wrapper-prefix \
    https://github.com/snakemake/snakemake-wrappers/raw/ \
    --notemp \
    --quiet \
    --use-conda \
    --no-hooks \
    --nolock \
    --mode \
    1 \
    --cores \
    1 \
    --allowed-rules \
    build_indices_rrna \
    --snakefile \
    /tmp/tmphaiub_o9/stgfb8379f6-a531-42bc-9005-7c7fafa5610b/Snakefile \
    vignette/index/yeast_rRNA.1.ht2 \
    vignette/index/yeast_rRNA.2.ht2 \
    vignette/index/yeast_rRNA.3.ht2 \
    vignette/index/yeast_rRNA.4.ht2 \
    vignette/index/yeast_rRNA.5.ht2 \
    vignette/index/yeast_rRNA.6.ht2 \
    vignette/index/yeast_rRNA.7.ht2 \
    vignette/index/yeast_rRNA.8.ht2
KeyError in line 7 of /tmp/tmphaiub_o9/stgfb8379f6-a531-42bc-9005-7c7fafa5610b/Snakefile:
'fq_files'
  File "/tmp/tmphaiub_o9/stgfb8379f6-a531-42bc-9005-7c7fafa5610b/Snakefile", line 7, in <module>
WARNING [job job-23] completed permanentFail
WARNING [step job-23] completed permanentFail
INFO [workflow ] completed permanentFail
{
    "job-0": null
}
WARNING Final process status is permanentFail
```

Cannot find configuration parameter.
