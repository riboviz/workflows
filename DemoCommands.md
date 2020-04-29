# Demonstration

## Snakemake

The following assumes you have:

* Installed [riboviz](https://github.com/riboviz/riboviz).
* Cloned https://github.com/riboviz/workflows
* Set up a `riboviz-snakemake` conda environment and have installed Snakemake, as described in the workflows repository's Snakemake [README.md](https://github.com/riboviz/workflows/blob/master/snakemake/README.md).

Edit `vignette/vignette_config.yaml` and comment out `WT3AT`.

```console
$ conda activate riboviz-snakemake
$ cd riboviz
$ cp ~/workflows/snakemake/Snakefile .
```

Dry-run:

```console
$ snakemake --configfile vignette/vignette_config.yaml -n
```

Dry run and bash commands:

```console
$ snakemake --configfile vignette/vignette_config.yaml --printshellcmds -n
```

Run:

```console
$ snakemake --configfile vignette/vignette_config.yaml
```

Edit `vignette/vignette_config.yaml` and uncomment `WT3AT`.

Incremental build, dry-run:

```console
$ snakemake --configfile vignette/vignette_config.yaml -n
$ snakemake --configfile vignette/vignette_config.yaml
```

Bash script:

```console
$ snakemake --configfile vignette/vignette_config.yaml --detailed-summary
```

Clean:

```console
$ snakemake clean --configfile vignette/vignette_config.yaml
```

---

## Nextflow

The following assumes you have:

* Installed [riboviz](https://github.com/riboviz/riboviz).
* Cloned https://github.com/riboviz/workflows
* Set up a `riboviz-nextflow` conda environment and have installed Nextflow, as described in the workflows repository's Nextflow [README.md](https://github.com/riboviz/workflows/blob/master/nextflow/README.md)

Edit `vignette/vignette_config.yaml` and comment out `WT3AT`.

```console
$ conda activate riboviz-nextflow
$ cd riboviz
```

Run:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false
```

Explore work directory, choose one for cutadapt:

```console
$ ls -lA work/NN/HASH/
$ cat work/NN/HASH/.command.sh 
$ head work/NN/HASH/.command.out 
$ head work/NN/HASH/.command.err 
$ cat work/NN/HASH/.exitcode
$ ls vignette/tmp/WTnone
```

Edit `vignette/vignette_config.yaml` and uncomment `WT3AT`.

Incremental build, resume:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false -resume
```

Reporting:

```console
$ nextflow run prep_riboviz.nf -params-file vignette/vignette_config.yaml -ansi-log false -with-report report.html -with-timeline timeline.html -with-dag workflow.svg -resume
```
