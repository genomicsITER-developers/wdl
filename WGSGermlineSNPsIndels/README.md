## Table of Contents

* [Prerequisites](#prerequisites)
* [Features](#features)
* [Usage](#usage)
    * [Set input parameters](#set-input-parameters)
    * [Build TSV for pair reads information](#build-tsv-for-pair-reads-information)
    * [Local](#local)
    * [HPC using Slurm scheduler](#hpc-using-slurm-scheduler)
* [Workflow](#workflow)
* [WDL Scripts](#wdl-scripts)

---

## Prerequisites

Basic software needed to run the pipeline:

* [GATK](https://software.broadinstitute.org/gatk/) (>= 4.0.0.0)
* [Cromwell](https://cromwell.readthedocs.io/en/stable/)
* [Womtool](https://cromwell.readthedocs.io/en/stable/WOMtool/) (Optional)
* [BWA](http://bio-bwa.sourceforge.net/)
* [Samtools, HTSlib](http://www.htslib.org/)
* [Picard](https://broadinstitute.github.io/picard/)
* [Python](https://www.python.org/)
* [Qualimap](http://qualimap.bioinfo.cipf.es/)

---

:warning:

---