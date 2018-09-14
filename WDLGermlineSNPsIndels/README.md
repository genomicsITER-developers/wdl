## Table of Contents

* [Prerequisites](#prerequisites) :warning:
* [Usage](#usage) :warning:
    * [Local](#local)
    * [HPC using Slurm scheduler](#hpc-using-slurm-scheduler)
* [Workflow](#workflow)
* [Features](#features)
* [WDL Scripts](#wdl-scripts) :warning:

---

## Prerequisites:

Basic software needed to run the pipeline:

* [GATK](https://software.broadinstitute.org/gatk/) (<= 4.0.0.0)
* [Cromwell](https://cromwell.readthedocs.io/en/stable/)
* [BWA](http://bio-bwa.sourceforge.net/)
* [Samtools, HTSlib](http://www.htslib.org/)
* [Picard](https://broadinstitute.github.io/picard/)
* [Python](https://www.python.org/)
* [Qualimap](http://qualimap.bioinfo.cipf.es/)

TODO :warning:

---

## Usage:

### Local

...

```
include required(classpath("application"))

backend {
  default="Local"
  providers {
    Local {
      config {
        filesystems {
          local {
            localization: [
              "soft-link", "copy", "hard-link"
            ]
          }
        }
      }
    }
  }
}
```

To run the workflow in your system using this configuration file:
```
java -Dconfig.file=./backends/backend-local.conf -jar cromwell-30.2.jar run ./pipeline/WholeGenomeSequencingGATK4.wdl -i ./inputs/WholeGenomeSequencingGATK4.inputs.json -m metadata.json
```

### HPC using Slurm scheduler

...

```
include required(classpath("application"))

# Backend for SLURM executions
backend {
  default = "Slurm"
  providers {
    # Slurm configuration
    Slurm {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
        String queue = "long"
        Int nodes = 1
        Int tasks_per_node = 16
        String runtime_minutes = "3-00:00:00"
        Int mem = 60000
        """

        submit = """
            sbatch \
            -J ${job_name} \
            -D ${cwd} \
            -o ${out} \
            -e ${err} \
            -p ${queue} \
            -N ${nodes} \
            --tasks-per-node ${tasks_per_node} \
            -t ${runtime_minutes} \
            --mem ${mem} \
            --wrap "/bin/bash ${script}"
        """

        # Int cpus = 8
        # Int requested_memory_mb_per_core = 1000
        # Other options
            #${"-n " + cpus} \
            #--mem-per-cpu=${requested_memory_mb_per_core} \

        kill = "scancel ${job_id}"
        check-alive = "squeue -j ${job_id}"
        job-id-regex = "Submitted batch job (\\d+).*"

        filesystems {
          local {
            localization: [
              "soft-link", "copy", "hard-link"
            ]
          }
        }
      }
    }
  }
}
```

To run the workflow in a system with Slurm scheduler (such as TeideHPC) using this configuration file:
```
screen -d -m -t wgs-wdl-pipeline -L sh ./slurm/run_pipeline.sh
```

TODO :warning:

---

## Workflow:

### A. Per sample, per lane

#### 1. Preprocessing
    1. FastqToSam....................................(step 1)
    2. SamToFastq....................................(step 2)
    3. Bwa Mem.......................................(step 3)
    4. MergeBamAlignments............................(step 4)
    5. MarkDuplicates................................(step 5)
    6. SortSam.......................................(step 6)
    7. SetNmMdAndUqTags..............................(step 7)

#### 2. Quality Control
    1. ValidateSam...................................(step 8)
    2. Qualimap......................................(step 9)
    3. CollectRawWgsMetrics (3-0)....................(step 10)
    4. CollectRawWgsMetrics (20-20)..................(step 11)
    5. CollectWgsMetricsWithNonZeroCoverage..........(step 12)
    6. ¿DepthOfCoverage?.............................(step 13)
    7. CollectMultipleMetrics........................(step 14)
    8. ¿CallableLoci?................................(step 15)
 
### B. Per sample
   
#### 1. Variant Calling
    1. MarkDuplicates (MergeBamsPerSample)............(step 16)
    2. BQSR
        1. BaseRecalibrator (by intervals)*...........(step 21)
        2. GatherBQSRReports..........................(step 22)
        3. ApplyBQSR (by intervals)*..................(step 23)
        4. GatherBamFiles.............................(step 24)
    3. HaplotypeCaller (by intervals)*................(step 25)
    4. MergeVcfs......................................(step 26)
 
#### 2. Quality Control (Con la salida del MergeBamsPerSample)
    1. ValidateSam....................................(step 17)
    2. ¿DepthOfCoverage?..............................(step 18)
    3. CollectMultipleMetrics.........................(step 19)
    4. ¿CallableLoci?.................................(step 20)
 
### C. Multi-sample
 
#### 1. Joing Genotyping
    if (useGenomicsDB == true)
        1. GenomicsDBImport (by intervals)*...........(step 21)
    else
        1. CombineGVCFs (by intervals)*...............(step 21)
    2. GenotypeGVCFs (by intervals)*..................(step 22)
    3. VariantFiltration (by intervals)*..............(step 23)
    4. MakeSitesOnlyVcf (by intervals)*...............(step 23)
    5. GatherVcfsCloud................................(step 24)
    6. VQSR
        1. SNPs
            1. VariantRecalibratorCreateModel.........(step 25)
        2. Indels
            1. VariantRecalibrator....................(step 28)
        3. ApplyRecalibration (by intervals)*.........(step 29)
        4. CollectVariantCallingMetrics...............(step 30)

---

## WDL Scripts:

#### WholeGenomeSequencingGATK4.wdl

#### SetWorkingDirectory.wdl

#### Preprocessing.wdl

#### QualityControl.wdl

#### VariantCallingPerSample.wdl

#### BaseQualityScoreRecalibration.wdl

#### QualityControlPerSample.wdl

#### JointGenotyping.wdl

#### VariantQualityScoreRecalibration.wdl

TODO :warning:

---

## Features:

- Possibility to run on a HPC infrastructure connecting the Cromwell engine and the SLURM scheduler.
- Starts from BCL data.
- Demultiplexing of samples pooled across the flowcell.
- Data processing both on a per-lane and a per-sample basis.
- Possibility to handle hg19 and hg38 reference genomes.
- Programmed to restart from every step in case of fail.

For benchmarking, we are following the guidelines of the Truth and Consistency precisionFDA challenges using Genome In A Bottle Consortium released genomes data.

---