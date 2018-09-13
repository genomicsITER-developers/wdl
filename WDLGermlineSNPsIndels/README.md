## Table of Contents

* [Prerequisites](#prerequisites)
* [Usage](#usage)
* [Workflow](#workflow)
* [Features](#features)

---

## Prerequisites:

TODO :warning:

---

## Usage:

TODO :warning:

---

## Workflow:

### [A. Per sample, per lane](#a-per-sample-per-lane)

#### [1. Preprocessing](#1-preprocessing)
    1. FastqToSam....................................(step 1)
    2. SamToFastq....................................(step 2)
    3. Bwa Mem.......................................(step 3)
    4. MergeBamAlignments............................(step 4)
    5. MarkDuplicates................................(step 5)
    6. SortSam.......................................(step 6)
    7. SetNmMdAndUqTags..............................(step 7)

#### 2. Quality Control (per sample, per lane)
    1. ValidateSam...................................(step 8)
    2. Qualimap......................................(step 9)
    3. CollectRawWgsMetrics (3-0)....................(step 10)
    4. CollectRawWgsMetrics (20-20)..................(step 11)
    5. CollectWgsMetricsWithNonZeroCoverage..........(step 12)
    6. ¿DepthOfCoverage?.............................(step 13)
    7. CollectMultipleMetrics........................(step 14)
    8. ¿CallableLoci?................................(step 15)
 
### [B. Per sample](#b-per-sample)
   
#### 1. Variant Calling (per sample)
    1. MarkDuplicates (MergeBamsPerSample)............(step 16)
    2. BQSR
        1. BaseRecalibrator (por intervalos)*...... ..(step 21)
        2. GatherBQSRReports..........................(step 22)
        3. ApplyBQSR (por intervalos)*........... ....(step 23)
        4. GatherBamFiles.............................(step 24)
    3. HaplotypeCaller (por intervalos)*......... ....(step 25)
    4. MergeVcfs......................................(step 26)
 
#### 2. Quality Control (Con la salida del MergeBamsPerSample)
    1. ValidateSam....................................(step 17)
    2. ¿DepthOfCoverage?..............................(step 18)
    3. CollectMultipleMetrics.........................(step 19)
    4. ¿CallableLoci?.................................(step 20)
 
### [C. Multi-sample](#c-multi-sample)
 
#### 1. Joing Genotyping (multi-sample)
    1. GenomicsDBImport (por intervalos)*.............(step 21)
    2. GenotypeGVCFs (por intervalos)*................(step 22)
    3. VariantFiltration (por intervalos)*............(step 23)
    4. MakeSitesOnlyVcf (por intervalos)*.............(step 23)
    5. GatherVcfsCloud................................(step 24)
    6. VQSR
        1. SNPs
            1. VariantRecalibratorCreateModel.........(step 25)
        2. Indels
            1. VariantRecalibrator....................(step 28)
        3. ApplyRecalibration (por intervalos)*... ...(step 29)
        4. CollectVariantCallingMetrics...............(step 30)

---

### A. Per sample, per lane

TODO :warning:

#### 1. Preprocessing

### B. Per sample

TODO :warning:

### C. Multi-sample

TODO :warning:

---

## Features:

- Possibility to run on a HPC infrastructure connecting the Cromwell engine and the SLURM scheduler.
- Starts from BCL data.
- Demultiplexing of samples pooled across the flowcell.
- Data processing both on a per-lane and a per-sample basis.
- Possibility to handle hg19 and hg38 reference genomes.
- Programmed to restart from every step in case of fail.

For benchmarking, we are following the guidelines of the Truth and Consistency precisionFDA challenges using Genome In A Bottle Consortium released genomes data

---