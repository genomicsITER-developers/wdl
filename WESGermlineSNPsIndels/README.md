## Table of Contents

* [Usage](#usage)
    * [Set input parameters](#set-input-parameters)
    * [Build TSV for pair reads information](#build-tsv-for-pair-reads-information)
    * [Local](#local)
    * [HPC using Slurm scheduler](#hpc-using-slurm-scheduler)
* [Workflow](#workflow)
* [WDL Scripts](#wdl-scripts)

---

## Usage

### Set input parameters

First of all, you need to configure all the inputs needed by the workflow. To do that, you need to change the values of the following parameters specified in the *inputs/WholeExomeSequencingGATK4.inputs.json*.

```
{
   "##_1": "REFERENCE GENOME DATA",
   "##_NOTE_1": "Set '...' to your own path",
   "WholeExomeSequencingGATK4WF.refBaseDir":                      ".../hg19/",
   "WholeExomeSequencingGATK4WF.refBaseName":                     "hg19_ref_genome",

   "WholeExomeSequencingGATK4WF.refFasta":                        ".../hg19/hg19_ref_genome.fasta",
   "WholeExomeSequencingGATK4WF.refAlt":                          ".../hg19/hg19_ref_genome.alt.fasta",
   "WholeExomeSequencingGATK4WF.refIndex":                        ".../hg19/hg19_ref_genome.fasta.fai",
   "WholeExomeSequencingGATK4WF.refDict":                         ".../hg19/hg19_ref_genome.dict",
   "WholeExomeSequencingGATK4WF.refAmb":                          ".../hg19/hg19_ref_genome.fasta.amb",
   "WholeExomeSequencingGATK4WF.refAnn":                          ".../hg19/hg19_ref_genome.fasta.ann",
   "WholeExomeSequencingGATK4WF.refBwt":                          ".../hg19/hg19_ref_genome.fasta.bwt",
   "WholeExomeSequencingGATK4WF.refPac":                          ".../hg19/hg19_ref_genome.fasta.pac",
   "WholeExomeSequencingGATK4WF.refSa":                           ".../hg19/hg19_ref_genome.fasta.sa",

   "WholeExomeSequencingGATK4WF.dbSnps":                          ".../hg19/dbs/dbsnp_138.hg19.vcf.gz",
   "WholeExomeSequencingGATK4WF.dbSnpsIdx":                       ".../hg19/dbs/dbsnp_138.hg19.vcf.gz.tbi",
   "WholeExomeSequencingGATK4WF.snp1000g":                        ".../hg19/dbs/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz",
   "WholeExomeSequencingGATK4WF.snp1000gIdx":                     ".../hg19/dbs/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz.tbi",
   "WholeExomeSequencingGATK4WF.hapMapResource":                  ".../hg19/dbs/hapmap_3.3.hg19.sites.vcf.gz",
   "WholeExomeSequencingGATK4WF.hapMapResourceIndex":             ".../hg19/dbs/hapmap_3.3.hg19.sites.vcf.gz.tbi",
   "WholeExomeSequencingGATK4WF.omniResource":                    ".../hg19/dbs/1000G_omni2.5.hg19.sites.vcf.gz",
   "WholeExomeSequencingGATK4WF.omniResourceIndex":               ".../hg19/dbs/1000G_omni2.5.hg19.sites.vcf.gz.tbi",
   "WholeExomeSequencingGATK4WF.millsResource":                   ".../hg19/dbs/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
   "WholeExomeSequencingGATK4WF.millsResourceIndex":              ".../hg19/dbs/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz.tbi",

   "WholeExomeSequencingGATK4WF.protectionTag":                        "''",

   "WholeExomeSequencingGATK4WF.useGenomicsDB": false,

   "##_ADDITIONAL_PARAMS": "tranches, annotations, filter levels...",
   "WholeExomeSequencingGATK4WF.snpRecalibrationTrancheValues":      ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"],
   "WholeExomeSequencingGATK4WF.snpRecalibrationAnnotationValues":   ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"],
   "WholeExomeSequencingGATK4WF.indelRecalibrationTrancheValues":    ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"],
   "WholeExomeSequencingGATK4WF.indelRecalibrationAnnotationValues": ["QD", "MQRankSum", "ReadPosRankSum", "FS", "SOR", "DP"],
   "WholeExomeSequencingGATK4WF.snpFilterLevel":                     99.7,
   "WholeExomeSequencingGATK4WF.indelFilterLevel":                   99.7,
   "WholeExomeSequencingGATK4WF.snpVQSRDownsampleFactor":            10,
   "WholeExomeSequencingGATK4WF.indelsMaxGaussians":                 2,
   "WholeExomeSequencingGATK4WF.snpsMaxGaussians":                   4,
   "WholeExomeSequencingGATK4WF.maxAltAlleles":                      3,
   "WholeExomeSequencingGATK4WF.intervalPadding":                    0,

   "WholeExomeSequencingGATK4WF.snpRecalibrationTrancheValues_SMALL_LIST": ["100.0", "99.9", "99.0", "90.0"],
   "WholeExomeSequencingGATK4WF.indelRecalibrationTrancheValues_SMALL_LIST": ["100.0", "99.9", "99.0", "95.0", "90.0"],


   "##_2": "READS",
   "##_NOTE_2": "Set '...' to reads path and libraryName, platform and sequencingCenter to your own values.",
   "WholeExomeSequencingGATK4WF.fastqReadsTSV": ".../data/read_pairs.tsv",
   "WholeExomeSequencingGATK4WF.libraryName": "LIB-01",
   "WholeExomeSequencingGATK4WF.platform": "illumina",
   "WholeExomeSequencingGATK4WF.sequencingCenter": "ITER-GENOMICA",


   "##_3": "INTERVALS FILES",
   "##_NOTE_3": "Set '...' to reads path.",
   "WholeExomeSequencingGATK4WF.bedFile":             ".../data/Illumina_TruSeq_Rapid_Exome/truseq-rapid-exome-targeted-regions-manifest-v1-2_6cols.bed",
   "WholeExomeSequencingGATK4WF.bedFilePadding":      ".../data/Illumina_TruSeq_Rapid_Exome/truseq-rapid-exome-targeted-regions-manifest-v1-2_6cols_padding100bp.bed",
   "WholeExomeSequencingGATK4WF.intervalList":        ".../data/Illumina_TruSeq_Rapid_Exome/truseq-rapid-exome-targeted-regions-manifest-v1-2_6cols.targets.interval_list",
   "WholeExomeSequencingGATK4WF.targets":             ".../data/Illumina_TruSeq_Rapid_Exome/truseq-rapid-exome-targeted-regions-manifest-v1-2_6cols.targets.interval_list",
   "WholeExomeSequencingGATK4WF.baits":               ".../data/Illumina_TruSeq_Rapid_Exome/truseq-rapid-exome-probes-manifest-v1-2_6cols.baits.interval_list",
   "WholeExomeSequencingGATK4WF.TSVIntervalsFile":    ".../data/sequence_grouping_small.txt",
   "WholeExomeSequencingGATK4WF.intervalContigsFile": ".../data/sequence_grouping_small.txt",
   "WholeExomeSequencingGATK4WF.intervalTargetsFile": ".../data/Nextera_Rapid_Capture_Expanded_Exome/nexterarapidcapture_expandedexome_targetedregions.bed",


   "##_4": "BWA",
   "##_NOTE_4": "Set '...' to BWA path and the -t value to your own configuration.",
   "WholeExomeSequencingGATK4WF.bwaPath": "bwa",
   "WholeExomeSequencingGATK4WF.bwaMemCommand": "mem -K 100000000 -p -v 3 -t 2 -Y $bashRefFasta $bashReadFastq",

   "WholeExomeSequencingGATK4WF.bwaMemCommand_pairedSplit": "mem -K 100000000 -v 3 -t 1 -Y $bashRefFasta $bashReadFastq1 $bashReadFastq2",


   "##_5": "SAMTOOLS",
   "##_NOTE_5": "Set '...' to Samtools path.",
   "WholeExomeSequencingGATK4WF.samtoolsPath": ".../samtools",


   "##_6": "GATK",
   "##_NOTE_6": "Set '...' to GATK4 path.",
   "WholeExomeSequencingGATK4WF.gatkPath": ".../gatk",


   "##_7": "PICARD",
   "##_NOTE_7": "Set '...' to Picard path.",
   "WholeExomeSequencingGATK4WF.picardPath": ".../picard.jar",


   "##_8": "PYTHON",
   "##_NOTE_8": "Set python path.",
   "WholeExomeSequencingGATK4WF.pythonPath": ".../python2",


   "##_9": "QUALIMAP",
   "##_NOTE_9": "Set '...' to Qualimap path.",
   "WholeExomeSequencingGATK4WF.qualimapPath": ".../qualimap",
   "WholeExomeSequencingGATK4WF.javaMemSize":  "8G",


   "##_12": "RESULTS DIR",
   "##_NOTE_10": "Set '...' to results path.",
   "WholeExomeSequencingGATK4WF.resultsDir": ".../data/results",


   "##_13": "CONFIGURATION PARAMETERS",
   "WholeExomeSequencingGATK4WF.wfPerSamplePerLane":    true,
   "WholeExomeSequencingGATK4WF.wfPreprocess":          true,
   "WholeExomeSequencingGATK4WF.wfQC_PerSamplePerLane": true,
   "WholeExomeSequencingGATK4WF.wfBQSR":                true,

   "WholeExomeSequencingGATK4WF.wfPerSample":           true,
   "WholeExomeSequencingGATK4WF.wfVariantCalling":      true,
   "WholeExomeSequencingGATK4WF.wfQC_PerSample":        true,

   "WholeExomeSequencingGATK4WF.wfMultiSample":         true,
   "WholeExomeSequencingGATK4WF.wfJointGenotyping":     true,
   "WholeExomeSequencingGATK4WF.wfVQSR":                true,
   "WholeExomeSequencingGATK4WF.wfQC_MultiSample":      true,
   "WholeExomeSequencingGATK4WF.wfAnnotation":          false,

   "WholeExomeSequencingGATK4WF.firstStep":             0,
   "WholeExomeSequencingGATK4WF.lastStep":              100,

   "##_14": "Java Opts: Select Xms16g for 30GB-RAM and 60GB-RAM nodes. Select Xmx28g for 30GB-RAM or Xmx56 for 60GB-RAM nodes",
   "WholeExomeSequencingGATK4WF.javaOpts": "-Xms8g -Xmx8g"
}
```

Optionally, you can create a clean JSON template with these parameters using Womtool with the following command:

```
java -jar womtool-<version>.jar inputs WholeExomeSequencingGATK4.wdl > WholeExomeSequencingGATK4.inputs.json
```

### Build TSV for pair reads information

In order to run the pipeline, you need to create a TSV file with the following columns for each pair of reads (per-sample and per-lane):

```
FLOWCELL  LANE  SAMPLE_ID  INDEX SAMPLE_PROJECT  SAMPLE_DIRECTORY READ1  READ2
```

For example:

```
H814YADXX 1 A CGATGT  A_Index  A_L001_001 /path/to/A_L001_R1.fastq.gz /path/to/A_L001_R2.fastq.gz
H814YADXX 2 A CGATGT  A_Index  A_L002_001 /path/to/A_L002_R1.fastq.gz /path/to/A_L002_R2.fastq.gz
H814YADXX 1 B TGACCA  B_Index  B_L001_001 /path/to/B_L001_R1.fastq.gz /path/to/B_L001_R2.fastq.gz
H814YADXX 2 B TGACCA  B_Index  B_L002_001 /path/to/B_L002_R1.fastq.gz /path/to/B_L002_R2.fastq.gz
```


### Local

You can run cromwell in local with its default settings or create a configuration file if you want to customize some settings. In the *backends/backend-local.conf* file, we set the **localization** option to create **soft-links** instead of create hard-links or copies of the input and output files.

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

Use the following command to run the workflow in your system using this configuration file:
```
java -Dconfig.file=./backends/backend-local.conf -jar cromwell-<version>.jar run ./pipelines/WholeExomeSequencingGATK4.wdl -i ./inputs/WholeExomeSequencingGATK4.inputs.json -m metadata.json
```

### HPC using Slurm scheduler

To run Cromwell in a HPC environment with Slurm scheduler, you need to change the Cromwell configuration to interact with a SLURM and dispatch jobs to it. We provide this configuration in the *backends/backend-slurm.conf* file:

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

In the file *slurm/run_pipeline.sh*, we provide a **Bash** script with the commands to run Cromwell using Slurm with the previous configuration file:

```
#!/bin/bash
#SBATCH -p long
#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH -t 3-00:00:00
#SBATCH --mem=60000
#SBATCH -o slurm-wgs-gatk-%j.out
#SBATCH -e slurm-wgs-gatk-%j.err

source /etc/profile.d/profile.modules.sh
module load java-jre/1.8.0_77
module load gatk/4.0.6.0
module load bwa/0.7.15
module load picard/2.10.10
module load python/2.7.14
module load python/3.4.5
module load proj/4.9.3/gcc gdal/2.2.4/gcc gcc/6.1.0 tiff/4.0.9/gcc pcre/3.38/gcc curl/7.40/gcc xz/5.2.2/gcc R/3.5-GA/gcc
module load samtools/1.7
module load htslib/1.3.1
module load java-jre/1.8.0_77 qualimap/2.2.1

java -Dconfig.file=./backends/backend-slurm.conf \
-jar .../cromwell-<version>.jar \
run ./pipelines/WholeExomeSequencingGATK4.wdl \
-i ./inputs/WholeExomeSequencingGATK4.inputs.json \
-m metadata.json
```

Use the following command to run the script in a system with Slurm scheduler (such as TeideHPC). Thanks to the **screen** tool, the workflow will continue running even if you log out of your session:
```
screen -d -m -t wgs-wdl-pipeline -L sh ./slurm/run_pipeline.sh
```

---

## Workflow

![](https://github.com/AdrianMBarrera/Presentations/blob/master/JBI-2018-images/wgs-gatk-pipeline-jbi.png?raw=true)

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
    3. CollectMultipleMetrics........................(step 10)
    4. CollectRawWgsMetrics..........................(step 11)
    5. CollectWgsMetrics.............................(step 12)
    6. CollectWgsMetricsWithNonZeroCoverage..........(step 13)
    7. CollectHsMetrics..............................(step 14)
    8. CollectOxoGMetrics............................(step 15)

#### 3. BaseQualityScoreRecalibration
    1. BaseRecalibrator (by intervals*)..............(step 16)
    2. GatherBQSRReports.............................(step 17)
    3. ApplyBQSR (by intervals*).....................(step 18)
    4. GatherBamFiles................................(step 19)

### B. Per sample

#### 1. Variant Calling
    1. MarkDuplicates (MergeBamsPerSample)...........(step 20)
    2. HaplotypeCaller (by intervals*)...............(step 21)
    3. MergeVcfs.....................................(step 22)

#### 2. Quality Control (Using MarkDuplicates (MergeBamsPerSample) output file)
    1. ValidateSam...................................(step 23)
    2. CollectMultipleMetrics........................(step 24)
    3. Qualimap......................................(step 25)

### C. Multi-sample

#### 1. Joing Genotyping
    if (useGenomicsDB == true)
        1. GenomicsDBImport (by intervals*)..........(step 26)───┐
    else                                                         |─ In the inputs JSON file, select if you are using GenomicsDB or not.
        1. CombineGVCFs (by intervals*)..............(step 26)───┘
    2. GenotypeGVCFs (by intervals*).................(step 27)
    3. VariantFiltration (by intervals*).............(step 28)
    4. MakeSitesOnlyVcf (by intervals*)..............(step 29)
    5. GatherVcfsCloud...............................(step 30)

#### 2. VariantQualityScoreRecalibration
    1. SNPs
        1. VariantRecalibratorCreateModel............(step 31)
    2. Indels
        1. VariantRecalibrator.......................(step 32)
    3. ApplyRecalibration (by intervals*)............(step 33)
    4. CollectVariantCallingMetrics..................(step 34)

#### 3. Quality Control
    1. CollectVariantCallingMetrics..................(step 35)


**\* Intervals:**
   * The -L argument directs the GATK engine to restrict processing to specific genomic intervals.
   * This argument allows the pipeline to run the Scatter-Gather functionality.

---

## WDL Scripts

Hierarchy tree:

```
WholeExomeSequencingGATK4.wdl
    |
    |
    └─── Utils.wdl
    |
    |
    └─── Preprocessing.wdl
    |    |
    |    |
    |    └─── QualityControl.wdl
    |    |
    |    |
    |    └─── BaseQualityScoreRecalibration.wdl
    |
    |
    └─── VariantCallingPerSample.wdl
    |    |
    |    |
    |    └─── QualityControlPerSample.wdl
    |
    |
    └─── JointGenotyping.wdl
         |
         |
         └─── VariantQualityScoreRecalibration.wdl
         |
         |
         └─── QualityControlMultiSample.wdl
```

---
