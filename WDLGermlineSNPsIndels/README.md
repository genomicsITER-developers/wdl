Workflow:
---------

### A. Per sample, per lane

1. Preprocessing
1.1 FastqToSam....................................(step 1)
1.2 SamToFastq....................................(step 2)
1.3 Bwa Mem.......................................(step 3)
1.4 MergeBamAlignments............................(step 4)
1.5 MarkDuplicates................................(step 5)
1.6 SortSam.......................................(step 6)
1.7 SetNmMdAndUqTags..............................(step 7)

2. Quality Control (per sample, per lane)
2.1 ValidateSam...................................(step 8)
2.2 Qualimap......................................(step 9)
2.3 CollectRawWgsMetrics (3-0)....................(step 10)
2.4 CollectRawWgsMetrics (20-20)..................(step 11)
2.5 CollectWgsMetricsWithNonZeroCoverage..........(step 12)
2.6 ¿DepthOfCoverage?.............................(step 13)
2.7 CollectMultipleMetrics........................(step 14)
2.8 ¿CallableLoci?................................(step 15)
 
### B. Per sample
   
1. Variant Calling (per sample)
1.1 MarkDuplicates (MergeBamsPerSample)...........(step 16)
1.2 BQSR
1.2.1 BaseRecalibrator (por intervalos)*........(step 21)
1.2.2 GatherBQSRReports.........................(step 22)
1.2.3 ApplyBQSR (por intervalos)*...............(step 23)
1.2.4 GatherBamFiles............................(step 24)
1.3 HaplotypeCaller (por intervalos)*.............(step 25)
1.4 MergeVcfs.....................................(step 26)
 
2. Quality Control (Con la salida del MergeBamsPerSample)
2.1 ValidateSam...................................(step 17)
2.2 ¿DepthOfCoverage?.............................(step 18)
2.3 CollectMultipleMetrics........................(step 19)
2.4 ¿CallableLoci?................................(step 20)
 
### C. Multi-sample
 
1. Joing Genotyping (multi-sample)
1.1 GenomicsDBImport (por intervalos)*............(step 21)
1.2 GenotypeGVCFs (por intervalos)*...............(step 22)
1.3 VariantFiltration (por intervalos)*...........(step 23)
1.4 MakeSitesOnlyVcf (por intervalos)*............(step 23)
1.5 GatherVcfsCloud...............................(step 24)
1.6 VQSR
1.6.1 SNPs
1.6.1.1 VariantRecalibratorCreateModel........(step 25)
1.6.2 Indels
1.6.2.1 VariantRecalibrator...................(step 28)
1.6.3 ApplyRecalibration (por intervalos)*......(step 29)
1.6.4 CollectVariantCallingMetrics..............(step 30)


Features:
- Possibility to run on a HPC infrastructure connecting the Cromwell engine and the SLURM scheduler.
- Starts from BCL data.
- Demultiplexing of samples pooled across the flowcell.
- Data processing both on a per-lane and a per-sample basis.
- Possibility to handle hg19 and hg38 reference genomes.
- Programmed to restart from every step in case of fail.

For benchmarking, we are following the guidelines of the Truth and Consistency precisionFDA challenges using Genome In A Bottle Consortium released genomes data
