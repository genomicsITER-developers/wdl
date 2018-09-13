Workflow:
---------

### A. Per sample, per lane

1. Preprocessing
··1.1 FastqToSam....................................(step 1)
··1.2 SamToFastq....................................(step 2)
··1.3 Bwa Mem.......................................(step 3)
··1.4 MergeBamAlignments............................(step 4)
··1.5 MarkDuplicates................................(step 5)
··1.6 SortSam.......................................(step 6)
··1.7 SetNmMdAndUqTags..............................(step 7)


---

Features:
- Possibility to run on a HPC infrastructure connecting the Cromwell engine and the SLURM scheduler.
- Starts from BCL data.
- Demultiplexing of samples pooled across the flowcell.
- Data processing both on a per-lane and a per-sample basis.
- Possibility to handle hg19 and hg38 reference genomes.
- Programmed to restart from every step in case of fail.

For benchmarking, we are following the guidelines of the Truth and Consistency precisionFDA challenges using Genome In A Bottle Consortium released genomes data
