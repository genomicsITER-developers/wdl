# WDL-based pipelines for whole genome and exome sequencing analysis
<!--
<p align="center">
   <img src="/home/administrador/Imágenes/JBI-poster/folder/ULL-recortada-small.png"/>
   <img src="/home/administrador/Imágenes/JBI-poster/folder/Taro-small.png"/>
   <img src="/home/administrador/Imágenes/JBI-poster/folder/ITER-small.png"/>
   <img src="/home/administrador/Imágenes/JBI-poster/folder/CabildodeTenerife-small.png"/>
</p>
-->

<a target="_blank" rel="noopener noreferrer" href="https://www.ull.es/">
<img src="/home/administrador/Imágenes/JBI-poster/folder/ULL-recortada-small.png" alt="ull" hspace="3"/>
</a>
<a target="_blank" rel="noopener noreferrer" href="https://www.ull.es/investigacion/grupos-investigacion/taro-ingenieria-del-software-y-bases-de-datos/">
<img src="/home/administrador/Imágenes/JBI-poster/folder/Taro-small.png" alt="taro" hspace="5"/>
</a>
<a target="_blank" rel="noopener noreferrer" href="http://www.iter.es/">
<img src="/home/administrador/Imágenes/JBI-poster/folder/ITER-small.png" alt="iter" hspace="5"/>
</a>
<a target="_blank" rel="noopener noreferrer" href="http://www.tenerife.es/portalcabtfe/es/">
<img src="/home/administrador/Imágenes/JBI-poster/folder/CabildodeTenerife-small.png" alt="cabildo" hspace="3"/>
</a>

<!--
<div class="row" style="clear:both; display:table;">
  <div class="column" style="float:left; padding-right:15px;">
    <a href="https://www.ull.es/">
      <img src="/home/administrador/Imágenes/JBI-poster/folder/ULL-recortada-small.png" alt="ull">
    </a>
  </div>
  <div class="column" style="float:left; padding-right:15px;">
    <a href="https://www.ull.es/investigacion/grupos-investigacion/taro-ingenieria-del-software-y-bases-de-datos/">
      <img src="/home/administrador/Imágenes/JBI-poster/folder/Taro-small.png" alt="taro">
    </a>
  </div>
  <div class="column" style="float:left; padding-right:15px;">
    <a href="http://www.iter.es/">
      <img src="/home/administrador/Imágenes/JBI-poster/folder/ITER-small.png" alt="iter">
    </a>
  </div>
  <div class="column" style="float:left;">
    <a href="http://www.tenerife.es/portalcabtfe/es/">
      <img src="/home/administrador/Imágenes/JBI-poster/folder/CabildodeTenerife-small.png" alt="cabildo">
    </a>
  </div>
</div>
-->

<!--
# ![](https://github.com/AdrianMBarrera/Presentations/blob/master/JBI-2018-images/Logos-1.png?raw=true)
-->

## Table of Contents

* [Description](#description)
* [Prerequisites](#prerequisites)
* [Features](#features)
* [WGSGermlineSNPsIndels](#wgsgermlinesnpsindels)
* [WESGermlineSNPsIndels](#wesgermlinesnpsindels)
* [Funding and Acknowledgement](#funding-and-acknowledgement)

---

## Description

Next Generation Sequencing data analysis comprises a series of computational tasks frequently based on the use of command line tools. These analyses are defined in workflows that group all the necessary tasks, improving data processing performance and results interpretation. Some Domain Specific Languages (DSLs), such as WDL and Nextflow, have been recently created to define and program complex pipelines, as well as to improve the parallelization, the scalability and the reusability. We have developed complete pipelines programmed in WDL via scripting and Rabix Composer based on the Broad Institute’s best practices and the Genome Analysis Toolkit (GATK4) to analyze whole-genome (WGS) and whole-exome (WES) data.

For benchmarking, we are following the guidelines of the Truth and Consistency precisionFDA challenges using Genome In A Bottle Consortium released genomes data. A full pipeline is currently running on TeideHPC to analyze WGS and WES germline data produced by an Illumina HiSeq4000 sequencing platform for research purposes.

We have developed two workflows based in GATK4 using WDL and Cromwell technologies, and both of them could run in local mode, over a HPC infrastructure or in a dockerized cluster.

![](https://github.com/AdrianMBarrera/Presentations/blob/master/JBI-2018-images/from-sequencing-to-execution.png?raw=true)

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

## Features

- Possibility to run on a HPC infrastructure connecting the Cromwell engine and the SLURM scheduler.
- Starts from BCL data.
- Demultiplexing of samples pooled across the flowcell.
- Data processing both on a per-lane and a per-sample basis.
- Possibility to handle hg19 and hg38 reference genomes.
- Programmed to restart from every step in case of fail.

For benchmarking, we are following the guidelines of the Truth and Consistency [precisionFDA challenges](https://precision.fda.gov/) using [Genome In A Bottle Consortium](http://jimb.stanford.edu/giab/) released genomes data.

---

## WGSGermlineSNPsIndels

Pipeline for whole genome and sequencing analysis.

* [Documentation](https://github.com/genomicsITER-developers/wdl/tree/master/WGSGermlineSNPsIndels)

---

## WESGermlineSNPsIndels

Pipeline for whole exome and sequencing analysis.

* [Documentation](https://github.com/genomicsITER-developers/wdl/tree/master/WESGermlineSNPsIndels)

---

## Funding and Acknowledgement

>Funded by Ministerio de Ciencia, Innovación y Universidades (RTC-2017-6471-1; MINECO/AEI/FEDER, UE). This work has been supported by the CEDeI program (Centro de Excelencia de Desarrollo e Innovación, Cabildo de Tenerife). The authors also thankfully acknowledge the computer resources and the technical support provided by TARO Research Group of the University of La Laguna.
>
>## ![](https://github.com/AdrianMBarrera/Presentations/blob/master/JBI-2018-images/Logos-2.png?raw=true)

For more information, see the following [poster](https://github.com/AdrianMBarrera/Presentations/raw/master/Poster_JBI_2018_WDL-based-pipelines-for-WGS-and-WES-analysis.pdf).
