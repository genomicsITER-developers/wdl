
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 25. Validate merged BAM with GATK (per-sample)
# 26. Collect filtered WGS Multiple Metrics with GATK (per-sample)
# 27. Compute DepthOfCoverage statistics with GATK (per-sample)
# 28. CallableLoci statistics with GATK (per-lane)
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils
import "SubWorkflows/QualityControl.wdl" as QC

workflow QualityControlPerSampleWF {

  File refFasta
  File? bamFile

  File bedFile
  File bedFilePadding

  String sampleName
  String resultsDir

  String gatkBaseCommand

  String qualimapPath
  String javaMemSize

  # ############################################## #
  # 25. Validate merged BAM with GATK (per-sample) #
  # ############################################## #

  call QC.ValidateSam as ValidateMergedBam {
    input:
      bamFile         = bamFile,
      outputBasename  = sampleName + ".ready.deduped.validation.summary",
      gatkBaseCommand = gatkBaseCommand
  }

  call utils.CopyResultsFilesToDir as copySummaryMergedBam {input: resultsDir = resultsDir, files = ValidateMergedBam.summary}

  # ################################################################ #
  # 26. Collect filtered WGS Multiple Metrics with GATK (per-sample) #
  # ################################################################ #

  call CollectMultipleMetrics {
    input:
      refFasta        = refFasta,
      bamFile         = bamFile,
      outputBasename  = sampleName + ".ready.deduped.multipleqc.metrics",
      gatkBaseCommand = gatkBaseCommand
  }

  call utils.CopyResultsFilesToDir as copyMultipleMetrics {input: resultsDir = resultsDir, files = CollectMultipleMetrics.collectedMetrics}

  # ############################################################# #
  # 27. Qualimap QC of sorted and duplicate-marked BAM (per lane) #
  # ############################################################# #

  call QC.Qualimap as QualimapZeroPaddingPerLane {
    input:
      bamFile      = bamFile,
      bedFile      = bedFile,
      javaMemSize  = javaMemSize,
      resultsDir   = resultsDir,
      qualimapDir  = "qualimapZeroPadding",
      outFile      = sampleName + ".aligned.merged.deduped.sorted.fixed.qualimapqc-zeropadding.pdf",
      qualimapPath = qualimapPath
  }

  call QC.Qualimap as QualimapWithPaddingPerLane {
    input:
      bamFile      = bamFile,
      bedFile      = bedFilePadding,
      javaMemSize  = javaMemSize,
      resultsDir   = resultsDir,
      qualimapDir  = "qualimapWithPadding",
      outFile      = sampleName + ".aligned.merged.deduped.sorted.fixed.qualimapqc-withpadding.pdf",
      qualimapPath = qualimapPath
  }

  # ############################################################# #
  # 28. Compute DepthOfCoverage statistics with GATK (per-sample) #
  # ############################################################# #

  #call DepthOfCoverage as depthOfCov {input:refFasta = refFasta,bamFile = bamFile,refSeq = refSeq,outputBasename = sampleName + ".ready.deduped.DephOfCoverage",gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  # ################################################ #
  # 29. CallableLoci statistics with GATK (per-lane) #
  # ################################################ #

  #call CallableLoci {input:refFasta = refFasta,bamFile = bamFile,outputBasename = sampleName + ".ready.deduped",gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  output {
    File? summary            = ValidateMergedBam.summary
    Array[File]? multMetrics = CollectMultipleMetrics.collectedMetrics

    # NOT IN GATK4
    # File summaryCoverage = depthOfCov.summaryCoverage
    # File callableLociBed = callableLoci.callableLociBed
  }
}


task CollectMultipleMetrics {

  File refFasta
  File? bamFile

  String outputBasename

  String gatkBaseCommand

  command {
    ${gatkBaseCommand} CollectMultipleMetrics \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --REFERENCE_SEQUENCE ${refFasta} \
      --ASSUME_SORTED true \
      --PROGRAM "null" \
      --PROGRAM "CollectBaseDistributionByCycle" \
      --PROGRAM "CollectInsertSizeMetrics" \
      --PROGRAM "MeanQualityByCycle" \
      --PROGRAM "QualityScoreDistribution" \
      --PROGRAM "CollectAlignmentSummaryMetrics" \
      --PROGRAM "CollectGcBiasMetrics" \
      --METRIC_ACCUMULATION_LEVEL "null" \
      --METRIC_ACCUMULATION_LEVEL "ALL_READS" \
      --METRIC_ACCUMULATION_LEVEL "SAMPLE" \
      --METRIC_ACCUMULATION_LEVEL "READ_GROUP"
  }

  output {
    Array[File] collectedMetrics = glob("${outputBasename}*")
  }

  meta {
    taskDescription: "Collect filtered WGS multiple metrics."
  }
}