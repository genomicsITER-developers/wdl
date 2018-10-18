
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
#  9. Validate merged sorted fixed coordinate indexed BAM (per-lane per-sample)
# 10. Qualimap QC of sorted and duplicate-marked BAM (per lane per sample)
# 11. Collect Multiple Metrics (per lane)
# 12. Collect WES Metrics with GATK (per-lane per-sample)
# 13. CollectHsMetrics statistics with GATK (per lane)
# 14. CollectOxoGMetrics statistics with GATK (per lane)
# 15. DepthOfCoverage statistics with GATK (per-lane per-sample)
# 16. CallableLoci statistics with GATK (per-lane per-sample)
# 17. FindCoveredIntervals statistics with GATK (per lane)
# 18. DiagnoseTargets statistics with GATK (per lane)
# 19. QualifyMissingIntervals statistics with GATK (per lane)
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow QualityControlWF {

  # CONFIGURATION PARAMETERS
  Int firstStep
  Int lastStep

  File refFasta
  File? bamFile

  File bedFile
  File bedFilePadding

  File intervalList
  File baits
  File targets

  String resultsDir

  String sampleName

  String gatkPath
  String qualimapPath

  String javaOpts
  String javaMemSize

  # Step 6 - Validate Sam File
  if ((firstStep <= 6) && (6 <= lastStep)) {

    # ############################################################################ #
    # 9. Validate merged sorted fixed coordinate indexed BAM (per-lane per-sample) #
    # ############################################################################ #

    call ValidateSam as ValidateReadyBam {
      input:
        bamFile        = bamFile,
        outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.summary",
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as copySummaryFile {input: resultsDir = resultsDir, files = ValidateReadyBam.summary}

  } # End step 6

  # Step 7 - Qualimap
  if ((firstStep <= 7) && (7 <= lastStep)) {

    # ######################################################################## #
    # 10. Qualimap QC of sorted and duplicate-marked BAM (per lane per sample) #
    # ######################################################################## #

    call Qualimap as QualimapZeroPadding {
      input:
        bamFile      = bamFile,
        bedFile      = bedFile,
        javaMemSize  = javaMemSize,
        resultsDir   = resultsDir,
        qualimapDir  = "qualimapZeroPadding",
        outFile      = sampleName + ".aligned.merged.deduped.sorted.fixed.qualimapqc-zeropadding.pdf",
        qualimapPath = qualimapPath
    }

    call Qualimap as QualimapWithPadding {
      input:
        bamFile      = bamFile,
        bedFile      = bedFilePadding,
        javaMemSize  = javaMemSize,
        resultsDir   = resultsDir,
        qualimapDir  = "qualimapWithPadding",
        outFile      = sampleName + ".aligned.merged.deduped.sorted.fixed.qualimapqc-withpadding.pdf",
        qualimapPath = qualimapPath
    }

  } # End step 7

  # Step 8 - Collect Metrics
  if ((firstStep <= 8) && (8 <= lastStep)) {

    # ####################################### #
    # 11. Collect Multiple Metrics (per lane) #
    # ####################################### #

    call CollectMultipleMetrics {
      input:
        refFasta       = refFasta,
        bamFile        = bamFile,
        intervalList   = intervalList,
        outputBasename = sampleName + ".clean.deduped.sorted.fixed.bam.qcmetrics",
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyMultipleMetrics {input: resultsDir = resultsDir, files = CollectMultipleMetrics.collectedMetrics}

    # ####################################################### #
    # 12. Collect WES Metrics with GATK (per-lane per-sample) #
    # ####################################################### #

    # Collect Raw metrics:
    call CollectRawWgsMetrics {
      input:
        refFasta       = refFasta,
        bamFile        = bamFile,
        outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.wes_metrics_raw",
        minBQ          = 3,
        minMQ          = 0,
        intervalList   = intervalList,
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    # Collect metrics:
    call CollectWgsMetrics {
      input:
        refFasta       = refFasta,
        bamFile        = bamFile,
        outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.wes_metrics_filtered",
        minBQ          = 20,
        minMQ          = 20,
        intervalList   = intervalList,
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    # Collect metrics ONLY in non zero coverage retions:
    call CollectWgsMetricsWithNonZeroCoverage {
      input:
        refFasta       = refFasta,
        bamFile        = bamFile,
        outputBasename = sampleName + ".clean.deduped.sorted.fixed.wes_metrics_nonzero_coverage",
        minBQ          = 2,
        minMQ          = 20,
        intervalList   = intervalList,
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyCollectMetrics {input: resultsDir = resultsDir, 
      files = [CollectRawWgsMetrics.collectedMetrics, CollectWgsMetrics.collectedMetrics, CollectWgsMetricsWithNonZeroCoverage.collectedMetrics, CollectWgsMetricsWithNonZeroCoverage.collectedMetricsPDF]}

    # #################################################### #
    # 13. CollectHsMetrics statistics with GATK (per lane) #
    # #################################################### #

    call CollectHsMetrics {
      input:
        refFasta       = refFasta,
        bamFile        = bamFile,
        outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.HsMetrics",
        baits          = baits,
        targets        = targets,
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyHsMetrics {input: resultsDir = resultsDir, files = [CollectHsMetrics.hsMetrics, CollectHsMetrics.perTargetCoverage]}

    # ###################################################### #
    # 14. CollectOxoGMetrics statistics with GATK (per lane) #
    # ###################################################### #

    call CollectOxoGMetrics {
      input:
        refFasta       = refFasta,
        bamFile        = bamFile,
        outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.oxoG_metrics",
        intervalList   = intervalList,
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyOxoGMetrics {input: resultsDir = resultsDir, files = CollectOxoGMetrics.oxoGMetrics}

  } # End step 8

  # ############################################################## #
  # 15. DepthOfCoverage statistics with GATK (per-lane per-sample) #
  # ############################################################## #

  #call DepthOfCoverage as depthOfCov {input:refFasta = refFasta,bamFile = bamFile,refSeq = refSeq,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.DepthOfCoverage",gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  # ########################################################### #
  # 16. CallableLoci statistics with GATK (per-lane per-sample) #
  # ########################################################### #

  #call CallableLoci {input:refFasta = refFasta,bamFile = bamFile,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed",gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  # ######################################################### #
  # 17. FindCoveredIntervals statistics with GATK (per lane)  #
  # ######################################################### #

  #call FindCoveredIntervals as FindCovered {input:refFasta = refFasta,bamFile = bamFile,regions = regions,intervalPadding = 0,coverageThreshold = 0,minBQ = 0,minMQ = 0,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.covered_regions.list",uncovered = "",gatkPath = gatkPath,javaOpts = javaOpts}

  #call FindCoveredIntervals as FindUncovered {input:refFasta = refFasta,bamFile = bamFile,regions = regions,intervalPadding = 0,coverageThreshold = 0,minBQ = 0,minMQ = 0,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.covered_regions.list",uncovered = "--uncovered"gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  # ################################################### #
  # 18. DiagnoseTargets statistics with GATK (per lane) #
  # ################################################### #

  #call DiagnoseTargets as DiagnoseTargetsCovered {input:refFasta = refFasta,bamFile = bamFile,regions = FindCovered.coveredList,missingIntervals = sampleName + ".aligned.merged.deduped.sorted.fixed.DiagnoseTargets.missing_intervals.list",intervalPadding = 0,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.DiagnoseTargets.vcf",gatkPath = gatkPath,javaOpts = javaOpts}

  #call DiagnoseTargets as DiagnoseTargetsUncovered {input:refFasta = refFasta,bamFile = bamFile,regions = FindUncovered.coveredList,missingIntervals = sampleName + ".aligned.merged.deduped.sorted.fixed.DiagnoseTargets.missing_intervals_uncovered.list",intervalPadding = 0,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.DiagnoseTargets.uncovered.vcf",gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  # ########################################################### #
  # 19. QualifyMissingIntervals statistics with GATK (per lane) #
  # ########################################################### #

  #call QualifyMissingIntervals {input:refFasta = refFasta,bamFile = bamFile,outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.QualifyMissingIntervalss.grp",missingIntervals = DiagnoseTargetsCovered.missingIntervalList,targets = targets,threads = 8,gatkPath = gatkPath,javaOpts = javaOpts}

  # TO DO: call utils.CopyResultsFilesToDir ...

  output {
    File? summary              = ValidateReadyBam.summary
    Array[File]? multMetrics   = CollectMultipleMetrics.collectedMetrics
    File? rawMetrics           = CollectRawWgsMetrics.collectedMetrics
    File? wgsMetrics           = CollectWgsMetrics.collectedMetrics
    File? wgsNonZeroMetrics    = CollectWgsMetricsWithNonZeroCoverage.collectedMetrics
    File? wgsNonZeroMetricsPDF = CollectWgsMetricsWithNonZeroCoverage.collectedMetricsPDF
    File? hsMetrics            = CollectHsMetrics.hsMetrics
    File? hsPerTargetMetrics   = CollectHsMetrics.perTargetCoverage
    File? oxoGMetrics          = CollectOxoGMetrics.oxoGMetrics

    # NOT IN GATK4
    # File summaryCoverage = depthOfCov.summaryCoverage
    # File callableLociBed = callableLoci.callableLociBed
    # File callableLociSummary = callableLoci.callableLociSummary
    # Add... FindCoveredIntervals
    # Add... DiagnoseTargets
    # Add... QualifyMissingIntervals
  }
}


task ValidateSam {

  File? bamFile

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" ValidateSamFile \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --MODE SUMMARY
  }

  output {
    File summary = "${outputBasename}"
  }

  meta {
    taskDescription: "Validate merged sorted fixed coordinate indexed BAM."
  }
}


task Qualimap {

  String resultsDir

  File? bamFile

  File bedFile

  String qualimapDir

  String javaMemSize

  String qualimapPath

  String outFile

  command <<<
    ${qualimapPath} bamqc \
      -bam ${bamFile} \
      --java-mem-size=${javaMemSize} \
      --paint-chromosome-limits \
      --feature-file ${bedFile} \
      --collect-overlap-pairs \
      --outdir ${resultsDir}/${qualimapDir}/html \
      -outfile ${outFile} \
      -outformat PDF:HTML \
      --outside-stats
      # \
    #&& \
    #mv ${resultsDir}/qualimap/html/${outFile} ${resultsDir}/qualimap/${outFile}
  >>>

  meta {
    taskDescription: "Qualimap QC of merged-sorted-fixed BAM."
  }
}


task CollectMultipleMetrics {

  File refFasta
  File? bamFile

  File intervalList

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectMultipleMetrics \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --REFERENCE_SEQUENCE ${refFasta} \
      --ASSUME_SORTED true \
      --PROGRAM "null" \
      --PROGRAM "CollectAlignmentSummaryMetrics" \
      --PROGRAM "CollectInsertSizeMetrics" \
      --PROGRAM "QualityScoreDistribution" \
      --PROGRAM "MeanQualityByCycle" \
      --PROGRAM "CollectBaseDistributionByCycle" \
      --PROGRAM "CollectGcBiasMetrics" \
      --PROGRAM "CollectSequencingArtifactMetrics" \
      --PROGRAM "CollectQualityYieldMetrics" \
      --METRIC_ACCUMULATION_LEVEL "null" \
      --METRIC_ACCUMULATION_LEVEL "SAMPLE" \
      --INTERVALS ${intervalList}
  }

  # Use different levels of metric aggregation:
  # METRIC_ACCUMULATION_LEVEL="ALL_READS" \
  # METRIC_ACCUMULATION_LEVEL="LIBRARY" \
  # METRIC_ACCUMULATION_LEVEL="READ_GROUP" \

  output {
    Array[File] collectedMetrics = glob("${outputBasename}*")
  }

  meta {
    taskDescription: "Collect filtered WGS multiple metrics."
  }
}


task DepthOfCoverage {

  File refFasta
  File refSeq

  File? bamFile

  File regions

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" DepthOfCoverage \
      -R ${refFasta} \
      -I ${bamFile} \
      -O ${outputBasename}
      --calculateCoverageOverGenes ${refSeq} \
      -L ${regions}
  }

  # --interval_padding ${intervalPadding}

  output {
    File summaryCoverage = "${outputBasename}"
  }

  meta {
    taskDescription: "Compute DepthOfCoverage statistics."
  }
}


task CollectRawWgsMetrics {

  File refFasta
  File? bamFile

  String outputBasename

  Int minBQ # Minimum base quality
  Int minMQ # Minimum mapping quality

  File intervalList

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectRawWgsMetrics \
      --REFERENCE_SEQUENCE ${refFasta} \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --INCLUDE_BQ_HISTOGRAM true \
      --COUNT_UNPAIRED true \
      --MINIMUM_BASE_QUALITY ${minBQ} \
      --MINIMUM_MAPPING_QUALITY ${minMQ} \
      --INTERVALS ${intervalList}
  }

  output {
    File collectedMetrics = "${outputBasename}"
  }

  meta {
    taskDescription: "Collect raw WGS multiple metrics."
  }
}


task CollectWgsMetrics {

  File refFasta
  File? bamFile

  String outputBasename

  Int minBQ # Minimum base quality
  Int minMQ # Minimum mapping quality

  File intervalList

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectWgsMetrics \
      --REFERENCE_SEQUENCE ${refFasta} \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --INCLUDE_BQ_HISTOGRAM true \
      --COUNT_UNPAIRED true \
      --MINIMUM_BASE_QUALITY ${minBQ} \
      --MINIMUM_MAPPING_QUALITY ${minMQ} \
      --INTERVALS ${intervalList}
  }

  output {
    File collectedMetrics = "${outputBasename}"
  }

  meta {
    taskDescription: "Collect filtered WGS multiple metrics."
  }
}


task CollectWgsMetricsWithNonZeroCoverage {

  File refFasta
  File? bamFile

  String outputBasename

  Int minBQ # Minimum base quality
  Int minMQ # Minimum mapping quality

  File intervalList

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectWgsMetricsWithNonZeroCoverage \
      --REFERENCE_SEQUENCE ${refFasta} \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --CHART_OUTPUT ${outputBasename}.pdf \
      --INCLUDE_BQ_HISTOGRAM true \
      --COUNT_UNPAIRED true \
      --MINIMUM_BASE_QUALITY ${minBQ} \
      --MINIMUM_MAPPING_QUALITY ${minMQ} \
      --INTERVALS ${intervalList}
  }

  output {
    File collectedMetrics    = "${outputBasename}"
    File collectedMetricsPDF = "${outputBasename}.pdf"
  }

  meta {
    taskDescription: "Collect WGS multiple metrics ONLY in non zero coverage regions."
  }
}


task CallableLoci {

  File refFasta
  File? bamFile

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CallableLoci \
      -R ${refFasta} \
      -I ${bamFile} \
      --summary ${outputBasename}.callable_loci_table.txt \
      -O ${outputBasename}.callable_loci.bed
  }

  output {
    File callableLociBed     = "${outputBasename}.callable_loci.bed"
    File callableLociSummary = "${outputBasename}.callable_loci_table.txt"
  }

  meta {
    taskDescription: "CallableLoci statistics with GATK."
  }
}

task FindCoveredIntervals {

  File refFasta

  File? bamFile

  File regions

  Int intervalPadding
  Int coverageThreshold
  Int minBQ
  Int minMQ

  String outputBasename

  String uncovered # "--uncovered"

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" FindCoveredIntervals \
      -R ${refFasta} \
      -I ${bamFile} \
      -L ${regions} \
      --interval_padding ${intervalPadding} \
      --coverage_threshold ${coverageThreshold} \
      --minBaseQuality ${minBQ} \
      --minMappingQuality ${minMQ} \
      -o ${outputBasename} \
      ${uncovered}
  }

  output {
    File coveredList = "${outputBasename}"
  }

  meta {
    taskDescription: "Find covered or uncovered intervals."
  }
}


task DiagnoseTargets {

  File refFasta
  File? bamFile

  File regions
  String missingIntervals

  Int intervalPadding

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" DiagnoseTargets \
      -R ${refFasta} \
      -I ${bamFile} \
      -L ${regions} \
      --interval_padding ${intervalPadding} \
      --missing_intervals ${missingIntervals} \
      -o ${outputBasename}
  }

  output {
    File diagnoseVCF         = "${outputBasename}"
    File missingIntervalList = "${missingIntervals}"
  }

  meta {
    taskDescription: "DiagnoseTargets statistics with GATK."
  }
}


task QualifyMissingIntervals {

  File refFasta
  File? bamFile

  String outputBasename

  File missingIntervals

  File targets

  Int threads

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" QualifyMissingIntervals \
      -R ${refFasta} \
      -I ${bamFile} \
      -o ${outputBasename} \
      -L ${missingIntervals} \
      -targets ${targets} \
      -nct ${threads}
  }

  output {
    File missedIntervals = "${outputBasename}"
  }

  meta {
    taskDescription: "QualifyMissingIntervals statistics with GATK."
  }
}


task CollectHsMetrics {

  File refFasta
  File? bamFile
  
  String outputBasename

  File baits
  File targets

  String gatkPath
  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectHsMetrics \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --REFERENCE_SEQUENCE ${refFasta} \
      --BAIT_INTERVALS ${baits} \
      --TARGET_INTERVALS ${targets} \
      --PER_TARGET_COVERAGE ${outputBasename}.PerTargetCoverage
  }

  output {
    File hsMetrics         = "${outputBasename}"
    File perTargetCoverage = "${outputBasename}.PerTargetCoverage"
  }

  meta {
    taskDescription: "CollectHsMetrics statistics."
  }
}


task CollectOxoGMetrics {

  File refFasta
  File? bamFile

  String outputBasename

  File intervalList

  String gatkPath
  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectOxoGMetrics \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --REFERENCE_SEQUENCE ${refFasta} \
      --INTERVALS ${intervalList} \
      --MINIMUM_INSERT_SIZE 0 \
      --MINIMUM_MAPPING_QUALITY 30 \
      --MINIMUM_QUALITY_SCORE 20
  }

  output {
    File oxoGMetrics = "${outputBasename}"
  }

  meta {
    taskDescription: "CollectOxoGMetrics statistics."
  }
}
