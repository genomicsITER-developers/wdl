#
# WORKFLOW: Preprocessing
# STEPS:
# 6. Quality Control (per-sample, per-lane)
#

# IMPORTS
import "SubWorkflows/SetWorkingDirectory.wdl" as workDir

workflow QualityControlWF {

  File refFasta
  File? bamFile

  String resultsDir

  String sampleName

  String gatkPath
  String qualimapPath

  String javaOpts

  call ValidateSam as ValidateReadyBam {
    input:
      bamFile = bamFile,
      outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.summary",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copySummaryFile {input: resultsDir = resultsDir, files = ValidateReadyBam.summary}

  call Qualimap {
    input:
      bamFile = bamFile,
      resultsDir = resultsDir,
      outFile = sampleName + ".aligned.merged.deduped.sorted.fixed.qualimapqc.pdf",
      qualimapPath = qualimapPath
  }

  call CollectWGSMetrics as collectMetrics3_0 {
    input:
      refFasta = refFasta,
      bamFile = bamFile,
      outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.raw-multipleqc_3-0.metrics",
      minBQ = 3,
      minMQ = 0,
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call CollectWGSMetrics as collectMetrics20_20 {
    input:
      refFasta = refFasta,
      bamFile = bamFile,
      outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.raw-metrics_20-20.filtered",
      minBQ = 20,
      minMQ = 20,
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call CollectWGSMetricsWithNonZeroCoverage as collectMetricsNonZero {
    input:
      refFasta = refFasta,
      bamFile = bamFile,
      outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.raw-metrics_20-20.filtered-non-zero-coverage",
      minBQ = 20,
      minMQ = 20,
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyCollectMetrics {input: resultsDir = resultsDir, 
    files = [collectMetrics3_0.collectMetrics, collectMetrics20_20.collectMetrics, collectMetricsNonZero.collectMetrics, collectMetricsNonZero.collectMetricsPDF]}

  #call DepthOfCoverage as depthOfCov {
  #  input:
  #    refFasta = refFasta,
  #    bamFile = bamFile,
  #    refSeq = refSeq,
  #    outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed.DepthOfCoverage",
  #    gatkPath = gatkPath,
  #    javaOpts = javaOpts
  #}

  call CollectMultipleMetrics as multipleMetrics {
    input:
      refFasta = refFasta,
      bamFile = bamFile,
      outputBasename = sampleName + ".aligned.duplicates_marked.sorted.fixed.multipleqc.metrics",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyMultipleMetrics {input: resultsDir = resultsDir, files = multipleMetrics.collectMetrics}

  #call CallableLociStatistics as callableLoci {
  #  input:
  #    refFasta = refFasta,
  #    bamFile = bamFile,
  #    outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed",
  #    gatkPath = gatkPath,
  #    javaOpts = javaOpts
  #}

  output {
    File summary = ValidateReadyBam.summary
    File metrics3_0 = collectMetrics3_0.collectMetrics
    File metrics20_20 = collectMetrics20_20.collectMetrics
    File metricsNonZero = collectMetricsNonZero.collectMetrics
    File metricsNonZeroPDF = collectMetricsNonZero.collectMetricsPDF
    #File summaryCoverage = depthOfCov.summaryCoverage
    Array[File] collectMultipleMetrics = multipleMetrics.collectMetrics
    #File callableLociBed = callableLoci.callableLociBed
    #File callableLociSummary = callableLoci.callableLociSummary
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

  String qualimapPath

  String outFile

  command <<<
    ${qualimapPath} bamqc \
      -bam ${bamFile} \
      --java-mem-size=56G \
      --paint-chromosome-limits \
      --collect-overlap-pairs \
      --outdir ${resultsDir}/qualimap/html \
      -outfile ${outFile} \
      -outformat PDF:HTML \
    && \
    mv ${resultsDir}/qualimap/html/${outFile} ${resultsDir}/qualimap/${outFile}
  >>>

  meta {
    taskDescription: "Qualimap QC of merged-sorted-fixed BAM."
  }
}


task CollectWGSMetrics {
  
  File refFasta
  File? bamFile

  String outputBasename

  Int minBQ # Minimum base quality
  Int minMQ # Minimum mapping quality

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectRawWgsMetrics \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --REFERENCE_SEQUENCE ${refFasta} \
      --VALIDATION_STRINGENCY SILENT \
      --COUNT_UNPAIRED true \
      --MINIMUM_BASE_QUALITY ${minBQ} \
      --MINIMUM_MAPPING_QUALITY ${minMQ} \
      --USE_FAST_ALGORITHM true \
      --INCLUDE_BQ_HISTOGRAM true
  }

  output {
    File collectMetrics = "${outputBasename}"
  }

  meta {
    taskDescription: "Collect raw WGS multiple metrics."
  }
}

# From GATK Docs: **BETA FEATURE - WORK IN PROGRESS**
task CollectWGSMetricsWithNonZeroCoverage {

  File refFasta
  File? bamFile

  String outputBasename

  Int minBQ # Minimum base quality
  Int minMQ # Minimum mapping quality

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectWgsMetricsWithNonZeroCoverage \
      -I ${bamFile} \
      -O ${outputBasename} \
      -R ${refFasta} \
      --CHART_OUTPUT ${outputBasename}.pdf \
      --COUNT_UNPAIRED true \
      --MINIMUM_BASE_QUALITY ${minBQ} \
      --MINIMUM_MAPPING_QUALITY ${minMQ} \
      --INCLUDE_BQ_HISTOGRAM true
  }

  output {
    File collectMetrics = "${outputBasename}"
    File collectMetricsPDF = "${outputBasename}.pdf"
  }

  meta {
    taskDescription: "Collect WGS metrics with non zero coverage."
  }
}


task DepthOfCoverage {

  File refFasta
  File refSeq
  File? bamFile

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" DepthOfCoverage \
      -R ${refFasta} \
      -I ${bamFile} \
      -O ${outputBasename}
      --calculateCoverageOverGenes ${refSeq} \
      --partitionType readgroup \
      --omitDepthOutputAtEachBase
  }

  # -nt ${cpuThreads} \
  # --summaryCoverageThreshold 5 \
  # --summaryCoverageThreshold 10 \
  # --summaryCoverageThreshold 15 \
  # --summaryCoverageThreshold 20 \
  # --summaryCoverageThreshold 25 \
  # --summaryCoverageThreshold 30 \
  # --summaryCoverageThreshold 40 \
  # --summaryCoverageThreshold 50 \
  # --summaryCoverageThreshold 100 \
  # --omitIntervalStatistics

  output {
    File summaryCoverage = "${outputBasename}"
  }

  meta {
    taskDescription: "Compute DepthOfCoverage statistics."
  }
}


task CollectMultipleMetrics {
  
  File refFasta
  File? bamFile

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectMultipleMetrics \
      --REFERENCE_SEQUENCE ${refFasta} \
      --INPUT ${bamFile} \
      --OUTPUT ${outputBasename} \
      --ASSUME_SORTED true \
      --PROGRAM CollectBaseDistributionByCycle \
      --PROGRAM CollectInsertSizeMetrics \
      --PROGRAM MeanQualityByCycle \
      --PROGRAM QualityScoreDistribution \
      --PROGRAM CollectAlignmentSummaryMetrics \
      --PROGRAM CollectGcBiasMetrics \
      --METRIC_ACCUMULATION_LEVEL null \
      --METRIC_ACCUMULATION_LEVEL READ_GROUP
  }

  # --METRIC_ACCUMULATION_LEVEL ALL_READS
  # --METRIC_ACCUMULATION_LEVEL SAMPLE

  output {
    Array[File] collectMetrics = glob("${outputBasename}*")
  }

  meta {
    taskDescription: "Collect filtered WGS multiple metrics."
  }
}


task CallableLociStatistics {

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
    File callableLociBed = "${outputBasename}.callable_loci.bed"
    File callableLociSummary = "${outputBasename}.callable_loci_table.txt"
  }

  meta {
    taskDescription: "CallableLoci statistics with GATK."
  }
}