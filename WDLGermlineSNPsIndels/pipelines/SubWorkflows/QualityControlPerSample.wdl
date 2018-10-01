#
# WORKFLOW: Variant Calling
# STEPS:
# ?. Quality Control (per-sample)
#

# IMPORTS
import "SubWorkflows/SetWorkingDirectory.wdl" as workDir
import "SubWorkflows/QualityControl.wdl" as QC

workflow QualityControlPerSampleWF {

  File refFasta
  File? bamFile

  String resultsDir

  String sampleName

  String gatkPath

  String javaOpts

  call QC.ValidateSam as ValidateMergeBam {
    input:
      bamFile = bamFile,
      outputBasename = sampleName + ".ready.deduped.validation.summary",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copySummaryFilePerSample {input: resultsDir = resultsDir, files = ValidateMergeBam.summary}

  #call DepthOfCoverage as depthOfCovPerSample {
  #  input:
  #    refFasta = refFasta,
  #    bamFile = bamFile,
  #    refSeq = refSeq,
  #    outputBasename = sampleName + ".ready.deduped.DepthOfCoverage",
  #    gatkPath = gatkPath,
  #    javaOpts = javaOpts
  #}

  # call workDir.CopyResultsFilesToDir...

  call QC.CollectMultipleMetrics as multipleMetricsPerSample {
    input:
      refFasta = refFasta,
      bamFile = bamFile,
      outputBasename = sampleName + ".ready.deduped.multipleqc.metrics",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyMultipleMetricsPerSample {input: resultsDir = resultsDir, files = multipleMetricsPerSample.collectMetrics}

  #call CallableLociStatistics as callableLociPerSample {
  #  input:
  #    refFasta = refFasta,
  #    bamFile = bamFile,
  #    outputBasename = sampleName + ".aligned.deduped",
  #    gatkPath = gatkPath,
  #    javaOpts = javaOpts
  #}

  # call workDir.CopyResultsFilesToDir...

  output {
    File summary = ValidateMergeBam.summary
    #File summaryCoverage = depthOfCov.summaryCoverage
    Array[File] collectMultipleMetrics = multipleMetricsPerSample.collectMetrics
    #File callableLociBed = callableLoci.callableLociBed
    #File callableLociSummary = callableLoci.callableLociSummary
  }

}