# IMPORTS
import "SubWorkflows/SetWorkingDirectory.wdl" as workDir

workflow BaseQualityScoreRecalibrationWF {

  Array[Array[String]] intervals

  File? inputBam
  File? inputBai

  File dbSnpsVcf
  File dbSnpsVcfIdx
  Array[File] knownIndelsSitesVcfs
  Array[File] knownIndelsSitesIndices

  File refFasta
  File refDict
  File refIndex

  String sampleName

  String gatkPath

  String resultsDir

  String javaOpts

  #scatter (interval in read_tsv(select_first([TSVIntervalsFile, CreateSequenceGroupingTSV.sequenceGrouping]))) {
  scatter (interval in intervals) {
    call BaseRecalibrator {
      input:
        #inputBam = if (firstStep < 9) then SortSam.sortedBam else resultsDir + "/" + sampleName + ".aligned.duplicates_marked.sorted.bam",
        #inputBai = if (firstStep < 9) then SortSam.sortedBamIndex else resultsDir + "/" + sampleName + ".aligned.duplicates_marked.sorted.bai",
        inputBam = inputBam,
        inputBai = inputBai,
        recalReportFilename = sampleName + "ready.deduped.parcial_recal_data.csv",
        sequenceGroupInterval = interval,
        dbSnpsVcf = dbSnpsVcf,
        dbSnpsVcfIdx = dbSnpsVcfIdx,
        knownIndelsSitesVcfs = knownIndelsSitesVcfs,
        knownIndelsSitesIndices = knownIndelsSitesIndices,
        refFasta = refFasta,
        refDict = refDict,
        refIndex = refIndex,
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }
  }

  call GatherBaseRecalReports {
    input:
      baseRecalReports = BaseRecalibrator.recalibrationReport,
      outputReportsFilename = sampleName + "ready.deduped.recal_data.csv",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyBqsrReports {input: resultsDir = resultsDir, files = GatherBaseRecalReports.bqsrReports}

  #scatter (interval in read_tsv(select_first([TSVIntervalsFile, CreateSequenceGroupingTSV.sequenceGroupingWithUnmapped]))) {
  scatter (interval in intervals) {
    call ApplyBQSR {
      input:
        #inputBam = if (firstStep < 9) then SortSam.sortedBam else resultsDir + "/" + sampleName + ".aligned.duplicates_marked.sorted.bam",
        #inputBai = if (firstStep < 9) then SortSam.sortedBamIndex else resultsDir + "/" + sampleName + ".aligned.duplicates_marked.sorted.bai",
        inputBam = inputBam,
        inputBai = inputBai,
        outputBasename = sampleName + ".ready.deduped.parcial_recalibrated",
        recalibrationReport = GatherBaseRecalReports.bqsrReports,
        sequenceGroupInterval = interval,
        refFasta = refFasta,
        refDict = refDict,
        refIndex = refIndex,
        gatkPath = gatkPath,
        javaOpts = javaOpts

    }
  }

  call GatherBamRecalibratedFiles {
    input:
      inputBams = ApplyBQSR.recalibratedBam,
      outputBasename = sampleName + ".ready.deduped.recalibrated",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyRecalibratedBam {input: resultsDir = resultsDir, 
    files = [GatherBamRecalibratedFiles.recalibratedBam, GatherBamRecalibratedFiles.recalibratedBamIndex, GatherBamRecalibratedFiles.recalibratedBamChecksum]}

  output {
  	File recalibrationReports = GatherBaseRecalReports.bqsrReports
    File recalibratedBam = GatherBamRecalibratedFiles.recalibratedBam
    File recalibratedBamIndex = GatherBamRecalibratedFiles.recalibratedBamIndex
    File recalibratedBamChecksum = GatherBamRecalibratedFiles.recalibratedBamChecksum
  }

}

task BaseRecalibrator {

  File? inputBam
  File? inputBai

  String recalReportFilename

  Array[String] sequenceGroupInterval

  File dbSnpsVcf
  File dbSnpsVcfIdx
  Array[File] knownIndelsSitesVcfs
  Array[File] knownIndelsSitesIndices

  File refFasta
  File refDict
  File refIndex

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" BaseRecalibrator \
      --reference ${refFasta} \
      --input ${inputBam} \
      --intervals ${sep=" -L " sequenceGroupInterval} \
      --use-original-qualities \
      --output ${recalReportFilename} \
      --known-sites ${dbSnpsVcf} \
      --known-sites ${sep=" --known-sites " knownIndelsSitesVcfs}
  }

  output {
    File recalibrationReport = "${recalReportFilename}"
  }

  meta {
    taskDescription: "Generates recalibration table for Base Quality Score Recalibration (BQSR) (per-sample)."
  }
}


task GatherBaseRecalReports {

  Array[File] baseRecalReports
  String outputReportsFilename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" GatherBQSRReports \
      -I ${sep=" -I " baseRecalReports} \
      -O ${outputReportsFilename}
  }

  output {
    File bqsrReports = "${outputReportsFilename}"
  }

  meta {
    taskDescription: "Gather scattered BQSR recalibration reports into a single file (per-sample)."
  }
}


task ApplyBQSR {

  File? inputBam
  File? inputBai

  String outputBasename

  File recalibrationReport

  Array[String] sequenceGroupInterval

  File refFasta
  File refDict
  File refIndex

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ${refFasta} \
      -I ${inputBam} \
      --use-original-qualities \
      -O ${outputBasename}.bam \
      -bqsr ${recalibrationReport} \
      --static-quantized-quals 10 \
      --static-quantized-quals 20 \
      --static-quantized-quals 30 \
      -L ${sep=" -L " sequenceGroupInterval}
  }

  output {
    File recalibratedBam = "${outputBasename}.bam"
    File recalibratedBamChecksum = "${outputBasename}.bam.md5"
  }

  meta {
    taskDescription: "Apply BQSR (per-sample)."
  }
}


task GatherBamRecalibratedFiles {

  Array[File] inputBams

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" GatherBamFiles \
      -I ${sep=" -I " inputBams} \
      -O ${outputBasename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }

  output {
    File recalibratedBam = "${outputBasename}.bam"
    File recalibratedBamIndex = "${outputBasename}.bai"
    File recalibratedBamChecksum = "${outputBasename}.bam.md5"
  }

  meta {
    taskDescription: "Concatenate efficiently BAM files that resulted from scattered parallel analysis (per-sample)."
  }
}