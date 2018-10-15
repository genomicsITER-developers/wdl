
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 20. Analyze patterns of covariation in the sequence dataset (per-sample per-lane)
# 21. Gather all the reports from BaseRecalibrator calls (per-sample per-lane)
# 22.
# 23.
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow BaseQualityScoreRecalibrationWF {

  File? inputBam
  File? inputBai

  File intervalTargetsFile # ¿Fichero con los intervalos de las regiones "targets"? ¿Haría falta agrupar intervalos para que no sean demasiados?
  Array[Array[String]] intervalTargetsList = read_tsv(intervalTargetsFile)

  File intervalContigsFile # ¿Fichero con los intervalos de "contigs" para el ApplyBQSR? ¿Agregar intervalos?
  Array[Array[String]] intervalContigsList = read_tsv(intervalContigsFile)

  File dbSnps
  File dbSnpsIdx
  File dbIndelsM1000g
  File dbIndelsM1000gIdx

  File refFasta
  File refDict
  File refIndex

  String sampleName

  # CONFIGURATION PARAMETERS
  Int firstStep
  Int lastStep

  # RESULTS DIR
  String resultsDir

  String gatkPath

  String javaOpts

  # Step 9 - Base recalibrator and gather the results
  if ((firstStep <= 9) && (9 <= lastStep)) {

    # ################################################################################# #
    # 20. Analyze patterns of covariation in the sequence dataset (per-sample per-lane) #
    # ################################################################################# #

    scatter (interval in intervalTargetsList) {
      call BaseRecalibrator {
        input:
          inputBam = inputBam,
          inputBai = inputBai,
          recalReportFilename = sampleName + "ready.deduped.parcial_recal_data.csv",
          sequenceGroupInterval = interval,
          dbSnps = dbSnps,
          dbSnpsIdx = dbSnpsIdx,
          dbIndelsM1000g = dbIndelsM1000g,
          dbIndelsM1000gIdx = dbIndelsM1000gIdx,
          refFasta = refFasta,
          refDict = refDict,
          refIndex = refIndex,
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }
    }

    # ############################################################################ #
    # 21. Gather all the reports from BaseRecalibrator calls (per-sample per-lane) #
    # ############################################################################ #

    call GatherBaseRecalReports {
      input:
        baseRecalReports = BaseRecalibrator.recalibrationReport,
        outputReportsFilename = sampleName + "ready.deduped.recal_data.csv",
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyBqsrReports {input: resultsDir = resultsDir, files = GatherBaseRecalReports.bqsrReports}

  }  # End step 9


  # Step 10 - Apply recalibration and gather the results
  if ((firstStep <= 10) && (10 <= lastStep)) {

    # ####################################################################### #
    # 22. Apply the recalibration to your sequence data (per-sample per-lane) #
    # ####################################################################### #

    scatter (interval in intervalContigsList) { # CAMBIAR ESTE INTERVAL-LIST POR EL QUE CORRESPONDA EN EL CASO DEL APPLYBQSR
      call ApplyBQSR {
        input:
          inputBam = inputBam,
          inputBai = inputBai,
          recalibrationReport = select_first([GatherBaseRecalReports.bqsrReports, resultsDir + "/" + sampleName + "ready.deduped.recal_data.csv"]),
          outputBasename = sampleName + ".ready.deduped.parcial_recalibrated",
          sequenceGroupInterval = interval,
          refFasta = refFasta,
          refDict = refDict,
          refIndex = refIndex,
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }
    }

    # ##################################################################### #
    # 23. Gather all the results from ApplyBQSR calls (per-sample per-lane) #
    # ##################################################################### #

    call GatherBamRecalibratedFiles {
      input:
        inputBams = ApplyBQSR.recalibratedBam,
        outputBasename = sampleName + ".ready.deduped.recalibrated",
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyRecalibratedBam {input: resultsDir = resultsDir, 
    files = [GatherBamRecalibratedFiles.recalibratedBam, GatherBamRecalibratedFiles.recalibratedBamIndex, GatherBamRecalibratedFiles.recalibratedBamChecksum]}

  }  # End step 10

  output {
    File? recalibrationReports = GatherBaseRecalReports.bqsrReports
    File? recalibratedBam = GatherBamRecalibratedFiles.recalibratedBam
    File? recalibratedBamIndex = GatherBamRecalibratedFiles.recalibratedBamIndex
    File? recalibratedBamChecksum = GatherBamRecalibratedFiles.recalibratedBamChecksum
  }

}


task BaseRecalibrator {

  File? inputBam
  File? inputBai

  String recalReportFilename

  Array[String] sequenceGroupInterval

  File dbSnps
  File dbSnpsIdx
  File dbIndelsM1000g
  File dbIndelsM1000gIdx

  File refFasta
  File refDict
  File refIndex

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" BaseRecalibrator \
      --reference ${refFasta} \
      --input ${inputBam} \
      --intervals ${sep=" --intervals " sequenceGroupInterval} \
      --use-original-qualities \
      --output ${recalReportFilename} \
      --known-sites ${dbSnps} \
      --known-sites ${dbIndelsM1000g}
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