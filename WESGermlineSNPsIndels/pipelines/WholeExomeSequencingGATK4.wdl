#
# WORKFLOW: WholeExomeSequencingGATK4.wdl
#
# Workflow:
# A. Per sample, per lane
# 
#   1. Preprocessing
#     1.1 FastqToSam....................................(step 1)
#     1.2 SamToFastq....................................(step 2)
#     1.3 Bwa Mem.......................................(step 3)
#     1.4 MergeBamAlignments............................(step 4)
#     1.5 MarkDuplicates................................(step 5)
#     1.6 SortSam.......................................(step 6)
#     1.7 SetNmMdAndUqTags..............................(step 7)
#
# TODO:
#  - ¿Añadir opción TMP_DIR para ficheros temporales?
#

# IMPORTS
import "SubWorkflows/Utils.wdl"                            as utils
import "SubWorkflows/Preprocessing.wdl"                    as Preprocessing
import "SubWorkflows/QualityControl.wdl"                   as QualityControl
import "SubWorkflows/BaseQualityScoreRecalibration.wdl"    as BQSR
import "SubWorkflows/VariantCalling.wdl"                   as VariantCalling

# WORKFLOW DEFINITION

workflow WholeExomeSequencingGATK4WF {

  # REFERENCE
  File refFasta
  String refBaseDir
  String refBaseName
  File? refAlt
  File refIndex
  File refDict
  File refAmb
  File refAnn
  File refBwt
  File refPac
  File refSa

  File dbSnps
  File dbSnpsIdx
  File dbIndelsM1000g
  File dbIndelsM1000gIdx

  # READS - TSV read pais for per-sample, per-lane stage
  File fastqReadsTSV
  Array[Array[File]] pairReads = read_tsv(fastqReadsTSV)
  String libraryName      # ¿Pasar al TSV?
  String platform         # ¿Pasar al TSV?
  String sequencingCenter # ¿Pasar al TSV?

  # Illumina TruSeq Rapid Exome Capture Kit manifest for QUALIMAP:
  File bedFile
  File bedFilePadding

  File intervalList
  File baits
  File targets

  # Interval files for BQSR
  File intervalTargetsFile
  File intervalContigsFile

  File TSVIntervalsFile
  Float? contamination
  Int maxAltAlleles
  Int intervalPadding

  # BWA
  String bwaPath
  String bwaMemCommand

  # SAMTOOLS
  String samtoolsPath

  # GATK
  String gatkPath

  # PICARD
  String picardPath

  # PYTHON
  String pythonPath

  # QUALIMAP
  String qualimapPath
  String javaMemSize

  # RESULTS DIR
  String resultsDir

  # CONFIGURATION PARAMETERS
  Boolean wfPerSamplePerLane
  Boolean wfPreprocess
  Boolean wfQC_PerSamplePerLane
  Boolean wfBQSR

  Boolean wfPerSample
  Boolean wfVariantCalling
  Boolean wfQC_PerSample

  # Fix this steps thing...
  Int firstStep
  Int lastStep

  # Java Opts
  String javaOpts

  # Fin declaración de variables


  # Step 0
  if ((firstStep <= 0) && (0 <= lastStep)) {
    call utils.ConfigureResultsDirectory as BaseDirectory {input: resultsDir = resultsDir}
  }
  # End step 0


  ######### PER-SAMPLE PER-LANE #########

  if (wfPerSamplePerLane == true) {

    scatter (pair in pairReads) {
                                     # Example:
      String flowcell      = pair[0] # H814YADXX
      String lane          = pair[1] # 1
      String sampleID      = pair[2] # U0a
      String index         = pair[3] # CGATGT
      String sampleProject = pair[4] # RM8398
      String sampleDir     = pair[5] # U0a_CGATGT_L001_001
      File readsFile1      = pair[6] # U0a_CGATGT_L001_R1_001.chr22.fastq.gz
      File readsFile2      = pair[7] # U0a_CGATGT_L001_R2_001.chr22.fastq.gz

      String readGroupName = flowcell + "." + lane + "." + index # H814YADXX.1.CGATGT
      String platformUnit  = flowcell + "." + lane + "." + index # H814YADXX.1.CGATGT

      String sampleName = sampleID + "_" + index

      String actualPairDir = resultsDir + "/" + sampleDir
      call utils.ConfigureResultsDirectory as ActualDirectory {input: resultsDir = actualPairDir}

      if (wfPreprocess == true) {
        call Preprocessing.PreprocessingWF as PreprocessingSubworkflow {
          input:
            firstStep = firstStep,
            lastStep = lastStep,
            refFasta = refFasta,
            refBaseDir = refBaseDir,
            refBaseName = refBaseName,
            refAlt = refAlt,
            refIndex = refIndex,
            refDict = refDict,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            readsFile1 = readsFile1,
            readsFile2 = readsFile2,
            sampleName = sampleName,
            readGroupName = readGroupName,
            libraryName = libraryName,
            platformUnit = platformUnit,
            platform = platform,
            sequencingCenter = sequencingCenter,
            gatkPath = gatkPath,
            bwaPath = bwaPath,
            bwaMemCommand = bwaMemCommand,
            resultsDir = ActualDirectory.directory,
            javaOpts = javaOpts
        }
      } # End Preprocessing

      if (wfQC_PerSamplePerLane == true) {
        call QualityControl.QualityControlWF as QCSubworkflow {
          input:
            firstStep = firstStep,
            lastStep = lastStep,
            refFasta = refFasta,
            bamFile = select_first([PreprocessingSubworkflow.fixedBam, actualPairDir + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam"]),
            bedFile = bedFile,
            bedFilePadding = bedFilePadding,
            intervalList = intervalList,
            baits = baits,
            targets = targets,
            resultsDir = actualPairDir,
            sampleName = sampleName,
            gatkPath = gatkPath,
            qualimapPath = qualimapPath,
            javaOpts = javaOpts,
            javaMemSize = javaMemSize
        }
      } # End QC

      if (wfBQSR == true) {
        call BQSR.BaseQualityScoreRecalibrationWF as BQSRSubworkflow {
          input:
            inputBam = select_first([PreprocessingSubworkflow.fixedBam, actualPairDir + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam"]),
            inputBai = select_first([PreprocessingSubworkflow.fixedBamIndex, actualPairDir + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bai"]),
            intervalTargetsFile = intervalTargetsFile, # CAMBIAR
            intervalContigsFile = intervalContigsFile,
            dbSnps = dbSnps,
            dbSnpsIdx = dbSnpsIdx,
            dbIndelsM1000g = dbIndelsM1000g,
            dbIndelsM1000gIdx = dbIndelsM1000gIdx,
            refFasta = refFasta,
            refDict = refDict,
            refIndex = refIndex,
            sampleName = sampleName,
            firstStep = firstStep,
            lastStep = lastStep,
            resultsDir = resultsDir,
            gatkPath = gatkPath,
            javaOpts = javaOpts
        }
      } # End BQSR

    }

  }

  ######### END PER-SAMPLE PER-LANE #########


  ######### PER-SAMPLE #########

  if (wfPerSample == true) {

    call utils.FindFilesInDir as findMD5Bams {
      input:
        dir = resultsDir,
        pattern1 = resultsDir + "/*/*.aligned.merged.deduped.sorted.fixed.bam.md5",
        pattern2 = resultsDir + "/*/*.aligned.merged.deduped.sorted.fixed.bai"
    }

    call utils.WritePerSampleTSV as writeTSVPerSample {
      input:
        fastqReadsTSV = fastqReadsTSV,
        resultsDir = resultsDir,
        bamSuffix = ".aligned.merged.deduped.sorted.fixed.bam",
        #md5Bams = select_first([PreprocessingSubworkflow.fixedBamMD5, actualPairDir + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam.md5"])
        md5Bams = if ((wfPreprocess == true) && (firstStep < 9)) then PreprocessingSubworkflow.fixedBamMD5 else findMD5Bams.files[0]
    }

    call utils.CopyResultsFilesToDir as CopyPerSampleTSV {input: resultsDir = resultsDir, files = writeTSVPerSample.perSampleArray}

    scatter (sample in read_tsv(writeTSVPerSample.perSampleArray)) {

      String dirName = sample[0]
      String inputBams = sample[1]

      String actualSampleDir = resultsDir + "/" + dirName

      if (wfVariantCalling == true) {
        call VariantCalling.VariantCallingWF as VCSubworkflow {
          input:
            actualSampleDir = actualSampleDir,
            inputBams = inputBams,
            dirName = dirName,
            refFasta = refFasta,
            refAlt = refAlt,
            refIndex = refIndex,
            refDict = refDict,
            TSVIntervalsFile = TSVIntervalsFile,
            contamination = contamination,
            maxAltAlleles = maxAltAlleles,
            intervalPadding = intervalPadding,
            firstStep = firstStep,
            lastStep = lastStep,
            wfQC_PerSample = wfQC_PerSample,
            gatkPath = gatkPath,
            javaOpts = javaOpts
        }
      } # End VC
    }
  }

  ######### END PER-SAMPLE #########


  ######### MULTI-SAMPLE (JointGenotypingWF.wdl) #########

  #if (wfMultiSample == true) {}

  ######### END MULTI-SAMPLE #########


  output {
    # PER-LANE, PER-SAMPLE
    Array[File?]? unmappedBam = PreprocessingSubworkflow.unmappedBam

    Array[File?]? unmappedTaggedBam = PreprocessingSubworkflow.unmappedTaggedBam
    Array[File?]? metricsIlluminaAdapters = PreprocessingSubworkflow.metricsIlluminaAdapters

    Array[File?]? convertedFastq = PreprocessingSubworkflow.convertedFastq
    Array[File?]? alignedBam = PreprocessingSubworkflow.alignedBam
    Array[File?]? alignedUnsortedBam = PreprocessingSubworkflow.alignedUnsortedBam

    Array[File?]? markedBam = PreprocessingSubworkflow.markedBam
    Array[File?]? markedBai = PreprocessingSubworkflow.markedBai
    Array[File?]? markedMD5 = PreprocessingSubworkflow.markedMD5
    Array[File?]? duplicateMetrics = PreprocessingSubworkflow.duplicateMetrics

    Array[File?]? sortedBam = PreprocessingSubworkflow.sortedBam
    Array[File?]? sortedBamIndex = PreprocessingSubworkflow.sortedBamIndex
    Array[File?]? sortedBamMD5 = PreprocessingSubworkflow.sortedBamMD5

    Array[File?]? fixedBam = PreprocessingSubworkflow.fixedBam
    Array[File?]? fixedBamIndex = PreprocessingSubworkflow.fixedBamIndex
    Array[File?]? fixedBamMD5 = PreprocessingSubworkflow.fixedBamMD5

    Array[File?]? summary = QCSubworkflow.summary
    Array[Array[File]?]? multMetrics = QCSubworkflow.multMetrics
    Array[File?]? rawMetrics = QCSubworkflow.rawMetrics
    Array[File?]? wgsMetrics = QCSubworkflow.wgsMetrics
    Array[File?]? wgsNonZeroMetrics = QCSubworkflow.wgsNonZeroMetrics
    Array[File?]? wgsNonZeroMetricsPDF = QCSubworkflow.wgsNonZeroMetricsPDF
    Array[File?]? hsMetrics = QCSubworkflow.hsMetrics
    Array[File?]? hsPerTargetMetrics = QCSubworkflow.hsPerTargetMetrics
    Array[File?]? oxoGMetrics = QCSubworkflow.oxoGMetrics

    Array[File?]? recalibrationReports = BQSRSubworkflow.recalibrationReports
    Array[File?]? recalibratedBam = BQSRSubworkflow.recalibratedBam
    Array[File?]? recalibratedBamIndex = BQSRSubworkflow.recalibratedBamIndex
    Array[File?]? recalibratedBamChecksum = BQSRSubworkflow.recalibratedBamChecksum

    # PER-SAMPLE
    Array[File?]? bamsPerSample = VCSubworkflow.bamsPerSample
    Array[File?]? baisPerSample = VCSubworkflow.baisPerSample
    Array[File?]? md5sPerSample = VCSubworkflow.md5sPerSample
    Array[File?]? dupMetricsPerSample = VCSubworkflow.dupMetricsPerSample

    Array[File?]? summaryVC = VCSubworkflow.summary
    Array[Array[File]?]? multMetricsVC = VCSubworkflow.multMetrics

    Array[File?]? mergedGVCFs = VCSubworkflow.mergedGVCFs
    Array[File?]? mergedGVCFsIndexes = VCSubworkflow.mergedGVCFsIndexes

    # MULTI-SAMPLE

  }
}
