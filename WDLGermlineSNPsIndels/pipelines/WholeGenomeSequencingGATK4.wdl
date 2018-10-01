#
# WORKFLOW: WholeGenomeSequencingGATK4.wdl
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
#   2. Quality Control (per sample, per lane)
#     2.1 ValidateSam...................................(step 8)
#     2.2 Qualimap......................................(step 9)
#     2.3 CollectRawWgsMetrics (3-0)....................(step 10)
#     2.4 CollectRawWgsMetrics (20-20)..................(step 11)
#     2.5 CollectWgsMetricsWithNonZeroCoverage..........(step 12)
#     2.6 ¿DepthOfCoverage?.............................(step 13)
#     2.7 CollectMultipleMetrics........................(step 14)
#     2.8 ¿CallableLoci?................................(step 15)
# 
# B. Per sample
#   
#   1. Variant Calling (per sample)
#     1.1 MarkDuplicates (MergeBamsPerSample)...........(step 16)
#     1.2 BQSR
#       1.2.1 BaseRecalibrator (por intervalos)*........(step 21)
#       1.2.2 GatherBQSRReports.........................(step 22)
#       1.2.3 ApplyBQSR (por intervalos)*...............(step 23)
#       1.2.4 GatherBamFiles............................(step 24)
#     1.3 HaplotypeCaller (por intervalos)*.............(step 25)
#     1.4 MergeVcfs.....................................(step 26)
# 
#   2. Quality Control (Con la salida del MergeBamsPerSample)
#     2.1 ValidateSam...................................(step 17)
#     2.2 ¿DepthOfCoverage?.............................(step 18)
#     2.3 CollectMultipleMetrics........................(step 19)
#     2.4 ¿CallableLoci?................................(step 20)
# 
# C. Multi-sample
# 
#   1. Joing Genotyping (multi-sample)
#     1.1 GenomicsDBImport (por intervalos)*............(step 21)
#     1.2 GenotypeGVCFs (por intervalos)*...............(step 22)
#     1.3 VariantFiltration (por intervalos)*...........(step 23)
#     1.4 MakeSitesOnlyVcf (por intervalos)*............(step 23)
#     1.5 GatherVcfsCloud...............................(step 24)
#     1.6 VQSR
#       1.6.1 SNPs
#         1.6.1.1 VariantRecalibratorCreateModel........(step 25)
#       1.6.2 Indels
#         1.6.2.1 VariantRecalibrator...................(step 28)
#       1.6.3 ApplyRecalibration (por intervalos)*......(step 29)
#       1.6.4 CollectVariantCallingMetrics..............(step 30)
#
# TODO:
#  - Añadir opciones de JAVA a las funciones
#  - ¿Añadir opción TMP_DIR para ficheros temporales?
#

# IMPORTS
import "SubWorkflows/SetWorkingDirectory.wdl"              as workDir
import "SubWorkflows/Preprocessing.wdl"                    as Preprocessing
import "SubWorkflows/QualityControl.wdl"                   as QC
import "SubWorkflows/VariantCallingPerSample.wdl"          as VCPerSample
import "SubWorkflows/QualityControlPerSample.wdl"          as QCPerSample
import "SubWorkflows/BaseQualityScoreRecalibration.wdl"    as BQSR
import "SubWorkflows/JointGenotyping.wdl"                  as JointGenotyping
#import "SubWorkflows/VariantQualityScoreRecalibration.wdl" as VQSR

# WORKFLOW DEFINITION

workflow WholeGenomeSequencingGATK4WF {

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
  File dbSnpsVcf
  File dbSnpsVcfIdx
  File hapMapResourceVcf
  File hapMapResourceVcfIndex
  File omniResourceVcf
  File omniResourceVcfIndex
  File oneThousandGenomesResourceVcf
  File oneThousandGenomesResourceVcfIndex
  File millsResourceVcf
  File millsResourceVcfIndex
  File axiomPolyResourceVcf
  File axiomPolyResourceVcfIndex
  Array[File] knownIndelsSitesVcfs
  Array[File] knownIndelsSitesIndices

  # READS
  # TSV read pais for per-sample, per-lane stage
  File fastqReadsTSV
  Array[Array[File]] pairReads = read_tsv(fastqReadsTSV)
  String libraryName      # ¿Pasar al TSV?
  String platform         # ¿Pasar al TSV?
  String sequencingCenter # ¿Pasar al TSV?

  # TSV read samples for per-sample stage
  File? readSamples

  # TSV Sample Name Map for multi-sample stage
  File? sampleNameMap

  # True for GenomicsDBImport, False for CombineGVCFs
  Boolean useGenomicsDB

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
  File? TSVIntervalsFile
  String protectionTag # Only for hg38

  # QUALIMAP
  String qualimapPath

  # RESULTS DIR
  String resultsDir

  # CONFIGURATION PARAMETERS
  Boolean wfPerSamplePerLane
  Boolean wfPreprocess
  Boolean wfQC_PerSamplePerLane

  Boolean wfPerSample
  Boolean wfVariantCalling
  Boolean wfBQSR
  Boolean wfQC_PerSample

  Boolean wfMultiSample
  Boolean wfJointGenotyping
  Boolean wfVQSR
  Boolean wfQC_MultiSample

  # Arreglar esto de los steps...
  Int firstStep
  Int lastStep

  # Java Opts
  String javaOpts

  # Interval List
  File intervalList
  Int scatterCount
  Int breakBandsAtMultiplesOf
  String intervalListOutDir

  # VQSR parameters
  Array[String] snpRecalibrationTrancheValues
  Array[String] snpRecalibrationAnnotationValues
  Array[String] indelRecalibrationTrancheValues
  Array[String] indelRecalibrationAnnotationValues
  Float snpFilterLevel
  Float indelFilterLevel
  Int snpVQSRDownsampleFactor
  Int indelsMaxGaussians
  Int snpsMaxGaussians

  # Fin declaración de variables


  # Step 0
  if ((firstStep <= 0) && (0 <= lastStep)) {
    call workDir.ConfigureResultsDirectory as BaseDirectory {input: resultsDir = resultsDir}
  }
  # End step 0

  ######### PER-SAMPLE PER-LANE #########

  if (wfPerSamplePerLane == true) {

    scatter (pair in pairReads) {

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
      call workDir.ConfigureResultsDirectory as ActualDirectory {input: resultsDir = actualPairDir}
      
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
        call QC.QualityControlWF as QCSubworkflow {
          input:
            refFasta = refFasta,
            bamFile = select_first([PreprocessingSubworkflow.fixedBam, actualPairDir + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam"]),
            resultsDir = actualPairDir,
            sampleName = sampleName,
            gatkPath = gatkPath,
            qualimapPath = qualimapPath,
            javaOpts = javaOpts
        }
      } # End QC
    }
  }

  ######### END PER-SAMPLE PER-LANE #########


  ######### PER-SAMPLE #########

  if (wfPerSample == true) {

    call workDir.FindFilesInDir as findMD5Bams {
      input:
        dir = resultsDir,
        pattern1 = resultsDir + "/*/*.aligned.merged.deduped.sorted.fixed.bam.md5",
        pattern2 = resultsDir + "/*/*.aligned.merged.deduped.sorted.fixed.bai"
    }

    #Array[File?]? md5Bams = if ((wfPreprocess == true) && (firstStep < 7)) then PreprocessingSubworkflow.fixedBamMD5 else findMD5Bams.files[0]

    call WritePerSampleTSV {
      input:
        fastqReadsTSV = fastqReadsTSV,
        resultsDir = resultsDir,
        bamSuffix = ".aligned.merged.deduped.sorted.fixed.bam",
        md5Bams = if ((wfPreprocess == true) && (firstStep < 7)) then PreprocessingSubworkflow.fixedBamMD5 else findMD5Bams.files[0]
    }

    call workDir.CopyResultsFilesToDir as CopyPerSampleTSV {input: resultsDir = resultsDir, files = WritePerSampleTSV.perSampleArray}

    scatter (sample in read_tsv(WritePerSampleTSV.perSampleArray)) {

      String dirName = sample[0]
      String inputBams = sample[1]

      String actualSampleDir = resultsDir + "/" + dirName

      if (wfVariantCalling == true) {
        call VCPerSample.VariantCallingPerSampleWF as VCPerSample {
          input:
            actualSampleDir = actualSampleDir,
            inputBams = inputBams,
            dirName = dirName,
            TSVIntervalsFile = TSVIntervalsFile,
            protectionTag = protectionTag,
            refFasta = refFasta,
            refAlt = refAlt,
            refIndex = refIndex,
            refDict = refDict,
            dbSnpsVcf = dbSnpsVcf,
            dbSnpsVcfIdx = dbSnpsVcfIdx,
            knownIndelsSitesVcfs = knownIndelsSitesVcfs,
            knownIndelsSitesIndices = knownIndelsSitesIndices,
            intervalList = intervalList,
            scatterCount = scatterCount,
            breakBandsAtMultiplesOf = breakBandsAtMultiplesOf,
            intervalListOutDir = intervalListOutDir,
            gatkPath = gatkPath,
            picardPath = picardPath,
            javaOpts = javaOpts
        }
      } # End Variant Calling

      if (wfQC_PerSample == true) {
        call QCPerSample.QualityControlPerSampleWF as QCPerSampleSubwf {
          input:
            refFasta = refFasta,
            bamFile = select_first([VCPerSample.bamsPerSample, actualSampleDir + "/" + dirName + ".ready.deduped.bam"]),
            resultsDir = actualSampleDir,
            sampleName = dirName,
            gatkPath = gatkPath,
            javaOpts = javaOpts
        }
      }

    }
  }

  ######### END PER-SAMPLE #########


  ######### MULTI-SAMPLE (JointGenotypingWF.wdl) #########

  if (wfMultiSample == true) {

    String multiSampleName = "multi-sample"
    String multiSampleDir = resultsDir + "/" + multiSampleName

    call workDir.ConfigureResultsDirectory as MultiSampleDirectory {input: resultsDir = multiSampleDir}

    if (!defined(sampleNameMap)) {

      # Create sampleNameMap file with the structure:
      #   HG00096  HG00096.g.vcf.gz
      #   NA19625 NA19625.g.vcf.gz
      #   HG00268 HG00268.g.vcf.gz

      call workDir.FindFilesInDir as findGVCFs {
        input:
          dir = resultsDir,
          pattern1 = resultsDir + "/*/*.ready.deduped.recalibrated.merged.g.vcf.gz.tbi",
          pattern2 = resultsDir + "/*/*.ready.deduped.recalibrated.merged.g.vcf.gz.md5"
      }

      File? perSampleTSV = select_first([WritePerSampleTSV.perSampleArray, resultsDir + "/read_samples.tsv"])

      call WriteSampleNameMapTSV {
        input:
          perSampleTSV = perSampleTSV,
          resultsDir = resultsDir,
          vcfSuffix = ".ready.deduped.recalibrated.merged.g.vcf.gz",
          gVcfsIndex = if ((wfPerSample == true) && (wfVariantCalling == true)) then VCPerSample.mergedGVCFsIndexes else findGVCFs.files[0]
      }

      call workDir.CopyResultsFilesToDir as CopySampleNameMapTSV {input: resultsDir = multiSampleDir, files = WriteSampleNameMapTSV.sampleNameMap}
    }

    # Multi-sample Variant Calling
    if (wfJointGenotyping == true) {

      call JointGenotyping.JointGenotypingWF as JGMultiSample {
        input:
          refFasta = refFasta,
          refIndex = refIndex,
          refDict = refDict,
          dbSnpsVcf = dbSnpsVcf,
          dbSnpsVcfIdx = dbSnpsVcfIdx,
          hapMapResourceVcf = hapMapResourceVcf,
          hapMapResourceVcfIndex = hapMapResourceVcfIndex,
          omniResourceVcf = omniResourceVcf,
          omniResourceVcfIndex = omniResourceVcfIndex,
          oneThousandGenomesResourceVcf = oneThousandGenomesResourceVcf,
          oneThousandGenomesResourceVcfIndex = oneThousandGenomesResourceVcfIndex,
          millsResourceVcf = millsResourceVcf,
          millsResourceVcfIndex = millsResourceVcfIndex,
          axiomPolyResourceVcf = axiomPolyResourceVcf,
          axiomPolyResourceVcfIndex = axiomPolyResourceVcfIndex,
          multiSampleName = multiSampleName,
          multiSampleDir = multiSampleDir,
          sampleNameMap = select_first([sampleNameMap, WriteSampleNameMapTSV.sampleNameMap]),
          intervalList = read_tsv(select_first([TSVIntervalsFile, VCPerSample.sequenceGrouping])),
          useGenomicsDB = useGenomicsDB,
          snpRecalibrationTrancheValues = snpRecalibrationTrancheValues,
          snpRecalibrationAnnotationValues = snpRecalibrationAnnotationValues,
          indelRecalibrationTrancheValues = indelRecalibrationTrancheValues,
          indelRecalibrationAnnotationValues = indelRecalibrationAnnotationValues,
          snpFilterLevel = snpFilterLevel,
          indelFilterLevel = indelFilterLevel,
          snpVQSRDownsampleFactor = snpVQSRDownsampleFactor,
          snpsMaxGaussians = snpsMaxGaussians,
          indelsMaxGaussians = indelsMaxGaussians,
          wfVQSR = wfVQSR,
          wfQC_MultiSample = wfQC_MultiSample,
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }

    } # End Variant Calling

  }

  ######### END MULTI-SAMPLE #########


  output {
    # PER-LANE, PER-SAMPLE
    Array[File?]? unmappedBam = PreprocessingSubworkflow.unmappedBam

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

    Array[File?]? sumaryQC = QCSubworkflow.summary
    Array[File?]? metrics3_0 = QCSubworkflow.metrics3_0
    Array[File?]? metrics20_20 = QCSubworkflow.metrics20_20
    Array[File?]? metricsNonZero = QCSubworkflow.metricsNonZero
    Array[File?]? metricsNonZeroPDF = QCSubworkflow.metricsNonZeroPDF
    Array[Array[File]?]? collectMultipleMetrics = QCSubworkflow.collectMultipleMetrics

    # PER-SAMPLE
    Array[File?]? bamsPerSample = VCPerSample.bamsPerSample
    Array[File?]? baisPerSample = VCPerSample.baisPerSample
    Array[File?]? md5sPerSample = VCPerSample.md5sPerSample
    Array[File?]? dupMetricsPerSample = VCPerSample.dupMetricsPerSample

    Array[File?]? sumaryQCPerSample = QCPerSampleSubwf.summary
    Array[Array[File]?]? collectMultipleMetricsPerSample = QCPerSampleSubwf.collectMultipleMetrics

    Array[File?]? sequenceGrouping = VCPerSample.sequenceGrouping
    Array[File?]? sequenceGroupingWithUnmapped = VCPerSample.sequenceGroupingWithUnmapped

    Array[File?]? recalibrationReports = VCPerSample.recalibrationReports
    Array[File?]? recalibratedBam = VCPerSample.recalibratedBam
    Array[File?]? recalibratedBamIndex = VCPerSample.recalibratedBamIndex
    Array[File?]? recalibratedBamChecksum = VCPerSample.recalibratedBamChecksum

    Array[File?]? mergedGVCFs = VCPerSample.mergedGVCFs
    Array[File?]? mergedGVCFsIndexes = VCPerSample.mergedGVCFsIndexes

    # MULTI-SAMPLE
    File? gatheredVCF = JGMultiSample.gatheredVCF
    File? gatheredVCFIndex = JGMultiSample.gatheredVCFIndex

    File? recalIndelsFile = JGMultiSample.recalIndelsFile
    File? recalIndexIndelsFile = JGMultiSample.recalIndexIndelsFile
    File? tranchesIndelsFile = JGMultiSample.tranchesIndelsFile
    File? rScriptIndels = JGMultiSample.rScriptIndels
    File? plotsIndelsFile = JGMultiSample.plotsIndelsFile

    File? recalSnpsFile = JGMultiSample.recalSnpsFile
    File? recalIndexSnpsFile = JGMultiSample.recalIndexSnpsFile
    File? tranchesSnpsFile = JGMultiSample.tranchesSnpsFile
    File? modelSnpsReport = JGMultiSample.modelSnpsReport
    File? rScriptSnps = JGMultiSample.rScriptSnps
    File? plotsSnpsFile = JGMultiSample.plotsSnpsFile

    File? recalibratedVCF = JGMultiSample.recalibratedVCF
  	File? recalibratedVCFIndex = JGMultiSample.recalibratedVCFIndex
  	File? recalibratedVCFMD5 = JGMultiSample.recalibratedVCFMD5

    File? detailMetrics = JGMultiSample.detailMetrics
    File? summaryMetrics = JGMultiSample.summaryMetrics
  }
}


task WritePerSampleTSV {

  File fastqReadsTSV
  String resultsDir
  String bamSuffix

  Array[File?]? md5Bams # Solo sirve para que espere por el proceso anterior

  command <<<
    python2 <<CODE
    import sys
    import csv
    import os
    with open("${fastqReadsTSV}", "r") as read_pairs_file, open("read_samples.tsv", "w") as bam_samples_file:
        tsvreader = csv.reader(read_pairs_file, delimiter="\t")
        d = {}
        for r in tsvreader:
            sample_name = r[2] + "_" + r[3]
            sub_dir = r[5]
            bam_file = sample_name + "${bamSuffix}"
            bam_file_path = os.path.join("${resultsDir}", sub_dir, bam_file)
            if os.path.isfile(bam_file_path):
                if not sample_name in d:
                    d[sample_name] = []
                d[sample_name].append(bam_file_path)

        tsvwriter = csv.writer(bam_samples_file, delimiter="\t")
        for elem in d:
            tsvwriter.writerow([elem] + [' '.join(d[elem])])
    CODE
  >>>

  output {
    File perSampleArray = "read_samples.tsv"
  }

  meta {
    taskDescription: "Create a new TSV file with per-sample information and files."
  }
}


task WriteSampleNameMapTSV {
  
  File? perSampleTSV

  String resultsDir
  String vcfSuffix

  Array[File?]? gVcfsIndex # Solo para que espere por el proceso anterior

  command <<<
    python2 <<CODE
    import sys
    import csv
    import os
    with open("${perSampleTSV}", "r") as per_sample_file, open("sample_name_map.tsv", "w") as sample_name_map:
        tsvreader = csv.reader(per_sample_file, delimiter="\t")
        d={}
        for r in tsvreader:
            sample_name = r[0]
            vcf_file = r[0] + "${vcfSuffix}"
            vcf_file_path = os.path.join("${resultsDir}", sample_name, vcf_file)
            if os.path.isfile(vcf_file_path):
                d[sample_name] = vcf_file_path

        tsvwriter = csv.writer(sample_name_map, delimiter="\t")
        for elem in d:
            tsvwriter.writerow([elem] + [d[elem]])
    CODE
  >>>

  output {
    File sampleNameMap = "sample_name_map.tsv"
  }

  meta {
    taskDescription: "Create a Sample Name Map TSV file for multi-sample stage with the next structure:\n\tHG00096\tHG00096.g.vcf.gz\n\tNA19625\tNA19625.g.vcf.gz\n\tHG00268\tHG00268.g.vcf.gz"
  }
}