#
# WORKFLOW: WholeExomeSequencingGATK4.wdl
#
# Workflow:
#
# A. wfPerSamplePerLane:
#
#  1. Preprocessing:
#
#     Step 1: Fastq to uBAM:
#        1.1 FastqsToUnmappedBam..............................................................1
#     Step 2: Mark Illumina Adapters:
#        2.1 MarkIlluminaAdapters.............................................................2
#
#     Step 3: SamToFastq, BWA-MEM, and MergeBamAlignment to generate a clean BAM:
#        3.1 UnmappedBamToFastq...............................................................3
#        3.2 BwaMem...........................................................................4
#        3.3 MergeBamAlignment................................................................5
#
#     Step 4: Mark duplicates:
#        4.1 MarkDuplicates...................................................................6
#
#     Step 5: SortSam:
#        5.1 SortSam..........................................................................7
#        5.2 FixBamTags.......................................................................8
#
#
#  2. QualityControl (Per Sample, Per Lane):
#
#     Step 6: Validate:
#        6.1 ValidateReadyBam.................................................................9
#
#     Step 7: Qualimap:
#        7.1 QualimapZeroPadding, QualimapWithPadding.........................................10
#
#     Step 8: Collect metrics:
#        8.1 CollectMultipleMetrics...........................................................11
#        8.2 CollectRawWgsMetrics, CollectWgsMetrics, CollectWgsMetricsWithNonZeroCoverage....12
#        8.3 CollectHsMetrics.................................................................13
#        8.4 CollectOxoGMetrics...............................................................14
#        8.5 DepthOfCoverage??................................................................15
#        8.6 CallableLoci??...................................................................16
#        8.7 FindCoveredIntervals??...........................................................17
#        8.8 DiagnoseTargets??................................................................18
#        8.9 QualifyMissingIntervals??........................................................19
#
#
#  3. BaseQualityScoreRecalibration:
#
#     Step 9: Base recalibrator and gather the results:
#        9.1 BaseRecalibrator.................................................................20
#        9.2 GatherBaseRecalReports...........................................................21
#
#     Step 10: Apply recalibration and gather the results:
#        10.1 ApplyBQSR.......................................................................22
#        10.2 GatherBamRecalibratedFiles......................................................23
#
#
# B. wfPerSample:
#
#  Step 11: Merge Bams per sample:
#     11.1 MarkDuplicates........................................................................24
#
#  Step 12: QualityControl (Per Sample):
#     12.1 ValidateMergedBam...............................................................25
#     12.2 CollectMultipleMetrics..........................................................26
#     12.3 QualimapZeroPadding, QualimapWithPadding........................................27
#     12.4 DepthOfCoverage??...............................................................28
#     12.5 CallableLoci??..................................................................29
#
#  4. VariantCalling:
#
#     Step 13:  Haplotype Caller and Merge results:
#        13.1 HaplotypeCaller.................................................................30
#        13.2 MergeGVCFs......................................................................31
#
#
# C.................................................... wfMultiSample:
#
#  5. JointGenotyping:
#
#     Step 14: Joint Genotyping:
#        IF genomicsDBImport:
#           14.1 ImportGVCFs..................................................................32
#           14.2 GenotypeGVCFsWithGenomicsDB..................................................33
#           14.3 HardFilterAndMakeSitesOnlyVcfWithGenomicsDB..................................34
#
#        ELSE IF CombineGVCFs:
#           14.4 CombineGVCFs.................................................................35
#           14.5 GenotypeGVCFsWithCohortGVCF..................................................36
#           14.6 HardFilterAndMakeSitesOnlyVcfWithCohortGVCF..................................37
#
#        14.7 Joint Genotyping: SitesOnlyGatherVcf............................................38
#
#
#  6. VQSR:
#     Step 15: VQSR for indels and snps:
#        15.1 IndelsVariantRecalibrator.......................................................39
#        15.2 SNPsVariantRecalibratorCreateModel..............................................40
#
#     Step 16: Apply VQSR:
#        16.1 ApplyRecalibration..............................................................41
#
#
#  7. QualityControl (Multi-sample):
#
#     Step 17: Collect VC metrics:
#        17.1 CollectVariantCallingMetrics....................................................42
#
#     Step 18: Variant callset evaluation:
#        18.1 VariantEval??...................................................................43
#
#
#  8. Annotation:
#
#     Step 19: Annotation with Annovar, SnpEff and GATK
#        19.1 table_annovar.pl................................................................44
#        19.2 snpeff..........................................................................45
#        19.3 VariantAnnotator................................................................46
#

# IMPORTS
import "SubWorkflows/Utils.wdl"                            as utils
import "SubWorkflows/Preprocessing.wdl"                    as Preprocessing
import "SubWorkflows/QualityControl.wdl"                   as QualityControl
import "SubWorkflows/BaseQualityScoreRecalibration.wdl"    as BQSR
import "SubWorkflows/VariantCalling.wdl"                   as VariantCalling
import "SubWorkflows/JointGenotyping.wdl"                  as JointGenotyping
import "SubWorkflows/VariantQualityScoreRecalibration.wdl" as VQSR
import "SubWorkflows/QualityControlMultiSample.wdl"        as QC_MultiSample
import "SubWorkflows/Annotation.wdl"                       as Annotation

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

  # FOR HG19
  File dbSnps
  File dbSnpsIdx
  File snp1000g
  File snp1000gIdx
  File hapMapResource
  File hapMapResourceIndex
  File omniResource
  File omniResourceIndex
  File millsResource
  File millsResourceIndex

  # FOR HG38
  #File dbSnps
  #File dbSnpsIdx
  #File dbIndelsM1000g
  #File dbIndelsM1000gIdx
  #File hapMapResource
  #File hapMapResourceIndex
  #File omniResource
  #File omniResourceIndex
  #File oneThousandGenomesResource
  #File oneThousandGenomesResourceIndex
  #File millsResource
  #File millsResourceIndex
  #File axiomPolyResource
  #File axiomPolyResourceIndex

  Array[String] indelRecalibrationTrancheValues
  Array[String] indelRecalibrationAnnotationValues
  Array[String] snpRecalibrationTrancheValues
  Array[String] snpRecalibrationAnnotationValues

  String protectionTag

  Float indelFilterLevel
  Float snpFilterLevel

  Int snpVQSRDownsampleFactor

  Int indelsMaxGaussians
  Int snpsMaxGaussians

  # READS - TSV read pais for per-sample, per-lane stage
  File fastqReadsTSV
  Array[Array[File]] pairReads = read_tsv(fastqReadsTSV)
  String libraryName      # Pasar al TSV??
  String platform         # Pasar al TSV??
  String sequencingCenter # Pasar al TSV??

  # True for GenomicsDBImport, False for CombineGVCFs
  Boolean useGenomicsDB

  # Illumina TruSeq Rapid Exome Capture Kit manifest for QUALIMAP:
  File bedFile
  File bedFilePadding

  File intervalList
  File baits
  File targets

  # Interval files for BQSR
  File? intervalTargetsFile
  File? intervalContigsFile

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
  String javaOpts
  String gatkBaseCommand = gatkPath + ' --java-options ' + '"' + javaOpts + '"' + ' '

  # PICARD
  String picardPath

  # PYTHON
  String pythonPath

  # QUALIMAP
  String qualimapPath
  String javaMemSize

  # ANNOVAR
  String annovarPath
  String humanDB

  # SNPEFF
  String snpEffPath

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

  Boolean wfMultiSample
  Boolean wfJointGenotyping
  Boolean wfVQSR
  Boolean wfQC_MultiSample
  Boolean wfAnnotation

  # Fix this steps thing...
  Int firstStep
  Int lastStep

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

      call utils.ConfigureResultsDirectory as ActualDirectory {input: resultsDir = resultsDir + "/" + sampleDir}

      if (wfPreprocess == true) {
        call Preprocessing.PreprocessingWF as PreprocessingSubworkflow {
          input:
            firstStep        = firstStep,
            lastStep         = lastStep,
            refFasta         = refFasta,
            refBaseDir       = refBaseDir,
            refBaseName      = refBaseName,
            refAlt           = refAlt,
            refIndex         = refIndex,
            refDict          = refDict,
            refAmb           = refAmb,
            refAnn           = refAnn,
            refBwt           = refBwt,
            refPac           = refPac,
            refSa            = refSa,
            readsFile1       = readsFile1,
            readsFile2       = readsFile2,
            sampleName       = sampleName,
            readGroupName    = readGroupName,
            libraryName      = libraryName,
            platformUnit     = platformUnit,
            platform         = platform,
            sequencingCenter = sequencingCenter,
            gatkBaseCommand  = gatkBaseCommand,
            bwaPath          = bwaPath,
            bwaMemCommand    = bwaMemCommand,
            resultsDir       = ActualDirectory.directory
        }
      } # End Preprocessing

      String qcDir = ActualDirectory.directory + "/qc"
      call utils.ConfigureResultsDirectory as QCDirectory {input: resultsDir = qcDir}

      if (wfQC_PerSamplePerLane == true) {
        call QualityControl.QualityControlWF as QCSubworkflow {
          input:
            firstStep       = firstStep,
            lastStep        = lastStep,
            refFasta        = refFasta,
            refDict         = refDict,
            refIndex        = refIndex,
            bamFile         = select_first([PreprocessingSubworkflow.fixedBam, ActualDirectory.directory + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam"]),
            bedFile         = bedFile,
            bedFilePadding  = bedFilePadding,
            intervalList    = intervalList,
            baits           = baits,
            targets         = targets,
            resultsDir      = qcDir,
            sampleName      = sampleName,
            gatkBaseCommand = gatkBaseCommand,
            qualimapPath    = qualimapPath,
            javaMemSize     = javaMemSize
        }
      } # End QC

      if (wfBQSR == true) {

        call utils.CreateSequenceGroupingTSV {input: refDict = refDict, protectionTag = protectionTag}

        call utils.CopyResultsFilesToDir as CopyTSVIntervals {input: resultsDir = ActualDirectory.directory,
          files = [CreateSequenceGroupingTSV.sequenceGrouping, CreateSequenceGroupingTSV.sequenceGroupingWithUnmapped]}

        call BQSR.BaseQualityScoreRecalibrationWF as BQSRSubworkflow {
          input:
            inputBam            = select_first([PreprocessingSubworkflow.fixedBam,      ActualDirectory.directory + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam"]),
            inputBai            = select_first([PreprocessingSubworkflow.fixedBamIndex, ActualDirectory.directory + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bai"]),
            targets             = targets,
            sequenceGroup       = select_first([intervalContigsFile, CreateSequenceGroupingTSV.sequenceGrouping]),
            dbSnps              = dbSnps,
            dbSnpsIdx           = dbSnpsIdx,
            millsResource       = millsResource,
            millsResourceIndex  = millsResourceIndex,
            refFasta            = refFasta,
            refDict             = refDict,
            refIndex            = refIndex,
            sampleName          = sampleName,
            firstStep           = firstStep,
            lastStep            = lastStep,
            resultsDir          = ActualDirectory.directory,
            gatkPath            = gatkPath,
            javaOpts            = javaOpts
        }
      } # End BQSR

    }

  }

  ######### END PER-SAMPLE PER-LANE #########


  ######### PER-SAMPLE #########

  if (wfPerSample == true) {

    call utils.FindFilesInDir as findMD5Bams {
      input:
        dir      = resultsDir,
        pattern1 = resultsDir + "/*/*.aligned.merged.deduped.sorted.fixed.bam.md5",
        pattern2 = resultsDir + "/*/*.aligned.merged.deduped.sorted.fixed.bai"
    }

    call utils.WritePerSampleTSV as writeTSVPerSample {
      input:
        fastqReadsTSV = fastqReadsTSV,
        resultsDir    = resultsDir,
        bamSuffix     = ".aligned.merged.deduped.sorted.fixed.bam",
        #md5Bams      = select_first([PreprocessingSubworkflow.fixedBamMD5, ActualDirectory.directory + "/" + sampleName + ".aligned.merged.deduped.sorted.fixed.bam.md5"])
        md5Bams       = if ((wfPreprocess == true) && (firstStep < 9)) then PreprocessingSubworkflow.fixedBamMD5 else findMD5Bams.files[0]
    }

    call utils.CopyResultsFilesToDir as CopyPerSampleTSV {input: resultsDir = resultsDir, files = writeTSVPerSample.perSampleArray}

    scatter (sample in read_tsv(writeTSVPerSample.perSampleArray)) {

      String dirName = sample[0]
      String inputBams = sample[1]

      String actualSampleDir = resultsDir + "/" + dirName

      if (wfVariantCalling == true) {
        call VariantCalling.VariantCallingWF as VCSubworkflow {
          input:
            actualSampleDir  = actualSampleDir,
            inputBams        = inputBams,
            dirName          = dirName,
            refFasta         = refFasta,
            refAlt           = refAlt,
            refIndex         = refIndex,
            refDict          = refDict,
            bedFile          = bedFile,
            bedFilePadding   = bedFilePadding,
            TSVIntervalsFile = TSVIntervalsFile,
            contamination    = contamination,
            maxAltAlleles    = maxAltAlleles,
            intervalPadding  = intervalPadding,
            firstStep        = firstStep,
            lastStep         = lastStep,
            wfQC_PerSample   = wfQC_PerSample,
            gatkPath         = gatkPath,
            qualimapPath     = qualimapPath,
            javaOpts         = javaOpts,
            javaMemSize      = javaMemSize
        }
      } # End VC
    }
  }

  ######### END PER-SAMPLE #########


  ######### MULTI-SAMPLE #########

  if (wfMultiSample == true) {

    String multiSampleName = "multi-sample"
    String multiSampleDir  = resultsDir + "/" + multiSampleName

    call utils.ConfigureResultsDirectory as MultiSampleDirectory {input: resultsDir = multiSampleDir}

    # Create sampleNameMap file with the structure:
    #   HG00096  HG00096.g.vcf.gz
    #   NA19625 NA19625.g.vcf.gz
    #   HG00268 HG00268.g.vcf.gz

    call utils.FindFilesInDir as findGVCFs {
      input:
        dir      = resultsDir,
        pattern1 = resultsDir + "/*/*.ready.deduped.recalibrated.merged.g.vcf.gz.tbi",
        pattern2 = resultsDir + "/*/*.ready.deduped.recalibrated.merged.g.vcf.gz.md5"
    }

    File? perSampleTSV = select_first([writeTSVPerSample.perSampleArray, resultsDir + "/read_samples.tsv"])

    call utils.WriteSampleNameMapTSV as writeSampleNameMap {
      input:
        perSampleTSV = perSampleTSV,
        resultsDir   = resultsDir,
        vcfSuffix    = ".ready.deduped.recalibrated.merged.g.vcf.gz",
        gVcfsIndex   = if ((wfPerSample == true) && (wfVariantCalling == true)) then VCSubworkflow.mergedGVCFsIndexes else findGVCFs.files[0]
    }

    call utils.CopyResultsFilesToDir as CopySampleNameMapTSV {input: resultsDir = multiSampleDir, files = writeSampleNameMap.sampleNameMap}

    if (wfJointGenotyping == true) {
      call JointGenotyping.JointGenotypingWF as JGSubworkflow {
        input:
          refFasta        = refFasta,
          refIndex        = refIndex,
          refDict         = refDict,
          dbSnpsVcf       = dbSnps,
          dbSnpsVcfIdx    = dbSnpsIdx,
          multiSampleName = multiSampleName,
          multiSampleDir  = multiSampleDir,
          sampleNameMap   = writeSampleNameMap.sampleNameMap,
          intervalList    = intervalList,
          useGenomicsDB   = useGenomicsDB,
          firstStep       = firstStep,
          lastStep        = lastStep,
          gatkPath        = gatkPath,
          javaOpts        = javaOpts
      }
    } # End JointGenotyping

#    if (wfVQSR == true) {
#      call VQSR.VariantQualityScoreRecalibrationWF as VQSRSubworkflow {
#        input:
#          refFasta                           = refFasta,
#          refIndex                           = refIndex,
#          refDict                            = refDict,
#          dbSnps                             = dbSnps,
#          dbSnpsIdx                          = dbSnpsIdx,
#          hapMapResource                     = hapMapResource,
#          hapMapResourceIndex                = hapMapResourceIndex,
#          omniResource                       = omniResource,
#          omniResourceIndex                  = omniResourceIndex,
#          oneThousandGenomesResource         = oneThousandGenomesResource,
#          oneThousandGenomesResourceIndex    = oneThousandGenomesResourceIndex,
#          millsResource                      = millsResource,
#          millsResourceIndex                 = millsResourceIndex,
#          axiomPolyResource                  = axiomPolyResource,
#          axiomPolyResourceIndex             = axiomPolyResourceIndex,
#          gatheredVCF                        = select_first([JGSubworkflow.gatheredVCF,      multiSampleDir + "/" + multiSampleName + ".cohort.genotyped.filtered.sites_only.gathered.vcf.gz"]),
#          gatheredVCFIndex                   = select_first([JGSubworkflow.gatheredVCFIndex, multiSampleDir + "/" + multiSampleName + ".cohort.genotyped.filtered.sites_only.gathered.vcf.gz.tbi"]),
#          indelRecalibrationTrancheValues    = indelRecalibrationTrancheValues,
#          indelRecalibrationAnnotationValues = indelRecalibrationAnnotationValues,
#          snpRecalibrationTrancheValues      = snpRecalibrationTrancheValues,
#          snpRecalibrationAnnotationValues   = snpRecalibrationAnnotationValues,
#          snpVQSRDownsampleFactor            = snpVQSRDownsampleFactor,
#          indelFilterLevel                   = indelFilterLevel,
#          snpFilterLevel                     = snpFilterLevel,
#          indelsMaxGaussians                 = indelsMaxGaussians,
#          snpsMaxGaussians                   = snpsMaxGaussians,
#          multiSampleName                    = multiSampleName,
#          multiSampleDir                     = multiSampleDir,
#          firstStep                          = firstStep,
#          lastStep                           = lastStep,
#          gatkPath                           = gatkPath,
#          javaOpts                           = javaOpts
#      }
#    } # End VQSR


#    if (wfQC_MultiSample == true) {
#      call QC_MultiSample.QualityControlMultiSampleWF as QCMultiSampleSubworkflow {
#        input:
#          refDict         = refDict,
#          recalVCF        = select_first([VQSRSubworkflow.recalibratedVCF, multiSampleDir + "/" + multiSampleName + ".SNP_INDEL.recalibrated.vcf.gz"]),
#          dbSnp           = dbSnps,
#          multiSampleName = multiSampleName,
#          multiSampleDir  = multiSampleDir,
#          firstStep       = firstStep,
#          lastStep        = lastStep,
#          gatkPath        = gatkPath,
#          javaOpts        = javaOpts
#      }
#    } # End Multisample QC


#    if (wfAnnotation == true) {
#      call Annotation.AnnotationWF as AnnotationSubworkflow {
#        input:
#          refFasta        = refFasta,
#          refIndex        = refIndex,
#          refDict         = refDict,
#          recalVCF        = select_first([VQSRSubworkflow.recalibratedVCF, multiSampleDir + "/" + multiSampleName + ".SNP_INDEL.recalibrated.vcf.gz"]),
#          multiSampleName = multiSampleName,
#          multiSampleDir  = multiSampleDir,
#          annovarPath     = annovarPath,
#          humanDB         = humanDB,
#          snpEffPath      = snpEffPath,
#          gatkPath        = gatkPath,
#          javaOpts        = javaOpts,
#          firstStep       = firstStep,
#          lastStep        = lastStep
#      }
#    } # End Annotation

  }
  ######### END MULTI-SAMPLE #########


  output {
    # PER-LANE, PER-SAMPLE
    Array[File?]? unmappedBam             = PreprocessingSubworkflow.unmappedBam

    Array[File?]? unmappedTaggedBam       = PreprocessingSubworkflow.unmappedTaggedBam
    Array[File?]? metricsIlluminaAdapters = PreprocessingSubworkflow.metricsIlluminaAdapters

    Array[File?]? convertedFastq          = PreprocessingSubworkflow.convertedFastq
    Array[File?]? alignedBam              = PreprocessingSubworkflow.alignedBam
    Array[File?]? alignedUnsortedBam      = PreprocessingSubworkflow.alignedUnsortedBam

    Array[File?]? markedBam               = PreprocessingSubworkflow.markedBam
    Array[File?]? markedBai               = PreprocessingSubworkflow.markedBai
    Array[File?]? markedMD5               = PreprocessingSubworkflow.markedMD5
    Array[File?]? duplicateMetrics        = PreprocessingSubworkflow.duplicateMetrics

    Array[File?]? sortedBam               = PreprocessingSubworkflow.sortedBam
    Array[File?]? sortedBamIndex          = PreprocessingSubworkflow.sortedBamIndex
    Array[File?]? sortedBamMD5            = PreprocessingSubworkflow.sortedBamMD5

    Array[File?]? fixedBam                = PreprocessingSubworkflow.fixedBam
    Array[File?]? fixedBamIndex           = PreprocessingSubworkflow.fixedBamIndex
    Array[File?]? fixedBamMD5             = PreprocessingSubworkflow.fixedBamMD5

    Array[File?]? summary                 = QCSubworkflow.summary
    Array[Array[File]?]? multMetrics      = QCSubworkflow.multMetrics
    Array[File?]? rawMetrics              = QCSubworkflow.rawMetrics
    Array[File?]? wgsMetrics              = QCSubworkflow.wgsMetrics
    Array[File?]? wgsNonZeroMetrics       = QCSubworkflow.wgsNonZeroMetrics
    Array[File?]? wgsNonZeroMetricsPDF    = QCSubworkflow.wgsNonZeroMetricsPDF
    Array[File?]? hsMetrics               = QCSubworkflow.hsMetrics
    Array[File?]? hsPerTargetMetrics      = QCSubworkflow.hsPerTargetMetrics
    Array[File?]? oxoGMetrics             = QCSubworkflow.oxoGMetrics

    Array[File?]? sequenceGroup           = CreateSequenceGroupingTSV.sequenceGrouping
    Array[File?]? sequenceGroupUnmapped   = CreateSequenceGroupingTSV.sequenceGroupingWithUnmapped

    Array[File?]? recalibrationReports    = BQSRSubworkflow.recalibrationReports
    Array[File?]? recalibratedBam         = BQSRSubworkflow.recalibratedBam
    Array[File?]? recalibratedBamIndex    = BQSRSubworkflow.recalibratedBamIndex
    Array[File?]? recalibratedBamChecksum = BQSRSubworkflow.recalibratedBamChecksum

    # PER-SAMPLE
#    Array[File?]? bamsPerSample           = VCSubworkflow.bamsPerSample
#    Array[File?]? baisPerSample           = VCSubworkflow.baisPerSample
#    Array[File?]? md5sPerSample           = VCSubworkflow.md5sPerSample
#    Array[File?]? dupMetricsPerSample     = VCSubworkflow.dupMetricsPerSample
#    Array[File?]? summaryVC               = VCSubworkflow.summary
#    Array[Array[File]?]? multMetricsVC    = VCSubworkflow.multMetrics
#    Array[File?]? mergedGVCFs             = VCSubworkflow.mergedGVCFs
#    Array[File?]? mergedGVCFsIndexes      = VCSubworkflow.mergedGVCFsIndexes

    # MULTI-SAMPLE
#    File? gatheredVCF                     = JGSubworkflow.gatheredVCF
#    File? gatheredVCFIndex                = JGSubworkflow.gatheredVCFIndex
#    File? recalIndelsFile                 = VQSRSubworkflow.recalIndelsFile
#    File? recalIndexIndelsFile            = VQSRSubworkflow.recalIndexIndelsFile
#    File? tranchesIndelsFile              = VQSRSubworkflow.tranchesIndelsFile
#    File? rScriptIndels                   = VQSRSubworkflow.rScriptIndels
#    File? plotsIndelsFile                 = VQSRSubworkflow.plotsIndelsFile
#    File? recalSnpsFile                   = VQSRSubworkflow.recalSnpsFile
#    File? recalIndexSnpsFile              = VQSRSubworkflow.recalIndexSnpsFile
#    File? tranchesSnpsFile                = VQSRSubworkflow.tranchesSnpsFile
#    File? modelSnpsReport                 = VQSRSubworkflow.modelSnpsReport
#    File? rScriptSnps                     = VQSRSubworkflow.rScriptSnps
#    File? plotsSnpsFile                   = VQSRSubworkflow.plotsSnpsFile
#    File? recalibratedVCF                 = VQSRSubworkflow.recalibratedVCF
#    File? recalibratedVCFIndex            = VQSRSubworkflow.recalibratedVCFIndex
#    File? recalibratedVCFMD5              = VQSRSubworkflow.recalibratedVCFMD5
#    File? VCMetrics                       = QCMultiSampleSubworkflow.VCMetrics
#    File? annovarVCF                      = AnnotationSubworkflow.annovarVCF
#    File? annovarTXT                      = AnnotationSubworkflow.annovarTXT
#    File? snpEffVCF                       = AnnotationSubworkflow.snpEffVCF
#    File? snpEffStats                     = AnnotationSubworkflow.snpEffStats
#    File? varAnnotatorVCF                 = AnnotationSubworkflow.varAnnotatorVCF
  }
}
