
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#  1. Convert FASTQ to uBAM and add Read Group information
#  2. Mark adapter sequences using MarkIlluminaAdapters
#  3. SamToFastq
#  4. BWA-MEM
#  5. MergeBamAlignment to generate a clean BAM
#  6. MarkDuplicates and Indexing of query-sorted merged BAM with PICARD
#  7. Sort merged BAM by coordinate order
#  8. Fix tag values for NM, MD, and UQ

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow PreprocessingWF {

  # CONFIGURATION PARAMETERS
  Int firstStep
  Int lastStep

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

  # READS
  String readGroupName
  String libraryName
  String platformUnit
  String platform
  String sequencingCenter

  String sampleName # pair[0]
  File readsFile1   # pair[1]
  File readsFile2   # pair[2]

  # GATK
  String gatkBaseCommand

  # BWA
  String bwaPath
  String bwaMemCommand

  # RESULTS DIR
  String resultsDir


  # PREPROCESSING WORKFLOW

  # Step 1 - Fastq to uBAM
  if ((firstStep <= 1) && (1 <= lastStep)) {

    # ############################################################################## #
    # 1. Convert FASTQ to uBAM and add Read Group information (per-lane per-sample)  #
    # ############################################################################## #

    call FastqsToUnmappedBam {
      input:
        readsFile1       = readsFile1,
        readsFile2       = readsFile1,
        sampleName       = sampleName,
        readGroupName    = readGroupName,
        libraryName      = libraryName,
        platformUnit     = platformUnit,
        platform         = platform,
        sequencingCenter = sequencingCenter,
        gatkBaseCommand  = gatkBaseCommand,
        outputBasename   = sampleName + ".unmapped"
    }

    call utils.CopyResultsFilesToDir as copyUnmappedBams {input: resultsDir = resultsDir, files = FastqsToUnmappedBam.unmappedBam}
  }

  # Step 2 - Mark Illumina Adapters
  if ((firstStep <= 2) && (2 <= lastStep)) {

    # ########################################################################### #
    # 2. Mark adapter sequences using MarkIlluminaAdapters (per-lane per-sample)  #
    # ########################################################################### #

    call MarkIlluminaAdapters {
      input:
        unmappedBam     = FastqsToUnmappedBam.unmappedBam,
        metricsFilename = sampleName + ".unmapped.markilluminaadapters.metrics",
        gatkBaseCommand = gatkBaseCommand,
        outputBasename  = sampleName + ".unmapped.markilluminaadapters"
    }

    call utils.CopyResultsFilesToDir as copyMarkedIlluminaAdapters  {input: resultsDir = resultsDir,
      files = [MarkIlluminaAdapters.unmappedTaggedBam, MarkIlluminaAdapters.metricsIlluminaAdapters]}

  } # End step 2


  # Step 3 - SamToFastq, BWA-MEM, and MergeBamAlignment to generate a clean BAM
  if ((firstStep <= 3) && (3 <= lastStep)) {

    # ############################################## #
    # 3. Convert uBAM to FASTQ (per-lane per-sample) #
    # ############################################## #

    call UnmappedBamToFastq {
      input:
        unmappedBam     = select_first([MarkIlluminaAdapters.unmappedTaggedBam, resultsDir + "/" + sampleName + ".unmapped.markilluminaadapters.bam"]),
        outputBasename  = sampleName,
        gatkBaseCommand = gatkBaseCommand
    }

    call utils.CopyResultsFilesToDir as copyConvertedFastq  {input: resultsDir = resultsDir, files = UnmappedBamToFastq.convertedFastq}

    # ################################ #
    # 4. BWA-MEM (per-lane per-sample) #
    # ################################ #

    call BwaMem {
      input:
        refFasta          = refFasta,
        refAlt            = refAlt,
        refIndex          = refIndex,
        refDict           = refDict,
        refAmb            = refAmb,
        refAnn            = refAnn,
        refBwt            = refBwt,
        refPac            = refPac,
        refSa             = refSa,
        readsFile         = UnmappedBamToFastq.convertedFastq,
        sampleName        = sampleName,
        bwaPath           = bwaPath,
        bwaMemCommand     = bwaMemCommand,
        outputBamBasename = sampleName + ".aligned"
    }

    call utils.CopyResultsFilesToDir as copyAlignedBams {input: resultsDir = resultsDir, files = BwaMem.alignedBam}

    # ################################################################## #
    # 5. MergeBamAlignment to generate a clean BAM (per-lane per-sample) #
    # ################################################################## #

    call MergeBamAlignment {
      input:
        refFasta          = refFasta,
        refIndex          = refIndex,
        refDict           = refDict,
        refAmb            = refAmb,
        refAnn            = refAnn,
        refBwt            = refBwt,
        refPac            = refPac,
        refSa             = refSa,
        sampleName        = sampleName,
        # Select unmappedBam from FastqsToUnmappedBam function out file or resultsDir file
        unmappedBam       = select_first([FastqsToUnmappedBam.unmappedBam, resultsDir + "/" + sampleName + ".unmapped.bam"]),
        alignedBam        = BwaMem.alignedBam,
        outputBamBasename = sampleName + ".aligned.merged",
        gatkBaseCommand   = gatkBaseCommand

    }

    call utils.CopyResultsFilesToDir as copyAlignedUnsortedBams {input: resultsDir = resultsDir, files = MergeBamAlignment.alignedUnsortedBam}

  } # End step 3


  # Step 4 - Mark duplicates
  if ((firstStep <= 4) && (4 <= lastStep)) {

    # ############################################################################### #
    # 6. MarkDuplicates and Indexing of query-sorted merged BAM (per-lane per-sample) #
    # ############################################################################### #

    call MarkDuplicates {
      input:
        bams              = select_first([MergeBamAlignment.alignedUnsortedBam, resultsDir + "/" + sampleName + ".aligned.merged.bam"]),
        outputBamBasename = sampleName + ".aligned.merged.deduped",
        metricsFilename   = sampleName + ".aligned.merged.deduped.MarkDuplicates.metrics.txt",
        sortOrder         = "queryname",
        gatkBaseCommand   = gatkBaseCommand
    }

    call utils.CopyResultsFilesToDir as copyMarkedFiles  {input: resultsDir = resultsDir, 
      files = [MarkDuplicates.markedBam, MarkDuplicates.markedMD5, MarkDuplicates.duplicateMetrics]}

  } # End step 4


  # Step 5 - Sort Sam
  if ((firstStep <= 5) && (5 <= lastStep)) {

    # ############################################################ #
    # 7. Sort merged BAM by coordinate order (per-lane per-sample) #
    # ############################################################ #

    call SortSam {
      input:
        bamToSort       = select_first([MarkDuplicates.markedBam, resultsDir + "/" + sampleName + ".aligned.merged.deduped.bam"]),
        outputBasename  = sampleName + ".aligned.merged.deduped.sorted",
        gatkBaseCommand = gatkBaseCommand

    }

    call utils.CopyResultsFilesToDir as copySortedFiles {input: resultsDir = resultsDir, files = [SortSam.sortedBam, SortSam.sortedBamIndex, SortSam.sortedBamMD5]}

    # ########################################################## #
    # 8. Fix tag values for NM, MD, and UQ (per-lane per-sample) #
    # ########################################################## #

    call FixBamTags {
      input:
        sortedBam      = SortSam.sortedBam,
        outputBasename = sampleName + ".aligned.merged.deduped.sorted.fixed",
        refFasta       = refFasta,
        gatkBaseCommand = gatkBaseCommand
    }

    call utils.CopyResultsFilesToDir as copySortedFixedFiles {input: resultsDir = resultsDir, 
      files = [FixBamTags.fixedBam, FixBamTags.fixedBamIndex, FixBamTags.fixedBamMD5]}

  } # End step 5


  output {
    File? unmappedBam             = FastqsToUnmappedBam.unmappedBam

    File? unmappedTaggedBam       = MarkIlluminaAdapters.unmappedTaggedBam
    File? metricsIlluminaAdapters = MarkIlluminaAdapters.metricsIlluminaAdapters

    File? convertedFastq          = UnmappedBamToFastq.convertedFastq
    File? alignedBam              = BwaMem.alignedBam
    File? alignedUnsortedBam      = MergeBamAlignment.alignedUnsortedBam

    File? markedBam               = MarkDuplicates.markedBam
    File? markedBai               = MarkDuplicates.markedBai
    File? markedMD5               = MarkDuplicates.markedMD5
    File? duplicateMetrics        = MarkDuplicates.duplicateMetrics

    File? sortedBam               = SortSam.sortedBam
    File? sortedBamIndex          = SortSam.sortedBamIndex
    File? sortedBamMD5            = SortSam.sortedBamMD5
    File? fixedBam                = FixBamTags.fixedBam
    File? fixedBamIndex           = FixBamTags.fixedBamIndex
    File? fixedBamMD5             = FixBamTags.fixedBamMD5
  }
}


task FastqsToUnmappedBam {

  File readsFile1
  File readsFile2

  String sampleName
  String readGroupName
  String libraryName
  String platformUnit
  String platform
  String sequencingCenter

  String gatkBaseCommand

  String outputBasename

  command {
    ${gatkBaseCommand} FastqToSam \
      --FASTQ ${readsFile1} \
      --FASTQ2 ${readsFile2} \
      --OUTPUT ${outputBasename}.bam \
      --READ_GROUP_NAME ${readGroupName} \
      --SAMPLE_NAME ${sampleName} \
      --LIBRARY_NAME ${libraryName} \
      --PLATFORM_UNIT ${platformUnit} \
      --PLATFORM ${platform} \
      --SEQUENCING_CENTER ${sequencingCenter}
  }

  output {
    File unmappedBam = "${outputBasename}.bam"
  }

  meta {
    taskDescription: "Convert a FASTQ file to an unmapped BAM (uBAM) and add Read Group Information (per-sample per-lane)."
  }
}


task MarkIlluminaAdapters {

  File unmappedBam

  String metricsFilename

  String gatkBaseCommand

  String outputBasename

  command {
    ${gatkBaseCommand} MarkIlluminaAdapters \
      --INPUT ${unmappedBam} \
      --METRICS ${metricsFilename} \
      --OUTPUT ${outputBasename}.bam \
      --ADAPTERS PAIRED_END
  }

  output {
    File unmappedTaggedBam       = "${outputBasename}.bam"
    File metricsIlluminaAdapters = "${metricsFilename}"
  }

  meta {
    taskDescription: "Mark adapter sequences using MarkIlluminaAdapters (per-sample per-lane)."
  }
}


task UnmappedBamToFastq {

  File unmappedBam

  String gatkBaseCommand

  String outputBasename

  command {
    ${gatkBaseCommand} SamToFastq \
      --INPUT ${unmappedBam} \
      --FASTQ ${outputBasename}.fastq \
      --CLIPPING_ATTRIBUTE XT \
      --CLIPPING_ACTION 2 \
      --INTERLEAVE true \
      --INCLUDE_NON_PF_READS true
  }

  output {
    File convertedFastq = "${outputBasename}.fastq"
  }

  meta {
    taskDescription: "Convert BAM to FASTQ and discount adapter sequences using SamToFastq (per-sample per-lane)."
  }
}


task BwaMem {

  File refFasta
  File? refAlt
  File refIndex
  File refDict
  File refAmb
  File refAnn
  File refBwt
  File refPac
  File refSa

  File readsFile
  String sampleName

  String bwaPath
  String bwaMemCommand

  String outputBamBasename

  command {

    set -eo pipefail

    bashRefFasta=${refFasta}
    bashReadFastq=${readsFile}

    # if refAlt exists and has data in it
    # if [ -s ${refAlt} ]; then

    ${bwaPath} ${bwaMemCommand} > ${outputBamBasename}.bam 2> >(tee ${outputBamBasename}.bwa.stderr.log >&2)

    # else refAlt is empty or not exists
    #else
    #  exit 1;
    #fi
  }

  output {
    File alignedBam   = "${outputBamBasename}.bam"
    File bwaStderrLog = "${outputBamBasename}.bwa.stderr.log"
  }

  meta {
    taskDescription: "Align reads and flag secondary hits using BWA-MEM (per-sample per-lane)."
  }
}


task MergeBamAlignment {

  File refFasta
  File refIndex
  File refDict
  File refAmb
  File refAnn
  File refBwt
  File refPac
  File refSa

  String sampleName

  File unmappedBam
  File alignedBam

  String gatkBaseCommand

  String outputBamBasename

  command {
    ${gatkBaseCommand} MergeBamAlignment \
      --ALIGNED_BAM ${alignedBam} \
      --UNMAPPED_BAM ${unmappedBam} \
      --OUTPUT ${outputBamBasename}.bam \
      --REFERENCE_SEQUENCE ${refFasta} \
      --SORT_ORDER "unsorted" \
      --ADD_MATE_CIGAR true \
      --CLIP_ADAPTERS false \
      --CLIP_OVERLAPPING_READS true \
      --INCLUDE_SECONDARY_ALIGNMENTS true \
      --MAX_INSERTIONS_OR_DELETIONS -1 \
      --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
      --ATTRIBUTES_TO_RETAIN XS
  }

  output {
    File alignedUnsortedBam = "${outputBamBasename}.bam"
  }

  meta {
    taskDescription: "Restore altered data and apply & adjust meta information using MergeBamAlignment (per-sample per-lane)."
  }
}


# Aggregate aligned+merged flowcell BAM files and mark duplicates
task MarkDuplicates {

  Array[File?] bams

  String gatkBaseCommand

  String outputBamBasename

  String metricsFilename

  String sortOrder # Preprocess -> queryname ; Variant Calling -> coordinate

  # The program default for READ_NAME_REGEX is appropriate in nearly every case.
  # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
  # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
  String? readNameRegex

  String? createIndex # Solo se añade esta opción en el MarkDuplicates del VariantCalling

  # IF --> Ya que no se pueden poner outputs opcionales... esto es un workaround para que cree un .bai vacío y que no de problemas por no encontrar el output markedBai
  command <<<
    ${gatkBaseCommand} MarkDuplicates \
      --INPUT ${sep=" -I " bams} \
      --OUTPUT ${outputBamBasename}.bam \
      ${"--CREATE_INDEX " + createIndex} \
      --CREATE_MD5_FILE true \
      --MAX_RECORDS_IN_RAM 500000 \
      --METRICS_FILE ${metricsFilename} \
      --REMOVE_DUPLICATES false \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --VALIDATION_STRINGENCY SILENT \
      --ASSUME_SORT_ORDER ${sortOrder} \
      ${"--READ_NAME_REGEX " + readNameRegex} \
      --CLEAR_DT false


    if [ ${default="false" createIndex} == "false" ]; then
      touch ${outputBamBasename}.bai
    fi
  >>>

  output {
    File markedBam        = "${outputBamBasename}.bam"
    File markedBai        = "${outputBamBasename}.bai"
    File markedMD5        = "${outputBamBasename}.bam.md5"
    File duplicateMetrics = "${metricsFilename}"
  }

  meta {
    taskDescription: "MarkDuplicates and Indexing of query-sorted merged BAM (per-sample per-lane).\nMerge and MarkDuplicates of BAM files from the same sample (per-sample)"
  }
}


task SortSam {

  File? bamToSort

  String outputBasename

  String gatkBaseCommand

  command {
    ${gatkBaseCommand} SortSam \
      --INPUT ${bamToSort} \
      --OUTPUT ${outputBasename}.bam \
      --SORT_ORDER coordinate \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }

  output {
    File sortedBam      = "${outputBasename}.bam"
    File sortedBamIndex = "${outputBasename}.bai"
    File sortedBamMD5   = "${outputBasename}.bam.md5"
  }

  meta {
    taskDescription: "Sort merged BAM by coordinate order (per-sample per-lane)."
  }
}


task FixBamTags {

  File? sortedBam

  File refFasta

  String outputBasename

  String gatkBaseCommand

  command {
    ${gatkBaseCommand} SetNmMdAndUqTags \
      --INPUT ${sortedBam} \
      --OUTPUT ${outputBasename}.bam \
      --REFERENCE_SEQUENCE ${refFasta} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  }

  output {
    File fixedBam      = "${outputBasename}.bam"
    File fixedBamIndex = "${outputBasename}.bai"
    File fixedBamMD5   = "${outputBasename}.bam.md5"
  }

  meta {
    taskDescription: "Fix tag values for NM, MD and UQ (per-sample per-lane)."
  }
}