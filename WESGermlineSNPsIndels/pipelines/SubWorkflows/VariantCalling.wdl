
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 24. Merge and MarkDuplicates of BAM files from the same sample with GATK (per-sample)
# 25-28. Quality Control (per-sample)
# 29. Call variants with GATK-HaplotypeCaller (per-sample)
# 30. Combines multiple variant files into a single variant file (per-sample)
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils
import "SubWorkflows/Preprocessing.wdl" as Preprocessing
import "SubWorkflows/QualityControlPerSample.wdl" as QCPerSample

workflow VariantCallingWF {

  String actualSampleDir

  String inputBams

  String dirName

  # REFERENCE
  File refFasta
  File? refAlt
  File refIndex
  File refDict

  File bedFile
  File bedFilePadding

  File? TSVIntervalsFile

  Float? contamination
  Int maxAltAlleles
  Int intervalPadding

  Int firstStep
  Int lastStep
  Boolean wfQC_PerSample

  # GATK
  String gatkPath

  # PICARD
  #String picardPath

  # Qualimap
  String qualimapPath

  # JAVA OPTS
  String javaOpts
  String javaMemSize

  call utils.ConfigureResultsDirectory as SampleDirectory {input: resultsDir = actualSampleDir}

  call utils.GetFilesFromString as GetBamsArray {input: inputBams = inputBams}

  # Step 11 - Merge BAMs per sample and validate
  if ((firstStep <= 11) && (11 <= lastStep)) {

    # ##################################################################################### #
    # 24. Merge and MarkDuplicates of BAM files from the same sample with GATK (per-sample) #
    # ##################################################################################### #

    call Preprocessing.MarkDuplicates as MergeBamsPerSample {
      input:
        bams              = GetBamsArray.bams,
        gatkPath          = gatkPath,
        outputBamBasename = dirName + ".ready.deduped",
        metricsFilename   = dirName + ".ready.deduped.metrics.txt",
        sortOrder         = "coordinate",
        createIndex       = "true",
        javaOpts          = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyMergedBams {input: resultsDir = SampleDirectory.directory,
      files = [MergeBamsPerSample.markedBam, MergeBamsPerSample.markedBai, MergeBamsPerSample.markedMD5, MergeBamsPerSample.duplicateMetrics]}

  } # End step 11

  if ((firstStep <= 12) && (12 <= lastStep) && (wfQC_PerSample == true)) {
    call QCPerSample.QualityControlPerSampleWF as QCPerSampleSubworkflow {
      input:
        refFasta       = refFasta,
        bamFile        = MergeBamsPerSample.markedBam,
        bedFile        = bedFile,
        bedFilePadding = bedFilePadding,
        sampleName     = dirName,
        resultsDir     = SampleDirectory.directory,
        gatkPath       = gatkPath,
        javaOpts       = javaOpts,
        qualimapPath   = qualimapPath,
        javaMemSize    = javaMemSize
    }
  } # End QC and Step 12

  # Step 13 - Haplotype Caller and merge the results
  if ((firstStep <= 13) && (13 <= lastStep)) {

    Pair[String?, String?] bam_bai = (select_first([MergeBamsPerSample.markedBam, actualSampleDir + "/" + dirName + ".ready.deduped.bam"]), 
                                      select_first([MergeBamsPerSample.markedBai, actualSampleDir + "/" + dirName + ".ready.deduped.bai"]))

    # ######################################################## #
    # 30. Call variants with GATK-HaplotypeCaller (per-sample) #
    # ######################################################## #

    #scatter (interval in read_tsv(select_first([TSVIntervalsFile, CreateSequenceGroupingTSV.sequenceGrouping]))) {
    scatter (interval in read_tsv(TSVIntervalsFile)) {
      call HaplotypeCaller {
        input:
          refFasta        = refFasta,
          refDict         = refDict,
          refIndex        = refIndex,
          inputBam        = bam_bai.left,
          inputBai        = bam_bai.right,
          contamination   = contamination,
          maxAltAlleles   = maxAltAlleles,
          intervalPadding = intervalPadding,
          intervalList    = interval,
          outputBasename  = dirName + ".ready.deduped",
          gatkPath        = gatkPath,
          javaOpts        = javaOpts
      }
    }

    # ########################################################################### #
    # 31. Combines multiple variant files into a single variant file (per-sample) #
    # ########################################################################### #

    call MergeVCFs {
      input:
        inputVCFs        = HaplotypeCaller.gVCF,
        inputVCFsIndexes = HaplotypeCaller.gVCFIndex,
        outputBasename   = dirName + ".ready.deduped.merged",
        gatkPath         = gatkPath,
        javaOpts         = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyMergedGVCFs {input: resultsDir = actualSampleDir,
      files = [MergeVCFs.mergedGVCFs, MergeVCFs.mergedGVCFsIndexes]}

  } # End step 13

  output {
    File? bamsPerSample       = MergeBamsPerSample.markedBam
    File? baisPerSample       = MergeBamsPerSample.markedBai
    File? md5sPerSample       = MergeBamsPerSample.markedMD5
    File? dupMetricsPerSample = MergeBamsPerSample.duplicateMetrics

    File? summary             = QCPerSampleSubworkflow.summary
    Array[File]? multMetrics  = QCPerSampleSubworkflow.multMetrics

    File? mergedGVCFs         = MergeVCFs.mergedGVCFs
    File? mergedGVCFsIndexes  = MergeVCFs.mergedGVCFsIndexes
  }

}


# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {

  File refDict

  String protectionTag

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python2 <<CODE
    with open("${refDict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ${protectionTag}
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]

    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt", "w") as tsv_file:
        tsv_file.write(tsv_string)
        tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt", "w") as tsv_file_with_unmapped:
        tsv_file_with_unmapped.write(tsv_string)
        tsv_file_with_unmapped.close()
    CODE
  >>>

  output {
    File sequenceGrouping             = "sequence_grouping.txt"
    File sequenceGroupingWithUnmapped = "sequence_grouping_with_unmapped.txt"
  }

  meta {
    taskDescription: "Generate sets of intervals for scatter-gathering over chromosomes."
  }
}


task HaplotypeCaller {

  File refFasta
  File refDict
  File refIndex

  File inputBam
  File inputBai

  Float? contamination
  Int maxAltAlleles
  Int intervalPadding

  Array[String] intervalList

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 ${javaOpts}" HaplotypeCaller \
      -R ${refFasta} \
      -I ${inputBam} \
      -L ${sep=" -L " intervalList} \
      -O ${outputBasename}.g.vcf.gz \
      --interval-padding ${intervalPadding} \
      ${"-contamination " + contamination} \
      --max-alternate-alleles ${maxAltAlleles} \
      -ERC GVCF \
      --annotation MappingQualityRankSumTest \
      --annotation QualByDepth \
      --annotation ReadPosRankSumTest \
      --annotation RMSMappingQuality \
      --annotation FisherStrand \
      --annotation Coverage
  }

  output {
    File gVCF      = "${outputBasename}.g.vcf.gz"
    File gVCFIndex = "${outputBasename}.g.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Call germline SNPs and indels via local re-assembly of haplotypes (per-sample)."
  }
}


task MergeVCFs {

  Array[File] inputVCFs
  Array[File] inputVCFsIndexes

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" MergeVcfs \
      -I ${sep=' -I ' inputVCFs} \
      -O ${outputBasename}.g.vcf.gz
  }

  output {
    File mergedGVCFs        = "${outputBasename}.g.vcf.gz"
    File mergedGVCFsIndexes = "${outputBasename}.g.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Combines multiple variant files into a single variant file (per-sample)."
  }
}