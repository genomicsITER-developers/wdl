#
# WORKFLOW: Variant Calling
# STEPS:
# 1. MarkDuplicates --> MergeBamsPerSample
# 2. BQSR (por intervalos --> TSVIntervalsFile)
# 3. HaplotypeCaller (por intervalos --> ScatterIntervalList)
# 4. MergeVCFs
#

# IMPORTS
import "SubWorkflows/SetWorkingDirectory.wdl" as workDir
import "SubWorkflows/Preprocessing.wdl" as Preprocessing
import "SubWorkflows/BaseQualityScoreRecalibration.wdl" as BQSR


workflow VariantCallingPerSampleWF {

  String actualSampleDir
  String inputBams

  String dirName

  File? TSVIntervalsFile
  String protectionTag # Only for hg38

  # REFERENCE
  File refFasta
  File? refAlt
  File refIndex
  File refDict
  File dbSnpsVcf
  File dbSnpsVcfIdx
  Array[File] knownIndelsSitesVcfs
  Array[File] knownIndelsSitesIndices

  # Interval List
  File intervalList
  Int scatterCount
  Int breakBandsAtMultiplesOf
  String intervalListOutDir

  # GATK
  String gatkPath

  # PICARD
  String picardPath

  # JAVA OPTS
  String javaOpts

  call workDir.ConfigureResultsDirectory as SampleDirectory {input: resultsDir = actualSampleDir}

  call workDir.GetFilesFromString as GetBamsArray {input: inputBams = inputBams}

  call Preprocessing.MarkDuplicates as MergeBamsPerSample {
    input:
      bams = GetBamsArray.bams,
      gatkPath = gatkPath,
      outputBamBasename = dirName + ".ready.deduped",
      metricsFilename = dirName + ".ready.MarkDuplicates.metrics.txt",
      sortOrder = "coordinate",
      createIndex = "true",
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyMergedBams {input: resultsDir = SampleDirectory.directory, 
    files = [MergeBamsPerSample.markedBam, MergeBamsPerSample.markedBai, MergeBamsPerSample.markedMD5, MergeBamsPerSample.duplicateMetrics]}

  # Crear el fichero de intervalos en caso de que no esté definido en el fichero de inputs
  if (!defined(TSVIntervalsFile)) {
    call CreateSequenceGroupingTSV {input: refDict = refDict, protectionTag = protectionTag}

    call workDir.CopyResultsFilesToDir as CopyTSVIntervals {input: resultsDir = SampleDirectory.directory, 
      files = [CreateSequenceGroupingTSV.sequenceGrouping, CreateSequenceGroupingTSV.sequenceGroupingWithUnmapped]}
  }

  # Descomentar toda esta parte SOLO si se añade la opción de que se pueda ejecutar el BQSR pero no lo anterior del Variant Calling

  #call workDir.FindFilesInDir as findMergedBams {
  #  input:
  #    dir = actualSampleDir,
  #    pattern1 = actualSampleDir + "/*.ready.deduped.bam",
  #    pattern2 = actualSampleDir + "/*.ready.deduped.bai"
  #}

  #Array[String?] aux_Bams = if (firstStep < 9) then MergeBamsPerSample.markedBam else findMergedBams.outBams[0]
  #Array[String?] aux_Bais = if (firstStep < 9) then MergeBamsPerSample.markedBai else findMergedBams.outBams[1]
  #Array[Pair[String?, String?]] bams = zip(aux_Bams, aux_Bais)

  # Base Quality Score Recalibration

  Pair[String?, String?] bam_bai = (MergeBamsPerSample.markedBam, MergeBamsPerSample.markedBai)

  # BQSR (Sub workflow)
  call BQSR.BaseQualityScoreRecalibrationWF as BaseQuality {
    input:
      intervals = read_tsv(select_first([TSVIntervalsFile, CreateSequenceGroupingTSV.sequenceGrouping])),
      inputBam = bam_bai.left,
      inputBai = bam_bai.right,
      dbSnpsVcf = dbSnpsVcf,
      dbSnpsVcfIdx = dbSnpsVcfIdx,
      knownIndelsSitesVcfs = knownIndelsSitesVcfs,
      knownIndelsSitesIndices = knownIndelsSitesIndices,
      refFasta = refFasta,
      refDict = refDict,
      refIndex = refIndex,
      sampleName = dirName,
      gatkPath = gatkPath,
      resultsDir = actualSampleDir,
      javaOpts = javaOpts
  }
  # End BQSR

  # Divide el fichero de intervalos "wgs_calling_regions.hg38.interval_list" en sub-intervalos (crea "scatterCount" sub-ficheros)
  #call ScatterIntervalList {
  #  input:
  #    intervalList = intervalList,
  #    scatterCount = scatterCount,
  #    breakBandsAtMultiplesOf = breakBandsAtMultiplesOf,
  #    intervalListOutDir = intervalListOutDir,
  #    picardPath = picardPath,
  #    javaOpts = javaOpts
  #}

  #scatter (interval in ScatterIntervalList.outIntervals) {
  scatter (interval in read_tsv(select_first([TSVIntervalsFile, CreateSequenceGroupingTSV.sequenceGrouping]))) {
    call HaplotypeCaller {
      input:
        inputBam = BaseQuality.recalibratedBam,
        inputBai = BaseQuality.recalibratedBamIndex,
        intervalList = interval,
        outputBasename = dirName + ".ready.deduped.recalibrated",
        refFasta = refFasta,
        refDict = refDict,
        refIndex = refIndex,
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }
  }

  call MergeVCFs {
  	input:
  	  inputVCFs = HaplotypeCaller.firstGVCF,
  	  inputVCFsIndexes = HaplotypeCaller.firstGVCFIndex,
  	  outputBasename = dirName + ".ready.deduped.recalibrated.merged",
  	  gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as copyMergedGVCFs {input: resultsDir = actualSampleDir, 
    files = [MergeVCFs.mergedGVCFs, MergeVCFs.mergedGVCFsIndexes]}

  output {
    File bamsPerSample = MergeBamsPerSample.markedBam
    File baisPerSample = MergeBamsPerSample.markedBai
    File md5sPerSample = MergeBamsPerSample.markedMD5
    File dupMetricsPerSample = MergeBamsPerSample.duplicateMetrics

    File? sequenceGrouping = CreateSequenceGroupingTSV.sequenceGrouping
    File? sequenceGroupingWithUnmapped = CreateSequenceGroupingTSV.sequenceGroupingWithUnmapped

    File recalibrationReports = BaseQuality.recalibrationReports
    File recalibratedBam = BaseQuality.recalibratedBam
    File recalibratedBamIndex = BaseQuality.recalibratedBamIndex
    File recalibratedBamChecksum = BaseQuality.recalibratedBamChecksum

    File mergedGVCFs = MergeVCFs.mergedGVCFs
    File mergedGVCFsIndexes = MergeVCFs.mergedGVCFsIndexes
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
    File sequenceGrouping = "sequence_grouping.txt"
    File sequenceGroupingWithUnmapped = "sequence_grouping_with_unmapped.txt"
  }

  meta {
    taskDescription: "Generate sets of intervals for scatter-gathering over chromosomes."
  }
}


# Código Python3: Cambia los nombres de los ficheros para que no se sobreescriban al hacerle el glob y cuenta el número de ficheros
task ScatterIntervalList {

  File intervalList

  Int scatterCount
  Int breakBandsAtMultiplesOf

  String intervalListOutDir

  String picardPath

  command <<<
    set -e

    mkdir ${intervalListOutDir}

    java -jar ${picardPath} IntervalListTools \
      INPUT=${intervalList} \
      OUTPUT=${intervalListOutDir} \
      SCATTER_COUNT=${scatterCount} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=${breakBandsAtMultiplesOf}

    python <<CODE
    import glob, os
    intervals = sorted(glob.glob("${intervalListOutDir}/*/*.interval_list"))
    for i, interval in enumerate(intervals):
        (dir, filename) = os.path.split(interval)
        newName = os.path.join(dir, str(i+1) + filename)
        os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>

  output {
    Array[File] outIntervals = glob("${intervalListOutDir}/*/*.interval_list")
    Int intervalCount = read_int(stdout())
  }

  meta {
    taskDescription: "Create sets of intervals from interval_list file."
  }
}


task HaplotypeCaller {

  File inputBam
  File inputBai

  Array[String] intervalList

  String outputBasename

  File refFasta
  File refDict
  File refIndex

  String gatkPath

  Float? contamination

  String javaOpts

  # --annotation HaplotypeScore \

  command {
    ${gatkPath} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 ${javaOpts}" \
      HaplotypeCaller \
      -R ${refFasta} \
      -I ${inputBam} \
      -L ${sep=" -L " intervalList} \
      -O ${outputBasename}.g.vcf.gz \
      --interval-padding 500 \
      ${"-contamination " + contamination} \
      --max-alternate-alleles 6 \
      -ERC GVCF \
      --annotation MappingQualityRankSumTest \
      --annotation QualByDepth \
      --annotation ReadPosRankSumTest \
      --annotation RMSMappingQuality \
      --annotation FisherStrand \
      --annotation Coverage
  }

  output {
  	File firstGVCF = "${outputBasename}.g.vcf.gz"
  	File firstGVCFIndex = "${outputBasename}.g.vcf.gz.tbi"
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
  	File mergedGVCFs = "${outputBasename}.g.vcf.gz"
    File mergedGVCFsIndexes = "${outputBasename}.g.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Combines multiple variant files into a single variant file (per-sample)."
  }
}
