import "SubWorkflows/SetWorkingDirectory.wdl" as workDir
#import "SubWorkflows/VariantQualityScoreRecalibration.wdl" as VQSR

workflow JointGenotypingWF {

  File refFasta
  File refIndex
  File refDict

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

  String multiSampleName
  String multiSampleDir

  File sampleNameMap

  Array[Array[String]] intervalList

  Boolean useGenomicsDB

  # Levels, thresholds, ...
  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme 
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excessHetThreshold = 54.69

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

  Boolean wfVQSR
  Boolean wfQC_MultiSample

  String gatkPath

  String javaOpts

  # ONLY FOR GENOMICSDBIMPORT
  if (useGenomicsDB == true) {
    # Aplanar lista para poder llamar a GenomicsDBImport con un único intervalo, solo si queremos usar GenomicsDBImport
    #call flatArray {
    #  input:
    #    intervalLists = intervalList
    #}

    #scatter (interval in flatArray.flattenedIntervals) {
    scatter (interval in intervalList) {
      call ImportGVCFs {
        input:
          sampleNameMap = sampleNameMap,
          interval = interval,
          workspaceDirName = "genomicsdb",
          batchSize = 50,
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }

      call GenotypeGVCFsWithGenomicsDB {
        input:
          workspaceTar = ImportGVCFs.outputGenomicsDB,
          interval = interval,
          outputBasename = multiSampleName + ".cohort.genotyped",
          refFasta = refFasta,
          refIndex = refIndex,
          refDict = refDict,
          dbSnpsVcf = dbSnpsVcf,
          dbSnpsVcfIdx = dbSnpsVcfIdx,
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }

      call HardFilterAndMakeSitesOnlyVcf as HardFilterAndMakeSitesOnlyVcfWithGenomicsDB {
        input:
          vcf = GenotypeGVCFsWithGenomicsDB.outputVCF,
          vcfIndex = GenotypeGVCFsWithGenomicsDB.outputVCFIndex,
          excessHetThreshold = excessHetThreshold,
          outputBasename = multiSampleName + ".cohort.genotyped.filtered",
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }
    }
  }


  # ONLY FOR COMBINEGVCFS
  if (useGenomicsDB == false) {
    call getGVcfFiles as getGVcfs {input: sampleNameMap = read_tsv(sampleNameMap)}

    call getGVcfFiles as getGVcfsIndex {input: sampleNameMap = read_tsv(sampleNameMap), suffix = ".tbi"}

    scatter (interval in intervalList) {
      call CombineGVCFs {
        input:
          refFasta = refFasta,
          refDict = refDict,
          refIndex = refIndex,
          gVCFs = getGVcfs.gVCFs,
          gVCFsIndex = getGVcfsIndex.gVCFs,
          interval = interval,
          outputBasename = multiSampleName + ".cohort",
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }

      call GenotypeGVCFsWithCohortGVCF {
        input:
          gVCF = CombineGVCFs.outputCohort,
          gVCFIndex = CombineGVCFs.outputCohortIndex,
          interval = interval,
          outputBasename = multiSampleName + ".cohort.genotyped",
          refFasta = refFasta,
          refIndex = refIndex,
          refDict = refDict,
          dbSnpsVcf = dbSnpsVcf,
          dbSnpsVcfIdx = dbSnpsVcfIdx,
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }

      call HardFilterAndMakeSitesOnlyVcf as HardFilterAndMakeSitesOnlyVcfWithCohortGVCF {
        input:
          vcf = GenotypeGVCFsWithCohortGVCF.outputVCF,
          vcfIndex = GenotypeGVCFsWithCohortGVCF.outputVCFIndex,
          excessHetThreshold = excessHetThreshold,
          outputBasename = multiSampleName + ".cohort.genotyped.filtered",
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }
    }
  }

  call GatherVCFs as SitesOnlyGatherVcf {
    input:
      inputVCFs = write_lines(select_first([HardFilterAndMakeSitesOnlyVcfWithGenomicsDB.sitesOnlyVcf, HardFilterAndMakeSitesOnlyVcfWithCohortGVCF.sitesOnlyVcf])),
      outputBasename = multiSampleName + ".cohort.genotyped.filtered.sites_only.gathered",
      gatkPath = gatkPath,
      javaOpts = javaOpts
  }

  call workDir.CopyResultsFilesToDir as CopyGatheredVCF {input: resultsDir = multiSampleDir, files = [SitesOnlyGatherVcf.gatheredVCF, SitesOnlyGatherVcf.gatheredVCFIndex]}

  # VQSR
  if (wfVQSR == true) {

    call IndelsVariantRecalibrator {
      input:
        gatheredVCF = SitesOnlyGatherVcf.gatheredVCF,
        gatheredVCFIndex = SitesOnlyGatherVcf.gatheredVCFIndex,
        recalTrancheValues = indelRecalibrationTrancheValues,
        recalAnnotationValues = indelRecalibrationAnnotationValues,
        millsResourceVcf = millsResourceVcf,
        millsResourceVcfIndex = millsResourceVcfIndex,
        axiomPolyResourceVcf = axiomPolyResourceVcf,
        axiomPolyResourceVcfIndex = axiomPolyResourceVcfIndex,
        dbSnpsVcf = dbSnpsVcf,
        dbSnpsVcfIdx = dbSnpsVcfIdx,
        indelsMaxGaussians = indelsMaxGaussians,
        outputBasename = multiSampleName + ".recalibrate_indels",
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }

    call workDir.CopyResultsFilesToDir as CopyIndelsRecalFiles {
      input:
        resultsDir = multiSampleDir, 
        files = [IndelsVariantRecalibrator.recalibrationFile, 
                 IndelsVariantRecalibrator.recalibrationIndex, 
                 IndelsVariantRecalibrator.tranchesFile, 
                 IndelsVariantRecalibrator.rScriptFile, 
                 IndelsVariantRecalibrator.plotsFile]
    }

    call SNPsVariantRecalibratorCreateModel {
      input:
        gatheredVCF = SitesOnlyGatherVcf.gatheredVCF,
        gatheredVCFIndex = SitesOnlyGatherVcf.gatheredVCFIndex,
        recalTrancheValues = snpRecalibrationTrancheValues,
        recalAnnotationValues = snpRecalibrationAnnotationValues,
        downsampleFactor = snpVQSRDownsampleFactor,
        hapMapResourceVcf = hapMapResourceVcf,
        hapMapResourceVcfIndex = hapMapResourceVcfIndex,
        omniResourceVcf = omniResourceVcf,
        omniResourceVcfIndex = omniResourceVcfIndex,
        oneThousandGenomesResourceVcf = oneThousandGenomesResourceVcf,
        oneThousandGenomesResourceVcfIndex = oneThousandGenomesResourceVcfIndex,
        dbSnpsVcf = dbSnpsVcf,
        dbSnpsVcfIdx = dbSnpsVcfIdx,
        snpsMaxGaussians = snpsMaxGaussians,
        outputBasename = multiSampleName + ".recalibrate_snps",
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }

    call workDir.CopyResultsFilesToDir as CopySnpsRecalCreateModelFiles {
      input:
        resultsDir = multiSampleDir, 
        files = [SNPsVariantRecalibratorCreateModel.recalibrationFile, 
                 SNPsVariantRecalibratorCreateModel.recalibrationIndex,
                 SNPsVariantRecalibratorCreateModel.tranchesFile, 
                 SNPsVariantRecalibratorCreateModel.modelReportFile, 
                 SNPsVariantRecalibratorCreateModel.rScriptFile, 
                 SNPsVariantRecalibratorCreateModel.plotsFile]
    }

    #scatter ... -> ¿SNPsVariantRecalibrator?

    #Int numGVCFs = length(read_lines(WriteSampleNameMapTSV.sampleNameMap))
    #Boolean isSmallMultisample = numGVCFs <= 1000

    #scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variantFilteredVcf))) {
    #  call ApplyRecalibration {
    #    input:
    #      inputVCF = HardFilterAndMakeSitesOnlyVcf.variantFilteredVcf[idx],
    #      inputVCFIndex = HardFilterAndMakeSitesOnlyVcf.variantFilteredVcfIndex[idx],
    #      indelsRecal = IndelsVariantRecalibrator.recalibrationFile,
    #      indelsRecalIndex = IndelsVariantRecalibrator.recalibrationIndex,
    #      indelsTranches = IndelsVariantRecalibrator.tranchesFile,
    #      snpsRecal = SNPsVariantRecalibratorCreateModel.recalibrationFile,
    #      snpsRecalIndex = SNPsVariantRecalibratorCreateModel.recalibrationIndex,
    #      snpsTranches = SNPsVariantRecalibratorCreateModel.tranchesFile,
    #      indelFilterLevel = indelFilterLevel,
    #      snpFilterLevel = snpFilterLevel,
    #      outputBasename = multiSampleName + ".SNP_INDEL.recalibrated." + idx,
    #      gatkPath = gatkPath
    #  }
    #}

    # Aplicar la recalibración sin el scatter
    call ApplyRecalibration {
      input:
        refFasta = refFasta,
        refIndex = refIndex,
        refDict = refDict,
        inputVCF = SitesOnlyGatherVcf.gatheredVCF,
        inputVCFIndex = SitesOnlyGatherVcf.gatheredVCFIndex,
        indelsRecal = IndelsVariantRecalibrator.recalibrationFile,
        indelsRecalIndex = IndelsVariantRecalibrator.recalibrationIndex,
        indelsTranches = IndelsVariantRecalibrator.tranchesFile,
        snpsRecal = SNPsVariantRecalibratorCreateModel.recalibrationFile,
        snpsRecalIndex = SNPsVariantRecalibratorCreateModel.recalibrationIndex,
        snpsTranches = SNPsVariantRecalibratorCreateModel.tranchesFile,
        indelFilterLevel = indelFilterLevel,
        snpFilterLevel = snpFilterLevel,
        outputBasename = multiSampleName + ".SNP_INDEL.recalibrated",
        gatkPath = gatkPath,
        javaOpts = javaOpts
    }

    call workDir.CopyResultsFilesToDir as CopyRecalibratedVCFs {
      input: resultsDir = multiSampleDir, files = [ApplyRecalibration.recalibratedVCF, ApplyRecalibration.recalibratedVCFIndex, ApplyRecalibration.recalibratedVCFMD5]
    }

    #if (isSmallMultisample) {
    #  call GatherVCFs as FinalGatherVcfs {
    #    inputVCFs = 
    #  }
    #}

    if (wfQC_MultiSample == true) {
      call CollectVariantCallingMetrics {
        input:
          refDict = refDict,
          recalVcf = ApplyRecalibration.recalibratedVCF,
          recalVcfIndex = ApplyRecalibration.recalibratedVCFIndex,
          dbSnpsVcf = dbSnpsVcf,
          dbSnpsVcfIdx = dbSnpsVcfIdx,
          outputBasename = multiSampleName + ".SNP_INDEL.recalibrated",
          gatkPath = gatkPath,
          javaOpts = javaOpts
      }

      call workDir.CopyResultsFilesToDir as CopyVariantCallingMetrics {
      	input: resultsDir = multiSampleDir, files = [CollectVariantCallingMetrics.detailMetrics, CollectVariantCallingMetrics.summaryMetrics]
      }
    }

  }

  output {
    File? gatheredVCF = SitesOnlyGatherVcf.gatheredVCF
    File? gatheredVCFIndex = SitesOnlyGatherVcf.gatheredVCFIndex

    File? recalIndelsFile = IndelsVariantRecalibrator.recalibrationFile
    File? recalIndexIndelsFile = IndelsVariantRecalibrator.recalibrationIndex
    File? tranchesIndelsFile = IndelsVariantRecalibrator.tranchesFile
    File? rScriptIndels = IndelsVariantRecalibrator.rScriptFile
    File? plotsIndelsFile = IndelsVariantRecalibrator.plotsFile

    File? recalSnpsFile = SNPsVariantRecalibratorCreateModel.recalibrationFile
    File? recalIndexSnpsFile = SNPsVariantRecalibratorCreateModel.recalibrationIndex
    File? tranchesSnpsFile = SNPsVariantRecalibratorCreateModel.tranchesFile
    File? modelSnpsReport = SNPsVariantRecalibratorCreateModel.modelReportFile
    File? rScriptSnps = SNPsVariantRecalibratorCreateModel.rScriptFile
    File? plotsSnpsFile = SNPsVariantRecalibratorCreateModel.plotsFile

    File? recalibratedVCF = ApplyRecalibration.recalibratedVCF
  	File? recalibratedVCFIndex = ApplyRecalibration.recalibratedVCFIndex
  	File? recalibratedVCFMD5 = ApplyRecalibration.recalibratedVCFMD5

  	File? detailMetrics = CollectVariantCallingMetrics.detailMetrics
  	File? summaryMetrics = CollectVariantCallingMetrics.summaryMetrics
  }

}



task flatArray {
  Array[Array[String]] intervalLists

  command <<<
    python2 <<CODE
    l = [${sep=',' intervalLists}]
    flattenedList = [item for interval in l for item in interval]
    for elem in flattenedList:
        print(elem)
    CODE
  >>>

  output {
    Array[String] flattenedIntervals = read_lines(stdout())
  }

  meta {
    taskDescription: "Flatten a list of sub-lists."
  }
}



task ImportGVCFs {

  File sampleNameMap

  Array[String] interval

  String workspaceDirName
  Int batchSize

  String gatkPath

  String javaOpts

  # En GATK 4.0.2.1 la función GenomicsDBImport admite un ÚNICO intervalo 
  command <<<
    set -e

    rm -rf ${workspaceDirName}

    ${gatkPath} --java-options "${javaOpts}" GenomicsDBImport \
      --genomicsdb-workspace-path ${workspaceDirName} \
      --batch-size ${batchSize} \
      -L ${sep=" -L " interval} \
      --sample-name-map ${sampleNameMap} \
      --reader-threads 5 \
      -ip 500

    tar -cf ${workspaceDirName}.tar ${workspaceDirName}
  >>> 

  output {
    File outputGenomicsDB = "${workspaceDirName}.tar"
  }

  meta {
    taskDescription: "Import VCFs from multiple samples to GenomicsDB (multi-sample)."
  }
}



task getGVcfFiles {

  Array[Array[String]] sampleNameMap

  String? suffix

  command <<<
    python2 <<CODE
    suf = "${default='' suffix}"
    for sample in [${sep="," sampleNameMap}]:
        print(sample[1] + suf)
    CODE
  >>>

  output {
    Array[File] gVCFs = read_lines(stdout())
  }

  meta {
    taskDescription: "Get a list of gVCF files from SampleNameMap tsv."
  }
}



task CombineGVCFs {

  File refFasta
  File refDict
  File refIndex

  Array[File] gVCFs
  Array[File] gVCFsIndex

  Array[String] interval

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CombineGVCFs \
      -R ${refFasta} \
      --variant ${sep=" --variant " gVCFs} \
      -L ${sep=" -L " interval} \
      -O ${outputBasename}.g.vcf.gz \
      --create-output-variant-index true
  }

  output {
    File outputCohort = "${outputBasename}.g.vcf.gz"
    File outputCohortIndex = "${outputBasename}.g.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropiate annotations"
  }
}


task GenotypeGVCFsWithGenomicsDB {

  File workspaceTar

  Array[String] interval

  String outputBasename

  File refFasta
  File refDict
  File refIndex
  File dbSnpsVcf
  File dbSnpsVcfIdx

  String gatkPath

  String javaOpts

  command <<<
    set -e

    tar -xf ${workspaceTar}
    WORKSPACE=$( basename ${workspaceTar} .tar)

    ${gatkPath} --java-options "${javaOpts}" GenotypeGVCFs \
      -R ${refFasta} \
      -O ${outputBasename}.vcf.gz \
      -D ${dbSnpsVcf} \
      -G StandardAnnotation \
      --only-output-calls-starting-in-intervals true \
      -new-qual true \
      -V gendb://$WORKSPACE \
      -L ${sep=" -L " interval}
  >>>

  output {
    File outputVCF = "${outputBasename}.vcf.gz"
    File outputVCFIndex = "${outputBasename}.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Perform joint genotyping after GenomicsDBImport on one or more samples pre-called with HaplotypeCaller (multi-sample)."
  }
}


task GenotypeGVCFsWithCohortGVCF {

  File gVCF
  File gVCFIndex

  Array[String] interval

  String outputBasename

  File refFasta
  File refDict
  File refIndex
  File dbSnpsVcf
  File dbSnpsVcfIdx

  String gatkPath

  String javaOpts

  command {

    ${gatkPath} --java-options "${javaOpts}" GenotypeGVCFs \
      -R ${refFasta} \
      -O ${outputBasename}.vcf.gz \
      -D ${dbSnpsVcf} \
      -G StandardAnnotation \
      --only-output-calls-starting-in-intervals true \
      -new-qual true \
      -V ${gVCF} \
      -L ${sep=" -L " interval}
  }

  output {
    File outputVCF = "${outputBasename}.vcf.gz"
    File outputVCFIndex = "${outputBasename}.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Perform joint genotyping after CombineGVCFs on one or more samples pre-called with HaplotypeCaller (multi-sample)."
  }
}


task HardFilterAndMakeSitesOnlyVcf {

  File vcf
  File vcfIndex

  Float excessHetThreshold

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    set -e

    ${gatkPath} --java-options "${javaOpts}" VariantFiltration \
      --filter-expression "ExcessHet > ${excessHetThreshold}" \
      --filter-name ExcessHet \
      -O ${outputBasename}.vcf.gz \
      -V ${vcf}

    ${gatkPath} --java-options "-Xmx3g -Xms3g" MakeSitesOnlyVcf \
      -I ${outputBasename}.vcf.gz \
      -O ${outputBasename}.sites_only.vcf.gz
  }

  output {
    File variantFilteredVcf = "${outputBasename}.vcf.gz"
    File variantFilteredVcfIndex = "${outputBasename}.vcf.gz.tbi"
    File sitesOnlyVcf = "${outputBasename}.sites_only.vcf.gz"
    File sitesOnlyVcfIndex = "${outputBasename}.sites_only.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Filter variant calls based on INFO and/or FORMAT annotations and remove all genotype information from it while retaining all site level information (multi-sample)."
  }
}


task GatherVCFs {

  File inputVCFs

  String outputBasename

  String gatkPath

  String javaOpts

  command <<<
    set -e

    mv ${inputVCFs} inputs.list

    ${gatkPath} --java-options "${javaOpts}" GatherVcfsCloud \
      --ignore-safety-checks true \
      --gather-type BLOCK \
      -I inputs.list \
      -O ${outputBasename}.vcf.gz

    tabix ${outputBasename}.vcf.gz
  >>>

  output {
    File gatheredVCF = "${outputBasename}.vcf.gz"
    File gatheredVCFIndex = "${outputBasename}.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Gathers multiple VCF files from a scatter operation into a single VCF file (multi-sample)."
  }
}


task IndelsVariantRecalibrator {

  File gatheredVCF
  File gatheredVCFIndex

  Array[String] recalTrancheValues
  Array[String] recalAnnotationValues

  File millsResourceVcf
  File millsResourceVcfIndex
  File axiomPolyResourceVcf
  File axiomPolyResourceVcfIndex
  File dbSnpsVcf
  File dbSnpsVcfIdx

  Int indelsMaxGaussians

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" VariantRecalibrator \
      -V ${gatheredVCF} \
      -O ${outputBasename}.recal \
      --tranches-file ${outputBasename}.tranches \
      --trust-all-polymorphic true \
      -tranche ${sep=" -tranche " recalTrancheValues} \
      -an ${sep=" -an " recalAnnotationValues} \
      -mode INDEL \
      --max-gaussians ${indelsMaxGaussians} \
      -resource mills,known=false,training=true,truth=true,prior=12:${millsResourceVcf} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPolyResourceVcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbSnpsVcf} \
      -rscript-file ${outputBasename}_plots.R
  }

  output {
    File recalibrationFile = "${outputBasename}.recal"
    File recalibrationIndex = "${outputBasename}.recal.idx"
    File tranchesFile = "${outputBasename}.tranches"
    File rScriptFile = "${outputBasename}_plots.R"
    File plotsFile = "${outputBasename}_plots.R.pdf"
  }

  meta {
    taskDescription: "Build a recalibration model to score variant (indels) quality for filtering purposes (VQSR) (multi-sample)."
  }
}


task SNPsVariantRecalibratorCreateModel {

  File gatheredVCF
  File gatheredVCFIndex

  Array[String] recalTrancheValues
  Array[String] recalAnnotationValues

  Int downsampleFactor

  File hapMapResourceVcf
  File hapMapResourceVcfIndex
  File omniResourceVcf
  File omniResourceVcfIndex
  File oneThousandGenomesResourceVcf
  File oneThousandGenomesResourceVcfIndex
  File dbSnpsVcf
  File dbSnpsVcfIdx

  Int snpsMaxGaussians

  String outputBasename

  String gatkPath

  String javaOpts

  # En el wdl de gatk tienen -Xmx100g y -Xms100g
  command {
    ${gatkPath} --java-options "${javaOpts}" VariantRecalibrator \
      -V ${gatheredVCF} \
      -O ${outputBasename}.recal \
      --tranches-file ${outputBasename}.tranches \
      --trust-all-polymorphic true \
      -tranche ${sep=" -tranche " recalTrancheValues} \
      -an ${sep=" -an " recalAnnotationValues} \
      -mode SNP \
      --output-model ${outputBasename}.model.report \
      --max-gaussians ${snpsMaxGaussians} \
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapMapResourceVcf} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omniResourceVcf} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${oneThousandGenomesResourceVcf} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbSnpsVcf} \
      -rscript-file ${outputBasename}_plots.R
  }

  #-sample-every ${downsampleFactor} \

  output {
    File recalibrationFile = "${outputBasename}.recal"
    File recalibrationIndex = "${outputBasename}.recal.idx"
    File tranchesFile = "${outputBasename}.tranches"
    File modelReportFile =  "${outputBasename}.model.report"
    File rScriptFile = "${outputBasename}_plots.R"
    File plotsFile = "${outputBasename}_plots.R.pdf"
  }

  meta {
    taskDescription: "Build a recalibration model to score variant (SNPs) quality for filtering purposes (VQSR) (multi-sample)."
  }
}


# ............


task ApplyRecalibration {

  File refFasta
  File refIndex
  File refDict

  File inputVCF
  File inputVCFIndex

  File indelsRecal
  File indelsRecalIndex
  File indelsTranches

  File snpsRecal
  File snpsRecalIndex
  File snpsTranches

  Float indelFilterLevel
  Float snpFilterLevel

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    set -e

    ${gatkPath} --java-options "${javaOpts}" ApplyVQSR \
      -R ${refFasta} \
      -O tmp.INDEL.recalibrated.vcf \
      -V ${inputVCF} \
      --recal-file ${indelsRecal} \
      --tranches-file ${indelsTranches} \
      --truth-sensitivity-filter-level ${indelFilterLevel} \
      --create-output-variant-index true \
      --create-output-variant-md5 true \
      -mode INDEL

    ${gatkPath} --java-options "${javaOpts}" ApplyVQSR \
      -R ${refFasta} \
      -O  ${outputBasename}.vcf.gz \
      -V tmp.INDEL.recalibrated.vcf \
      --recal-file ${snpsRecal} \
      --tranches-file ${snpsTranches} \
      --truth-sensitivity-filter-level ${snpFilterLevel} \
      --create-output-variant-index true \
      --create-output-variant-md5 true \
      -mode SNP
  }

  output {
    File recalibratedVCF = "${outputBasename}.vcf.gz"
    File recalibratedVCFIndex = "${outputBasename}.vcf.gz.tbi"
    File recalibratedVCFMD5 = "${outputBasename}.vcf.gz.md5"
  }

  meta {
    taskDescription: "Apply a score cutoff to filter variants based on a recalibration table (multi-sample)."
  }
}


task CollectVariantCallingMetrics {

  File refDict

  File recalVcf
  File recalVcfIndex

  File dbSnpsVcf
  File dbSnpsVcfIdx

  File? intervalList

  String outputBasename

  String gatkPath

  String javaOpts

  command {
    ${gatkPath} --java-options "${javaOpts}" CollectVariantCallingMetrics \
      -I ${recalVcf} \
      --DBSNP ${dbSnpsVcf} \
      --SEQUENCE_DICTIONARY ${refDict} \
      -O ${outputBasename} \
      --THREAD_COUNT 8
  }
  #--TARGET_INTERVALS ${intervalList}

  output {
    File detailMetrics = "${outputBasename}.variant_calling_detail_metrics"
    File summaryMetrics = "${outputBasename}.variant_calling_summary_metrics"
  }

  meta {
    taskDescription: "Collects per-sample and aggregate (spanning all samples) metrics from provided VCF file (multi-sample)."
  }
}