
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 38. Build a recalibration model to score variant quality for filtering purposes on INDELs
# 39. Build a recalibration model to score variant quality for filtering purposes on SNPs
# 40. Apply the desired level of recalibration
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow VariantQualityScoreRecalibrationWF {

  File refFasta
  File refIndex
  File refDict

  File dbSnps
  File dbSnpsIdx
  File hapMapResource
  File hapMapResourceIndex
  File omniResource
  File omniResourceIndex
  File oneThousandGenomesResource
  File oneThousandGenomesResourceIndex
  File millsResource
  File millsResourceIndex
  File axiomPolyResource
  File axiomPolyResourceIndex

  File gatheredVCF
  File gatheredVCFIndex

  Array[String] indelRecalibrationTrancheValues
  Array[String] indelRecalibrationAnnotationValues
  Array[String] snpRecalibrationTrancheValues
  Array[String] snpRecalibrationAnnotationValues
  Float indelFilterLevel
  Float snpFilterLevel

  Int indelsMaxGaussians
  Int snpsMaxGaussians

  String multiSampleName
  String multiSampleDir

  Int firstStep
  Int lastStep

  String gatkPath
  String javaOpts

  # Step 15 - VQSR for indels and snps
  if ((firstStep <= 15) && (15 <= lastStep)) {

    # ######################################################################################### #
    # 38. Build a recalibration model to score variant quality for filtering purposes on INDELs #
    # ######################################################################################### #

    call IndelsVariantRecalibrator {
      input:
        gatheredVCF               = gatheredVCF,
        gatheredVCFIndex          = gatheredVCFIndex,
        recalTrancheValues        = indelRecalibrationTrancheValues,
        recalAnnotationValues     = indelRecalibrationAnnotationValues,
        millsResource             = millsResource,
        millsResourceIndex        = millsResourceIndex,
        axiomPolyResource         = axiomPolyResource,
        axiomPolyResourceIndex    = axiomPolyResourceIndex,
        dbSnps                    = dbSnps,
        dbSnpsIdx                 = dbSnpsIdx,
        indelsMaxGaussians        = indelsMaxGaussians,
        outputBasename            = multiSampleName + ".recalibrate_indels",
        gatkPath                  = gatkPath,
        javaOpts                  = javaOpts
    }

    call utils.CopyResultsFilesToDir as CopyIndelsRecalFiles {
      input:
        resultsDir = multiSampleDir, 
        files = [IndelsVariantRecalibrator.recalibrationFile, 
                 IndelsVariantRecalibrator.recalibrationIndex, 
                 IndelsVariantRecalibrator.tranchesFile, 
                 IndelsVariantRecalibrator.rScriptFile, 
                 IndelsVariantRecalibrator.plotsFile]
    }

    # ####################################################################################### #
    # 39. Build a recalibration model to score variant quality for filtering purposes on SNPs #
    # ####################################################################################### #

    call SNPsVariantRecalibrator {
      input:
        gatheredVCF                        = gatheredVCF,
        gatheredVCFIndex                   = gatheredVCFIndex,
        recalTrancheValues                 = snpRecalibrationTrancheValues,
        recalAnnotationValues              = snpRecalibrationAnnotationValues,
        #downsampleFactor                   = snpVQSRDownsampleFactor,
        hapMapResource                     = hapMapResource,
        hapMapResourceIndex                = hapMapResourceIndex,
        omniResource                       = omniResource,
        omniResourceIndex                  = omniResourceIndex,
        oneThousandGenomesResource         = oneThousandGenomesResource,
        oneThousandGenomesResourceIndex    = oneThousandGenomesResourceIndex,
        dbSnps                             = dbSnps,
        dbSnpsIdx                          = dbSnpsIdx,
        snpsMaxGaussians                   = snpsMaxGaussians,
        outputBasename                     = multiSampleName + ".recalibrate_snps",
        gatkPath                           = gatkPath,
        javaOpts                           = javaOpts
    }

    call utils.CopyResultsFilesToDir as CopySnpsRecalCreateModelFiles {
      input:
        resultsDir = multiSampleDir, 
        files = [SNPsVariantRecalibrator.recalibrationFile, 
                 SNPsVariantRecalibrator.recalibrationIndex,
                 SNPsVariantRecalibrator.tranchesFile, 
                 SNPsVariantRecalibrator.modelReportFile, 
                 SNPsVariantRecalibrator.rScriptFile, 
                 SNPsVariantRecalibrator.plotsFile]
    }
  } # End step 15

  # Step 16 - Apply recalibration
  if ((firstStep <= 16) && (16 <= lastStep)) {

    # ############################################ #
    # 40. Apply the desired level of recalibration #
    # ############################################ #

    call ApplyRecalibration {
      input:
        refFasta         = refFasta,
        refIndex         = refIndex,
        refDict          = refDict,
        inputVCF         = gatheredVCF,
        inputVCFIndex    = gatheredVCFIndex,
        indelsRecal      = IndelsVariantRecalibrator.recalibrationFile,
        indelsRecalIndex = IndelsVariantRecalibrator.recalibrationIndex,
        indelsTranches   = IndelsVariantRecalibrator.tranchesFile,
        snpsRecal        = SNPsVariantRecalibrator.recalibrationFile,
        snpsRecalIndex   = SNPsVariantRecalibrator.recalibrationIndex,
        snpsTranches     = SNPsVariantRecalibrator.tranchesFile,
        indelFilterLevel = indelFilterLevel,
        snpFilterLevel   = snpFilterLevel,
        outputBasename   = multiSampleName + ".SNP_INDEL.recalibrated",
        gatkPath         = gatkPath,
        javaOpts         = javaOpts
    }

    call utils.CopyResultsFilesToDir as CopyRecalibratedVCFs {
      input: resultsDir = multiSampleDir, files = [ApplyRecalibration.recalibratedVCF, ApplyRecalibration.recalibratedVCFIndex, ApplyRecalibration.recalibratedVCFMD5]
    }

  } # End step 16

  output {
    File? recalIndelsFile      = IndelsVariantRecalibrator.recalibrationFile
    File? recalIndexIndelsFile = IndelsVariantRecalibrator.recalibrationIndex
    File? tranchesIndelsFile   = IndelsVariantRecalibrator.tranchesFile
    File? rScriptIndels        = IndelsVariantRecalibrator.rScriptFile
    File? plotsIndelsFile      = IndelsVariantRecalibrator.plotsFile

    File? recalSnpsFile        = SNPsVariantRecalibrator.recalibrationFile
    File? recalIndexSnpsFile   = SNPsVariantRecalibrator.recalibrationIndex
    File? tranchesSnpsFile     = SNPsVariantRecalibrator.tranchesFile
    File? modelSnpsReport      = SNPsVariantRecalibrator.modelReportFile
    File? rScriptSnps          = SNPsVariantRecalibrator.rScriptFile
    File? plotsSnpsFile        = SNPsVariantRecalibrator.plotsFile

    File? recalibratedVCF      = ApplyRecalibration.recalibratedVCF
    File? recalibratedVCFIndex = ApplyRecalibration.recalibratedVCFIndex
    File? recalibratedVCFMD5   = ApplyRecalibration.recalibratedVCFMD5
  }

}


task IndelsVariantRecalibrator {

  File gatheredVCF
  File gatheredVCFIndex

  Array[String] recalTrancheValues
  Array[String] recalAnnotationValues

  File millsResource
  File millsResourceIndex
  File axiomPolyResource
  File axiomPolyResourceIndex
  File dbSnps
  File dbSnpsIdx

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
      -resource mills,known=false,training=true,truth=true,prior=12:${millsResource} \
      -resource axiomPoly,known=false,training=true,truth=false,prior=10:${axiomPolyResource} \
      -resource dbsnp,known=true,training=false,truth=false,prior=2:${dbSnps} \
      -rscript-file ${outputBasename}_plots.R
  }

  output {
    File recalibrationFile  = "${outputBasename}.recal"
    File recalibrationIndex = "${outputBasename}.recal.idx"
    File tranchesFile       = "${outputBasename}.tranches"
    File rScriptFile        = "${outputBasename}_plots.R"
    File plotsFile          = "${outputBasename}_plots.R.pdf"
  }

  meta {
    taskDescription: "Build a recalibration model to score variant (indels) quality for filtering purposes (VQSR) (multi-sample)."
  }
}


task SNPsVariantRecalibrator {

  File gatheredVCF
  File gatheredVCFIndex

  Array[String] recalTrancheValues
  Array[String] recalAnnotationValues

  Int downsampleFactor

  File hapMapResource
  File hapMapResourceIndex
  File omniResource
  File omniResourceIndex
  File oneThousandGenomesResource
  File oneThousandGenomesResourceIndex
  File dbSnps
  File dbSnpsIdx

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
      -resource hapmap,known=false,training=true,truth=true,prior=15:${hapMapResource} \
      -resource omni,known=false,training=true,truth=true,prior=12:${omniResource} \
      -resource 1000G,known=false,training=true,truth=false,prior=10:${oneThousandGenomesResource} \
      -resource dbsnp,known=true,training=false,truth=false,prior=7:${dbSnps} \
      -rscript-file ${outputBasename}_plots.R
  }

  #-sample-every ${downsampleFactor} \

  output {
    File recalibrationFile  = "${outputBasename}.recal"
    File recalibrationIndex = "${outputBasename}.recal.idx"
    File tranchesFile       = "${outputBasename}.tranches"
    File modelReportFile    = "${outputBasename}.model.report"
    File rScriptFile        = "${outputBasename}_plots.R"
    File plotsFile          = "${outputBasename}_plots.R.pdf"
  }

  meta {
    taskDescription: "Build a recalibration model to score variant (SNPs) quality for filtering purposes (VQSR) (multi-sample)."
  }
}


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
    File recalibratedVCF      = "${outputBasename}.vcf.gz"
    File recalibratedVCFIndex = "${outputBasename}.vcf.gz.tbi"
    File recalibratedVCFMD5   = "${outputBasename}.vcf.gz.md5"
  }

  meta {
    taskDescription: "Apply a score cutoff to filter variants based on a recalibration table (multi-sample)."
  }
}