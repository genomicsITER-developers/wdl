
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 42. Collect Variant Calling Metrics with GATK
# 43. Evaluation of callset with GATK VariantEval
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow QualityControlMultiSampleWF {

  File refDict

  File recalVCF

  File dbSnp

  String multiSampleName
  String multiSampleDir

  Int firstStep
  Int lastStep

  String gatkPath
  String javaOpts
  String gatkBaseCommand = gatkPath + ' --java-options ' + '"' + javaOpts + '"' + ' '

  # Step 17 - Collect VC metrics
  if ((firstStep <= 17) && (17 <= lastStep)) {

    # ############################################# #
    # 42. Collect Variant Calling Metrics with GATK #
    # ############################################# #

    call CollectVariantCallingMetrics {
      input:
        refDict         = refDict,
        recalVCF        = recalVCF,
        dbSnp           = dbSnp,
        threads         = 1,
        outputBasename  = multiSampleName + ".SNP_INDEL.recalibrated.metrics",
        gatkBaseCommand = gatkBaseCommand
    }

    call utils.CopyResultsFilesToDir as copyVCMetrics {input: resultsDir = multiSampleDir, files = CollectVariantCallingMetrics.collectedMetrics}
  } # End step 17


  # Step 18 - Variant Evaluation
  #if ((firstStep <= 18) && (18 <= lastStep)) {

    # ############################################### #
    # 43. Evaluation of callset with GATK VariantEval #
    # ############################################### #

    #call VariantEval {input:...}
    #call utils.CopyResultsFilesToDir ...
  #} # End step 18

  output {
    File? VCMetrics = CollectVariantCallingMetrics.collectedMetrics
  }
}


task CollectVariantCallingMetrics {

  File refDict

  File recalVCF
  File dbSnp

  Int threads

  String outputBasename

  String gatkBaseCommand

  command {
    ${gatkBaseCommand} CollectVariantCallingMetrics \
      --INPUT ${recalVCF} \
      --DBSNP ${dbSnp} \
      --SEQUENCE_DICTIONARY ${refDict} \
      --THREAD_COUNT ${threads} \
      --GVCF_INPUT true \
      --OUTPUT ${outputBasename}
  }

  # Omitimos de momento la lista de intervalos
  # --TARGET_INTERVALS ${interval_list} \

  output {
    File collectedMetrics = "${outputBasename}"
  }

  meta {
    taskDescription: "Collects per-sample and aggregate (spanning all samples) metrics from the provided VCF file."
  }
}


#task VariantEvaluation {}