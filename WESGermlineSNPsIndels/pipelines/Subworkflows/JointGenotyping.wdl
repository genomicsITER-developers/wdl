
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 31. 
# 32. 
# 33. 
# 34. Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
# 35. Perform joint genotyping on gVCF files produced by HaplotypeCaller
# 36. 
#

# IMPORTS
import "SubWorkflows/SetWorkingDirectory.wdl" as workDir

workflow JointGenotypingWF {

  File refFasta
  File refIndex
  File refDict

  Array[Array[String]] intervalList

  Boolean useGenomicsDB

  String gatkPath

  String javaOpts

  # Step 14 - Joint Genotyping with genomicsDBImport or CombineGVCFs
  if ((firstStep <= 14) && (14 <= lastStep)) {

    # ONLY FOR GENOMICSDBIMPORT
    if (useGenomicsDB == true) {

      scatter (interval in intervalList) {

        # ########################################################################################### #
        # 31.  #
        # ########################################################################################### #

        call ImportGVCFs {}

        # ########################################################################################### #
        # 32.  #
        # ########################################################################################### #

        call GenotypeGVCFsWithGenomicsDB {}

        # ########################################################################################### #
        # 33.  #
        # ########################################################################################### #

        call HardFilterAndMakeSitesOnlyVcf as HardFilterAndMakeSitesOnlyVcfWithGenomicsDB {}
      }
    }

    # ONLY FOR COMBINEGVCFS
    if (useGenomicsDB == false) {
      scatter (interval in intervalList) {

        # ########################################################################################### #
        # 34. Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file #
        # ########################################################################################### #

        call CombineGVCFs {}

        # ###################################################################### #
        # 35. Perform joint genotyping on gVCF files produced by HaplotypeCaller #
        # ###################################################################### #

        call GenotypeGVCFsWithCohortGVCF {}


        # ########################################################################################### #
        # 36.  #
        # ########################################################################################### #

        call HardFilterAndMakeSitesOnlyVcf as HardFilterAndMakeSitesOnlyVcfWithCohortGVCF {}
      }
    }
  } # End step 14

}

task ... {

}