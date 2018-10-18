
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 32. Import VCFs into GenomicsDB before joint genotyping
# 33. Perform joint genotyping on gVCF files produced by HaplotypeCaller using GenomicsDB
# 34. Filter variant calls and removes all genotype information from a VCF/VCF.gz/BCF file
# 35. Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
# 36. Perform joint genotyping on gVCF files produced by HaplotypeCaller
# 37. Filter variant calls and removes all genotype information from a VCF/VCF.gz/BCF file
# 38. Merge all VCFs from the previous scatter
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow JointGenotypingWF {

  File refFasta
  File refIndex
  File refDict

  File dbSnpsVcf
  File dbSnpsVcfIdx

  String multiSampleName
  String multiSampleDir

  File sampleNameMap

  Array[Array[String]] intervalList

  Boolean useGenomicsDB

  # Levels, thresholds, ...
  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  Float excessHetThreshold = 54.69

  Int firstStep
  Int lastStep

  String gatkPath

  String javaOpts

  # Step 14 - Joint Genotyping with genomicsDBImport or CombineGVCFs
  if ((firstStep <= 14) && (14 <= lastStep)) {

    # ONLY FOR GENOMICSDBIMPORT
    if (useGenomicsDB == true) {

      scatter (interval in intervalList) { # CHANGE THIS INTERVALLIST FILE

        # ####################################################### #
        # 32. Import VCFs into GenomicsDB before joint genotyping #
        # ####################################################### #

        call ImportGVCFs {
          input:
            sampleNameMap    = sampleNameMap,
            interval         = interval,
            workspaceDirName = "genomicsdb",
            batchSize        = 50,
            gatkPath         = gatkPath,
            javaOpts         = javaOpts
        }

        # ####################################################################################### #
        # 33. Perform joint genotyping on gVCF files produced by HaplotypeCaller using GenomicsDB #
        # ####################################################################################### #

        call GenotypeGVCFsWithGenomicsDB {
          input:
            workspaceTar   = ImportGVCFs.outputGenomicsDB,
            interval       = interval,
            outputBasename = multiSampleName + ".cohort.genotyped",
            refFasta       = refFasta,
            refIndex       = refIndex,
            refDict        = refDict,
            dbSnpsVcf      = dbSnpsVcf,
            dbSnpsVcfIdx   = dbSnpsVcfIdx,
            gatkPath       = gatkPath,
            javaOpts       = javaOpts
        }

        # ######################################################################################## #
        # 34. Filter variant calls and removes all genotype information from a VCF/VCF.gz/BCF file #
        # ######################################################################################## #

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

      scatter (interval in intervalList) { # CHANGE THIS INTERVALLIST FILE

        # ########################################################################################### #
        # 35. Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file #
        # ########################################################################################### #

        call CombineGVCFs {
          input:
            refFasta       = refFasta,
            refDict        = refDict,
            refIndex       = refIndex,
            gVCFs          = getGVcfs.gVCFs,
            gVCFsIndex     = getGVcfsIndex.gVCFs,
            interval       = interval,
            outputBasename = multiSampleName + ".cohort",
            gatkPath       = gatkPath,
            javaOpts       = javaOpts
        }

        # ###################################################################### #
        # 36. Perform joint genotyping on gVCF files produced by HaplotypeCaller #
        # ###################################################################### #

        call GenotypeGVCFsWithCohortGVCF {
          input:
            gVCF           = CombineGVCFs.outputCohort,
            gVCFIndex      = CombineGVCFs.outputCohortIndex,
            interval       = interval,
            outputBasename = multiSampleName + ".cohort.genotyped",
            refFasta       = refFasta,
            refIndex       = refIndex,
            refDict        = refDict,
            dbSnpsVcf      = dbSnpsVcf,
            dbSnpsVcfIdx   = dbSnpsVcfIdx,
            gatkPath       = gatkPath,
            javaOpts       = javaOpts
        }

        # ######################################################################################## #
        # 37. Filter variant calls and removes all genotype information from a VCF/VCF.gz/BCF file #
        # ######################################################################################## #

        call HardFilterAndMakeSitesOnlyVcf as HardFilterAndMakeSitesOnlyVcfWithCohortGVCF {
          input:
            vcf                = GenotypeGVCFsWithCohortGVCF.outputVCF,
            vcfIndex           = GenotypeGVCFsWithCohortGVCF.outputVCFIndex,
            excessHetThreshold = excessHetThreshold,
            outputBasename     = multiSampleName + ".cohort.genotyped.filtered",
            gatkPath           = gatkPath,
            javaOpts           = javaOpts
        }
      }
    }

    # ############################################ #
    # 38. Merge all VCFs from the previous scatter #
    # ############################################ #

    call GatherVCFs as SitesOnlyGatherVcf {
      input:
        #inputVCFs = write_lines(select_first([HardFilterAndMakeSitesOnlyVcfWithGenomicsDB.sitesOnlyVcf, HardFilterAndMakeSitesOnlyVcfWithCohortGVCF.sitesOnlyVcf])),
        inputVCFs      = write_lines(HardFilterAndMakeSitesOnlyVcfWithCohortGVCF.sitesOnlyVcf),
        outputBasename = multiSampleName + ".cohort.genotyped.filtered.sites_only.gathered",
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as CopyGatheredVCF {input: resultsDir = multiSampleDir, files = [SitesOnlyGatherVcf.gatheredVCF, SitesOnlyGatherVcf.gatheredVCFIndex]}

  } # End step 14

  output {
    File? gatheredVCF      = SitesOnlyGatherVcf.gatheredVCF
    File? gatheredVCFIndex = SitesOnlyGatherVcf.gatheredVCFIndex
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


# -----------------------------------------------------
# - Tasks for Joint Genotyping using GenomicsDBImport -
# -----------------------------------------------------

task ImportGVCFs {

  File sampleNameMap

  Array[String] interval

  String workspaceDirName
  Int batchSize

  String gatkPath

  String javaOpts

  # In GATK 4.0.2.1 version, GenomicsDBImport allows ONLY one interval
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



# -------------------------------------------------
# - Tasks for Joint Genotyping using CombineGVCFs -
# -------------------------------------------------

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
    File outputCohort      = "${outputBasename}.g.vcf.gz"
    File outputCohortIndex = "${outputBasename}.g.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Merges one or more HaplotypeCaller GVCF files into a single GVCF with appropiate annotations"
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
      -new-qual true \
      -V ${gVCF} \
      -L ${sep=" -L " interval}
  }

  # --only-output-calls-starting-in-intervals true \

  output {
    File outputVCF      = "${outputBasename}.vcf.gz"
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
    File variantFilteredVcf      = "${outputBasename}.vcf.gz"
    File variantFilteredVcfIndex = "${outputBasename}.vcf.gz.tbi"
    File sitesOnlyVcf            = "${outputBasename}.sites_only.vcf.gz"
    File sitesOnlyVcfIndex       = "${outputBasename}.sites_only.vcf.gz.tbi"
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

  # --gather-type:
  #    * BLOCK: Faster but onlyworks when you have both bgzipped inputs and outputs.
  #    * CONVENTIONAL: Much slower but should work on all vcf files.
  #    * AUTOMATIC: Chooses BLOCK if possible and CONVENTIONAL otherwise.

  output {
    File gatheredVCF      = "${outputBasename}.vcf.gz"
    File gatheredVCFIndex = "${outputBasename}.vcf.gz.tbi"
  }

  meta {
    taskDescription: "Gathers multiple VCF files from a scatter operation into a single VCF file (multi-sample)."
  }
}