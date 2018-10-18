
# ################# #
# Pipeline Workflow #
# ################# #

# All steps per-lane per-sample:
#
# 44. ANNOVAR
# 45. SNPEFF
# 46. GATK Variant Annotator
#

# IMPORTS
import "SubWorkflows/Utils.wdl" as utils

workflow AnnotationWF {

  File refFasta
  File refIndex
  File refDict

  File recalVCF

  String multiSampleName
  String multiSampleDir

  String annovarPath
  String humanDB

  String snpEffPath

  String gatkPath
  String javaOpts

  Int firstStep
  Int lastStep

  # Step 19 - Annovar, snpEff and GATK
  if ((firstStep <= 19) && (19 <= lastStep)) {

    # ########### #
    # 44. ANNOVAR #
    # ########### #

    call Annovar {
      input:
        outputBasename = multiSampleName + ".SNP_INDEL.recalibrated.ANNOVAR.hg38",
        annovarPath    = annovarPath,
        humanDB        = humanDB,
        recalVCF       = recalVCF
    }

    call utils.CopyResultsFilesToDir as copyAnnovarFiles {input: resultsDir = multiSampleDir, files = [Annovar.annovarVCF, Annovar.annovarTXT]}

    # ########## #
    # 45. SNPEFF #
    # ########## #

    call SnpEff {
      input:
        annovarVCF     = Annovar.annovarVCF,
        outputBasename = multiSampleName + ".SNP_INDEL.recalibrated.ANNOVAR.snpEff.GATK-style.hg38",
        snpEffDb       = "/data/d2/jlorsal/snpEff", #?????
        snpEffRef      = "hg38",
        snpEffConfig   = "/data/d2/jlorsal/snpEff/snpEff.config", #?????
        snpEffPath     = snpEffPath,
        javaOpts       = javaOpts,
    }

    call utils.CopyResultsFilesToDir as copySnpEffFiles {input: resultsDir = multiSampleDir, files = [SnpEff.snpEffVCF, SnpEff.snpEffStats]}

    # ########################## #
    # 46. GATK Variant Annotator #
    # ########################## #

    call VariantAnnotator {
      input:
        refFasta       = refFasta,
        refIndex       = refIndex,
        refDict        = refDict,
        annovarVCF     = Annovar.annovarVCF,
        snpEffVCF      = SnpEff.snpEffVCF,
        outputBasename = multiSampleName + ".SNP_INDEL.recalibrated.ANNOVAR.snpEff.GATK.annotations.hg38",
        gatkPath       = gatkPath,
        javaOpts       = javaOpts
    }

    call utils.CopyResultsFilesToDir as copyVariantAnnotatorFiles {input: resultsDir = multiSampleDir, files = VariantAnnotator.varAnnotatorVCF}

  } # End step 19

  output {
    File? annovarVCF      = Annovar.annovarVCF
    File? annovarTXT      = Annovar.annovarTXT

    File? snpEffVCF       = SnpEff.snpEffVCF
    File? snpEffStats     = SnpEff.snpEffStats

    File? varAnnotatorVCF = VariantAnnotator.varAnnotatorVCF
  }
}


task Annovar {

  String outputBasename

  String annovarPath

  File recalVCF

  String humanDB

  String annovarref = "hg38+" # ?????
  String refgene    = "refGene"
  String cytoband   = "cytoBand"
  String exac       = "exac03"
  String avsnp      = "avsnp150"
  String bravo      = "bravo_dbsnp_all"
  String dbnsfp     = "dbnsfp33a"
  String clinvar    = "clinvar_20180603"
  String gnomad     = "gnomad_exome"
  String intervar   = "intervar_20180118"

  command {
    ${annovarPath}/table_annovar.pl ${recalVCF} ${humanDB} \
      -buildver ${annovarref} \
      -out ${outputBasename} \
      -remove \
      -protocol ${refgene}, ${cytoband}, ${exac}, ${avsnp}, ${bravo}, ${dbnsfp}, ${clinvar}, ${gnomad}, ${intervar} \
      -operation g,r,f,f,f,f,f,f,f \
      -nastring . \
      -vcfinput
  }

  output {
    File annovarVCF = "${outputBasename}.vcf"
    File annovarTXT = "${outputBasename}.txt"
  }

  meta {
    taskDescription: "ANNOVAR is a rapid, efficient tool to annotate functional consequences of genetic variation from high-throughput sequencing data."
  }
}


task SnpEff {

  File annovarVCF

  String outputBasename

  String snpEffDb
  String snpEffRef
  String snpEffConfig

  String snpEffPath
  String javaOpts

  command {
    java ${javaOpts} -jar ${snpEffPath} eff ${snpEffRef} \
      -c ${snpEffConfig} \
      -v \
      -noLog \
      ${annovarVCF} \
      -s ${outputBasename}.stats \
      -noShiftHgvs \
      -o gatk > ${outputBasename}.vcf
  }

  output {
    File snpEffVCF   = "${outputBasename}.vcf"
    File snpEffStats = "${outputBasename}.stats"
  }

  meta {
    taskDescription: "Genomic variant annotations and functional effect prediction toolbox."
  }
}


task VariantAnnotator {

  File refFasta
  File refIndex
  File refDict

  File annovarVCF
  File snpEffVCF

  String outputBasename

  String gatkPath
  String javaOpts

  command {
    ${gatkPath} --java-opts ${javaOpts} VariantAnnotator \
      --reference ${refFasta} \
      --annotation SnpEff \ # ?????
      --variant ${annovarVCF} \
      --snpEffFile ${snpEffVCF} \
      --intervals ${snpEffVCF} \ #?????
      --output ${outputBasename}
  }

  output {
    File varAnnotatorVCF = "${outputBasename}.vcf"
  }

  meta {
    taskDescription: "Adding annotations to VCF files using GATK's Variant Annotator tool."
  }
}