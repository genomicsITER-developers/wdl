#

# ConfigureResultsDirectory
task ConfigureResultsDirectory {

  String resultsDir

  command {
    set -eo pipefail

    if [ ! -d ${resultsDir} ]; then
      mkdir -p ${resultsDir}
    fi
  }

  runtime {
    backend: "Local"
  }

  output {
    String directory = "${resultsDir}"
  }

  meta {
    taskDescription: "Create the directory where result files are going copied."
  }
}

# 
task CopyResultsFilesToDir {
  
  String resultsDir

  Array[File?] files

  # cp $file ${resultsDir}
  command {
    set -eo pipefail

    for file in ${sep=" " files}
    do
      filename=$(basename $file)
      ln -sf $file ${resultsDir}/$filename
    done
  }

  runtime {
    backend: "Local"
  }

  meta {
    taskDescription:  "Copy the output files from a previous function in the results directory."
  }
}

# 
task FindFilesInDir {

  String dir
  String pattern1
  String pattern2

  String dollar = "$"

  command <<<
    files1=(${dollar}(ls ${pattern1}))
    files2=(${dollar}(ls ${pattern2}))
    echo -e "${dollar}{files1[@]}\n${dollar}{files2[@]}" | tr [:blank:] \\t
  >>>

  runtime {
    backend: "Local"
  }

  output {
    Array[Array[String]] files = read_tsv(stdout())
  }

  meta {
    taskDescription: "Search files based in one or two patterns in a directory."
  }
}


task GetFilesFromString {

  String inputBams

  command <<<
    set -eo pipefail
    for b in ${inputBams}
    do
      echo $b
    done
  >>>

  runtime {
    backend: "Local"
  }

  output {
    Array[String] bams = read_lines(stdout())
  }

  meta {
    taskDescription: "Get full path of files."
  }
}


task WritePerSampleTSV {

  File fastqReadsTSV
  String resultsDir
  String bamSuffix

  Array[File?]? md5Bams # Solo sirve para que espere por el proceso anterior

  command <<<
    python2 <<CODE
    import sys
    import csv
    import os
    with open("${fastqReadsTSV}", "r") as read_pairs_file, open("read_samples.tsv", "w") as bam_samples_file:
        tsvreader = csv.reader(read_pairs_file, delimiter="\t")
        d = {}
        for r in tsvreader:
            sample_name = r[2] + "_" + r[3]
            sub_dir = r[5]
            bam_file = sample_name + "${bamSuffix}"
            bam_file_path = os.path.join("${resultsDir}", sub_dir, bam_file)
            if os.path.isfile(bam_file_path):
                if not sample_name in d:
                    d[sample_name] = []
                d[sample_name].append(bam_file_path)

        tsvwriter = csv.writer(bam_samples_file, delimiter="\t")
        for elem in d:
            tsvwriter.writerow([elem] + [' '.join(d[elem])])
    CODE
  >>>

  output {
    File perSampleArray = "read_samples.tsv"
  }

  meta {
    taskDescription: "Create a new TSV file with per-sample information and files."
  }
}


task WriteSampleNameMapTSV {
  
  File? perSampleTSV

  String resultsDir
  String vcfSuffix

  Array[File?]? gVcfsIndex # Solo para que espere por el proceso anterior

  command <<<
    python2 <<CODE
    import sys
    import csv
    import os
    with open("${perSampleTSV}", "r") as per_sample_file, open("sample_name_map.tsv", "w") as sample_name_map:
        tsvreader = csv.reader(per_sample_file, delimiter="\t")
        d={}
        for r in tsvreader:
            sample_name = r[0]
            vcf_file = r[0] + "${vcfSuffix}"
            vcf_file_path = os.path.join("${resultsDir}", sample_name, vcf_file)
            if os.path.isfile(vcf_file_path):
                d[sample_name] = vcf_file_path

        tsvwriter = csv.writer(sample_name_map, delimiter="\t")
        for elem in d:
            tsvwriter.writerow([elem] + [d[elem]])
    CODE
  >>>

  output {
    File sampleNameMap = "sample_name_map.tsv"
  }

  meta {
    taskDescription: "Create a Sample Name Map TSV file for multi-sample stage with the next structure:\n\tHG00096\tHG00096.g.vcf.gz\n\tNA19625\tNA19625.g.vcf.gz\n\tHG00268\tHG00268.g.vcf.gz"
  }
}