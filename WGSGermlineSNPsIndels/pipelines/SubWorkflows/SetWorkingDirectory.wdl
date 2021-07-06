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
