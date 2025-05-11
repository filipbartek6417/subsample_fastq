version 1.0

task stream_and_sample {
  input {
    String fastq_url
    Float sampling_fraction
  }

  command {
    set -e

    # Install seqtk
    apt-get update && apt-get install -y apt-utils curl seqtk

    # Stream the FASTQ file and subsample using seqtk
    curl -s ~{fastq_url} | seqtk sample -s100 - ~{sampling_fraction} | gzip > subsampled.fastq.gz
  }

  output {
    File subsampled_fastq = "subsampled.fastq.gz"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "64G"
    disks: "local-disk 150 SSD"
    cpu: 16
  }
}

workflow sample_fastq {
  input {
    String fastq_url = "https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_PoreC/Dorado_v4/HG002_1_Dorado_v4_R1041_PoreC_400bps_sup.fastq.gz"
    Float sampling_fraction = 0.2
  }

  call stream_and_sample {
    input:
      fastq_url = fastq_url,
      sampling_fraction = sampling_fraction
  }

  output {
    File subsampled_fastq = stream_and_sample.subsampled_fastq
  }
}
