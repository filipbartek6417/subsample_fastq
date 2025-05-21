version 1.0

task stream_and_sample {
  input {
    File fasta
    Float sampling_fraction
  }

  command {
    set -e
    apt-get update && apt-get install -y apt-utils seqtk
    seqtk sample -s100 ~{fasta} ~{sampling_fraction} | gzip > subsampled.fastq.gz
  }

  output {
    File subsampled_fastq = "subsampled.fastq.gz"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "64G"
    disks: "local-disk 2000 SSD"
    cpu: 32
  }
}

workflow sample_fastq {
  input {
    String fasta = "gs://fc-c3eed389-0be2-4bbc-8c32-1a40b8696969/submissions/bd394ffa-b4db-4031-81d6-8bf319b60390/porec_qc/3aa79ef5-c067-4d32-a29a-32fa5c9c378e/call-minimap2_align/porec.fa.gz"
    Float sampling_fraction = 0.2
  }

  call stream_and_sample {
    input:
      fasta = fasta,
      sampling_fraction = sampling_fraction
  }

  output {
    File subsampled_fastq = stream_and_sample.subsampled_fastq
  }
}
