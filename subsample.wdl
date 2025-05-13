version 1.0

task stream_and_sample {
  input {
    String fastq_url
    String reference_path
    Float sampling_fraction
  }

  command {
    set -e

    # Install seqtk
    apt-get update && apt-get install -y apt-utils curl seqtk python3 pip bwa coreutils
    pip install pairtools
    curl -s ~{fastq_url} | seqtk sample -s100 - ~{sampling_fraction} | gzip > subsampled.fastq.gz

    curl -s ~{reference_path} | gzip -d > hs1.fa
    bwa index hs1.fa
    timeout 20m bwa mem -5SP -T0 -t16 hs1.fa subsampled.fastq.gz -o aligned.sam
    awk 'NR > 1 && length($10) != length($11) {print NR, $0}' aligned.sam
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hs1.genome aligned.sam > parsed.pairsam
    mkdir ebs
    mkdir ebs/temp
    pairtools sort --nproc 16 --tmpdir=./ebs/temp/ parsed.pairsam > sorted.pairsam
    pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
  }

  output {
    File stats = "stats.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "64G"
    disks: "local-disk 200 SSD"
    cpu: 16
  }
}

workflow sample_fastq {
  input {
    String fastq_url = "https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/5b73fa0e-658a-4248-b2b8-cd16155bc157--UCSC_GIAB_R1041_nanopore/HG002_R1041_PoreC/Dorado_v4/HG002_1_Dorado_v4_R1041_PoreC_400bps_sup.fastq.gz"
    String reference_path = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
    Float sampling_fraction = 0.2
  }

  call stream_and_sample {
    input:
      fastq_url = fastq_url,
      sampling_fraction = sampling_fraction,
      reference_path = reference_path
  }

  output {
    File stats = stream_and_sample.stats
  }
}
