version 1.0

task stream_and_sample {
  input {
    String r1
    String r2
    String reference_path
    Float sampling_fraction
  }

  command {
    set -e

    # Install seqtk
    apt-get update && apt-get install -y apt-utils curl seqtk python3 pip bwa coreutils
    pip install pairtools

    curl -o r1.fa.gz ~{r1}
    curl -o r2.fa.gz ~{r2}
    curl -s ~{reference_path} | gzip -d > hs1.fa
    bwa index hs1.fa
    bwa mem -5SP -T0 -t16 hs1.fa r1.fa.gz r2.fa.gz -o aligned.sam
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 8 --nproc-out 8 --chroms-path hs1.genome aligned.sam > parsed.pairsam
    mkdir ebs
    mkdir ebs/temp
    pairtools sort --nproc 32 --tmpdir=./ebs/temp/ parsed.pairsam > sorted.pairsam
    pairtools dedup --nproc-in 8 --nproc-out 8 --mark-dups --output-stats stats.txt --output dedup.pairsam sorted.pairsam
  }

  output {
    File stats = "stats.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: "64G"
    disks: "local-disk 300 SSD"
    cpu: 32
  }
}

workflow sample_fastq {
  input {
    String r1 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S2_R1_001.fastq.gz"
    String r2 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S2_R2_001.fastq.gz"
    String reference_path = "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"
    Float sampling_fraction = 0.2
  }

  call stream_and_sample {
    input:
      r1 = r1,
      r2 = r2,
      sampling_fraction = sampling_fraction,
      reference_path = reference_path
  }

  output {
    File stats = stream_and_sample.stats
  }
}
