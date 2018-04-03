workflow bam_to_signals {
        Array[File] input_bams
        File chrom_sizes
        String strandedness
        String bamroot

    scatter (i in range(length(input_bams))) {
        call bam_to_signals { input:
            input_bam = input_bams[i],
            chrom_sizes = chrom_sizes,
            strandedness = strandedness,
            bamroot = "bamroot_"i"_genome",
        }
    }
}


task  bam_to_signals {
    File input_bam
    File chrom_sizes
    String strandedness
    String bamroot

    command {
        python3 $(which bam_to_signals.py) \
            --bamfile ${input_bam} \
            --chrom_sizes ${chrom_sizes} \
            --strandedness ${strandedness} \
            --bamroot ${bamroot}
    }

    output {
        Array[File] unique = glob("*niq.bw")
        Array[File] all = glob("*ll.bw")
    }

    runtime {
        docker : "quay.io/encode-dcc/rna-seq-pipeline:latest"
    }
}
