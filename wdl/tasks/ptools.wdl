version 1.0

task makepbam_genome {
    input {
        File bam
        File reference_fasta
        Int ncpus
        Int ramGB
        String disks
    }

    String bam_prefix = basename(bam, ".bam")
    String out = bam_prefix + ".sorted.p.bam"

    command {
        makepBAM_genome.sh \
            ~{bam} \
            ~{reference_fasta}
    }

    output {
        File pbam = out
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
    }
}


task makepbam_transcriptome {
    input {
        File bam
        File genome_fasta
        File transcriptome_fasta
        File annotation_gtf
        Int ncpus
        Int ramGB
        String disks
    }

    String bam_prefix = basename(bam, ".bam")
    String out = bam_prefix + ".p.bam"

    command {
        $(which makepBAM_transcriptome.sh) \
            ~{bam} \
            ~{genome_fasta} \
            ~{transcriptome_fasta} \
            ~{annotation_gtf}
    }

    output {
        File pbam = out
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
    }
}
