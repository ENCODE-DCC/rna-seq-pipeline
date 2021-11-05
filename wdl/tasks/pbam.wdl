version 1.0

import "../structs/runtime.wdl"

task make_genome_pbam {
    input {
        File bam
        File reference_genome
        Int ncpus
        Int ramGB
        String disks
        RuntimeEnvironment runtime_environment
    }

    String bam_base = basename(bam, ".bam")

    command {
        makepBAM_genome.sh \
            ~{bam} \
            ~{reference_genome}
    }

    output {
        File out = "~{bam_base}.sorted.p.bam"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}

task make_transcriptome_pbam {
    input {
        File bam
        File reference_genome
        File reference_transcriptome
        File reference_annotation
        Int ncpus
        Int ramGB
        String disks
        RuntimeEnvironment runtime_environment
    }

    String bam_base = basename(bam, ".bam")

    command {
        makepBAM_transcriptome.sh \
            ~{bam} \
            ~{reference_genome} \
            ~{reference_transcriptome} \
            ~{reference_annotation}
    }

    output {
        File out = "~{bam_base}.p.bam"
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
