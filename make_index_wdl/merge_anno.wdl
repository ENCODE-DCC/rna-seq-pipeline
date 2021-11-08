version 1.0

#ENCODE DCC RNA-Seq pipeline merge-annotation


workflow merge_anno {
    meta {
        author: "Otto Jolanki"
        version: "1.2.4"
        caper_docker: "encodedcc/rna-seq-pipeline:1.2.4"
        caper_singularity:"docker://encodedcc/rna-seq-pipeline:1.2.4"
    }

    input {
        File annotation
        File tRNA
        File spikeins
        String output_filename
        Int? cpu
        Int? memGB
        String? disks
        String docker = "encodedcc/rna-seq-pipeline:1.2.4"
        String singularity = "docker://encodedcc/rna-seq-pipeline:1.2.4"
    }

    call merge_annotation { input :
        annotation=annotation,
        tRNA=tRNA,
        spikeins=spikeins,
        output_filename=output_filename,
        runtime_environment=runtime_environment,
    }
}

task merge_annotation {
    input {
        File annotation
        File tRNA
        File spikeins
        String output_filename
        Int? cpu
        Int? memGB
        String? disks
        RuntimeEnvironment runtime_environment
    }

    command {
        python3 $(which merge_annotation.py) \
            ~{"--annotation " + annotation} \
            ~{"--tRNA " + tRNA} \
            ~{"--spikeins " + spikeins} \
            ~{"--output_filename " + output_filename}
    }
    output {
        File merged_annotation = output_filename
    }
    runtime {
        cpu : select_first([cpu,2])
        memory : "~{select_first([memGB,'8'])} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
