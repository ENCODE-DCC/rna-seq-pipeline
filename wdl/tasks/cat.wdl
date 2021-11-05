version 1.0

import "../structs/runtime.wdl"

task cat {
    input {
        Array[File] files
        String output_filename = "concatenated"
        Int ncpus
        Int ramGB
        String disks
        RuntimeEnvironment runtime_environment
    }

    command {
        cat \
            ~{sep=" " files} \
            > ~{output_filename}
    }

    output {
        File concatenated = output_filename
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
