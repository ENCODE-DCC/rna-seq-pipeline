version 1.0

import "../structs/runtime.wdl"

task decompress {
    input {
        File input_file
        String output_filename = "out"
        Int ncpus
        Int ramGB
        String disks
        RuntimeEnvironment runtime_environment
    }

    command {
        gzip \
            -d \
            -c \
            ~{input_file} \
            > ~{output_filename}
    }

    output {
        File out = output_filename
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
        docker: runtime_environment.docker
        singularity: runtime_environment.singularity
    }
}
