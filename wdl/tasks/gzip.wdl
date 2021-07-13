version 1.0


task decompress {
    input {
        File input_file
        String output_filename = "out"
        Int ncpus
        Int ramGB
        String disks
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
    }
