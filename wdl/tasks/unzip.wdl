version 1.0


task unzip {
    input {
        File input_file
        String output_filename = "out"
        Int cpu
        Int memory_gb
        String disk
    }

    command {
        gzip \
            -cd \
            ~{input_file} \
            > ~{output_filename}
    }

    output {
        File out = output_filename
    }

    runtime {
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: disk
    }
}
