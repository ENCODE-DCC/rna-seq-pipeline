version 1.0


task cat {
    input {
        Array[File] files
        String output_filename = "concatenated"
        Int ncpus
        Int ramGB
        String disks
    }

    command {
        cat \
            ~{sep=" " files} \
            > ~{out}
    }

    output {
        File concatenated = out
    }

    runtime {
        cpu: ncpus
        memory: "~{ramGB} GB"
        disks: disks
    }
}
