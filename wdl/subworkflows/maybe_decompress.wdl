version 1.0


import "../tasks/unzip.wdl"


workflow maybe_decompress {
    input {
        File? input_plain
        File? input_gz
        Int cpu
        Int memory_gb
        String disk
    }

    if (defined(input_gz)) {
        call unzip.unzip {
            input:
                input_file=select_first([input_gz]),
                cpu=cpu,
                memory_gb=memory_gb,
                disk=disk,
        }
    }

    output {
        File out = select_first([input_plain, unzip.out])
    }
}
