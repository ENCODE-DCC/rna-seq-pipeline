#ENCODE DCC RNA-Seq pipeline merge-annotation
#Maintainer: Otto Jolanki

#CAPER docker quay.io/encode-dcc/rna-seq-pipeline:v1.1
#CAPER singularity docker://quay.io/encode-dcc/rna-seq-pipeline:v1.1

workflow merge_anno {
    # input filenames
    File annotation
    File tRNA
    File spikeins
    # output filename
    String output_filename
    Int? cpu
    Int? memGB
    String? disks

    call merge_annotation { input :
        annotation = annotation,
        tRNA = tRNA,
        spikeins = spikeins,
        output_filename = output_filename,
    }
}

task merge_annotation {
    File annotation
    File tRNA
    File spikeins
    String output_filename
    Int? cpu
    Int? memGB
    String? disks

    command {
        python3 $(which merge_annotation.py) \
            ${"--annotation " + annotation} \
            ${"--tRNA " + tRNA} \
            ${"--spikeins " + spikeins} \
            ${"--output_filename " + output_filename}
    }
    output {
        File merged_annotation = glob("${output_filename}")[0]
    }
    runtime {
        cpu : select_first([cpu,2])
        memory : "${select_first([memGB,'8'])} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
    }
}
