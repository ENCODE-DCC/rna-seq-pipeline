# ENCODE DCC RNA-Seq pipeline build_genome_index
# Maintainer: Otto Jolanki

#CAPER docker quay.io/encode-dcc/rna-seq-pipeline:v1.1
#CAPER singularity docker://quay.io/encode-dcc/rna-seq-pipeline:v1.1

workflow build_index {
    # Inputs
    # reference genome or transcriptome (in prep_kallisto mode)in gzipped fasta
    File reference_sequence
    # spikeins in gzipped fasta
    File? spikeins
    # annotation in gzipped gtf
    File? annotation
    # annotation version (e.g 'v24')
    String? anno_version
    # genome (e.g 'GRCh38')
    String? genome
    # Flavor of the index that gets built
    # available options:
    # prep_rsem, prep_srna, prep_star, prep_kallisto
    String index_type
    Int ncpu = 8
    Int? memGB
    String? disks

    call make_index { input:
        reference_sequence = reference_sequence,
        spikeins = spikeins,
        annotation = annotation,
        anno_version = anno_version,
        genome = genome,
        index_type = index_type,
        ncpu = ncpu,
        memGB = memGB,
        disks = disks,
    }
}

task make_index {
    File reference_sequence
    File? spikeins
    File? annotation
    String? anno_version
    String? genome
    String index_type
    Int ncpu
    Int? memGB
    String? disks

    command {
        $(which ${index_type + ".sh"}) \
            ${reference_sequence} \
            ${spikeins} \
            ${annotation} \
            ${anno_version} \
            ${genome} \
            ${ncpu}
    }

    output {
        File index = if (index_type == "prep_kallisto") then glob("*.idx")[0] else glob("*.tgz")[0]
    }

    runtime {
        cpu : ncpu
        memory : "${select_first([memGB,'8'])} GB"
        disks : select_first([disks,"local-disk 100 SSD"])
    }
}
