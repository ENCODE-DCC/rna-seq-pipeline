#ENCODE DCC RNA-Seq pipeline prep-tophat
#Maintainer: Otto Jolanki

workflow build_index {
    # Inputs
    # reference genome in gzipped fasta
    File reference_genome
    # spikeins in gzipped fasta 
    File? spikeins
    # annotation in gzipped gtf
    File annotation
    # fake gzipped fastq to trick tophat
    File? tiny_fq
    # annotation version (e.g 'v24')
    String anno_version
    # genome (e.g 'GRCh38')
    String genome

    call make_index { input:
        reference_genome = reference_genome,
        spikeins = spikeins,
        annotation = annotation,
        tiny_fq = tiny_fq,
        anno_version = anno_version,
        genome = genome,
    }
}

task make_index {
    File reference_genome
    File spikeins
    File annotation
    File tiny_fq
    String anno_version
    String genome

    command {
        $(which prep_tophat.sh) \
            ${reference_genome} \
            ${spikeins} \
            ${annotation} \
            ${tiny_fq} \
            ${anno_version} \
            ${genome}
    }
    output {
        File index = glob("*.tgz")[0]
    }
}