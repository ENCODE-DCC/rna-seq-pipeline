#ENCODE DCC RNA-Seq pipeline prep-STAR
#Maintainer: Otto Jolanki

workflow build_index {
    # Inputs
    # reference genome in gzipped fasta
    File reference_genome
    # spikeins in gzipped fasta 
    File? spikeins
    # annotation in gzipped gtf
    File annotation
    File? tiny_fq
    # annotation version (e.g 'v24')
    String anno_version
    # genome (e.g 'GRCh38')
    String genome

    call make_index { input:
        reference_genome = reference_genome,
        spikeins = spikeins,
        annotation = annotation,
        anno_version = anno_version,
        genome = genome,
    }
}

task make_index {
    File reference_genome
    File spikeins
    File annotation
    String anno_version
    String genome

    command {
        $(which prep_star.sh) \
            ${reference_genome} \
            ${spikeins} \
            ${annotation} \
            ${anno_version} \
            ${genome}
    }
    output {
        File index = glob("*.tgz")[0]
    }
}