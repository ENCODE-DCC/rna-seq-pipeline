# ENCODE DCC RNA-Seq pipeline build_genome_index
# Maintainer: Otto Jolanki

workflow build_index {
    # Inputs
    # reference genome in gzipped fasta
    File reference_genome
    # spikeins in gzipped fasta 
    File? spikeins
    # annotation in gzipped gtf
    File annotation
    # Tiny fastq use to fake an alignment, required for complete tophat build.
    # Available in   
    File? tiny_fq
    # annotation version (e.g 'v24')
    String anno_version
    # genome (e.g 'GRCh38')
    String genome
    # Flavor of the index that gets built
    # available options:
    # prep_rsem, prep_srna, prep_star, prep_tophat
    String index_type


    call make_index { input:
        reference_genome = reference_genome,
        spikeins = spikeins,
        annotation = annotation,
        tiny_fq = tiny_fq,
        anno_version = anno_version,
        genome = genome,
        index_type = index_type,
    }
}

task make_index {
    File reference_genome
    File? spikeins
    File annotation
    File? tiny_fq
    String anno_version
    String genome
    String index_type
    
    command {
        $(which ${index_type + ".sh"}) \
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