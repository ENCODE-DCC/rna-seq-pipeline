version 1.0


import "../tasks/ptools.wdl"


workflow make_transcriptome_pbam {
    input {
        File bam
        File genome_fasta
        File transcriptome_fasta
        File annotation_gtf
        Int ncpus
        Int ramGB
        String disks
    }

    call ptools.makepbam_transcriptome {
        input:
            bam=bam,
            genome_fasta=genome_fasta,
            transcriptome_fasta=transcriptome_fasta,
            annotation_gtf=annotation_gtf,
            ncpus=ncpus,
            ramGB=ramGB,
            disks=disks,
    }
}
