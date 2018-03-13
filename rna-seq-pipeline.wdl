# ENCODE DCC RNA-seq pipeline
# Maintainer: Otto Jolanki

workflow rna {
    # endedness: paired or single
    String endedness
    # fastqs_R1: fastq.gz files for Read1 (only one if single-ended)
    Array[File] fastqs_R1
    # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
    # corresponding to fastqs_R1
    Array[File]? fastqs_R2 
    # aligner: star for now, more added if/when needed
    String aligner
    # index: aligner index (tar.gz)
    File index
    # indexdir: where to extract the index, relative to cwd
    String? indexdir
    # libraryid: identifier which will be added to bam headers
    String? libraryid
    # bamroot: root name for output bams. For example foo_bar will
    # create foo_bar_genome.bam and foo_bar_anno.bam
    String? bamroot

    call align { input:
        endedness = endedness,
        fastq_R1 = fastqs_R1[0],
        fastq_R2 = fastqs_R2[0],
        endedness = endedness,
        aligner = aligner,
        indexdir = indexdir,
        libraryid = libraryid,
        bamroot = bamroot,
    }


    ## tasks
    task align {
        File fastq_R1
        File? fastq_R2
        String endedness
        String aligner
        File index
        String? indexdir
        String? libraryid
        String? bamroot

        command {
            python3 $(which aligner.py) \
                --fastqs ${fastq_R1} ${fastq_R2} \
                --endedness ${endedness} \
                --aligner ${aligner} \
                --index ${index} \
                ${"--indexdir " + indexdir} \
                ${"--libraryid " + libraryid} \
                ${"--bamroot " + bamroot}
        }

        output{
            File genomebam = glob("*_genome.bam")[0]
            File annobam = glob("*_anno.bam")[0]
            File log = glob("*_Log.final.out")[0]
        }
    }
}