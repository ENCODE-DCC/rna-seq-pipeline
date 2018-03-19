# ENCODE DCC RNA-seq pipeline
# Maintainer: Otto Jolanki

workflow rna {
    # endedness: paired or single
    String endedness
    # fastqs_R1: fastq.gz files for Read1 (only one if single-ended)
    Array[File] fastqs_R1
    # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
    # corresponding to fastqs_R1
    Array[File] fastqs_R2 = [] 
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
    String bamroot = ""
    # strandedness: is the library strand specific (stranded or unstranded)
    String strandedness = "unstranded"
    # chrom_sizes: chromosome sizes file
    File chrom_sizes = ""

    Int? align_ncpus

    Int? align_ramGB

    Array[Array[File]] fastqs_ = if length(fastqs_R2)>0 then transpose([fastqs_R1, fastqs_R2]) else transpose([fastqs_R1])

    scatter (i in range(length(fastqs_))) {
        call align { input:
            endedness = endedness,
            fastqs = fastqs_[i],
            index = index,
            aligner = aligner,
            indexdir = indexdir,
            libraryid = libraryid,
            bamroot = "rep"+(i+1)+bamroot,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
        }
        call bam_to_signals as genome_signal { input:
            input_bam = align.genomebam,
            chrom_sizes = chrom_sizes,
            strandedness = strandedness,
            bamroot = "rep"+(i+1)+bamroot,
        }

        call bam_to_signals as anno_signal { input:
            input_bam = align.annobam,
            chrom_sizes = chrom_sizes,
            strandedness = strandedness,
            bamroot = "rep"+(i+1)+bamroot,
        }
    }
}


    ## tasks
    task align {
        Array[File] fastqs
        String endedness
        String aligner
        File index
        String? indexdir
        String? libraryid
        String bamroot
        Int? ncpus
        Int? ramGB

        command {
            python3 $(which align.py) \
                ${if length(fastqs)<2 then "--fastqs " + fastqs[0] else "--fastqs " + fastqs[0] + " " + fastqs[1]} \
                --endedness ${endedness} \
                --aligner ${aligner} \
                --index ${index} \
                ${"--indexdir " + indexdir} \
                ${"--libraryid " + libraryid} \
                ${"--bamroot " + bamroot} \
                ${"--ncpus " + ncpus} \
                ${"--ramGB " + ramGB}
        }

        output {
            File genomebam = glob("*_genome.bam")[0]
            File annobam = glob("*_anno.bam")[0]
            File genome_flagstat = glob("*_genome_flagstat.txt")[0]
            File anno_flagstat = glob("*_anno_flagstat.txt")[0]
            File log = glob("*_Log.final.out")[0]
        }

        runtime {
        docker : "quay.io/encode-dcc/rna-seq-pipeline:latest"
        dx_instance_type : "mem3_ssd1_x16"
        }
    }

    task  bam_to_signals {
        File input_bam
        File chrom_sizes
        File strandedness
        String bamroot

        command {
            python3 $(which bam_to_signals.py) \
                --bamfile ${input_bam} \
                --chrom_sizes ${chrom_sizes}
                --strandedness ${strandedness}
                --bamroot ${bamroot}
        }

        output {
            Array[File] unique = glob("*niq.bw")
            Array[File] all = glob("*ll.bw")
        }

        runtime {
            docker : "quay.io/encode-dcc/rna-seq-pipeline:latest"
            dx_instance_type : "mem3_ssd1_x16"
        }
    }
