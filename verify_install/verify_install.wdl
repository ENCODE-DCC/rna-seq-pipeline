# ENCODE DCC Install verifier workflow
# Maintainer: Otto Jolanki

workflow verify_install {
    # endedness: paired or single
    String endedness
    # fastqs_R1: fastq.gz files for Read1 (only these if single-ended)
    Array[File] fastqs_R1
    # fastqs_R2: fastq.gz files for Read2 (omit if single-ended) in order
    # corresponding to fastqs_R1
    Array[File] fastqs_R2 = [] 
    # aligner: star for now, more added if/when needed
    String aligner
    # index: aligner index (tar.gz)
    File index
    # bamroot: root name for output bams. For example foo_bar will
    # create foo_bar_genome.bam and foo_bar_anno.bam
    String bamroot 
    # strandedness: is the library strand specific (stranded or unstranded)
    String strandedness 
    # strandedness_direction (forward, reverse, unstranded)
    String strandedness_direction
    # chrom_sizes: chromosome sizes file
    File chrom_sizes 
    # rsem_index: location of the RSEM index archive (tar.gz)
    File rsem_index
    # rnd_seed: random seed used for rsem
    Int rnd_seed = 12345
    # indexdir: where to extract the star index, relative to cwd
    String? indexdir
    # libraryid: identifier which will be added to bam headers
    String? libraryid

    Int align_ncpus = align_ncpus

    Int align_ramGB = align_ramGB

    # correct md5 sums go here
    File comparison_results_json

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
            bamroot = "rep"+(i+1)+bamroot+"_genome",
        }

        call rsem_quant { input:
            rsem_index = rsem_index,
            rnd_seed = rnd_seed,
            anno_bam = align.annobam,
            endedness = endedness,
            read_strand = strandedness_direction,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
        }

        call compare_md5 { input:
            inputs = [align.anno_flagstat, align.genome_flagstat, bam_to_signals.unique[0], 
                      bam_to_signals.all[0], rsem_quant.genes_results, rsem_quant.isoforms_results]

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
        Int ncpus
        Int ramGB

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
        }
    }

    task  bam_to_signals {
        File input_bam
        File chrom_sizes
        String strandedness
        String bamroot

        command {
            python3 $(which bam_to_signals.py) \
                --bamfile ${input_bam} \
                --chrom_sizes ${chrom_sizes} \
                --strandedness ${strandedness} \
                --bamroot ${bamroot}
        }

        output {
            Array[File] unique = glob("*niq.bw")
            Array[File] all = glob("*ll.bw")
        }

        runtime {
            docker : "quay.io/encode-dcc/rna-seq-pipeline:latest"
        }
    }

    task rsem_quant {
        File rsem_index
        File anno_bam
        String endedness
        String read_strand
        Int rnd_seed
        Int ncpus
        Int ramGB

        command {
            python3 $(which rsem_quant.py) \
                --rsem_index ${rsem_index} \
                --anno_bam ${anno_bam} \
                --endedness ${endedness} \
                --read_strand ${read_strand} \
                --rnd_seed ${rnd_seed} \
                --ncpus ${ncpus} \
                --ramGB ${ramGB}
        }

        output {
            File genes_results = glob("*.genes.results")[0]
            File isoforms_results = glob("*.isoforms.results")[0]
        }

        runtime {
            docker : "quay.io/encode-dcc/rna-seq-pipeline:latest"
        }
    }

    task compare_md5 {
    Array[File] inputs
    File reference_json
        command {
            python3 $(which compare_md5.py) \
                --input_files ${sep=' ' inputs} \
                --reference_json ${reference_json} \
                --outfile comparison_result.json
        }

        output {
            File comparison_result = glob("comparison_result.json")[0]
        }

        runtime {
            docker: "quay.io/encode-dcc/rna-seq-pipeline:latest" 
        }

    }