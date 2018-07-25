# ENCODE DCC Install verifier workflow
# Maintainer: Otto Jolanki

import "./rna-seq-pipeline.wdl" as rna

workflow test_wf {
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

    String? disks

    # correct md5 sums go here
    File comparison_results_json

    Array[Array[File]] fastqs_ = if length(fastqs_R2)>0 then transpose([fastqs_R1, fastqs_R2]) else transpose([fastqs_R1])

    scatter (i in range(length(fastqs_))) {
        call rna.align as align { input:
            endedness = endedness,
            fastqs = fastqs_[i],
            index = index,
            aligner = aligner,
            indexdir = indexdir,
            libraryid = libraryid,
            bamroot = "rep"+(i+1)+bamroot,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
            disks = disks,
        }

        call rna.bam_to_signals as genome_signal { input:
            input_bam = align.genomebam,
            chrom_sizes = chrom_sizes,
            strandedness = strandedness,
            bamroot = "rep"+(i+1)+bamroot+"_genome",
            ncpus = align_ncpus,
            ramGB = align_ramGB,
            disks = disks,
        }

        call rna.rsem_quant as rsem_quant{ input:
            rsem_index = rsem_index,
            rnd_seed = rnd_seed,
            anno_bam = align.annobam,
            endedness = endedness,
            read_strand = strandedness_direction,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
            disks = disks,
        }

        call rna.compare_md5 { input:
            inputs = [align.anno_flagstat, align.genome_flagstat, genome_signal.unique[0], 
                      genome_signal.all[0]],
            reference_json = comparison_results_json
        }
    }
}
