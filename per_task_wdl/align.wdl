workflow align{
    String endedness
    Array[File] fastqs_R1
    Array[File] fastqs_R2 = []
    String aligner = "star"
    File index
    String bamroot
    String? indexdir
    String? libraryid
    Int ncpus
    Int ramGB

    Array[Array[File]] fastqs_ = if length(fastqs_R2)>0 then transpose([fastqs_R1, fastqs_R2]) else transpose([fastqs_R1])

    scatter (i in range(length(fastqs_))) {
        call align_ { input:
            endedness = endedness,
            fastqs = fastqs_[i],
            index = index,
            aligner = aligner,
            indexdir = indexdir,
            libraryid = libraryid,
            bamroot = "rep"+(i+1)+bamroot,
            ncpus = ncpus,
            ramGB = ramGB,
        }
    }
}

task align_ {
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
          dx_instance_type : "mem3_ssd1_x16"
        }
    }