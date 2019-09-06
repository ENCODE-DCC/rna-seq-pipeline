import "../../rna-seq-pipeline.wdl" as rna

workflow test_align_PE {
    String endedness
    Array[Array[File]] fastqs_R1
    Array[Array[File]] fastqs_R2 = []
    String bamroot 
    File align_index
    Int align_ncpus
    Int align_ramGB
    String? indexdir
    String? align_disk
}