import "../../rna-seq-pipeline.wdl" as rna

workflow test_align {
    String endedness
    Array[Array[File]] fastqs_R1
    Array[Array[File]] fastqs_R2 = []
    String bamroot
    File align_index
    Int align_ncpus
    Int align_ramGB
    String? align_disk

    Array[Array[File]] fastqs_R2_ = if (endedness == "single") then fastqs_R1 else fastqs_R2

    scatter (i in range(length(fastqs_R1))) {
        call rna.align { input:
            endedness = endedness,
            fastqs_R1 = fastqs_R1[i],
            fastqs_R2 = fastqs_R2_[i],
            index = align_index,
            bamroot = "rep"+(i+1)+bamroot,
            ncpus = align_ncpus,
            ramGB = align_ramGB,
            disks = align_disk,

        }
    }
}
