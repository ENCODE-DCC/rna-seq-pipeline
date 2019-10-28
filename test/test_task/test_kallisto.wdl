import "../../rna-seq-pipeline.wdl" as rna

workflow test_kallisto {
    Array[Array[File]] fastqs_R1
    Array[Array[File]] fastqs_R2 = []
    File kallisto_index
    String endedness
    String strandedness_direction
    Int kallisto_number_of_threads
    Int kallisto_ramGB
    String out_prefix
    String kallisto_disk
    Int? kallisto_fragment_length
    Float? kallisto_sd_of_fragment_length

    Array[Array[File]] fastqs_R2_ = if (endedness == "single") then fastqs_R1 else fastqs_R2

    scatter (i in range(length(fastqs_R1))) {
        call rna.kallisto { input:
            fastqs_R1 = fastqs_R1[i],
            fastqs_R2 = fastqs_R2_[i],
            endedness = endedness,
            strandedness_direction = strandedness_direction,
            kallisto_index = kallisto_index,
            number_of_threads = kallisto_number_of_threads,
            ramGB = kallisto_ramGB,
            fragment_length = kallisto_fragment_length,
            sd_of_fragment_length = kallisto_sd_of_fragment_length,
            disks = kallisto_disk,
            out_prefix = "rep"+(i+1)+out_prefix,
        }
    }
}
