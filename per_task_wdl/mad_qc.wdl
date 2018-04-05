workflow mad_qc {
    File quants1
    File quants2

    call mad_qc_ { input:
        quants1 = quants1,
        quants2 = quants2,
    }
}

task mad_qc_ {
    File quants1
    File quants2

    command {
        python3 $(which mad_qc.py) \
            --quants1 ${quants1} \
            --quants2 ${quants2} \
            --MAD_R_path $(which MAD.R)
    }

    output {
        File madQCplot = glob("*_mad_plot.png")[0]
        File madQCmetrics = glob("*_mad_qc_metrics.json")[0]
    }

    runtime {
        docker : "quay.io/encode-dcc/rna-seq-pipeline:latest"
    }
}