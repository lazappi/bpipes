fastqc = {
    doc "Quality control of fastq files using FASTQC"

    requires FQC_DIR   : "Directory to store FastQC reports"
    requires NTHREADS  : "Number of threads to use"
    requires FASTQ_DIR : "Directory with FASTQ files"

    def shortPath = {path ->
        filename = path.split("/")[-1]
        "$FASTQ_DIR/$filename"
    }

    def files = "${inputs}".split(" ")
    files = files.collect{file -> shortPath(file)}
    files = files.join(" ")

    output.dir = FQC_DIR

    uses(threads:NTHREADS) {
        transform(".fastq.gz") to ("_fastqc.zip") {
            exec """
                fastqc
                    --threads   $NTHREADS
                    --noextract
                    -o          $output.dir
                    $files
            """
        }
        forward inputs
    }
}

fastq_screen = {
    doc "Screening of fastq files using fastq_screen"

    requires FQS_DIR   : "Directory for fastq_screen output"
    requires NTHREADS  : "Number of threads to use"
    requires FASTQ_DIR : "Directory containing FASTQ files"

    def shortPath = {path ->
        filename = path.split("/")[-1]
        "$FASTQ_DIR/$filename"
    }

    def files = "${inputs}".split(" ")
    files = files.collect{file -> shortPath(file)}
    files = files.join(" ")

    output.dir = FQS_DIR

    uses(threads:NTHREADS) {
        transform(".fastq.gz") to ("_screen.txt") {
            exec """
                fastq_screen
                    --aligner bowtie2
                    --force
                    --outdir $output.dir
                    --subset 100000
                    --threads $NTHREADS
                    $files
            """
        }
        forward inputs
    }
}

multiqc = {
    doc "Build a report using MultiQC"

    requires MQC_DIR : "Directory for MultiQC output"

    output.dir = MQC_DIR

    uses(threads:1) {
        produce("multiqc_report.html") {
            exec """
                multiqc
                    --force
                    --outdir $output.dir
                    .
            """
        }
    }
}
