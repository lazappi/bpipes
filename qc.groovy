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

