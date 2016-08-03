count_GTF = {
    doc "Count reads aligned to genes in GTF using featureCounts"

    requires NTHREADS  : "Number of threads to use"
    requires GTF       : "Reference annotation GTF file"
    requires COUNT_DIR : "Directory for count output"

    output.dir = COUNT_DIR

    produce("counts.txt") {
        uses(threads:NTHREADS) {
            exec """
                featureCounts
                    --primary
                    -p
                    -T        $NTHREADS
                    -a        $GTF
                    -o        $output
                    $inputs.bam
            """, "count_GTF"
        }
    }

    forward inputs
}

sort_bam = {
    doc "Sort BAM file"

    requires BAM_DIR  : "Directory for BAM output"
    requires BAM_N    : "Number of threads for processing individual BAM files"

    output.dir = BAM_DIR

    uses(threads:BAM_N) {
        filter("sorted") {
            exec """
                samtools sort
                    -@ $BAM_N
                    -o $output.bam
                    $input.bam
            """, "sort_bam"
        }
    }
}

index_bam = {
    doc = "Index BAM file"

    requires BAM_DIR : "Directory for BAM output"

    output.dir = BAM_DIR

    uses(threads:1) {
        transform("bam") to ("bam.bai") {
            exec "samtools index $input.bam", "index_bam"
        }
        forward inputs
    }
}

align_stats = {
    doc = "Get alignment stats from BAM files"

    requires STAT_DIR : "Directory for stats output"
    requires GTF      : "Reference annotation GTF file"
    requires NTHREADS : "Number of threads to use"

    output.dir = STAT_DIR

    int parallel = "$NTHREADS".toInteger() - 1

    produce("align_stats.txt") {
        uses(threads:NTHREADS) {
            exec """
                alignStats
                    -o $output
                    -g $GTF
                    -i $STAT_DIR/gtf.stat.index
                    -t bam
                    -p ${parallel}
                    $inputs.bam
            """, "align_stats"
        }
    }

    forward inputs
}

sra_download = {
    doc = "Download files from the SRA"

    requires SRADB   : "Path to SRA database"
    requires SRAIDS  : "SRA ids to download"
    requires SRA_DIR : "Directory to store SRA files"

    output.dir = SRA_DIR

    produce("$output.dir/*.sra"){
        exec """
            sraDownload
                -d $SRADB
                -o $output.dir
                -t sra
                $SRAIDS
        """, "sra_download"
    }
}

fastq_dump = {
    doc = "Extract FASTQ from SRA files"

    requires FASTQ_DIR : "Directory to store dumped FASTQ files"

    output.dir = FASTQ_DIR

    transform(".sra") to ("_1.fastq.gz", "_2.fastq.gz") {
        uses(threads:1) {
            exec """
                fastq-dump
                    --split-files
                    --gzip
                    --outdir $output.dir
                    $input.sra
            """, "fastq_dump"
        }
    }
}

cleanup = {
    doc = "Cleanup RNA-seq directory after pipeline"

    requires LOG_DIR  : "Directory for log files"
    requires SRA_DIR  : "Directory to store SRA files"
    requires STAR_DIR : "Directory for STAR output"

    exec """
        echo Processing complete;
        echo The directory will now be cleaned up;
        echo .;
        echo The following actions will be performed:;
        echo You have 20 seconds to kill the pipeline if necessary;
        echo .;
        echo 1. Log files will be moved to the log directory;
        echo 2. STAR directories will be deleted;
        echo 3. Downloaded SRA files will be deleted;
        echo .;
        sleep 20;

        echo Creating log directory...;
        mkdir -p $LOG_DIR;
        echo Moving log files...;
        mv *.log $LOG_DIR;
        mv Log.* $LOG_DIR;
        mv *log.txt $LOG_DIR;
        echo Removing STAR directories...;
        rm -rf _STARtmp;
        rm -rf $STAR_DIR;
        echo Removing SRA directory...;
        rm -rf $SRA_DIR;
        echo Done!;
    """, "cleanup"
}
