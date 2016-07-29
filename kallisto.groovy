kallisto_multi = {
    doc "Quantify transcript abundances of multiple samples using kallisto"

    requires KAL_DIR   : "kallisto output directory"
    requires KAL_IDX   : "Path to kallisto index"
    requires NTHREADS  : "Number of threads to use"
    requires FASTQ_DIR : "Directory for FASTQ files"

    def shortPath = {path ->
        filename = path.split("/")[-1]
        "$FASTQ_DIR/$filename"
    }

    String in_files = "${inputs}";
    // Divide input files into R1 and R2
    def all_files = in_files.split(" "); // Split file list on spaces
    all_files.sort()                     // Sort alphabetically
    all_files = all_files.collect{file -> shortPath(file)}

    def R1_files = []
    def R2_files = []
    all_files.eachWithIndex{item, index ->
        if(index % 2 == 0) {         // All R1 files
            R1_files.add(item)
        } else {
            R2_files.add(item)
        }
    }

    // Join file lists with ','
    R1_files = R1_files.join(" ")
    R2_files = R2_files.join(" ")

    output.dir = KAL_DIR

    produce("kallistoMulti.log") {
        uses(thread:NTHREADS) {
            exec """
                kallistoMulti
                    --index             $KAL_IDX
                    --output-dir        $KAL_DIR
                    --read1             ${R1_files}
                    --read2             ${R2_files}
                    --bootstrap-samples 20
                    --threads           $NTHREADS
            """
        }
    }

    forward inputs
}

kallisto_merge = {
    doc "Merge kallisto output"

    requires KAL_DIR   : "kallisto output directory"
    requires TX2GENE   : "File mapping transcripts to genes"
    requires COUNT_DIR : "Directory for count output"

    output.dir = COUNT_DIR

    produce("kallisto_gene_counts.tsv") {
        uses(thread:1) {
            exec """
                kallistoMerge
                    --directory $KAL_DIR
                    --tx2gene   $TX2GENE
                    --outpath   $COUNT_DIR
            """
        }
    }

    forward inputs
}
