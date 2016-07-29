star_1pass_SE = {
    doc "Map single-end reads using the STAR aligner: 1st pass"

    requires MAP_DIR  : "Missing directory for mapping output"
    requires NTHREADS : "Missing number of threads to use"
    requires GEN_DIR  : "Missing genome directory"
    requires GTF      : "Missing GTF file"

    String files = "${inputs}";
    files        = files.replaceAll(' ', ',');
    output.dir   = MAP_DIR + "/1pass"

    exec """
        STAR
            --runThreadN        $NTHREADS
            --genomeDir         $GEN_DIR
            --readFilesIn       ${files}
            --readFilesCommand  zcat
            --sjdbGTFfile       $GTF
            --sjdbOverhang      100
            --outSAMtype        None
            --outFileNamePrefix ${output.dir}/
    """
}

star_2pass_genome = {
    doc "Generate STAR genome for 2nd pass mapping"

    requires STAR_DIR  : "STAR output directory"
    requires GEN_FASTA : "Reference genome FASTA file"
    requires NTHREADS  : "Number of threads to use"

    output.dir = STAR_DIR + "/genome_2pass"

    produce("Genome") {
        uses(threads:NTHREADS) {
            exec """
                STAR
                    --runMode             genomeGenerate
                    --runThreadN          $NTHREADS
                    --genomeDir           ${output.dir}
                    --genomeFastaFiles    $GEN_FASTA
                    --sjdbGTFfile         $GTF
                    --sjdbOverhang        100
                    --limitSjdbInsertNsj  2000000
                    --sjdbFileChrStartEnd ${STAR_DIR}/1pass/SJ.out.tab
                    --outFileNamePrefix   ${output.dir}/
            """
        }
        forward inputs
    }
}

star_2pass_load = {
    doc "Load STAR genome into shared memory for 2nd pass mapping"

    requires STAR_DIR : "Missing directory for mapping output"

    exec """
        STAR
            --genomeDir  ${STAR_DIR}/genome_2pass
            --genomeLoad LoadAndExit
    """

    forward inputs
}

star_2pass_remove = {
    doc "Remove STAR genome from shared memory"

    requires STAR_DIR : "STAR output directory"

    exec """
        STAR
            --genomeDir  ${STAR_DIR}/genome_2pass
            --genomeLoad Remove
    """

    forward inputs
}

star_2pass_SE = {
    doc "Map paired-end reads using the STAR aligner: 2nd pass"

    requires MAP_DIR  : "Missing directory for mapping output"
    requires GEN_DIR  : "Missing genome directory"
    requires GTF      : "Missing GTF file"

    output.dir = MAP_DIR + "/2pass"

    // Assumes file names
    transform("(.*).fastq.gz") to ("\$1.Aligned.out.bam") {
        exec """

            STAR
                --genomeDir         mapped/genome_2pass
                --genomeLoad        LoadAndRemove
                --readFilesIn       ${input}
                --readFilesCommand  zcat
                --outSAMtype        BAM Unsorted SortedByCoordinate
                --outFileNamePrefix ${output.prefix.prefix.prefix}.
        """
    }
}

star_1pass_PE = {
    doc "Map paired-end reads using the STAR aligner: 1st pass"

    requires STAR_DIR  : "STAR output directory"
    requires NTHREADS  : "Number of threads to use"
    requires GEN_DIR   : "Directory containing STAR index"
    requires GTF       : "GTF file containing reference annotation"
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
    R1_files = R1_files.join(",")
    R2_files = R2_files.join(",")

    output.dir   = STAR_DIR + "/1pass"

    produce("SJ.out.tab") {
        uses(thread:NTHREADS) {
            exec """
                STAR
                    --runThreadN          $NTHREADS
                    --genomeDir           $GEN_DIR
                    --readFilesIn         ${R1_files} ${R2_files}
                    --readFilesCommand    zcat
                    --sjdbGTFfile         $GTF
                    --sjdbOverhang        100
                    --limitOutSJcollapsed 5000000
                    --limitIObufferSize   300000000
                    --outSAMtype          None
                    --outFileNamePrefix   ${output.dir}/
            """
        }
    }

    forward inputs
}

star_2pass_PE_GTF = {
    doc "Map paired-end reads using the STAR aligner: 2nd pass"

    requires STAR_DIR : "STAR output directory"
    requires GEN_DIR  : "Directory containing STAR index"
    requires GTF      : "GTF file containing reference genome"

    output.dir = STAR_DIR + "/2pass"

    // Assumes file names
    transform("(.*)_1.fastq.gz", "(.*)_2.fastq.gz") to ("\$1.Aligned.out.bam") {
        exec """
            STAR
                --genomeDir           $GEN_DIR
                --readFilesIn         $input1.gz $input2.gz
                --readFilesCommand    zcat
                --sjdbGTFfile         $GTF
                --sjdbOverhang        100
                --sjdbFileChrStartEnd ${STAR_DIR}/1pass/SJ.out.tab
                --outSAMtype          BAM Unsorted SortedByCoordinate
                --outFileNamePrefix   ${output.prefix.prefix.prefix}.
        """
    }
}
star_2pass_PE = {
    doc "Map paired-end reads using the STAR aligner: 2nd pass"

    requires STAR_DIR : "STAR output directory"
    requires STAR_N   : "Number of threads for individual STAR 2pass"

    output.dir = STAR_DIR + "/2pass"

    // Assumes file names
    transform("(.*)_1.fastq.gz", "(.*)_2.fastq.gz") to ("\$1.Aligned.out.bam") {
        uses(threads:STAR_N) {
            exec """
                STAR
                    --runThreadN          $STAR_N
                    --genomeDir           ${STAR_DIR}/genome_2pass
                    --genomeLoad          LoadAndKeep
                    --readFilesIn         $input1.gz $input2.gz
                    --readFilesCommand    zcat
                    --outSAMtype          BAM Unsorted
                    --outSAMunmapped      Within KeepPairs
                    --outFileNamePrefix   ${output.prefix.prefix.prefix}.
                    --limitBAMsortRAM     10000000000
            """
        }
    }
}

star_stats = {
    doc = "Get stats from STAR log files"

    requires STAT_DIR : "Stats output directory"

    String files  = "${inputs}";
    def log_files = files.replaceAll('Aligned.out.bam', 'Log.final.out');
    output.dir    = STAT_DIR

    produce("star_stats.txt") {
        exec """starStats -o $output ${log_files}"""
    }

    forward inputs
}
