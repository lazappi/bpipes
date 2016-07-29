describe = {
    doc "Describe current pipeline"

    println "The following variables have been set:"

    sleep(1000)

    println "You will have 20 seconds to kill the pipeline if necessary"
    println()

    sleep(2000)

    def out = new StringBuffer()

    binding.variables.each{
        def key = it.key
        def val = it.value

        if(val && !(key == "BPIPE_NO_EXTERNAL_STAGES")) {
            if(!val.toString().contains("_run_closure")) {
                out << key.padRight(15) << val << "\n"
            }
        }
    }

    println out.toString()

    sleep(20000)

    print "Starting pipeline now..."

    sleep(1000)
}

