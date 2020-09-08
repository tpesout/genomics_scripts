version 1.0

workflow runHappy {
	call happy
}

task happy {
    input {
        String sampleIdentifier
        File truthVcf
        File queryVcf
        File highConfBed
        File referenceFasta
        File referenceFastaIdx
        Int threadCount
        Int memoryGBThreadRatio = 2
        String dockerImage="pkrusche/hap.py"
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace


        mkdir output_happy
        ln -s ~{truthVcf}
        ln -s ~{queryVcf}
        ln -s ~{referenceFasta}
        ln -s ~{referenceFastaIdx}

        /opt/hap.py/bin/hap.py \
            `basename ~{truthVcf}` \
            `basename ~{queryVcf}` \
            -f ~{highConfBed} \
            -r `basename ~{referenceFasta}` \
            -o output_happy/~{sampleIdentifier} \
            --pass-only \
            --engine=vcfeval \
            --threads=~{threadCount}

        cp output_happy/~{sampleIdentifier}.summary.csv .
        tar czvf ~{sampleIdentifier}.happy.tar.gz output_happy/*

    >>>

    output {
        File pepperDeepVariantVCF = sampleIdentifier + ".summary.csv"
        File fullOutput = sampleIdentifier + ".happy.tar.gz"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: threadCount * memoryGBThreadRatio + "GB"
    }
}
