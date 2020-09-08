version 1.0

workflow runMarginVcfPhase {
	call marginVcfPhase
}

task marginVcfPhase {
    input {
        String sampleIdentifier
        File alignmentBam
        File alignmentBamIdx
        File phasingVariants
        File referenceFasta
        File parameters
        Int threadCount
        Int? memoryGBThreadRatio = 2
        String dockerImage
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

        mkdir output_mvp
        ln -s ~{alignmentBam}
        ln -s ~{alignmentBamIdx}
        cat ~{parameters} | sed 's#"include" : ".*"#"include" : "/opt/MarginPolish/params/base_params.json"#' >output_mvp/params.json

        marginVcfPhase \
            `basename ~{alignmentBam}` \
            ~{referenceFasta} \
            ~{phasingVariants} \
            output/params.json \
            -t ~{threadCount} \
            -a info \
            -o output_mvp/~{sampleIdentifier}.mvp \
            2>&1 | tee output_mvp/~{sampleIdentifier}.mvp.log

        mv output_mvp/~{sampleIdentifier}.mvp.haplotagged.bam .
        samtools index -@ ~{threadCount} ~{sampleIdentifier}.mvp.haplotagged.bam
        tar czvf ~{sampleIdentifier}.mvp.tar.gz output_mvp/*

    >>>

    output {
        File haplotaggedBam = sampleIdentifier + ".mvp.haplotagged.bam"
        File haplotaggedBamIdx = sampleIdentifier + ".mvp.haplotagged.bam.bai"
        File accessoryInformation = sampleIdentifier + ".mvp.tar.gz"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: threadCount * memoryGBThreadRatio + "GB"
    }
}
