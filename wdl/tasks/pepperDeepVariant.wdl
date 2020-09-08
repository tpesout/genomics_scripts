version 1.0

workflow runPepperDeepVariant {
	call pepperDeepVariant
}

task pepperDeepVariant {
    input {
        String sampleIdentifier
        File haplotaggedBam
        File haplotaggedBamIdx
        File referenceFasta
        String regionName
        Int threadCount
        Int? memoryGBThreadRatio = 2
        String dockerImage="kishwars/pepper_deepvariant_cpu:latest"
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


        mkdir output_pdv
        ln -s ~{haplotaggedBam}
        ln -s ~{haplotaggedBamIdx}

        /bin/bash /opt/run_pepper_deepvariant_only_hp.sh \
            -b `basename ~{haplotaggedBam}` \
            -f ~{referenceFasta} \
            -t ~{threadCount} \
            -s ~{sampleIdentifier} \
            -r ~{regionName} \
            -o output_pdv/

        cp output_pdv/PEPPER_HP_DEEPVARIANT_FINAL_OUTPUT.vcf.gz .
        gunzip PEPPER_HP_DEEPVARIANT_FINAL_OUTPUT.vcf.gz
        cp PEPPER_HP_DEEPVARIANT_FINAL_OUTPUT.vcf.gz ~{sampleIdentifier}.PEPPER_HP_DEEPVARIANT_FINAL_OUTPUT.vcf.gz
        tar czvf ~{sampleIdentifier}.pdv.tar.gz output_pdv/*

    >>>

    output {
        File pepperDeepVariantVCF = sampleIdentifier + ".PEPPER_HP_DEEPVARIANT_FINAL_OUTPUT.vcf.gz"
        File accessoryInformation = sampleIdentifier + ".pdv.tar.gz"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: threadCount * memoryGBThreadRatio + "GB"
    }
}
