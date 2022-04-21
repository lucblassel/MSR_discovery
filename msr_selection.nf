ssrDir = file("data/SSRs")
resDir = file("results/SSR_eval")

coverage = 0.1 // PBSim setting
number_of_reads = 35000 // NanoSim setting

// Input Data
reference = file("data/whole_human_genome.fa")
pbsim_model = file("data/models/P6C4.model")
nanosim_model = file("data/models/human_NA12878_DNA_FAB49712_guppy_flipflop/training")

process make_MSRs {
    publishDir "$ssrDir", mode: "copy"
    output:
    file "*.json" into msrs

    script:
    """
    makeAllUniqueFunctions.py -o .
    """
}

process nanosim {
    container "docker://lucblassel/nanosim:latest"

    output:
    tuple val("nanosim"), file("nanosim.fa") into nanosim

    memory 20.GB
    cpus 16

    script:
    """
    nanosim.py genome \
        -rg $reference \
        -n $number_of_reads \
        -c $nanosim_model \
        -o simulated \
        -t 16 \
        -dna_type linear

    convertNanoSimNames.py \
        --fasta simulated_aligned_reads.fasta \
        --errorProfile simulated_aligned_error_profile \
        --output nanosim.fa
    """
}

process pbsim {
    output:
    tuple val("pbsim"), file("pbsim.fa") into pbsim

    memory 20.GB
    cpus 16

    script:
    """
    samtools faidx $reference -o ref.fai
    pbsim --depth $coverage --hmm_model $pbsim_model $reference
    paftoolsCustom.js pbsim2fq ref.fai *maf > pbsim.fa
    """
}

msrs
    .flatMap()
    .into{ msrs_reference; msrs_reads }

process reduce_reference {
    input:
    file msr from msrs_reference

    output:
    tuple val("${msr.baseName}"), file("reduced.fa") into reduced_reference_sequence
    tuple val("${msr.baseName}"), file("offsets.json") into reduced_reference_offset

    memory 30.GB
    cpus 16

    script:
    """
    echo "reducing reference with ${msr.baseName}"

    if [[ ${msr.baseName} == "raw" ]]; then
        cp $reference reduced.fa
        touch offsets.json
    else
        reduce_sequences \
            -reduction ${msr} \
            -sequences $reference \
            -offsets offsets.json \
            -output reduced.fa \
            -threads 8
    fi
    """
}

simulators = Channel.from("pbsim", "nanosim")
reduced_reference_sequence
    .combine(simulators)
    .set{ reference_to_index }

process index_reference {
    input:
    tuple val(msr), file(reduced_reference), val(simulator) from reference_to_index

    output:
    tuple val(simulator), val(msr), file("reference.idx") into reference_index

    memory 30.GB
    cpus 16

    errorStrategy {task.exitStatus == 143 ? "retry": "ignore"}
    maxRetries 10

    script:
    k = (simulator == "pbsim") ? 19 : 15
    """
    echo "${simulator}, ${msr}"
    minimap2 -k ${k} -t 16 -d reference.idx ${reduced_reference}
    """
}

nanosim.into{ nanosim_count; nanosim_map}
pbsim.into{ pbsim_count; pbsim_map}

nanosim_map
    .concat(pbsim_map)
    .splitFasta(by: 2000, file: true, elem: 1)
    .set{ read_chunks }

msrs_reads
    .map{ it -> [it.baseName, it] }
    .combine(reduced_reference_offset, by: 0)
    .combine(read_chunks)
    .map{ it -> [it[3], it[0], it[1], it[2], it[4]]}
    .set{ reads_to_map }

reference_index
    .combine(reads_to_map, by:[0,1])
    .set{ mapper_input }

process map_reads {
    input:
    tuple val(simulator), val(msr), file(index), file(msr_file), file(offset), file(reads) from mapper_input

    output:
    tuple val(msr), val(simulator), file("chunk.paf.gz") into mappings

    memory 20.GB
    cpus 16
    time "20m"

    errorStrategy {task.exitStatus == 143 ? "retry": "ignore"}

    script:
    k = (simulator == "pbsim") ? 19 : 15
    """
    set -o pipefail
    echo "${msr}, ${simulator}"

    if [[ ${msr} != "raw" ]]; then
        # Reduce reads
        reduce_sequences \
            -reduction ${msr_file} \
            -sequences ${reads} \
            -output reduced.fa \
            -threads 16
        rename_sequences \
            -sequences reduced.fa \
            -offsetsPath ${offset} \
            -output renamed.fa
        rm reduced.fa
    else
        cp ${reads} renamed.fa
    fi

    minimap2 -k ${k} -t 16 -c ${index} -2 renamed.fa | pigz -9 -p 16 > chunk.paf.gz
    rm renamed.fa
    """
}

process eval_mapping {
    input:
    tuple val(msr), val(simulator), file("chunk") from mappings.groupTuple(by: [0, 1])

    publishDir "$resDir/${msr}", mode: "copy"
    output:
    file "${simulator}.paf.gz"
    file "${simulator}.acc" into mapevals

    script:
    """
    cat chunk* > mapping.paf.gz
    echo "${msr}, ${simulator}, minimap" > ${simulator}.acc
    paftoolsCustom.js mapeval mapping.paf.gz >> ${simulator}.acc
    mv mapping.paf.gz ${simulator}.paf.gz
    """
}

process select_MSRs {

    conda "pandas tqdm", useMamba: true

    input:
    file "eval" from mapevals
    file "nanosim.fa" from nanosim_count
    file "pbsim.fa" from pbsim_count

    publishDir "$resDir", mode: "copy"
    output:
    file "evaluation.csv"

    cpus 1
    memory 5.GB

    script:
    """
    N_COUNT=\$(fastatools count -i nanosim.fa)
    P_COUNT=\$fastatools count -i pbsim.fa)

    echo "{\"nanosim\":\$N_COUNT, \"pbsim\":\$P_COUNT}" > nreads.json

    mapevalGatherer.py --files eval* --nreads nreads.json --output evaluation.csv --unzipped
    """

}
