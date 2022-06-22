/////////////////////////////////////////////////////////////////////////////////////////////////
//////////                         GLOBAL VARIABLES + CONSTANTS                        //////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Constants
meanNanoSimLength = 8000
resDir = file("./results/MSR_eval/")
fastaSplitSize = 2000
mappers = Channel.from("minimap2", "winnowmap")
simulators = Channel.from("pbsim", "nanosim")

// Data
repeats = file("./repeats/chm13.repeats.bigBed")
MSRs = Channel.fromPath("./data/SSRs/MSRs/*2out*.json")

// Models
nanosimModel = file("./data/models/human_NA12878_DNA_FAB49712_guppy_flipflop/training")
pbsimModel = file("./data/models/P6C4.model")

// Genomes
// [organism, length, coverage, FASTA]
genomes = Channel.from(
    ["human", 3500000000, 1.5, file("./data/whole_human_genome.fa")],
    ["centromere", 1014904, 100, file("./data/tandemtools.ref.fa")],
    ["drosophila", 137547960, 1.5, file("./data/whole_drosophila_genome.fa")],
    ["ecoli", 4639675, 50, file("./data/whole_ecoli_genome.fa")]
)

// Returns number of reads to simulate from coverage
// genome length and average read length
def getRN(len, cov, rlen) {
    return (
        ((cov * len) / rlen).intValue()
    )
}

genomes.into{ genome_pbsim; genome_nanosim_prepare; genome_MSRs_prepare }

genome_nanosim_prepare
    // [organism, readNumber, file]
    .map{[it[0], getRN(it[1], it[2], meanNanoSimLength), it[3]]}
    .set{genome_nanosim}

genome_MSRs_prepare
    // [organims, file]
    .map{[it[0], it[3]]}
    .set{genome_MSRs}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////                               READ SIMULATION                               //////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// Simulate reads with NanoSim
 process nanosim {
    conda '/pasteur/zeus/projets/p01/Evolbioinfo/users/lblassel/miniconda3/envs/nanosim'

    input:
    tuple val(organism), val(nreads), file(reference) from genome_nanosim

    publishDir "$resDir/${organism}/reads/", mode: "copy"
    output:
    tuple val("nanosim"), val(organism), file("nanosim.fa") into nanosim

    memory 20.GB
    cpus 16

    script:
    """
    nanosim.py genome -rg ${reference} -o simulated -n $nreads -t 16 -dna_type linear -c $nanosimModel
    convertNanoSimNames.py --fasta simulated_aligned_reads.fasta --errorProfile simulated_aligned_error_profile --output nanosim.fa
    """
 }

// Simulate reads with Pbsim
process pbsim {
    input:
    tuple val(organism), val(_), val(coverage), file(reference) from genome_pbsim

    publishDir "$resDir/${organism}/reads/", mode: "copy"
    output:
    tuple val("pbsim"), val(organism), file("pbsim.fa") into pbsim

     memory 20.GB
     cpus 16

    script:
    """
    pbsim --depth ${coverage} --hmm_model $pbsimModel ${reference}
    samtools faidx ${reference} -o ref.fai
    paftools.js pbsim2fq ref.fai *maf > pbsim.fa
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//////////                            REFERENCE MANIPULATION                           //////////
/////////////////////////////////////////////////////////////////////////////////////////////////

MSRs.into{ MSRs_ref; MSRs_reads }

// Reduce reference sequence
process reduceRef {
    input:
    tuple val(organism), file(reference), file(msr) from genome_MSRs.combine(MSRs_ref)

    publishDir "$resDir/${organism}/eval/${msr.baseName}", mode: "copy"
    output:
    tuple val(organism), val("${msr.baseName}"), file("ref.fa") into reduced_ref
    tuple val(organism), val("${msr.baseName}"), file("offsets.json") into offsets_ref

     memory 30.GB
     cpus 8

    script:
    """
    if [ ${msr.baseName} == "raw" ]; then
        cp ${reference} ref.fa
        touch offsets.raw.json
    else
        reduceDelete_linux_amd64 -reduction ${msr} -sequences ${reference} -output ref.fa -offsets offsets.json -threads 8
    fi
    """
}

 reduced_ref
    .combine(simulators) // red,ref_f,sim
    .into{ ref_index; ref_count; ref_winnowmap }


// Index reference sequence with minimap for given k value
process indexRef {
    input:
    tuple val(organism), val(msr), file(reference), val(simulator) from ref_index

    output:
    tuple val(organism), val(simulator), val(msr), file("ref.idx") into indexes

     memory 30.GB
     cpus 16

    script:
    k_i = (simulator == "pbsim") ? 19 : 15
    """
    minimap2 -k ${k_i} -t 16 -d ref.idx ${reference}
    """
}

// Count repetitive k-mers with meryl
process countRef {
    input:
    tuple val(organism), val(msr), file(reference), val(simulator) from ref_count

    output:
    tuple val(organism), val(simulator), val(msr), file("repetitive.kmers.txt") into counts

    memory 30.GB
    cpus 16

    script:
    k_c = (simulator == "pbsim") ? 19 : 15
    """
    meryl count k=${k_c} threads=16 memory=25 output merylDB ${reference}
    meryl print greater-than distinct=0.9998 merylDB > repetitive.kmers.txt
    rm -rf merylDB
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////                                CHANNEL MAGIC                                //////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// For minimap:
// FILES: reads, msr, offsets, ref_index
// VALUES: simulator, msr

// For winnowmap:
// FILES: reads, msr, offsets, ref_fa, repetitive_kmers
// VALUES: simulator, msr

// Common Vals / files
// FILES: reads, msr, offsets
// VALUES: simulator, msr

// Split Reads into chunks
 nanosim
     .concat(pbsim)
     .splitFasta(by: fastaSplitSize, file: true, elem: 2) //sim,org,readchunk_f
     .map{ it -> [it[0], it[2], it[1]] } //sim,readchunk_f,org
     .set{ read_chunks }

// Get msrs into shape
MSRs_reads
    .map{ it -> [it, it.baseName] } //msr,msr_f
    .combine(offsets_ref, by: 1)    //msr,msr_f,org,offset_f
    .combine(read_chunks, by: 2) //org,msr,msr_f,offset_f,sim,readchunk_f
    .map{ it -> [it[0], it[4], it[1], it[2], it[3], it[5]] } //org,sim,msr,msr_f,offset_f,readchunk_f
    .into{ reads_minimap; reads_winnowmap }

indexes //org,sim,msr,index_f
    .combine(reads_minimap, by:[0,1,2]) //org,sim,msr,index_f,msr_f,offset_f,readchunk_f
    .set{ minimap_input }

ref_winnowmap //org,msr,ref_f,sim
    .map{ it -> [it[0], it[3], it[1], it[2]]} //org,sim,msr,ref_f
    .combine(counts, by:[0,1,2]) //org,sim,msr,ref_f,counts_f
    .combine(reads_winnowmap, by:[0,1,2]) //org,sim,msr,msr_f,counts_f,msr_f,offset_f,readchunk_f
    .set{ winnowmap_input }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////                                 READ MAPPING                                //////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// map reduced reads to reduced ref with minimap2
process minimap {
    input:
    tuple val(organism), val(simulator), val(msr), file(index), file(msr_f), file(offsets), file(reads) from minimap_input

    output:
    tuple val(organism), val(msr), val(simulator), val("minimap"), file("chunk.paf.gz") into minimap

    memory 30.GB
    cpus 16
    time "20m"

    script:
    k = (simulator == "pbsim") ? 19 : 15
    """
    set -o pipefail

    if [ "${msr}" != "raw" ]; then
        # Reduce Reads
        reduceDelete_linux_amd64 -msr ${msr_f} -sequences ${reads} -output r.fa -threads 16

        # Rename Reads
        applyOffsets_linux_amd64 -sequences r.fa -offsetsPath ${offsets} -output renamed.fa
        rm r.fa
    else
        cp ${reads} renamed.fa
    fi

    # Map reads
    minimap2 -k ${k} -t 16 -c ${index} -2 renamed.fa | pigz -9 -p 16 > chunk.paf.gz
    rm renamed.fa
    """
}

// map reduced reads to reduced ref with winnowmap2
process winnowmap {
    input:
    tuple val(organism), val(simulator), val(msr), file(ref), file(counts), file(msr_f), file(offsets), file(reads) from winnowmap_input

    output:
    tuple val(organism), val(msr), val(simulator), val("winnowmap"), file("chunk.paf.gz") into winnowmap

    memory 20.GB
    cpus 16
    time "20m"

    script:
    k = (simulator == "pbsim") ? 19 : 15
    """
    set -o pipefail

    if [ "${msr}" != "raw" ]; then
        # Reduce Reads
        reduceDelete_linux_amd64 -msr ${msr_f} -sequences ${reads} -output r.fa -threads 16

        # Rename Reads
        applyOffsets_linux_amd64 -sequences r.fa -offsetsPath ${offsets} -output renamed.fa
        rm r.fa
    else
        cp ${reads} renamed.fa
    fi

    # Map reads
    winnowmap -t 16 -W ${counts} -c -k ${k} -2 ${ref} renamed.fa | pigz -9 -p 16 > chunk.paf.gz
    rm renamed.fa
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////
//////////                              MAPPING EVALUATION                             //////////
/////////////////////////////////////////////////////////////////////////////////////////////////

minimap //org,msr,sim,mapper,mapping_f
    .concat(winnowmap) //org,msr,sim,mapper,mapping_f
    .groupTuple(by:[0,1,2,3])
    .into{ mapeval_full; mapeval_repeat }

// Evaluate mappings
process mapeval {
    input:
    tuple val(organism), val(msr), val(simulator), val(mapper), file("chunk") from mapeval_full

    publishDir "$resDir/${organism}/eval/${msr}", mode: "copy"
    output:
    file "${simulator}.${mapper}.paf.gz"
    file "${simulator}.${mapper}.acc.gz"

    script:
    """
    cat chunk* > mapping.paf.gz

    echo "${organism}, ${msr}, ${simulator}, ${mapper}" > ${simulator}.${mapper}.acc
    paftoolsCustom.js mapeval mapping.paf.gz >> ${simulator}.${mapper}.acc
    gzip -9 ${simulator}.${mapper}.acc
    mv mapping.paf.gz ${simulator}.${mapper}.paf.gz
    """
}

mapeval_repeat //org,msr,sim,mapper,[chunks_f]
    .filter{ it[0] == "human" }
    .set{ repeats }

// Evaluate mappings of reads in repeated regions of human genome
process mapeval_repeats {
    input:
    tuple val(organism), val(msr), val(simulator), val(mapper), file("chunk") from repeats

    publishDir "$resDir/${organism}/eval/${msr}/", mode: "copy"
    output:
    file "${simulator}.${mapper}.repeats.acc.gz"

    script:
    """
    # convert paf to bed
    cat chunk* > mapping.paf.gz
    paftools.js splice2bed mapping.paf.gz > mapping.bed

    # get reads with repeated region overlaps > 50%
    bedtools intersect -wa -u -f 0.50 -a mapping.bed -b $repeats > intersection.bed
    cut -f 4 intersection.bed | uniq > read.names

    # subset paf to keep reads in repeated regions
    gunzip mapping.paf.gz
    awk 'BEGIN{
        while (getline < "'"read.names"'") {
            select[\$0] = 1
        }
        close("'"read.names"'")
    }{
        if (select[\$1] == 1) {
            print \$0
        }
    }' mapping.paf > selected.paf
    echo "${organism}, ${msr}, ${simulator}, ${mapper}" > ${simulator}.${mapper}.repeats.acc
    paftoolsCustom.js mapeval selected.paf >> ${simulator}.${mapper}.repeats.acc
    gzip -9 ${simulator}.${mapper}.repeats.acc

    # cleaning up
    rm mapping.bed mapping.paf intersection.bed read.names selected.paf
    """
}
