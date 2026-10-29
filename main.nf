nextflow.enable.dsl=2

include { FASTQC as FASTQC_RAW    } from './modules/fastqc'
include { FASTP                   } from './modules/fastp'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { SPADES                  } from './modules/spades'
include { JELLYFISH               } from './modules/jellyfish'
include { SUMMARIZE_JELLYFISH } from './modules/summarize_jellyfish'
include { TIARA                   } from './modules/tiara'
include { REPEATMASKER            } from './modules/repeatmasker'
include { AUGUSTUS                } from './modules/augustus'
include { GFF_TO_GBK              } from './modules/gff_to_gbk'
include { ANTISMASH               } from './modules/antismash'
include { BUSCO                   } from './modules/busco'
include { SUMMARIZE_BUSCO         } from './modules/summarize_busco'
include { PROTEIN_FASTA           } from './modules/protein_fasta'
include { EGGNOG_MAPPER           } from './modules/eggnog_mapper'
include { INTERPROSCAN            } from './modules/interproscan'
include { DEEPSIG                 } from './modules/deepsig'
include { EFFECTORP               } from './modules/effectorp'
include { ORTHOFINDER             } from './modules/orthofinder'
include { MAFFT_ALIGN             } from './modules/mafft_align'
include { TRIMAL                  } from './modules/trimal'
include { ALIGNMENT_CONCATENATION } from './modules/alignment_concatenation'
include { RAXML                   } from './modules/raxml'

workflow {
    input_ch = Channel
        .fromPath(params.samples_table)
        .splitCsv(header: true)
        .map { row -> 
            [ row.sample, row.type, row.reads_1 ? file(row.reads_1) : null, row.reads_2 ? file(row.reads_2) : null, row.fasta ? file(row.fasta) : null ]
        }

    fastq_raw_ch = input_ch.filter { it[1] == 'fastq' }.map { id, type, r1, r2, f -> [ id, [ r1, r2 ] ] }
    fastq_run_ch = input_ch.filter { it[1] == 'fastq' }.map { id, type, r1, r2, f -> [ id, r1, r2 ] }
    fasta_input_ch = input_ch.filter { it[1] == 'fasta' }.map { id, type, r1, r2, f -> [ id, f ] }

    FASTQC_RAW(fastq_raw_ch)
    
    FASTP(fastq_run_ch)
    
    FASTQC_TRIM(FASTP.out.reads.map { id, r1, r2 -> [ id, [ r1, r2 ] ] })
    
    SPADES(FASTP.out.reads)

    raw_genomes_ch = SPADES.out.contigs.mix(fasta_input_ch)
    
    TIARA(raw_genomes_ch)

    JELLYFISH(TIARA.out.cleaned_fasta)
    
    BUSCO(TIARA.out.cleaned_fasta)
    
    SUMMARIZE_BUSCO(BUSCO.out.report.collect())
    
    SUMMARIZE_JELLYFISH(
        JELLYFISH.out.histo
            .map { sample_id, histo -> histo }
            .collect()
    )

    REPEATMASKER(TIARA.out.cleaned_fasta)
    
    AUGUSTUS(REPEATMASKER.out.masked_fasta)
    
    GFF_TO_GBK(AUGUSTUS.out.gff_and_fasta)
    
    ANTISMASH(GFF_TO_GBK.out.gbk)

    PROTEIN_FASTA(AUGUSTUS.out.gff_and_fasta)

    EGGNOG_MAPPER(PROTEIN_FASTA.out.proteins)
    
    INTERPROSCAN(PROTEIN_FASTA.out.proteins)
    
    DEEPSIG(PROTEIN_FASTA.out.proteins)
    
    EFFECTORP(PROTEIN_FASTA.out.proteins)

    all_proteins_ch = PROTEIN_FASTA.out.proteins.map { id, prot -> prot }.collect()
    
    ORTHOFINDER(all_proteins_ch)
    
    single_copy_ch = ORTHOFINDER.out.results
        .map { results_dir -> file("$results_dir/Single_Copy_Orthologue_Sequences/*.fa") }
        .flatten()

    MAFFT_ALIGN(single_copy_ch)

    TRIMAL(MAFFT_ALIGN.out.aligned_fa)

    ALIGNMENT_CONCATENATION(TRIMAL.out.trimmed_fa.collect())

    RAXML(ALIGNMENT_CONCATENATION.out.alignment)
}
