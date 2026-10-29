process EGGNOG_MAPPER {
    tag "$sample_id"
    publishDir "${params.outdir}/eggnog_mapper", mode: 'copy'

    input:
    tuple val(sample_id), path(protein_fasta)

    output:
    tuple val(sample_id), path("${sample_id}.emapper.annotations"), emit: annotations

    script:
    """
    set -euo pipefail

    echo "Running emapper.py for ${sample_id} using ${task.cpus} cores..."

    emapper.py -i ${protein_fasta} \
                 --itype proteins \
                 -o ${sample_id} \
                 --cpu ${task.cpus} \
                 --override \
                 --sensmode fast \
                 --data_dir /databases/eggnog
    """
}
