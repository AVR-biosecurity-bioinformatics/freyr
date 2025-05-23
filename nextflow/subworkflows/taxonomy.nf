/*
 *  Taxonomically assign ASVs
 */


//// modules to import
include { TAX_IDTAXA                                } from '../modules/tax_idtaxa'
include { TAX_BLAST                                 } from '../modules/tax_blast'
include { JOINT_TAX                                 } from '../modules/joint_tax'
include { MERGE_TAX                                 } from '../modules/merge_tax'
include { ASSIGNMENT_PLOT                           } from '../modules/assignment_plot'
include { TAX_SUMMARY                               } from '../modules/tax_summary'
include { TAX_SUMMARY_MERGE                         } from '../modules/tax_summary_merge'

workflow TAXONOMY {

    take:

    ch_seqtab
    ch_loci_params
    ch_idtaxa_db_new

    main:

    //// combine loci_params to seqtab
    ch_seqtab
        .combine ( ch_loci_params, by: 0 )
        .map { pcr_primers, fcid, seqtab, fasta, loci_params -> 
            [ pcr_primers, fcid, loci_params, fasta ] }
        .set { ch_seqs_params }

    //// use newly trained IDTAXA model, if it exists
    if ( params.train_idtaxa ) {
        ch_seqs_params
            .combine ( ch_idtaxa_db_new, by: 0 ) // join by pcr_primers
            .map { pcr_primers, fcid, loci_params, fasta, new_idtaxa ->
                [ pcr_primers, fcid, loci_params + [ idtaxa_db: new_idtaxa ], fasta ] }
            .set { ch_seqs_params }
    }

    //// split .fasta into chunks for taxonomic assignment
    ch_seqs_params
        .splitFasta (
            by: params.chunk_taxassign,
            file: true,
            elem: 3
        )
        .set { ch_taxassign_input }

    //// use IDTAXA to assign taxonomy
    TAX_IDTAXA ( 
        ch_taxassign_input
    )

    //// group IDTAXA assignment outputs by pcr_primers, fcid and loci_params
    TAX_IDTAXA.out.tax
        .groupTuple( by: [0,1,2] )
        .set { ch_idtaxa_tax }

    TAX_IDTAXA.out.ids
        .groupTuple( by: [0,1,2] )
        .set { ch_idtaxa_ids }

    //// use blastn to assign taxonomy
    TAX_BLAST ( 
        ch_taxassign_input
    )

    //// group BLAST assignment outputs by pcr_primers, fcid and loci_params
    TAX_BLAST.out.blast
        .groupTuple( by: [0,1,2] )
        .set { ch_blast_tax }

    //// group BLAST "low stringency" assignment outputa by pcr_primers and fcid
    TAX_BLAST.out.blast_assignment
        .groupTuple (by: [0,1] )
        .set { ch_blast_low }

    //// merge tax assignment outputs and filtered seqtab (pre-assignment)
    ch_idtaxa_tax // pcr_primers, fcid, loci_params, tax
        .join ( ch_blast_tax, by: [0,1,2] ) 
        .set { ch_joint_tax_input } // pcr_primers, fcid, loci_params, tax, blast

    //// aggregate taxonomic assignment
    JOINT_TAX ( ch_joint_tax_input )

    //// group taxtab across flowcells (per locus)
    JOINT_TAX.out.joint
        .map { pcr_primers, fcid, tax_tibble -> [ pcr_primers, tax_tibble ] }
        .groupTuple ( by: 0 ) // group into tuples using pcr_primers
        .set { ch_mergetax_input }

    //// merge tax tables across flowcells
    MERGE_TAX ( ch_mergetax_input )

    ch_mergetax_output = MERGE_TAX.out.merged_tax

    //// create assignment_plot input merging filtered fasta, taxtab, and blast output
    /// channel has one item per fcid x pcr_primer combo
    ch_seqs_params // pcr_primers, fcid, loci_params, fasta
        .join ( ch_blast_low, by: [0,1] ) // combine by pcr_primers, fcid 
        .join ( JOINT_TAX.out.joint, by: [0,1] ) // combine by pcr_primers, fcid
        .set { ch_assignment_plot_input }

    //// do assignment plot
    ASSIGNMENT_PLOT ( 
        ch_assignment_plot_input 
    )

    /// generate taxonomic assignment summary per locus
    ch_seqs_params // pcr_primers, fcid, loci_params, fasta
        .join ( ch_idtaxa_tax, by: [0,1,2] ) // + tax_csv (list)
        .join ( ch_idtaxa_ids, by: [0,1,2] ) // + "*_idtaxa_ids.rds" (list)
        .join ( ASSIGNMENT_PLOT.out.joint, by: [0,1,2] ) // + "*_joint.rds"
        .set { ch_tax_summary_input }

    //// create taxonomic assignment summaries per locus x flowcell
    TAX_SUMMARY ( 
        ch_tax_summary_input 
    )

    // create channel containing a single list of all TAX_SUMMARY outputs
    TAX_SUMMARY.out.csv
        .map { pcr_primers, fcid, loci_params, tax_summary ->
            [ tax_summary ] } 
        .collect()
        .set { ch_tax_summaries } 

    //// merge TAX_SUMMARY outputs together across loci and flow cells
    TAX_SUMMARY_MERGE ( 
        ch_tax_summaries 
    )


    emit:

    ch_tax_summaries
    ch_mergetax_output

}