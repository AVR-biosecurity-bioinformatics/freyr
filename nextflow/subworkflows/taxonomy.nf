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

    ch_seqtab_meta
    ch_loci_params


    main:

    //// use IDTAXA to assign taxonomy
    TAX_IDTAXA ( ch_seqtab_meta )
    
    ch_tax_idtaxa_tax = 
        TAX_IDTAXA.out.tax
        .map { pcr_primers, fcid, meta, tax -> // remove meta
            [ pcr_primers, fcid, tax ] }

    ch_tax_idtaxa_ids = 
        TAX_IDTAXA.out.ids 
        .map { pcr_primers, fcid, meta, ids -> // remove meta
            [ pcr_primers, fcid, ids ] }

    ch_seqtab = 
        ch_seqtab_meta
        .map { pcr_primers, fcid, meta, seqtab -> // remove meta
            [ pcr_primers, fcid, seqtab ] }

    //// use blastn to assign taxonomy
    TAX_BLAST ( ch_seqtab_meta )

    ch_tax_blast = 
        TAX_BLAST.out.blast
        .map { pcr_primers, fcid, meta, blast -> // remove meta
            [ pcr_primers, fcid, blast ] }

    //// merge tax assignment outputs and filtered seqtab (pre-assignment)
    ch_tax_idtaxa_tax
        .combine ( ch_tax_blast, by: [0,1] ) 
        .combine ( ch_seqtab, by: [0,1] )
        .combine ( ch_loci_params, by: 0 ) // adds map of loci_params
        .set { ch_joint_tax_input } // pcr_primers, fcid, tax, blast, seqtab, loci_params

    // ch_joint_tax_input.view()

    //// aggregate taxonomic assignment
    JOINT_TAX ( ch_joint_tax_input )

    //// group taxtab across flowcells (per locus)
    JOINT_TAX.out.taxtab
        .groupTuple ( by: 0 ) // group into tuples using pcr_primers
        .set { ch_mergetax_input }

    //// merge tax tables across flowcells
    MERGE_TAX ( ch_mergetax_input )

    //// create assignment_plot input merging filtered seqtab, taxtab, and blast output
    /// channel has one item per fcid x pcr_primer combo
    ch_seqtab
        .combine ( TAX_BLAST.out.blast_assignment, by: [0,1] ) // combine by pcr_primers, fcid 
        .combine ( JOINT_TAX.out.taxtab, by: [0,1] ) // combine by pcr_primers, fcid
        .combine ( ch_loci_params, by: 0 ) // add loci info (TODO: use ch_loci_params [stripped down] instead)
        .set { ch_assignment_plot_input }

    //// do assignment plot
    ASSIGNMENT_PLOT ( ch_assignment_plot_input )

    /// generate taxonomic assignment summary per locus (also hash seq)
    ch_tax_idtaxa_tax // pcr_primers, fcid, "*_idtaxa_tax.rds"
        .combine ( ch_tax_idtaxa_ids, by: [0,1] ) // + "*_idtaxa_ids.rds"
        .combine ( ASSIGNMENT_PLOT.out.joint, by: [0,1] ) // + "*_joint.rds", map(loci_params)
        .set { ch_tax_summary_input }

    //// create taxonomic assignment summaries per locus x flowcell
    TAX_SUMMARY ( ch_tax_summary_input )

    // create channel containing a single list of all TAX_SUMMARY outputs
    TAX_SUMMARY.out.rds
        .map { pcr_primers, fcid, loci_params, tax_summary ->
            [ tax_summary ] } 
        .collect()
        .set { ch_tax_summaries } 

    //// merge TAX_SUMMARY outputs together across loci and flow cells
    TAX_SUMMARY_MERGE ( ch_tax_summaries )


    emit:

    ch_tax_summaries

}