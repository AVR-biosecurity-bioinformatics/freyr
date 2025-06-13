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
    ch_primer_params

    main:

    // //// just keep fasta files for input
    // ch_seqtab
    //     .map { primers, read_group, seqtab, fasta -> 
    //         [ primers, read_group, fasta ] }
    //     .combine ( ch_primer_params, by: 0 )
    //     .map { primers, read_group, fasta, primer_params -> 
    //         [ primers, read_group, primer_params, fasta ] }
    //     .set { ch_seqs_params }

    //// just keep fasta files for input
    ch_seqtab
        .map { primers, read_group, seqtab, fasta -> 
            [ primers, read_group, fasta ] }
        .set { ch_fasta }

    //// split .fasta into chunks for taxonomic assignment
    ch_fasta
        .splitFasta (
            by: params.chunk_taxassign,
            file: true,
            elem: 2
        )
        .set { ch_taxassign_input }


    //// combine TAX_IDTAXA process params to input channel
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = primer_params.subMap('idtaxa_confidence', 'idtaxa_db') 
            [ primers, process_params ] }
        .set { ch_tax_idtaxa_params }
    
    ch_taxassign_input
        .combine ( ch_tax_idtaxa_params, by: 0 )
        .set { ch_tax_idtaxa_input }

    //// use IDTAXA to assign taxonomy
    TAX_IDTAXA ( 
        ch_tax_idtaxa_input
    )

    //// group IDTAXA assignment outputs by primers, read_group
    TAX_IDTAXA.out.tax
        .groupTuple( by: [0,1] )
        .set { ch_idtaxa_tax }

    TAX_IDTAXA.out.ids
        .groupTuple( by: [0,1] )
        .set { ch_idtaxa_ids }


    //// combine TAX_BLAST process params to input channel 
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = primer_params.subMap('ref_fasta', 'blast_min_identity', 'blast_min_coverage', 'run_blast') 
            [ primers, process_params ] }
        .set { ch_tax_blast_params }
    
    ch_taxassign_input
        .combine ( ch_tax_blast_params, by: 0 )
        .set { ch_tax_blast_input }

    //// use blastn to assign taxonomy
    TAX_BLAST ( 
        ch_tax_blast_input
    )

    //// group BLAST assignment outputs by primers, read_group
    TAX_BLAST.out.blast
        .groupTuple( by: [0,1] )
        .set { ch_blast_tax }

    //// group BLAST "low stringency" assignment outputa by primers
    TAX_BLAST.out.blast_assignment
        .groupTuple (by: [0,1] )
        .set { ch_blast_low }


    //// merge tax assignment outputs and filtered seqtab (pre-assignment)
    ch_idtaxa_tax // primers, read_group, tax
        .join ( ch_blast_tax, by: [0,1] ) 
        .set { ch_joint_tax_input } // primers, read_group, tax, blast

    //// aggregate taxonomic assignment
    JOINT_TAX ( 
        ch_joint_tax_input
    )

    //// group taxtab across flowcells (per locus)
    JOINT_TAX.out.joint
        .map { primers, read_group, tax_tibble -> [ primers, tax_tibble ] }
        .groupTuple ( by: 0 ) // group into tuples using primers
        .set { ch_mergetax_input }

    //// merge tax tables across flowcells
    MERGE_TAX ( 
        ch_mergetax_input
    )

    ch_mergetax_output = MERGE_TAX.out.merged_tax

    //// get ASSIGNMENT_PLOT process params 
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = primer_params.subMap('idtaxa_db','ref_fasta') 
            [ primers, process_params ] }
        .set { ch_assignment_plot_params }
    
    //// create assignment_plot input merging filtered fasta, taxtab, blast output and process_params
    /// channel has one item per read_group x pcr_primer combo
    ch_fasta // primers, read_group, fasta
        .join ( ch_blast_low, by: [0,1] ) // combine by primers, read_group 
        .join ( JOINT_TAX.out.joint, by: [0,1] ) // combine by primers, read_group
        .combine ( ch_assignment_plot_params, by: 0 )
        .set { ch_assignment_plot_input }

    //// do assignment plot
    ASSIGNMENT_PLOT ( 
        ch_assignment_plot_input 
    )


    //// get TAX_SUMMARY process params 
    ch_primer_params
        .map { primers, primer_params ->
            def process_params = primer_params.subMap('locus','idtaxa_db','ref_fasta') 
            [ primers, process_params ] }
        .set { ch_tax_summary_params }

    /// generate taxonomic assignment summary per locus
    ch_fasta // primers, read_group, fasta
        .join ( ch_idtaxa_tax, by: [0,1] ) // + tax_csv (list)
        .join ( ch_idtaxa_ids, by: [0,1] ) // + "*_idtaxa_ids.rds" (list)
        .join ( ASSIGNMENT_PLOT.out.joint, by: [0,1] ) // + "*_joint.rds"
        .combine ( ch_tax_summary_params, by: 0 )
        .set { ch_tax_summary_input }

    //// create taxonomic assignment summaries per locus x flowcell
    TAX_SUMMARY ( 
        ch_tax_summary_input 
    )

    // create channel containing a single list of all TAX_SUMMARY outputs
    TAX_SUMMARY.out.csv
        .map { primers, read_group, tax_summary ->
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