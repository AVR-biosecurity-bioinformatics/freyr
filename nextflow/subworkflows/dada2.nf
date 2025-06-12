/*
 *  Infer ASVs from processed reads
 */


//// modules to import
include { ERROR_MODEL as ERROR_MODEL_F              } from '../modules/error_model'
include { ERROR_MODEL as ERROR_MODEL_R              } from '../modules/error_model'
include { ERROR_MODEL as ERROR_MODEL_S              } from '../modules/error_model'
include { DENOISE as DENOISE1_F                     } from '../modules/denoise'
include { DENOISE as DENOISE1_R                     } from '../modules/denoise'
include { DENOISE as DENOISE1_S                     } from '../modules/denoise'
include { PRIORS as PRIORS_F                        } from '../modules/priors'
include { PRIORS as PRIORS_R                        } from '../modules/priors'
include { PRIORS as PRIORS_S                        } from '../modules/priors'
include { DENOISE as DENOISE2_F                     } from '../modules/denoise'
include { DENOISE as DENOISE2_R                     } from '../modules/denoise'
include { DENOISE as DENOISE2_S                     } from '../modules/denoise'
include { MAKE_SEQTAB_PAIRED                        } from '../modules/make_seqtab_paired'
include { MAKE_SEQTAB_SINGLE                        } from '../modules/make_seqtab_single'

workflow DADA2 {

    take:
    processed_reads
    ch_primer_params

    main:

    //// create empty channels
    ch_read_tracker_grouped = Channel.empty()

    //// populate read channels with branched processed_reads input channel
    processed_reads
        .branch { direction, primers, read_group, sample, sample_primers, reads ->
            forward: direction == "forward"
            reverse: direction == "reverse"
            single: direction == "single"
            }
        .set { ch_processed_reads }

    //// get primer_params 'concat_unmerged' parameter
    ch_primer_params
        .map { primers, primer_params ->  [ primers, primer_params.concat_unmerged ] }
        .set { ch_concat_unmerged } 

    //// infer ASVs depending on type of sequencing
    if ( params.paired == true && params.seq_type == "illumina" ) {

        //// check both forward and reverse read channels are non-empty
        if ( !ch_processed_reads.forward || !ch_processed_reads.reverse ) {
            error ( "Data are asserted to be paired, but forward and/or reverse read channel is empty" )
        }

        //// group reads into list for error modeling
        ch_processed_reads.forward
            .map { direction, primers, read_group, sample, sample_primers, reads ->
                [ direction, primers, read_group, reads ] }
            .groupTuple( by: [0,1,2] )
            .set { ch_error_input_forward }

        ch_processed_reads.reverse
            .map { direction, primers, read_group, sample, sample_primers, reads ->
                [ direction, primers, read_group, reads ]  }
            .groupTuple( by: [0,1,2] )
            .set { ch_error_input_reverse }
        
        //// error model on forward reads
        ERROR_MODEL_F ( 
            ch_error_input_forward
        )

        //// error model on reverse reads
        ERROR_MODEL_R ( 
            ch_error_input_reverse
        )

        //// input channel for denoising forward reads
        ch_processed_reads.forward
            .combine ( ERROR_MODEL_F.out.errormodel, by: [0,1,2] ) // combine with error model
            .map { direction, primers, read_group, sample, sample_primers, reads, errormodel -> // add empty file path for priors
                [ direction, primers, read_group, sample, sample_primers, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
            .set { ch_denoise1_input_forward }

        //// input channel for denoising reverse reads
        ch_processed_reads.reverse
            .combine ( ERROR_MODEL_R.out.errormodel, by: [0,1,2] ) // combine with error model
            .map { direction, primers, read_group, sample, sample_primers, reads, errormodel -> // add empty file path for priors
                [ direction, primers, read_group, sample, sample_primers, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
            .set { ch_denoise1_input_reverse }

        //// first pass of denoising per flowcell, primer and sample
        DENOISE1_F (
            ch_denoise1_input_forward, 
            "first" 
        )
        
        //// denoise forward reads per flowcell, primer and sample
        DENOISE1_R ( 
            ch_denoise1_input_reverse, 
            "first" 
        )

        // high sensitivity mode condition
        if ( params.high_sensitivity ) { // run prior extraction and second pass denoising
            //// group priors for each read file
            
            /// forward reads
            DENOISE1_F.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, priors ->
                    [ direction, primers, read_group, priors ] }
                .groupTuple ( by: [0,1,2] )
                .set { ch_priors_f_pre }

            //// get priors for forward reads
            PRIORS_F ( ch_priors_f_pre )

            /// combine with forward read data channel
            ch_denoise1_input_forward
                .map { direction, primers, read_group, sample, sample_primers, reads, errormodel, priors -> // remove null priors
                    [ direction, primers, read_group, sample, sample_primers, reads, errormodel ] }
                .combine ( PRIORS_F.out.priors, by: [0,1,2] )
                .set { ch_denoise2_input_forward }

            //// run pseudo-pooled 2nd denoising with priors on forward reads
            DENOISE2_F ( 
                ch_denoise2_input_forward, 
                "second" 
            )

            /// remove direction
            DENOISE2_F.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, seqs ->
                    [ primers, read_group, sample, sample_primers, reads, seqs ] }
                .set { ch_denoised_f }


            /// reverse reads
            DENOISE1_R.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, priors ->
                    [ direction, primers, read_group, priors ] }
                .groupTuple ( by: [0,1,2] )
                .set { ch_priors_r_pre }

            //// get priors for reverse reads
            PRIORS_R ( ch_priors_r_pre )

            /// combine with reverse read data channel
            ch_denoise1_input_reverse
                .map { direction, primers, read_group, sample, sample_primers, reads, errormodel, priors -> // remove null priors
                    [ direction, primers, read_group, sample, sample_primers, reads, errormodel ] }
                .combine ( PRIORS_R.out.priors, by: [0,1,2] )
                .set { ch_denoise2_input_reverse }

            //// run pseudo-pooled 2nd denoising with priors on reverse reads
            DENOISE2_R ( 
                ch_denoise2_input_reverse,
                "second"
            )

            /// remove direction
            DENOISE2_R.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, seqs ->
                    [ primers, read_group, sample, sample_primers, reads, seqs ] }
                .set { ch_denoised_r }  


        } else { // don't run second denoising step with priors
            /// join F and R DENOISE1 outputs
            // prepare forward reads
            DENOISE1_F.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, seqs ->
                    [ primers, read_group, sample, sample_primers, reads, seqs ] }
                .set { ch_denoised_f }

            // prepare reverse reads
            DENOISE1_R.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, seqs ->
                    [ primers, read_group, sample, sample_primers, reads, seqs ] }
                .set { ch_denoised_r }

        }

        //// join F and R denoised outputs
        ch_denoised_f
            .join ( ch_denoised_r, by: [0,1,2,3] ) // join by primers, read_group, sample, sample_primers
            .map { primers, read_group, sample, sample_primers, readsF, seqF, readsR, seqR -> 
                    [ primers, read_group, sample, sample_primers, readsF, readsR, seqF, seqR ] } 
            .groupTuple ( by: [0,1] ) // group by primers, read_group
            .combine ( ch_concat_unmerged, by: 0 )
            .set { ch_seq_combined }

        //// merge paired-end reads per flowcell x locus combo
        MAKE_SEQTAB_PAIRED ( 
            ch_seq_combined 
        )

        ch_read_tracker_grouped = 
            ch_read_tracker_grouped.concat(MAKE_SEQTAB_PAIRED.out.read_tracking)

        ch_seqtab = MAKE_SEQTAB_PAIRED.out.seqtab


    } else if ( params.paired == false && params.seq_type == "nanopore" )  {

        //// check both forward and reverse read channels are non-empty
        if ( !ch_processed_reads.single ) {
            error ( "Data are asserted to be single-end, but single read channel is empty" )
        }

        //// group reads into list for error modeling
        ch_processed_reads.single
            .map { direction, primers, read_group, sample, sample_primers, reads ->
                [ direction, primers, read_group, reads ]  }
            .groupTuple( by: [0,1,2] )
            .set { ch_error_input_single }
        
        //// error model on single reads
        ERROR_MODEL_S ( 
            ch_error_input_single 
        )

        //// input channel for denoising single reads
        ch_processed_reads.single
            .combine ( ERROR_MODEL_S.out.errormodel, by: [0,1,2] ) // combine with error model
            .map { direction, primers, read_group, sample, sample_primers, reads, errormodel -> // add empty file path for priors
                    [ direction, primers, read_group, sample, sample_primers, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
            .set { ch_denoise1_input_single }

        //// first pass of denoising per flowcell, primer and sample
        DENOISE1_S ( 
            ch_denoise1_input_single, 
            "first" 
        )

        // high sensitivity mode condition
        if ( params.high_sensitivity ) { // run prior extraction and second pass denoising
            //// group priors for each read file
            /// single reads
            DENOISE1_S.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, priors ->
                        [ direction, primers, read_group, priors ] }
                .groupTuple ( by: [0,1,2] )
                .set { ch_priors_s_pre }

            //// get priors for single reads
            PRIORS_S ( 
                ch_priors_s_pre 
            )
            
            /// combine with single read data channel
            ch_denoise1_input_single
                .map { direction, primers, read_group, sample, sample_primers, reads, errormodel, priors -> // remove null priors
                        [ direction, primers, read_group, sample, sample_primers, reads, errormodel ] }
                .combine ( PRIORS_S.out.priors, by: [0,1,2] )
                .set { ch_denoise2_input_single }

            //// run pseudo-pooled 2nd denoising with priors on forward reads
            DENOISE2_S ( 
                ch_denoise2_input_single, 
                "second" 
            )

            // prepare single reads
            DENOISE2_S.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, seqs ->
                        [ primers, read_group, sample, sample_primers, reads, seqs ] }
                .groupTuple ( by: [0,1] ) // group by primers, read_group
                .combine ( ch_concat_unmerged, by: 0 )
                .set { ch_seq_single }

        } else { // don't run second denoising step with priors
            DENOISE1_S.out.seq
                .map { direction, primers, read_group, sample, sample_primers, reads, seqs ->
                        [ primers, read_group, sample, sample_primers, reads, seqs ] }
                .groupTuple ( by: [0,1] ) // group by primers, read_group
                .combine ( ch_concat_unmerged, by: 0 )
                .set { ch_seq_single }
        }

        //// merge paired-end reads per flowcell x locus combo
        MAKE_SEQTAB_SINGLE ( 
            ch_seq_single 
        )

        ch_read_tracker_grouped = 
            ch_read_tracker_grouped.concat(MAKE_SEQTAB_SINGLE.out.read_tracking)

        ch_seqtab = MAKE_SEQTAB_SINGLE.out.seqtab


    } else {
        error ("Disallowed 'params.paired' and 'params.seq_type' combination")
    }


    emit:

    ch_seqtab
    ch_read_tracker_grouped


}