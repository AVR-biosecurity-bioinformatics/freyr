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

    main:

    //// create empty channels
    ch_processed_reads_forward = 
        Channel.empty()
    ch_processed_reads_reverse = 
        Channel.empty()
    ch_processed_reads_single = 
        Channel.empty()
    ch_read_tracker_grouped = Channel.empty()

    //// populate read channels with branched processed_reads input channel
    processed_reads
        .branch { direction, pcr_primers, fcid, meta, reads ->
            forward: direction == "forward"
            reverse: direction == "reverse"
            single: direction == "single"
            }
        .set { ch_processed_reads }

    ch_processed_reads_forward
        .concat ( ch_processed_reads.forward )
        .set { ch_processed_reads_forward }
    
    ch_processed_reads_reverse
        .concat ( ch_processed_reads.reverse )
        .set { ch_processed_reads_reverse }

    ch_processed_reads_single
        .concat ( ch_processed_reads.single )
        .set { ch_processed_reads_single }

    //// infer ASVs depending on type of sequencing
    if ( params.paired == true && params.seq_type == "illumina" ) {

        //// check both forward and reverse read channels are non-empty
        if ( !ch_processed_reads_forward || !ch_processed_reads_reverse ) {
            error ( "Data are asserted to be paired, but forward and/or reverse read channel is empty" )
        }

        //// group reads into list for error modeling
        ch_processed_reads_forward
            .map { direction, pcr_primers, fcid, meta, reads ->
            [ direction, pcr_primers, fcid, reads ]  }
            .groupTuple( by: [0,1,2] )
            .set { ch_error_input_forward }

        ch_processed_reads_reverse
            .map { direction, pcr_primers, fcid, meta, reads ->
            [ direction, pcr_primers, fcid, reads ]  }
            .groupTuple( by: [0,1,2] )
            .set { ch_error_input_reverse }
        
        //// error model on forward reads
        ERROR_MODEL_F ( ch_error_input_forward )

        //// error model on reverse reads
        ERROR_MODEL_R ( ch_error_input_reverse )

        //// input channel for denoising forward reads
        ch_processed_reads_forward
            .combine ( ERROR_MODEL_F.out.errormodel, by: [0,1,2] ) // combine with error model
            .map { direction, pcr_primers, fcid, meta, reads, errormodel -> // add empty file path for priors
                    [ direction, pcr_primers, fcid, meta, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
            .set { ch_denoise_input_forward }

        //// input channel for denoising reverse reads
        ch_processed_reads_reverse
            .combine ( ERROR_MODEL_R.out.errormodel, by: [0,1,2] ) // combine with error model
            .map { direction, pcr_primers, fcid, meta, reads, errormodel -> // add empty file path for priors
                    [ direction, pcr_primers, fcid, meta, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
            .set { ch_denoise_input_reverse }


        //// first pass of denoising per flowcell, primer and sample
        DENOISE1_F ( ch_denoise_input_forward, "first" )
        
        //// denoise forward reads per flowcell, primer and sample
        DENOISE1_R ( ch_denoise_input_reverse, "first" )

        // high sensitivity mode condition
        if ( params.high_sensitivity ) { // run prior extraction and second pass denoising
            //// group priors for each read file
            /// forward reads
            DENOISE1_F.out.seq
                .map { direction, pcr_primers, fcid, meta, reads, priors ->
                        [ direction, pcr_primers, fcid, priors ] }
                .groupTuple ( by: [0,1,2] )
                .set { ch_priors_f_pre }

            /// reverse reads
            DENOISE1_R.out.seq
                .map { direction, pcr_primers, fcid, meta, reads, priors ->
                        [ direction, pcr_primers, fcid, priors ] }
                .groupTuple ( by: [0,1,2] )
                .set { ch_priors_r_pre }

            //// get priors for forward reads
            PRIORS_F ( ch_priors_f_pre )
            
            /// combine with forward read data channel
            ch_denoise_input_forward
                .map { direction, pcr_primers, fcid, meta, reads, errormodel, priors -> // remove null priors
                        [ direction, pcr_primers, fcid, meta, reads, errormodel ] }
                .combine ( PRIORS_F.out.priors, by: [0,1,2] )
                .set { ch_denoise2_input_forward }

            //// get priors for reverse reads
            PRIORS_R ( ch_priors_r_pre )

            /// combine with reverse read data channel
            ch_denoise_input_reverse
                .map { direction, pcr_primers, fcid, meta, reads, errormodel, priors -> // remove null priors
                        [ direction, pcr_primers, fcid, meta, reads, errormodel ] }
                .combine ( PRIORS_R.out.priors, by: [0,1,2] )
                .set { ch_denoise2_input_reverse }

            //// run pseudo-pooled 2nd denoising with priors on forward reads
            DENOISE2_F ( ch_denoise2_input_forward, "second" )

            //// run pseudo-pooled 2nd denoising with priors on reverse reads
            DENOISE2_R ( ch_denoise2_input_reverse, "second" )

            /// join F and R denoise2 outputs
            // prepare forward reads
            DENOISE2_F.out.seq
                .map { direction, pcr_primers, fcid, meta, readsF, seqF ->
                        [ meta.sample_id, pcr_primers, fcid, meta, readsF, seqF ] }
                .set { ch_seq_forward }

            // prepare reverse reads
            DENOISE2_R.out.seq
                .map { direction, pcr_primers, fcid, meta, readsR, seqR ->
                        [ meta.sample_id, pcr_primers, fcid, meta, readsR, seqR ] }
                .set { ch_seq_reverse }

            // join
            ch_seq_forward
                .combine ( ch_seq_reverse, by: [0,1,2,3] ) // combine by sample_id
                .map { sample_id, pcr_primers, fcid, meta, readsF, seqF, readsR, seqR -> // remove sample_id and meta
                        [ pcr_primers, fcid, meta.concat_unmerged, meta,
                        file(readsF, checkIfExists: true),
                        file(readsR, checkIfExists: true), 
                        file(seqF, checkIfExists: true),
                        file(seqR, checkIfExists: true) ] } 
                .groupTuple ( by: [0,1,2] ) // assumes concat_unmerged is the same for all samples, which it should be
                .set { ch_seq_combined }

        } else { // don't run second denoising step with priors
            /// join F and R DENOISE1 outputs
            // prepare forward reads
            DENOISE1_F.out.seq
                .map { direction, pcr_primers, fcid, meta, readsF, seqF ->
                        [ meta.sample_id, pcr_primers, fcid, meta, readsF, seqF ] }
                .set { ch_seq_forward }

            // prepare reverse reads
            DENOISE1_R.out.seq
                .map { direction, pcr_primers, fcid, meta, readsR, seqR ->
                        [ meta.sample_id, pcr_primers, fcid, meta, readsR, seqR ] }
                .set { ch_seq_reverse }

            // join
            ch_seq_forward
                .combine ( ch_seq_reverse, by: [0,1,2,3] ) // combine by sample_id
                .map { sample_id, pcr_primers, fcid, meta, readsF, seqF, readsR, seqR -> // remove sample_id and meta
                        [ pcr_primers, fcid, meta.concat_unmerged, meta,
                        file(readsF, checkIfExists: true),
                        file(readsR, checkIfExists: true), 
                        file(seqF, checkIfExists: true),
                        file(seqR, checkIfExists: true) ] } 
                .groupTuple ( by: [0,1,2] ) // assumes concat_unmerged is the same for all samples
                .set { ch_seq_combined }
        }

        //// merge paired-end reads per flowcell x locus combo
        MAKE_SEQTAB_PAIRED ( 
            ch_seq_combined 
        )

        ch_read_tracker_grouped = 
            ch_read_tracker_grouped.concat(MAKE_SEQTAB_PAIRED.out.read_tracking)

        ch_seqtab = MAKE_SEQTAB_PAIRED.out.seqtab


    } else if ( params.paired == false && params.seq_type == "nanopore" )  {

        //// check both forward and reverse read channels are non-empty
        if ( !ch_processed_reads_single ) {
            error ( "Data are asserted to be single-end, but single read channel is empty" )
        }

        //// group reads into list for error modeling
        ch_processed_reads_single
            .map { direction, pcr_primers, fcid, meta, reads ->
            [ direction, pcr_primers, fcid, reads ]  }
            .groupTuple( by: [0,1,2] )
            .set { ch_error_input_single }
        
        //// error model on single reads
        ERROR_MODEL_S ( ch_error_input_single )

        //// input channel for denoising single reads
        ch_processed_reads_single
            .combine ( ERROR_MODEL_S.out.errormodel, by: [0,1,2] ) // combine with error model
            .map { direction, pcr_primers, fcid, meta, reads, errormodel -> // add empty file path for priors
                    [ direction, pcr_primers, fcid, meta, reads, errormodel, "$projectDir/assets/NO_FILE" ] }
            .set { ch_denoise_input_single }

        //// first pass of denoising per flowcell, primer and sample
        DENOISE1_S ( ch_denoise_input_single, "first" )

        // high sensitivity mode condition
        if ( params.high_sensitivity ) { // run prior extraction and second pass denoising
            //// group priors for each read file
            /// single reads
            DENOISE1_S.out.seq
                .map { direction, pcr_primers, fcid, meta, reads, priors ->
                        [ direction, pcr_primers, fcid, priors ] }
                .groupTuple ( by: [0,1,2] )
                .set { ch_priors_s_pre }

            //// get priors for single reads
            PRIORS_S ( ch_priors_s_pre )
            
            /// combine with single read data channel
            ch_denoise_input_single
                .map { direction, pcr_primers, fcid, meta, reads, errormodel, priors -> // remove null priors
                        [ direction, pcr_primers, fcid, meta, reads, errormodel ] }
                .combine ( PRIORS_S.out.priors, by: [0,1,2] )
                .set { ch_denoise2_input_single }

            //// run pseudo-pooled 2nd denoising with priors on forward reads
            DENOISE2_S ( ch_denoise2_input_single, "second" )

            /// create sequence table from dada object
            // prepare single reads
            DENOISE2_S.out.seq
                .map { direction, pcr_primers, fcid, meta, reads, seqs ->
                        [ pcr_primers, fcid, meta.concat_unmerged, meta, file(reads, checkIfExists: true), file(seqs, checkIfExists: true) ] }
                .groupTuple ( by: [0,1,2] ) // assumes concat_unmerged is the same for all samples
                .set { ch_seq_single }

        } else { // don't run second denoising step with priors
            DENOISE1_S.out.seq
                .map { direction, pcr_primers, fcid, meta, reads, seqs ->
                        [ pcr_primers, fcid, meta.concat_unmerged, meta, file(reads, checkIfExists: true), file(seqs, checkIfExists: true) ] }
                .groupTuple ( by: [0,1,2] ) // assumes concat_unmerged is the same for all samples
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