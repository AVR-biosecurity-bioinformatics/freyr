/*
 *  Create result summaries, using phyloseq
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
include { FILTER_SEQTAB                             } from '../modules/filter_seqtab'

workflow RESULT_SUMMARIES {




}