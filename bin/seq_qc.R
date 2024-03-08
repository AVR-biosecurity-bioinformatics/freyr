#!/usr/bin/env Rscript

print(paste0(projectDir,"/",data_loc,"/",flowcell_id))

step_seq_qc(flowcell_id, write_all=T)

step_switching_calc(flowcell_id)

quit(status = 0)

# tar_target(seq_qc, {
#             process <- temp_samdf1_grouped %>%
#             dplyr::group_by(fcid) %>%
#             tidyr::nest() %>%
#             dplyr::mutate(seq_qc = purrr::map(fcid, step_seq_qc, quiet=FALSE, write_all=FALSE))
#             out <- paste0("output/logs/", unique(process$fcid),"/",unique(process$fcid),"_flowcell_qc.pdf")
#             if(is.na(process$seq_qc[[1]]$reads_pf)){
#             pdf(file=out, paper="A4")
#             plot.new()
#             text(x=.5, y=.5, "ERROR: InterOp folder or RunInfo.xml not present") 
#             try(dev.off(), silent=TRUE)
#             }
#             return(out)
#             }, pattern = map(temp_samdf1_grouped), format="file",  iteration = "vector"),

# tar_target(switching_qc,{
#             process <- temp_samdf1_grouped %>%
#             dplyr::group_by(fcid) %>%
#             tidyr::nest() %>%
#             dplyr::mutate(switching_qc = purrr::map(fcid, step_switching_calc, quiet=TRUE))
#             out <- paste0("output/logs/",unique(process$fcid),"/",unique(process$fcid),"_index_switching.pdf")
#             if(is.na(process$switching_qc[[1]]$switch_rate)){
#             pdf(file=out, paper="A4")
#             plot.new()
#             text(x=.5, y=.5, "ERROR: Undetermined reads file not found") 
#             try(dev.off(), silent=TRUE)
#             }
#             return(out)
#             },pattern = map(temp_samdf1_grouped), format="file", iteration = "vector"),
            