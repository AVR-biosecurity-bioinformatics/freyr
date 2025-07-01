#!/usr/bin/env Rscript
### load only required packages
process_packages <- c(
    "Biostrings",
    "dplyr",
    "phyloseq",
    "readr",
    "seqateurs",
    "stringr",
    "tibble",
    "data.table",
    "tidyr",
    NULL
)
invisible(lapply(head(process_packages,-1), library, character.only = TRUE, warn.conflicts = FALSE))

### check Nextflow environment variables
nf_vars <- c(
    "projectDir",
    "primers",
    "seqtab_file",
    "taxtab_file",
    "samdf_file",
    "raw_file",
    "summary_file",
    "ps_file",
    "clusters_file",
    "merge_clusters"
)
lapply(nf_vars, nf_var_check)

## check and define variables
ps <- readRDS(ps_file)

seqtab_tibble <- readr::read_csv(seqtab_file)

taxtab_tibble <- readr::read_csv(taxtab_file)

samdf_tibble <- readr::read_csv(samdf_file)

raw_tibble <- readr::read_csv(raw_file)

summary_tibble <- readr::read_csv(summary_file)

clusters_tibble <- readr::read_csv(clusters_file)

if (!merge_clusters %in% c("central","lca","frequency","rank","abundance")){
    stop(paste0("'--merge_clusters' value must be one of: '",stringr::str_flatten(c("central","lca","frequency","rank","abundance"), " / "),"'!"))
}

### run R code


# check clusters are defined (ie. integers not NA)
if ( !any(clusters_tibble$cluster %>% is.na) ){
      
    summary_clusters <- 
        summary_tibble %>%
        # remove cluster
        dplyr::select(-cluster) %>%
        # add marked cluster
        dplyr::left_join(., clusters_tibble, by = "seq_name") %>%
        # mark central sequences
        dplyr::mutate(
            central = dplyr::if_else(cluster < 0, TRUE, FALSE),
            cluster = abs(cluster)
        ) %>%
        dplyr::relocate(cluster, central, .after = sequence) %>%
        # convert to long format 
        tidyr::pivot_longer(
            cols = !c(seq_name, sequence, cluster, central, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species),
            names_to = "sample_primers",
            values_to = "abundance"
        ) %>%
        # sum abundance within cluster for each sample
        dplyr::group_by(cluster, sample_primers) %>%
        dplyr::mutate(abundance = sum(abundance)) %>%
        dplyr::ungroup()

    # summary without abundance
    seq_clusters <- 
        summary_clusters %>%
        dplyr::select(-c(sample_primers, abundance)) %>%
        dplyr::distinct()

    # keep central sequence, with taxonomy determined by parameter
    # merge_clusters <- "central" # taxonomic assignment of central sequence
    # merge_clusters <- "lca" # lowest common ancestor
    # merge_clusters <- "frequency" # most common assignment, with tie-breaking by assignment with "lowest" hash name
    # merge_clusters <- "rank" # lowest rank assignment, with tie-breaking by frequency then by hash name
    # merge_clusters <- "abundance" # choose taxonomy of sequence with the highest read count across all samples

    if ( merge_clusters == "central" ){
    
        ## keep taxonomic assignment of central sequence
        summary_reps <- 
            summary_clusters %>%
            dplyr::filter(central == TRUE)

    } else if ( merge_clusters == "lca" ){

        ## assign taxonomy using the lowest common rank of all cluster members
        summary_reps <- 
            seq_clusters %>%
            dplyr::group_by(cluster) %>%
            # convert rank to NA if more than one value is present within cluster
            dplyr::mutate(
                dplyr::across(
                    Root:Species,
                    ~ ifelse(length(unique(.)) > 1, NA, .)
                )
            ) %>% 
            dplyr::ungroup() %>%
            # convert NA to double-underscore format
            dplyr::mutate(
                # unsure Root is always "Root"
                Root = "Root",
                # determine the lowest assigned rank, choosing the first true case
                lowest_assignment = dplyr::case_when(
                    is.na(Kingdom)  ~ paste0("R__",Root),
                    is.na(Phylum)   ~ paste0("K__",Kingdom),
                    is.na(Class)    ~ paste0("P__",Phylum),
                    is.na(Order)    ~ paste0("C__",Class),
                    is.na(Family)   ~ paste0("O__",Order),
                    is.na(Genus)    ~ paste0("F__",Family),
                    is.na(Species)  ~ paste0("G__",Genus),
                    .default        = NA
                ),
                # replace NA ranks with the "G__Genus" format ('propagate')
                dplyr::across(
                    Kingdom:Species, 
                    ~dplyr::case_when(
                        is.na(.) & !is.na(lowest_assignment) ~ lowest_assignment, 
                        .default = .
                    )
                )
            ) %>%
            dplyr::select(-lowest_assignment) %>% 
            dplyr::filter(central == TRUE) %>%
            # join back to long abundance data
            dplyr::left_join(
                ., 
                summary_clusters %>% dplyr::select(seq_name, cluster, sample_primers, abundance), 
                by = c("seq_name", "cluster")
            )

    } else if ( merge_clusters == "frequency" ) {

        # assign taxonomy using most common taxonomic assignment, breaking ties with sequence hash 
        summary_reps <-  
            seq_clusters %>%
            # replace rank columns with lineage string column
            tidyr::unite(Root:Species, col = "lineage_string", sep = ";", remove = T) %>%
            # count lineage strings within each cluster
            dplyr::group_by(cluster, lineage_string) %>%
            dplyr::mutate(count = n()) %>%
            # replace lineage string of cluster with the lineage string of the highest count + lowest hash
            dplyr::group_by(cluster) %>%
            dplyr::arrange(desc(count), seq_name) %>% 
            dplyr::mutate(lineage_string = first(lineage_string)) %>% 
            dplyr::ungroup() %>%
            # filter to central only
            dplyr::select(-count) %>%
            dplyr::filter(central == TRUE) %>%
            # regenerate rank columns from lineage string
            tidyr::separate(
                col = lineage_string, 
                into = c("Root", "Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                sep = ";", 
                extra = "merge",
                remove = T
            ) %>% 
            # join back to long abundance data
            dplyr::left_join(
                ., 
                summary_clusters %>% dplyr::select(seq_name, cluster, sample_primers, abundance), 
                by = c("seq_name", "cluster")
            )
    
    } else if ( merge_clusters == "rank" ){
    
        # assign taxonomy using the lowest assignment level (eg. species over genus), breaking ties with frequency then sequence hash
        summary_reps <- 
            seq_clusters %>%
            # get lowest assigned rank
            dplyr::mutate(
                lowest_assignment = dplyr::case_when(
                    stringr::str_detect(Kingdom, "^R__")  ~ "8_Root",
                    stringr::str_detect(Phylum, "^K__")  ~ "7_Kingdom",
                    stringr::str_detect(Class, "^P__")  ~ "6_Phylum",
                    stringr::str_detect(Order, "^C__")  ~ "5_Class",
                    stringr::str_detect(Family, "^O__")  ~ "4_Order",
                    stringr::str_detect(Genus, "^F__")  ~ "3_Family",
                    stringr::str_detect(Species, "/")  ~ "2_Genus",
                    .default  = "1_Species"
                )
            ) %>% 
            # replace rank columns with lineage string column
            tidyr::unite(Root:Species, col = "lineage_string", sep = ";", remove = T) %>%
            # count lineage strings within each cluster
            dplyr::group_by(cluster, lineage_string) %>%
            dplyr::mutate(count = n()) %>%
            # replace lineage string of cluster with the lineage string of the lowest assignment, then highest count, then lowest hash
            dplyr::group_by(cluster) %>%
            dplyr::arrange(lowest_assignment, desc(count), seq_name) %>% 
            dplyr::mutate(lineage_string = first(lineage_string)) %>% 
            dplyr::ungroup() %>%
            # filter to central only
            dplyr::select(-count, -lowest_assignment) %>%
            dplyr::filter(central == TRUE) %>%
            # regenerate rank columns from lineage string
            tidyr::separate(
                col = lineage_string, 
                into = c("Root", "Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                sep = ";", 
                extra = "merge",
                remove = T
            ) %>%
            # join back to long abundance data
            dplyr::left_join(
                ., 
                summary_clusters %>% dplyr::select(seq_name, cluster, sample_primers, abundance), 
                by = c("seq_name", "cluster")
            )

    } else if ( merge_clusters == "abundance" ) {
    
        # assign taxonomy of sequence with the highest read count across all samples
        # get total abundance per sequence
        seqs_totalab <- 
            summary_tibble %>%
            # pivot longer
            tidyr::pivot_longer(
                cols = !c(seq_name, sequence, cluster, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species),
                names_to = "sample_primers",
                values_to = "abundance"
            ) %>%
            dplyr::group_by(seq_name) %>%
            dplyr::mutate(abundance = sum(abundance)) %>% 
            dplyr::ungroup() %>%
            dplyr::select(seq_name, cluster, abundance) %>%
            dplyr::distinct() 
        
        summary_reps <- 
            seq_clusters %>%
            # add total abundance to clusters with taxonomy
            dplyr::left_join(., seqs_totalab, by = c("seq_name", "cluster")) %>%
            # replace rank columns with lineage string column
            tidyr::unite(Root:Species, col = "lineage_string", sep = ";", remove = T) %>%
            # arrange by highest count within cluster (and hash of sequence to break ties) and replace lineage string
            dplyr::group_by(cluster) %>%
            dplyr::arrange(desc(abundance), seq_name) %>%
            dplyr::mutate(lineage_string = first(lineage_string)) %>%
            dplyr::ungroup() %>%
            # filter to central only
            dplyr::select(-abundance) %>%
            dplyr::filter(central == TRUE) %>%
            # regenerate rank columns from lineage string
            tidyr::separate(
                col = lineage_string, 
                into = c("Root", "Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                sep = ";", 
                extra = "merge",
                remove = T
            ) %>%
            # join back to long abundance data
            dplyr::left_join(
                ., 
                summary_clusters %>% dplyr::select(seq_name, cluster, sample_primers, abundance), 
                by = c("seq_name", "cluster")
            ) 

    } else {
        stop(paste0("Invalid value of '--merge_clusters': '",merge_clusters,"'"))
    }

    # regenerate summary tibble
    summary_clustered <- 
        summary_reps %>%
        dplyr::select(-central) %>%
        tidyr::pivot_wider(
            names_from = sample_primers,
            values_from = abundance,
            values_fill = 0
        )
    
    # make clustered seqtab
    seqtab_clustered <- 
        summary_clustered %>%
        dplyr::select(-c(sequence, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species))

    # make clustered taxtab
    taxtab_clustered <- 
        summary_clustered %>%
        dplyr::select(seq_name, cluster, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species)

    # make clustered raw
    raw_meta <- 
        raw_tibble %>%
        dplyr::select(-c(seq_name, Abundance, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species, sequence, cluster)) %>%
        dplyr::distinct()

    raw_clustered <- 
        summary_clustered %>%
        tidyr::pivot_longer(
            cols = !c(seq_name, sequence, cluster, Root, Kingdom, Phylum, Class, Order, Family, Genus, Species),
            names_to = "sample_primers",
            values_to = "Abundance"
        ) %>%
        dplyr::filter(Abundance > 0) %>%
        dplyr::left_join(., raw_meta, by = "sample_primers") %>%
        dplyr::relocate(seq_name, sample_primers, Abundance, sample, read_group, primers, any_of(c("fwd", "rev", "single"))) %>%
        dplyr::relocate(Root, Kingdom, Phylum, Class, Order, Family, Genus, Species, sequence, cluster, .after = last_col())

    # create clustered phyloseq object
    otutab_clustered_ps <- 
        seqtab_clustered %>%
        dplyr::select(-cluster) %>%
        tibble::column_to_rownames(var = "seq_name") %>%
        as.matrix() %>%
        phyloseq::otu_table(taxa_are_rows = TRUE)

    taxtab_clustered_ps <- 
        taxtab_clustered %>%
        dplyr::select(-cluster) %>%
        tibble::column_to_rownames(var = "seq_name") %>%
        as.matrix() %>%
        phyloseq::tax_table()

    ps_clustered <- 
        phyloseq::phyloseq(
            otutab_clustered_ps,
            taxtab_clustered_ps,
            phyloseq::sample_data(ps),
            phyloseq::refseq(ps)[names(phyloseq::refseq(ps)) %in% seqtab_clustered$seq_name]
        )

    # save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
    seqs_output <- phyloseq::refseq(ps_clustered)

    seq_names_new <- 
        phyloseq::tax_table(ps_clustered) %>%
        as("matrix") %>%
        tibble::as_tibble(rownames = "seq_name") %>%
        # ensure sequence name order is same as the DSS object
        dplyr::arrange(factor(seq_name, levels = names(seqs_output))) %>%
        # add primers
        dplyr::mutate(primers = primers, .after = seq_name) %>%
        # unite columns into a single header string per sequence
        tidyr::unite(col = "lineage", Root:Species, sep = ";") %>%
        tidyr::unite(col = "id", c(seq_name, primers), sep = "|") %>%
        tidyr::unite(col = "header", c(id, lineage), sep = ";") %>%
        dplyr::pull(header)

    names(seqs_output) <- seq_names_new

    write_fasta(seqs_output, paste0("asvs_clustered_", primers, ".fasta"))  

    # Write out
    readr::write_csv(raw_clustered, paste0("raw_clustered_",primers,".csv"))
    readr::write_csv(summary_clustered, paste0("summary_clustered_",primers,".csv"))
    readr::write_csv(seqtab_clustered, paste0("seqtab_clustered_",primers,".csv"))
    readr::write_csv(taxtab_clustered, paste0("taxtab_clustered_",primers,".csv"))
    readr::write_csv(samdf_tibble, paste0("samdf_clustered",primers,".csv"))
    saveRDS(ps_clustered, paste0("ps_clustered_",primers,".rds"))
    
  
} else if (all(clusters_tibble$cluster %>% is.na)){
    
    ## save inputs as outputs

    # save sequences as .fasta file (with taxonomy in header, in format "seq_name|primers;Root;Kingdom;Phylum;Class;Order;Family;Genus;Species")
    seqs_output <- phyloseq::refseq(ps)

    seq_names_new <- 
        phyloseq::tax_table(ps) %>%
        as("matrix") %>%
        tibble::as_tibble(rownames = "seq_name") %>%
        # ensure sequence name order is same as the DSS object
        dplyr::arrange(factor(seq_name, levels = names(seqs_output))) %>%
        # add primers
        dplyr::mutate(primers = primers, .after = seq_name) %>%
        # unite columns into a single header string per sequence
        tidyr::unite(col = "lineage", Root:Species, sep = ";") %>%
        tidyr::unite(col = "id", c(seq_name, primers), sep = "|") %>%
        tidyr::unite(col = "header", c(id, lineage), sep = ";") %>%
        dplyr::pull(header)

    names(seqs_output) <- seq_names_new

    write_fasta(seqs_output, paste0("asvs_clustered_", primers, ".fasta"))  
    
    # write unmerged outputs
    readr::write_csv(raw_tibble, paste0("raw_clustered_",primers,".csv"))
    readr::write_csv(summary_tibble, paste0("summary_clustered_",primers,".csv"))
    readr::write_csv(seqtab_tibble, paste0("seqtab_clustered_",primers,".csv"))
    readr::write_csv(taxtab_tibble, paste0("taxtab_clustered_",primers,".csv"))
    readr::write_csv(samdf_tibble, paste0("samdf_clustered",primers,".csv"))
    saveRDS(ps, paste0("ps_clustered_",primers,".rds"))
    
} else {
    stop(paste0("Clusters are a mixture of 'NA' and other values!"))
}

# stop(" *** stopped manually *** ") ##########################################
