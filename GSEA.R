# GSEA Preranked [CCBR] (7ae8d90f-c3a2-4684-aa30-00751f0d2c56): v179
GSEAPreranked <- function(DEGAnalysis, msigdb_v6_2_with_orthologs) {
  
  # This function calculates pre-ranked GSEA for multiple contrasts
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(dplyr); library(fgsea); library(grid); library(gridExtra); library(gtable); library(patchwork); library(data.table)
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  ## datasets
  DEGAnalysis = read.csv("pptx1_deg_results.csv")[,-1]
  msigdb_v6_2_with_orthologs = read.table("/rstudio-files/ccbr-data/data/msigdb_v6_2_with_orthologs.txt", sep = "\t", header = T)
  
  # Primary inputs
  deg_table = DEGAnalysis
  pathways_database = msigdb_v6_2_with_orthologs
  
  # Basic parameters
  gene_names_column = "Gene"
  species = "Human"
  collections_to_include = c("H: hallmark gene sets")
  
  # GSEA parameters
  fdr_correction_mode = "within each collection"
  minimum_gene_set_size = 15
  maximum_gene_set_size = 500
  number_of_permutations = 5000
  random_seed = 246642
  gene_scores_column_s_suffix = "_tstat"
  gene_score_alternative = c()
  collapse_pathway_redundancy = FALSE
  
  # Advanced parameters
  contrasts_filter = 'none'
  contrasts = c()    
  sort_output_by = c("pval")
  sort_output_in_decreasing_order = FALSE
  number_of_processing_units = 0
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  #1# Begin run.gsea() function: will be applied to dplyr::group_by(contrast)
  run.gsea <-
    function(dx,
             mode,
             collections,
             db,
             minimum_size,
             maximum_size,
             number_perms,
             organism,
             Np,
             randomSeed) {
      # compute gsea stats
      
      ranked = dx$genescores
      names(ranked) = dx$gene_id
      db$inPathway = sapply(db$gene_symbol, function(x)
        paste(sort(x[x %in% names(ranked)]), collapse = ","))
      
      if (mode == "over all collections") {
        set.seed(randomSeed)
        gsea <-
          fgsea(
            pathways = collections,
            stats = ranked,
            minSize = minimum_size,
            maxSize = maximum_size,
            nperm = number_perms,
            nproc = Np
          )
        gsea$size_leadingEdge <- sapply(gsea$leadingEdge, length)
        gsea$fraction_leadingEdge <- gsea$size_leadingEdge / gsea$size
        gsea$leadingEdge <-
          sapply(gsea$leadingEdge, function(x)
            paste(x, collapse = ","))
        gsea <-
          dplyr::inner_join(
            gsea,
            select(db, collection, gene_set_name, inPathway),
            by = c("pathway" = "gene_set_name")
          ) %>% dplyr::select(collection, dplyr::everything())
        
      } else {
        included_collections <-
          setNames(unique(db$collection), unique(db$collection))
        gsea <- lapply(included_collections, function(x) {
          set.seed(randomSeed)
          gsea_collection = fgsea(
            pathways = collections[names(collections) %in% dplyr::filter(db, collection == x)$gene_set_name],
            stats = ranked,
            minSize = minimum_size,
            maxSize = maximum_size,
            nperm = number_perms,
            nproc = Np
          )
          gsea_collection$size_leadingEdge <-
            sapply(gsea_collection$leadingEdge, length)
          gsea_collection$fraction_leadingEdge <-
            gsea_collection$size_leadingEdge / gsea_collection$size
          gsea_collection$leadingEdge <-
            sapply(gsea_collection$leadingEdge, function(x)
              paste(x, collapse = ","))
          return(
            dplyr::inner_join(
              gsea_collection,
              select(db, collection, gene_set_name, inPathway) %>% filter(collection == x),
              by = c("pathway" = "gene_set_name")
            ) %>% dplyr::select(collection, dplyr::everything())
          )
        }) %>% dplyr::bind_rows()
      }
      gsea$species = organism
      return(gsea)
    } # run.gsea() function
  
  #2# Begin edit fgsea::collapsePathways() that is add nproc argument and set.seed() for fgsea runs
  collapsePathways <-
    function (fgseaRes,
              pathways,
              stats,
              pval.threshold = 0.05,
              nperm = 10 / pval.threshold,
              Nproc,
              gseaParam = 1,
              rSeed){
      universe <- names(stats)
      pathways <- pathways[fgseaRes$pathway]
      pathways <- lapply(pathways, intersect, universe)
      parentPathways <-
        setNames(rep(NA, length(pathways)), names(pathways))
      for (i in seq_along(pathways)) {
        p <- names(pathways)[i]
        if (!is.na(parentPathways[p])) {
          next
        }
        pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))),
                                   p)
        if (length(pathwaysToCheck) == 0) {
          break
        }
        minPval <- setNames(rep(1, length(pathwaysToCheck)),
                            pathwaysToCheck)
        u1 <- setdiff(universe, pathways[[p]])
        
        set.seed(rSeed)
        fgseaRes1 <- fgsea(
          pathways = pathways[pathwaysToCheck],
          stats = stats[u1],
          nperm = nperm,
          maxSize = length(u1) -
            1,
          nproc = Nproc,
          gseaParam = gseaParam
        )
        minPval[fgseaRes1$pathway] <- pmin(minPval[fgseaRes1$pathway],
                                           fgseaRes1$pval)
        u2 <- pathways[[p]]
        
        set.seed(rSeed)
        fgseaRes2 <- fgsea(
          pathways = pathways[pathwaysToCheck],
          stats = stats[u2],
          nperm = nperm,
          maxSize = length(u2) -
            1,
          nproc = Nproc,
          gseaParam = gseaParam
        )
        minPval[fgseaRes2$pathway] <- pmin(minPval[fgseaRes2$pathway],
                                           fgseaRes2$pval)
        parentPathways[names(which(minPval > pval.threshold))] <- p
      }
      return(list(
        mainPathways = names(which(is.na(parentPathways))),
        parentPathways = parentPathways
      ))
    } # End collapsePathways() edit
  
  #3# Begin collapse.gsea() function (from Matt Angel's code)
  collapse.gsea <- function(grp, dx, collections, Np, randomSeed) {
    # filter ranked variable
    temp = dx %>% filter(contrast %in% grp$contrast)
    ranked = temp$genescores
    names(ranked) = temp$gene_id
    
    # collapse function
    run.collapse <-
      function(cp.input,
               pvalue,
               collections,
               ranked,
               Nprocs,
               rS) {
        collapsedPathways <-
          collapsePathways(
            as.data.table(cp.input),
            collections,
            ranked,
            pval.threshold = pvalue,
            Nproc = Nprocs,
            rSeed = rS
          ) # requires the data.table library
        return(data.frame(pathway = collapsedPathways$mainPathways))
      }
    filter_gsea = grp %>% dplyr::filter(pval < 0.05) %>% dplyr::arrange(pval) %>% dplyr::group_by(collection)
    collapsedResults = dplyr::group_modify(
      filter_gsea ,
      ~ run.collapse(
        .,
        pvalue = 0.05,
        collections = collections,
        ranked = ranked,
        Nprocs = Np,
        rS = randomSeed
      )
    ) %>% dplyr::ungroup()
    out <-
      grp %>% dplyr::inner_join(collapsedResults,
                                by = c('pathway' = 'pathway', 'collection' = 'collection'))
    
    return(out)
    
  } # End collapse.gsea() function
  
  #4# Begin table.pvalue() pvalue cutoffs table
  table.pvalue <- function(gsea) {
    cuts <- c(-Inf, 1e-04, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 1)
    cutsLab <- paste("<", cuts[-1], sep = "")
    p = cumsum(table(
      cut(
        gsea$pval,
        breaks = cuts,
        labels = cutsLab,
        include.lowest = FALSE,
        right = TRUE
      )
    ))
    q = cumsum(table(
      cut(
        gsea$padj,
        breaks = cuts,
        labels = cutsLab,
        include.lowest = FALSE,
        right = TRUE
      )
    ))
    tab <- data.frame(cutsLab, p, q)
    colnames(tab) <- c("alpha", "p-value", "*adjusted\np-value")
    rownames(tab) <- NULL
    return(tab)
  } # End table.pvalue() function
  
  #5# Begin plot.table() pvalue cutoffs table
  plot.table <- function(dtab, score) {
    title <-
      textGrob(paste0(unique(dtab$contrast), score), gp = gpar(fontsize = 10))
    tab <- as.data.frame.matrix(dtab %>%
                                  dplyr::select(-contrast))
    table <-
      tableGrob(tab, theme = ttheme_default(
        core = list(fg_params = list(cex = 0.9)),
        colhead = list(fg_params = list(cex = 0.9, parse = FALSE)),
        rowhead = list(fg_params = list(cex = 0.6))
      ))
    table <-
      gtable_add_rows(table,
                      heights = grobHeight(title) + unit(2, "line"),
                      pos = 0)
    table <-
      gtable_add_grob(
        table,
        list(title),
        t = c(1),
        l = c(1),
        r = ncol(table)
      )
    wrap_elements(table)
  }  # End plot.table() function
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  ## INPUT HANDLING AND FILTERING ====
  library(dplyr)
  
  ## pathway collection
  pathways_database <- pathways_database %>%
    filter(species == !!species, collection %in% !!collections_to_include)
  
  db_unique <- pathways_database %>%
    select(gene_set_name) %>%
    distinct()
  
  db_selected <- pathways_database %>%
    select(collection, gene_set_name) %>%
    distinct()
  
  # determine if there are duplicated gene set names overall (i.e. more rows in db_selected than unique names)
  db_isDuplicated <- nrow(db_selected) > nrow(db_unique)
  
  if (db_isDuplicated) {
    # bring the selected pairs into memory (no-op if already local)
    db_selected_local <- collect(db_selected)
    
    # count gene_set_name duplicates within each collection
    within_collection <- db_selected_local %>%
      group_by(collection, gene_set_name) %>%
      filter(n() > 1) %>%             # rows where same gene_set_name appears >1 time in same collection
      ungroup() %>%
      distinct(collection, gene_set_name) %>%
      nrow()
    
    if (within_collection == 0) {
      stop(
        "ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter"
      )
      
    } else if (within_collection > 0) {
      # count gene_set_names that appear in the intersection across collections
      between_collection <-
        length(Reduce(intersect, split(db_selected_local$gene_set_name, db_selected_local$collection)))
      
      if (between_collection == 0) {
        stop(
          "ERROR: duplicated gene set names found in the 'Gene set database' due to not unique gene set names within a collection"
        )
        
      } else if (between_collection > 0) {
        # kept the original behavior (returned string). If you want to stop, replace with stop(...)
        "ERROR: duplicated gene set names found in the 'Gene set database' due to overlapping collections selected by the 'Collections to include' parameter and not unique gene set names within a collection"
      }
    }
    
  } else if (!db_isDuplicated) {
    # bring the filtered pathways_database into memory (no-op if already local)
    pathways_database <- collect(pathways_database)
  }
  pathways_database <-
    pathways_database %>% dplyr::group_by(collection, gene_set_name) %>% dplyr::summarize(gene_symbol = as.list(strsplit(paste0(
      unique(gene_symbol), collapse = " "
    ), " "))) %>% dplyr::ungroup()
  geneset_list = pathways_database$gene_symbol
  names(geneset_list) = pathways_database$gene_set_name
  
  ## ranking
  
  if (!is.null(gene_score_alternative)) {
    if (gsub("", "", gene_score_alternative) == "") {
      stop(
        "'ERROR: Gene score alternative' parameter is empty - remove the entry or specify it correctly"
      )
    }
    gene_scores_column_s_suffix <- gene_score_alternative
  }
  rank_columns = colnames(deg_table)[grepl(paste0("\\Q", gene_scores_column_s_suffix, "\\E$"),
                                           colnames(deg_table))]
  rank_contrasts = unlist(strsplit(rank_columns, gene_scores_column_s_suffix))
  
  if (contrasts_filter == "remove") {
    if (!is.null(contrasts)) {
      if (gsub("", "", contrasts) == "") {
        stop("'ERROR: Contrasts' parameter is empty - remove the entry or specify it correctly")
      }
      
      all_contrasts = rank_contrasts
      index = match(contrasts, rank_contrasts)
      rank_columns = rank_columns[-index]
      rank_contrasts = rank_contrasts[-index]
      removed = setdiff(all_contrasts, rank_contrasts)
      if (length(removed) < 1) {
        cat(
          sprintf(
            'WARNING:contrast(s) to remove (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n',
            paste(contrasts, collapse = ", "),
            paste(rank_contrasts, collapse = ", ")
          )
        )
        
      } else {
        cat(sprintf(
          "Removed contrast(s): %s\n",
          paste(removed, collapse = ", ")
        ))
        cat(sprintf(
          "Kept contrast(s): %s\n",
          paste(rank_contrasts, collapse = ", ")
        ))
      }
      
    } else if (is.null(contrasts)) {
      cat(
        sprintf(
          'WARNING:contrast(s) to remove (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n',
          paste(contrasts, collapse = ", "),
          paste(rank_contrasts, collapse = ", ")
        )
      )
    }
    
  } else if (contrasts_filter == "keep") {
    if (!is.null(contrasts)) {
      if (gsub("", "", contrasts) == "") {
        stop("'ERROR: Contrasts' parameter is empty - remove the entry or specify it correctly")
      }
      
      all_contrasts = rank_contrasts
      index = match(contrasts, rank_contrasts)
      rank_columns = rank_columns[index]
      rank_contrasts = rank_contrasts[index]
      removed = setdiff(all_contrasts, rank_contrasts)
      if (length(rank_contrasts) < 1) {
        cat(
          sprintf(
            'WARNING:contrast(s) to keep (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n',
            paste(contrasts, collapse = ", "),
            paste(rank_contrasts, collapse = ", ")
          )
        )
        
      } else {
        cat(sprintf(
          "Removed contrast(s): %s\n",
          paste(removed, collapse = ", ")
        ))
        cat(sprintf(
          "Kept contrast(s): %s\n",
          paste(rank_contrasts, collapse = ", ")
        ))
      }
      
    } else if (is.null(contrasts)) {
      cat(
        sprintf(
          'WARNING:contrast(s) to keep (%s) not found; filter not applied\nIdentified contrast(s) used: %s\n',
          paste(contrasts, collapse = ", "),
          paste(rank_contrasts, collapse = ", ")
        )
      )
    }
    
  } else if (contrasts_filter == "none") {
    if (!is.null(contrasts)) {
      cat(
        sprintf(
          'WARNING:contrast filter not specified correctly; filter not applied\nIdentified contrast(s) used: %s\n',
          paste(rank_contrasts, collapse = ", ")
        )
      )
    } else {
      cat(sprintf(
        'Filter contrast ("none"); Identified contrast(s) used: %s\n',
        paste(rank_contrasts, collapse = ", ")
      ))
    }
  }
  
  deg_table <-
    deg_table %>%
    dplyr::select(gene_names_column, rank_columns) %>% tidyr::pivot_longer(
      !gene_names_column,
      names_to = "contrast",
      values_to = "genescores",
      values_drop_na = TRUE
    ) %>% 
    dplyr::rename("gene_id" = gene_names_column) %>% dplyr::mutate(contrast = sub(gene_scores_column_s_suffix, "", contrast))
  duplicates <-
    dplyr::group_by(deg_table, contrast, gene_id) %>% dplyr::filter(dplyr::n() > 1)
  if (nrow(duplicates) > 0) {
    genescore_grouped <-
      deg_table %>% dplyr::group_by(contrast, gene_id) %>% 
      dplyr::summarize(genescores = mean(genescores, na.rm = TRUE))
    cat(
      sprintf(
        "WARNING: duplicated gene names found of %g gene(s), duplicated values of gene scores were averaged per gene",
        length(unique(duplicates$gene_id))
      )
    )
  } else {
    genescore_grouped <- dplyr::group_by(deg_table, contrast)
  }
  
  # ANALYSIS ====
  library(fgsea)
  ## GSEA
  gsea <-
    dplyr::group_modify(
      genescore_grouped,
      ~ run.gsea(
        .,
        db = pathways_database,
        collections = geneset_list,
        minimum_size = minimum_gene_set_size,
        maximum_size = maximum_gene_set_size,
        number_perms = number_of_permutations,
        mode = fdr_correction_mode,
        organism = species,
        Np = number_of_processing_units,
        randomSeed = random_seed
      )
    )
  # OUTPUT ====
  
  ## visualization
  library(patchwork)
  tab <-
    dplyr::group_modify(gsea, ~ table.pvalue(.x)) %>% dplyr::ungroup()
  ltab <-
    split(tab, tab$contrast) %>% lapply(function(x)
      plot.table(x, score = gene_scores_column_s_suffix)) %>% wrap_plots()
  print(
    ltab + plot_annotation(
      title = "Cumulative number of significant calls (GSEA)",
      subtitle = sprintf(
        "*p value adjusted %s by the method of Benjamini and Hochberg (1995)",
        fdr_correction_mode
      ),
      tag_levels = 'A',
      theme = theme(
        plot.title = element_text(
          size = 20,
          face = 'bold',
          hjust = 0.5,
          margin = margin(t = 0)
        ),
        plot.subtitle =  element_text(
          size = 11,
          face = 'italic',
          hjust = 0.5,
          margin = margin(t = 10, b = 20)
        )
      )
    )
  )
  ## logs
  cat("\nThe number of tested gene sets per each collection and contrast\n")
  N <- dplyr::count(gsea, collection)
  print(N, n = nrow(N))
  cat(
    sprintf(
      "\nCumulative number of significant calls\np-value adjusted for the false discovery rate %s by the method of Benjamini and Hochberg (1995)\n",
      fdr_correction_mode
    )
  )
  tab %>% print(n = nrow(tab))
  
  ## collapse redundant?
  
  if (collapse_pathway_redundancy == TRUE) {
    gsea_grouped <-
      gsea %>% dplyr::mutate(group_contrast = contrast) %>% 
      dplyr::group_by(group_contrast)
    gsea <-
      dplyr::group_modify(
        gsea_grouped,
        ~ collapse.gsea(
          .,
          dx = deg_table,
          collections = geneset_list,
          Np = number_of_processing_units,
          randomSeed = random_seed
        )
      ) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(-group_contrast) %>%
      dplyr::group_by(contrast)
    cat("\nThe number of non-redundant gene sets per each collection and contrast\n")
    N <- dplyr::count(gsea, collection)
    print(N, n = nrow(N))
  }
  
  ## return dataset
  
  if (sort_output_in_decreasing_order) {
    sort_output_by = sapply(sort_output_by, function(x)
      sprintf("desc(%s)", x))
  }
  
  gsea <-
    gsea %>% 
    dplyr::arrange_(.dots = sort_output_by) %>% 
    tibble::add_column(geneScore = gene_scores_column_s_suffix, .after = "contrast") %>% 
    tibble::add_column(fdr_correction_mode = fdr_correction_mode, .after = 'geneScore')
  
  #return(gsea)
  write.csv(gsea,"ccbr1321_gsea_preranked.csv")
  
}

# GSEA Filters [CCBR] (d3c1d012-4ecd-4ac3-9e85-fd68cda28ba0): v108
GSEA_Filtered <- function(GSEAPreranked, msigdb_v6_2_with_orthologs ) {
  
  ## This function filters GSEA Table
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(dplyr)
  library(ggplot2)
  library(plotly)
  library(RColorBrewer)
  library(tidyverse)    
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  # Primary inputs
  gsea_table <- read.csv("ccbr1321_gsea_preranked.csv")
  pathways_database <- msigdb_v6_2_with_orthologs
  
  # Basic parameters
  species <- "Human"
  p_value_filter = 'adjusted p-value'
  p_value_threshold = 0.05
  enrichment_score_filter = 'NES (Normalized Enrichment Score)'
  enrichment_score_threshold = 0
  enrichment_score_sign = "+/-"
  size_filter = "Pathway size"
  size_cutoff = 0
  top_rank_filter = "all"
  
  # Advanced parameters
  columns_to_sort_output_by = c()
  sort_output_in_decreasing_order = FALSE
  collections_to_include =  c()
  pathways_to_include = c("HALLMARK_MTORC1_SIGNALING","HALLMARK_GLYCOLYSIS")
  gene_filter_universe = "Leading Edge (LE)"
  genes_to_include = c()
  contrast_filter = "none"
  contrasts = c()
  
  #Visualization parameters
  bubble_color_variable = "collection" # collection
  bubble_color_opacity = 0.95
  bubble_maximal_size = 2
  x_axis_minimum = c()
  x_axis_maximum = c()
  y_axis_minimum = c()
  y_axis_maximum = c()
  
  # Legacy parameters
  tested_contrast = c()
  gene_score = c()
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  ## --------- ##
  ## Functions ##
  ## --------- #
  
  filter.message <-
    function(filter, warn = FALSE, condition, output) {
      n = length(unique(output$pathway))
      if (!warn) {
        if (n == 0) {
          cat(
            sprintf(
              "ERROR: Filter by %s (%s) returned %g unique pathway(s)\n",
              filter,
              condition,
              n
            )
          )
          stop("Filter condition error\n")
        } else {
          cat(
            sprintf(
              "OK: Filter by %s (%s) returned %g unique pathway(s)\n",
              filter,
              condition,
              n
            )
          )
        }
      } else {
        cat(
          sprintf(
            "WARNING: Filter by %s (%s) not specified correctly; this filter is not applied\n",
            filter,
            condition
          )
        )
      }
    }
  
  # adjust old GSEA table (v67 or lower)
  adjust.v67 <-
    function(input,
             gS,
             tC,
             sp,
             db,
             required_columns = c(
               "contrast",
               "geneScore",
               "fdr_correction_mode",
               "collection",
               "pathway",
               "pval",
               "padj",
               "ES",
               "NES",
               "nMoreExtreme",
               "size",
               "leadingEdge",
               "size_leadingEdge",
               "inPathway"
             )) {
      missing_columns = setdiff(required_columns, colnames(input))
      
      if (length(missing_columns) > 0) {
        cat("WARNING: Output from outdated 'Preranked GSEA [CCBR]' detected.\n")
        cat(
          sprintf(
            "\nWARNING: Missing columns are added (%s)\n\t 'inPathway' column includes all genes annotated to a gene set - run the latest released version of 'Preranked GSEA [CCBR]'\n\t to return only genes mapped in the dataset.\n",
            paste(missing_columns, collapse = ", ")
          )
        )
        
        if (is.null(gS)) {
          gS = "NULL"
          cat(
            "\nWARNING: 'Gene score' parameter not provided therefore 'geneScore' column is assigned the default NULL value\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it will fail due to this specification of the 'Gene score' parameter.\n"
          )
        } else {
          if (length(gS) > 1) {
            gS = gS[1]
            cat(
              sprintf(
                "\nWARNING: too many values for 'Gene score' provided; only the first one used, '%s'\n",
                gS
              )
            )
          }
          has.underscore = substring(gS, 1, 1) == "_"
          if (!has.underscore) {
            gS = paste0("_", gS[1])
            cat(
              sprintf(
                "\nWARNING: 'Gene score' should start with underscore to match column naming convention in a DEG table used for GSEA run (e.g. '_tstat');\n\t column 'geneScore' is assigned the provided value with '_' added ('%s');\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it may fail if this is not an adequate specification of the 'Gene score' parameter.\n",
                gS
              )
            )
          }
        }
        
        if (is.null(tC)) {
          tC = "NULL"
          cat(
            "\nWARNING: 'Tested contrast' parameter not provided therefore 'contrast' column is assigned the default NULL value;\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it will fail due to this specification of the 'Tested contrast' parameter.\n"
          )
        } else {
          if (length(tC) > 1) {
            tC = tC[1]
            cat(
              sprintf(
                "\nWARNING: too many values for 'Tested contrast' provided; only the first one used, '%s'\n",
                tC
              )
            )
          }
          is.contrast = any(grepl("-", tC))
          if (!is.contrast) {
            cat(
              sprintf(
                "\nWARNING: 'Tested contrast' should be specified exactly the same way as when Preranked GSEA was run (e.g. treated-control);\n\t column 'contrast' is assigned the provided '%s' value;\n\t if downstream 'GSEA Running Score Diagram & Leading Edge Heatmap [Bulk] [CCBR]' will be linked to this output,\n\t it may fail if this is not an adequate specification of the 'Tested contrast' parameter.\n",
                tC[1]
              )
            )
          }
        }
        
        input <-
          input %>% dplyr::mutate(
            contrast = tC,
            geneScore = gS,
            fdr_correction_mode = "over all collections",
            size_leadingEdge = sapply(strsplit(leadingEdge, ","), length)
          )
        input <-
          input[, match(required_columns[required_columns != 'inPathway'], colnames(input))]
        db <-
          db %>% SparkR::filter(db[["species"]] == sp) %>% SparkR::filter(SparkR::`%in%`(db[["collection"]] , unique(input$collection))) %>% SparkR::filter(SparkR::`%in%`(db[["gene_set_name"]], unique(input$pathway))) %>% SparkR::collect()
        db <-
          db %>% dplyr::group_by(collection, gene_set_name, species) %>% dplyr::summarize(inPathway = paste(unique(gene_symbol), collapse = ",")) %>% dplyr::ungroup() %>% data.frame()
        input <-
          input %>% dplyr::left_join(db,
                                     by = c("collection" = "collection", "pathway" = "gene_set_name"))
        input <-
          input[, na.omit(match(c(required_columns, "species"), colnames(input)))]
        return(input)
        
      } else {
        cat("Filtering Preranked GSEA table.\n")
        return(input)
      }
    }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # translate filters to column names
  filterInput_byScore = switch(
    enrichment_score_filter,
    'NES (Normalized Enrichment Score)' = 'NES',
    'ES (Enrichment Score)' = 'ES'
  )
  filterInput_byPvalue = switch(p_value_filter,
                                'p-value' = 'pval',
                                'adjusted p-value' = 'padj')
  
  ## adjust old Preranked GSEA output (release v67 or lower)
  gsea_table <-
    adjust.v67(
      input = gsea_table,
      gS = gene_score,
      tC = tested_contrast,
      sp = species,
      db = pathways_database
    )
  
  ## apply filters
  cat("\n\nFiltering steps\n\n")
  
  gsea_filtered <- gsea_table %>%
    dplyr::filter(get(filterInput_byPvalue) <= p_value_threshold)
  filter.message(filter = p_value_filter,
                 condition = p_value_threshold,
                 output = gsea_filtered)
  
  if (enrichment_score_sign == "+") {
    gsea_filtered <- gsea_filtered %>%
      dplyr::filter(get(filterInput_byScore) >= enrichment_score_threshold)
    filter.message(
      filter = "value of GSEA score",
      condition = sprintf("%s > %g", filterInput_byScore, enrichment_score_threshold),
      output = gsea_filtered
    )
    
  } else if (enrichment_score_sign == "-") {
    gsea_filtered <- gsea_filtered %>%
      dplyr::filter(get(filterInput_byScore) <= -enrichment_score_threshold)
    filter.message(
      filter = "value of GSEA score",
      condition = sprintf(
        "%s < %s%g",
        filterInput_byScore,
        ifelse(enrichment_score_threshold == 0, "", "-"),
        enrichment_score_threshold
      ),
      output = gsea_filtered
    )
    
  } else {
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(abs(get(filterInput_byScore)) >= enrichment_score_threshold)
    filter.message(
      filter = "value of GSEA score",
      condition = sprintf("|%s| > %g", filterInput_byScore, enrichment_score_threshold),
      output = gsea_filtered
    )
  }
  
  if (size_cutoff > 0) {
    if (size_filter == "Pathway size") {
      gsea_filtered <-
        gsea_filtered %>% dplyr::filter(size >= size_cutoff)
      filter.message(filter = size_filter,
                     condition = size_cutoff,
                     output = gsea_filtered)
      
    } else {
      gsea_filtered <-
        gsea_filtered %>% dplyr::filter(size_leadingEdge >= size_cutoff)
      filter.message(filter = size_filter,
                     condition = size_cutoff,
                     output = gsea_filtered)
    }
    
  } else {
    filter.message(
      filter = size_filter,
      condition = paste0(">", size_cutoff),
      output = gsea_filtered
    )
  }
  
  if (!is.null(collections_to_include)) {
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(collection %in% collections_to_include)
    filter.message(filter = "Collection filter ON",
                   condition = "Collections to include",
                   output = gsea_filtered)
    cat(sprintf("    Found collections: %s\n", paste(
      unique(gsea_filtered$collection), collapse = ", "
    )))
  }
  
  if (!is.null(pathways_to_include)) {
    gsea_filtered <-
      gsea_filtered %>% dplyr::filter(pathway %in% pathways_to_include)
    filter.message(filter = "Pathway filter ON",
                   condition = "Pathways to include",
                   output = gsea_filtered)
    cat(sprintf("    Found pathways: %s\n", paste(
      unique(gsea_filtered$pathway), collapse = ", "
    )))
  }
  
  if (!is.null(genes_to_include)) {
    if (gene_filter_universe == 'Leading Edge (LE)') {
      index_pathway = lapply(genes_to_include, function(x)
        grep(paste0("\\b\\Q", x, "\\E\\b"), gsea_filtered$leadingEdge))
      found = sapply(index_pathway, function(x) {
        length(x) > 0
      })
      index_pathway = Reduce(union, index_pathway)
      gsea_filtered <- gsea_filtered %>% dplyr::slice(index_pathway)
      
      filter.message(filter = "Gene filter ON",
                     condition = "Leading Edge",
                     output = gsea_filtered)
      cat(sprintf(
        "    Found genes: %s\n",
        paste(genes_to_include[which(found)], collapse = ", ")
      ))
      if (sum(!found) > 0) {
        cat(sprintf(
          "    Missing genes: %s\n",
          paste(genes_to_include[which(!found)], collapse = ", ")
        ))
      }
      
    } else if (gene_filter_universe == "Pathway") {
      index_pathway = lapply(genes_to_include, function(x)
        grep(paste0("\\b\\Q", x, "\\E\\b"), gsea_filtered$inPathway))
      found = sapply(index_pathway, function(x) {
        length(x) > 0
      })
      index_pathway = Reduce(union, index_pathway)
      gsea_filtered <- gsea_filtered %>% dplyr::slice(index_pathway)
      
      filter.message(filter = "Gene filter ON",
                     condition = "in Pathway",
                     output = gsea_filtered)
      cat(sprintf(
        "    Found genes: %s\n",
        paste(genes_to_include[which(found)], collapse = ", ")
      ))
      if (sum(!found) > 0) {
        cat(sprintf(
          "    Missing genes: %s\n",
          paste(genes_to_include[which(!found)], collapse = ", ")
        ))
      }
      
    }
  }
  
  if (contrast_filter == "remove") {
    if (!is.null(contrasts)) {
      all_contrasts = unique(gsea_filtered$contrast)
      gsea_filtered <-
        gsea_filtered %>% dplyr::filter(!contrast %in% contrasts)
      removed = setdiff(all_contrasts, unique(gsea_filtered$contrast))
      if (length(removed) < 1) {
        filter.message(
          warn = TRUE,
          filter = "Contrast filter",
          condition = sprintf("remove %s missing", paste(contrasts, collapse = ", ")),
          output = gsea_filtered
        )
      } else {
        filter.message(filter = "Contrast filter",
                       condition = contrast_filter,
                       output = gsea_filtered)
        cat(sprintf(
          "    Removed contrast(s): %s\n",
          paste(removed, collapse = ", ")
        ))
        cat(sprintf("    Keep contrast(s): %s\n", paste(
          unique(gsea_filtered$contrast), collapse = ", "
        )))
      }
      
    } else if (is.null(contrasts)) {
      filter.message(
        warn = TRUE,
        filter = "Contrast filter",
        condition = class(contrasts),
        output = gsea_filtered
      )
    }
    
  } else if (contrast_filter == "keep") {
    if (!is.null(contrasts)) {
      all_contrasts = unique(gsea_filtered$contrast)
      gsea_filtered <-
        gsea_filtered %>% dplyr::filter(contrast %in% contrasts)
      kept = intersect(all_contrasts, unique(gsea_filtered$contrast))
      removed = setdiff(all_contrasts, unique(gsea_filtered$contrast))
      if (length(kept) < 1) {
        filter.message(
          warn = TRUE,
          filter = "Contrast filter",
          condition = sprintf("keep %s missing", paste(contrasts, collapse = ", ")),
          output = gsea_filtered
        )
      } else {
        filter.message(filter = "Contrast filter",
                       condition = contrast_filter,
                       output = gsea_filtered)
        cat(sprintf(
          "    Removed contrast(s): %s\n",
          paste(removed, collapse = ", ")
        ))
        cat(sprintf("    Kept contrast(s): %s\n", paste(
          unique(gsea_filtered$contrast), collapse = ", "
        )))
      }
      
    } else if (is.null(contrasts)) {
      filter.message(
        warn = TRUE,
        filter = "Contrast filter",
        condition = class(contrasts),
        output = gsea_filtered
      )
    }
    
  } else if (contrast_filter == "none") {
    if (!is.null(contrasts)) {
      filter.message(
        warn = TRUE,
        filter = "Contrast filter",
        condition = sprintf("none; %s", paste(contrasts, collapse = ", ")),
        output = gsea_filtered
      )
    }
  }
  
  top_rank_filter = tolower(top_rank_filter)
  
  if (top_rank_filter == "all") {
    filterInput_keep = paste("Inf (ALL)")
    top_rank_filter = Inf
    
  } else {
    top_rank_filter = as.numeric(top_rank_filter)
    
    if (is.na(top_rank_filter)) {
      stop("ERROR in Top rank filter; enter ALL (case insensitive) or a numeric rank.\n")
    }
    else if (top_rank_filter <= 0) {
      top_rank_filter = 1
      cat("WARNING: 'Top rank filter' cannot be 0 or less; its value was changed to 1.\n")
    }
    filterInput_keep = top_rank_filter
  }
  
  gsea_filtered <-
    gsea_filtered %>% dplyr::group_by(contrast, collection) %>% dplyr::mutate(p_rank = rank(pval, ties.method =
                                                                                              "min")) %>% dplyr::filter(p_rank <= top_rank_filter) %>% dplyr::select(-p_rank)
  filter.message(
    filter = "Top significant filter",
    condition = sprintf(
      "up to p-value rank of %s per contrast and collection",
      filterInput_keep
    ),
    output = gsea_filtered
  )
  
  # OUTPUT ====
  
  ## sort output
  if (sort_output_in_decreasing_order) {
    columns_to_sort_output_by = sapply(columns_to_sort_output_by, function(x)
      sprintf("desc(%s)", x))
  }
  gsea_filtered <-
    gsea_filtered %>% dplyr::arrange_(.dots = columns_to_sort_output_by)
  
  ## do plot and return dataset
  
  if (nrow(gsea_filtered) == 0) {
    stop("ERROR: filtering returned 0 pathways")
    
  } else {
    cat("\n\nFiltered pathways\n")
    tab = table(gsea_filtered$collection, gsea_filtered$contrast) %>% addmargins(margin =
                                                                                   c(1, 2))
    print(tab)
    
    cat("\n\nGSEA statistics (filtered pathways)\n\n")
    print(tibble(gsea_filtered))
    
    df <-
      gsea_filtered %>% dplyr::select(-leadingEdge, -inPathway) %>% dplyr::mutate(textContrast =
                                                                                    sprintf(
                                                                                      "%s: NES = %g, P-value = %g",
                                                                                      contrast,
                                                                                      signif(NES, 2),
                                                                                      signif(pval, 1)
                                                                                    )) %>% dplyr::group_by(pathway, collection) %>%
      dplyr::summarize(
        mean_pval = mean(pval),
        n_contrast = length(contrast),
        Pathway_size = mean(size),
        mean_NES = mean(NES),
        individual_values = paste(textContrast, collapse = "\n")
      ) %>%
      dplyr::mutate(
        textPathway = sprintf(
          "%s<br>%s<br>mean NES = %g, mean P-value = %g<br>Pathway size = %g<br><br>N contrasts = %g\n%s",
          collection,
          pathway,
          signif(mean_NES, 2),
          signif(mean_pval, 1),
          Pathway_size,
          n_contrast,
          individual_values
        )
      ) %>%
      dplyr::select(
        collection,
        pathway,
        Pathway_size,
        n_contrast,
        mean_NES,
        mean_pval,
        individual_values,
        textPathway
      )
    
    find_sort = grepl(paste(columns_to_sort_output_by, collapse = "|"),
                      colnames(df))
    if (sum(find_sort) == 0) {
      sort_by = c("n_contrast", "collection", "mean_pval")
    } else {
      sort_by = colnames(df)[which(find_sort)]
    }
    
    if (sort_output_in_decreasing_order) {
      sort_by = sapply(sort_by, function(x)
        sprintf("desc(%s)", x))
    }
    df <- df %>% dplyr::arrange_(.dots = sort_by)
    
    #cat("\n\nGSEA contrast summary (filtered pathways)\n\n")
    #print(tibble(df %>% dplyr::select(-textPathway) %>% dplyr::rename("mean_size"="Pathway_size")))
    
    if (bubble_color_variable == "pathway size") {
      ggp <-
        plot_ly(
          data = df,
          x = ~ mean_NES,
          y = ~ -log10(mean_pval),
          color = ~ Pathway_size,
          size = ~ n_contrast,
          text = ~ textPathway,
          hoverinfo = "text",
          opacity = bubble_color_opacity,
          marker = list(sizeref = ~ 2 * max(n_contrast) / bubble_maximal_size ** 2)
        ) %>%
        
        plotly::layout(
          xaxis = list(
            title = "Normalized Enrichment Score (mean)",
            tickfont = list(size = 15),
            titlefont = list(size = 15),
            showgrid = TRUE
          ),
          yaxis = list(
            title = "-log10 P-value (mean)",
            tickfont = list(size = 15),
            titlefont = list(size = 15),
            showgrid = TRUE
          ),
          legend = list(itemsizing = 'constant'),
          title = list(
            text = "Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway",
            x = 0,
            font = list(size = 15)
          ),
          showlegend = TRUE
        )
      
    } else {
      ggp <-
        plot_ly(
          data = df,
          x = ~ mean_NES,
          y = ~ -log10(mean_pval),
          color = ~ collection,
          size = ~ n_contrast,
          text = ~ textPathway,
          hoverinfo = "text",
          opacity = bubble_color_opacity,
          marker = list(sizeref = ~ 2 * max(n_contrast) / bubble_maximal_size ** 2)
        ) %>%
        
        plotly::layout(
          xaxis = list(
            title = "Normalized Enrichment Score (mean)",
            tickfont = list(size = 15),
            titlefont = list(size = 15),
            showgrid = TRUE
          ),
          yaxis = list(
            title = "-log10 P-value (mean)",
            tickfont = list(size = 15),
            titlefont = list(size = 15),
            showgrid = TRUE
          ),
          legend = list(itemsizing = 'constant'),
          title = list(
            text = "Filtered pathways; bubble area proportional to the number of filtered contrasts per pathway",
            x = 0,
            font = list(size = 15)
          ),
          showlegend = TRUE
        )
    }
    # custom axis range?
    
    if (!(is.null(x_axis_minimum) & is.null(x_axis_maximum))) {
      if (is.null(x_axis_minimum))
        x_axis_minimum = floor(min(df$mean_NES))
      if (is.null(x_axis_maximum))
        x_axis_maximum = ceiling(max(df$mean_NES))
      ggp <-
        ggp %>% plotly::layout(xaxis = list(range = list(x_axis_minimum, x_axis_maximum)))
    }
    
    if (!(is.null(y_axis_minimum) & is.null(y_axis_maximum))) {
      if (is.null(y_axis_minimum))
        y_axis_minimum = floor(min(-log10(df$mean_pval)))
      if (is.null(y_axis_maximum))
        y_axis_maximum = ceiling(max(-log10(df$mean_pval)))
      ggp <-
        ggp %>% plotly::layout(yaxis = list(range = list(y_axis_minimum, y_axis_maximum)))
    }
    
    print(ggp)
    #return(gsea_filtered)
    write.csv(gsea_filtered, "ccbr1321_gsea_filtered.csv")
    
  }
  
}

# GSEA Visualization [CCBR] (20896e36-771a-419f-a987-f90726837351): v398
GSEA_Visualization <- function(DEGAnalysis, GSEA_Filtered, NormalizedCounts, Ccbr1321_metadata, msigdb_v6_2_with_orthologs) {
  
  graphics.off()
  
  ## LIBRARIES ====
  library(tibble); library(dplyr); library(tidyr)
  library(ggplot2); library(ComplexHeatmap); library(colorspace); library(RColorBrewer)
  library(fgsea); library(patchwork)
  
  ## PARAMETERS
  
  which_plot = "ES+RNK+LE"
  plot_limit = "20"
  topline = 'coordinate'
  contrasts = c("Y_Y-Y_H","Y_Y-hYP7_CAR_T","Y_Y-CD19")
  pathways = c("HALLMARK_MTORC1_SIGNALING","HALLMARK_GLYCOLYSIS")
  geneid = "Gene"
  geneid_gex = "Gene"
  
  species = "Human"
  
  genescore_df <- read.csv("ccbr1432_volc_tstat4gsea.csv")[,-1]
  colnames(genescore_df) <- gsub("\\.","-",colnames(genescore_df))
  
  gex_df <- read.csv("ccbr1321_normalized_counts.csv")[,-1]
  metadata_df <- read.csv("Ccbr1321 metadata.csv")
  gsea_df <- read.csv("ccbr1321_gsea_filtered.csv")
  gsea_df$contrast <- gsub("\\.","-",gsea_df$contrast)
  
  pathway_db <- msigdb_v6_2_with_orthologs
  
  gex_transformation = 'median centering'
  reference_phenotype = FALSE
  drop_ref = FALSE
  filter_samples = FALSE
  
  pdfType = "common PDF"
  pageWidth=11
  pageHeight=8
  
  geneScore = "geneScore"
  testedContrast = c()
  
  
  #..GSEA
  
  ## adjust old Preranked GSEA output (release v67 or lower)
  gsea_df <- adjust.v67(input=gsea_df, gS=geneScore, tC=testedContrast, sp=species, db=pathway_db)   
  
  # gene score name
  genescore = unique(gsea_df$geneScore)
  
  # pathway set
  if ( !is.null(pathways) ) { gsea_df <- gsea_df %>% dplyr::filter(pathway %in% pathways) }
  
  # gsea stats
  gsea_df <- gsea_df %>% dplyr::select( colnames(gsea_df)[colnames(gsea_df) %in% c("contrast", "collection", "pathway", "ES", "NES", "pval", "padj", "leadingEdge","size_leadingEdge", "inPathway")] )
  if ( !is.null(contrasts) ) { 
    gsea_df <- gsea_df %>% dplyr::filter(contrast %in% contrasts)
  } 
  
  #..GSEA ES  
  
  if (grepl("ES", which_plot)) {
    
    # gene scores    
    rank_columns = colnames(genescore_df)[grepl(paste0("\\Q", genescore, "\\E$"), colnames(genescore_df))]
    rank_contrasts = unlist(strsplit(rank_columns, genescore))
    if ( length(rank_columns)==0 ) stop("ERROR: 'Gene score' not specified correctly")
    if ( length(rank_contrasts)==0 ) stop("ERROR: 'Tested contrasts' not specified correctly")
    
    
    if (!is.null(contrasts)) {  
      
      index = match(contrasts, rank_contrasts)
      rank_columns = rank_columns[index]
      rank_contrasts = rank_contrasts[index]
      groups_from_contrasts = unique(unlist(strsplit(rank_contrasts,"-"))) 
      
    } else if (is.null(contrasts)) {         
      rank_contrasts = unique(gsea_df$contrast)
    }
    
    genescore_df <- genescore_df %>% dplyr::select(geneid, rank_columns) %>% tidyr::pivot_longer(!geneid, names_to="contrast", values_to="genescores", values_drop_na=TRUE) %>% dplyr::rename("geneid"=geneid) %>% dplyr::mutate(contrast=sub(genescore, "", contrast)) %>% tidyr::drop_na() %>% dplyr::arrange(desc(genescores))
    
  }
  
  #..LE HEATMAP
  
  if (grepl("LE", which_plot)) {
    
    # samples in gene expression dataset
    samples_to_include = setdiff(c("X1_CD19_314_T_resent","X2_Y_Y_319_T_resent","X3_Y_Y_304_T_resent","X4_Y_Y_729_T_resent","X5_Y_Y_727_T_resent","X6_CD19_738_T_resent","X7_CD19_271_T_resent","X8_Y_H_305_T_resent","X9_CD19_324_T_resent","X10_hYP7_731_T_resent","X11_Y_H_314_T_resent","X12_hYP7_740_T_resent","X13_Y_H_737_T_resent","X14_hYP7_732_T_resent","X15_Y_H_730_T_resent","X16_hYP7_303_T_resent"), geneid)
    
    # gene expression
    le_genes <- gsea_df %>% dplyr::select(leadingEdge) %>% tidyr::separate_rows(leadingEdge, sep=",") %>% dplyr::distinct()
    gex_df = gex_df %>% dplyr::select(geneid_gex, samples_to_include) %>% dplyr::rename("Gene"=geneid_gex) %>% inner_join(le_genes, by=c("Gene"="leadingEdge"))
  }    
  
  # DO PLOT ====
  
  cat(sprintf("Saving files in the workbook-output:\n\nGSEA-RunningES.csv\n\n"))
  
  # retain only top significant pathways if top rank requested
  
  plot_limit = tolower(plot_limit)
  
  if ( plot_limit != "all" ) {    
    
    n_input = nrow(gsea_df)
    
    plot_limit = as.numeric(plot_limit)
    
    if (is.na(plot_limit)) { 
      stop("ERROR in Top rank filter; enter ALL (case insensitive) or a numeric rank.\n")
    } else if (plot_limit <= 0 ) { 
      plot_limit = 1 
      cat("WARNING: 'Top rank filter' cannot be 0 or less; its value was changed to 1.\n")
    }
    
    gsea_df <- gsea_df %>% dplyr::group_by(contrast, collection) %>% dplyr::mutate(p_rank = rank(pval, ties.method="min")) %>% filter(p_rank <= plot_limit) %>% dplyr::select(-p_rank) %>% dplyr::ungroup()
    
    if(n_input > nrow(gsea_df)) {
      cat(sprintf("WARNING: Preparing top-ranked plots from each contrast and collection (%g out of %g) based on max P-value rank of %g\n\t Change 'Top rank filter'parameter if you intend to generate more plots\n", nrow(gsea_df), n_input, plot_limit))
    } else {
      cat(sprintf("Preparing top-ranked plots from each contrast and collection (%g out of %g) based on max P-value rank of %g\n", nrow(gsea_df), n_input, plot_limit))
    } } else {
    
    cat(sprintf("Preparing all available plots (%g)\n", nrow(gsea_df)))
  }
  
  gsea_list <- dplyr::group_split(gsea_df, contrast, pathway) 
  class(gsea_list) <- "list"
  names(gsea_list) = sapply(gsea_list, function(x) paste(x$pathway, x$contrast, sep="_"))
  
  if ( reference_phenotype ) {
    catch_reference = strsplit(unique(gsea_df$contrast), "-")
    names(catch_reference) = unique(gsea_df$contrast)
    catch_reference = sapply(names(catch_reference), function(x) paste(x,catch_reference[[x]][2], sep=": "))
    cat(sprintf("\nReference group for a contrast in LE heatmap (gene %s):", gex_transformation))
    cat(sprintf("\n%s",catch_reference),"\n\n")
  }
  
  header=c('contrast', 'pathway', 'geneRank', 'runningES', 'gene', 'leadingEdge', 'geneScore')
  write.table(t(header), "gsea_vis.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
  
  for ( i in 1:length(gsea_list) ) {
    
    gsea = gsea_list[[i]]
    txt <- sprintf( "%s, %s, ES=%g, NES=%g, pval=%g, padj=%g", gsea$collection, gsea$contrast,  round(gsea$ES,2), round(gsea$NES,2), signif(gsea$pval,2), signif(gsea$padj,2) )
    name = gsea$pathway
    
    fontsize_row = 0
    fontsize_col = 0
    
    counter1 = seq(0,length(gsea_list),100)
    counter2 = seq(0,length(gsea_list),25)
    if (i %in% counter1 | i == length(gsea_list)) { cat(i,"\n") } else if (i %in% counter2) { cat(i, " ") } else { cat(".")}
    
    plotES = gg.plotES(ranks=genescore_df, gsea=gsea, ntop=10, add_top=TRUE, cex_top=3, image_size='reduced', gset=gset, line_top=topline, ES_colo = 'ES sign')
    gsea_out <- plotES[[3]]
    write.table(gsea_out, "gsea_out.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',', append=TRUE)
    
    
    if (which_plot == "ES+RNK+LE") {
      
      # ES with heatbar and RNK barplot
      pRunes <- plotES[[1]]
      pRank <- plotES[[2]]        
      
      # gex transformation
      le = unlist(strsplit(gsea$leadingEdge, ","))
      gex_le = gex_df %>% dplyr::filter(Gene %in% le) %>% tibble::column_to_rownames("Gene") 
      metadata = prep.metadata(meta=metadata_df, gsea=gsea, samples=colnames(gex_le), Sample="Sample", Group="Group", transformation=gex_transformation, dropREF=drop_ref, reference=reference_phenotype, own_palette=c(), label_palette='Accent', filterSamples=filter_samples)
      gex_le = gex_le[, match(metadata$Sample, colnames(gex_le))]
      gex_trans = transform.gex(gex_le, meta=metadata, transformation = gex_transformation )
      if (drop_ref) {
        if(grepl("reference", gex_transformation)) {
          gex_trans <- gex_trans[, metadata$Condition != 'Reference']
          metadata <- metadata[metadata$Condition != 'Reference', ]
        }
      }
      
      # clustering        
      clustOrder <- prep.clustOrder( df=gex_trans, linkage='complete', distance='Euclidean', way='rows and columns' )
      rowv = clustOrder$rowv
      colv = clustOrder$colv
      
      # font size
      if (fontsize_row == 0) {
        fontsize_row = find.fontsize(max_size=8, min_size=3, max_n=c(500, gsea$size_leadingEdge)[which.max(c(500, gsea$size_leadingEdge))], min_n=25, stepdown=0.6, n=gsea$size_leadingEdge) }
      if (fontsize_col == 0) {
        fontsize_col = find.fontsize(max_size=4, min_size=3, max_n=c(500, nrow(metadata))[which.max(c(500, nrow(metadata)))], min_n=25, stepdown=0.6, n=nrow(metadata)) }
      
      # heatmap
      pLEdge = plot.heatmap( gex=gex_trans, meta=metadata, heat_colors=NULL, limit=NULL, rowv=rowv, colv=colv, transformation=gex_transformation, show_rownames = TRUE, show_colnames = FALSE, show_coldend = FALSE, show_rowdend = TRUE, row_size=fontsize_row, col_size=fontsize_col, heatmap_leg = TRUE, sample_leg = TRUE)
      
      # plot layout     
      lay <- c(
        patchwork::area(t = 1, l = 1, b = 2, r = 2), 
        patchwork::area(t = 3, l = 1, b = 4, r = 2),
        patchwork::area(t = 1, l = 3, b = 4, r = 3))
      
      # plot
      patch = pRunes + pRank + patchwork::wrap_elements(pLEdge) + plot_annotation(title = name, subtitle = txt) + plot_layout( design=lay) + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(size = 14))
      
    } else if (which_plot == "ES+RNK") {
      
      # ES with heatbar and RNK barplot
      pRunes <- plotES[[1]]
      pRank <- plotES[[2]]
      
      # plot layout
      lay <- c(
        patchwork::area(t = 1, l = 1, b = 2, r = 2), 
        patchwork::area(t = 3, l = 1, b = 4, r = 2))
      
      # plot
      patch = pRunes + pRank + plot_annotation(title = name, subtitle = txt) + plot_layout( design=lay) + plot_annotation(tag_levels='A') & theme(plot.tag = element_text(size = 14))
      
    } else if (which_plot == "ES") {
      
      # ES with heatbar and RNK barplot
      pRunes <- plotES[[1]]
      
      # plot
      patch = pRunes + plot_annotation(title = name, subtitle = txt)
      
    } else if (which_plot == "LE") {
      
      # gex transformation
      le = unlist(strsplit(gsea$leadingEdge, ","))
      gex_le = gex_df %>% dplyr::filter(Gene %in% le) %>% tibble::column_to_rownames("Gene") 
      metadata = prep.metadata(meta=metadata_df, gsea=gsea, samples=colnames(gex_le), Sample="Sample", Group="Group", transformation=gex_transformation, dropREF=drop_ref, reference=reference_phenotype, own_palette=c(), label_palette='Accent', filterSamples=filter_samples)
      gex_le = gex_le[, match(metadata$Sample, colnames(gex_le))]
      gex_trans = transform.gex(gex_le, meta=metadata, transformation = gex_transformation )
      if (drop_ref) {
        if(grepl("reference", gex_transformation)) {
          gex_trans <- gex_trans[, metadata$Condition != 'Reference']
          metadata <- metadata[metadata$Condition != 'Reference', ]
        }
      }
      
      # clustering        
      clustOrder <- prep.clustOrder( df=gex_trans, linkage='complete', distance='Euclidean', way='rows and columns' )
      rowv = clustOrder$rowv
      colv = clustOrder$colv
      
      # font size
      if (fontsize_row == 0) {
        fontsize_row = find.fontsize(max_size=8, min_size=3, max_n=c(500, gsea$size_leadingEdge)[which.max(c(500, gsea$size_leadingEdge))], min_n=25, stepdown=0.6, n=gsea$size_leadingEdge)
      }
      if (fontsize_col == 0) {
        fontsize_col = find.fontsize(max_size=4, min_size=3, max_n=c(500, nrow(metadata))[which.max(c(500, nrow(metadata)))], min_n=25, stepdown=0.6, n=nrow(metadata))
      }
      
      # heatmap
      pLEdge = plot.heatmap( gex=gex_trans, meta=metadata, heat_colors=NULL, limit=NULL, rowv=rowv, colv=colv, transformation=gex_transformation, show_rownames = TRUE, show_colnames = FALSE, show_coldend = FALSE, show_rowdend = TRUE, row_size=fontsize_row, col_size=fontsize_col, heatmap_leg = TRUE, sample_leg = TRUE)
      
      # plot
      patch = wrap_elements(pLEdge) + plot_annotation(title = name, subtitle = txt)
      
    }
    
    if (pdfType == 'common PDF') {
      
      gsea_list[[i]] = patch
      if (i == 1) { preview = patch }
      
    } else {
      
      fileName =  make.names(sprintf("%s.pdf", names(gsea_list)[i]))
      pdf(output_fs$get_path(fileName, 'w'), height=pageHeight, width=pageWidth)
      print(patch)
      #dev.off()
      if (i == 1) { preview = patch }
      
    }
  }
  
  if (pdfType == 'common PDF') {
    
    fileName=sprintf("GSEA-Plot_%s.pdf", make.names(which_plot))
    cat(sprintf("\n%s\n\n", fileName))          
    pdf(output_fs$get_path(fileName, 'w'), height=pageHeight, width=pageWidth)
    lapply(gsea_list, function(patch) {
      print(patch)
      text_contrast = sapply(strsplit(patch$patches$annotation$subtitle, ", "), function(x) x[2])
      text_pathway = patch$patches$annotation$title
      cat(sprintf("%s (%s)\n", text_pathway, text_contrast))
    })
    #dev.off()
    
  } else {
    cat(sprintf("\n%s.pdf", make.names(names(gsea_list))))
  }
  
  png(filename=graphicsFile, width=(pageWidth+5)*300, height=(pageHeight+5)*300, units="px", pointsize=4, bg="white", res=300, type="cairo") 
  print(preview)
  
}

#################################################
## Global imports and functions included below ##
#################################################

# Functions defined here will be available to call in
# the code for any table.

gg.plotES <- function(ranks, gsea, ntop, add_top, image_size, gset, ES_colo, line_top, cex_top) {
  
  # get all ranks
  ranks <- ranks %>% dplyr::filter(contrast %in% gsea$contrast)
  rnk=ranks$genescores
  names(rnk) = ranks$geneid
  le = gsea %>% dplyr::select(leadingEdge) %>% tidyr::separate_rows(leadingEdge) 
  lep=is.element(names(rnk), le$leadingEdge) & rnk>0
  len=is.element(names(rnk), le$leadingEdge) & rnk<0
  zero = sum(rnk>0)
  yrange = diff(range(rnk))
  ymax = max(rnk)
  ymin = min(rnk)
  
  # get gene set ranks and running scores
  ES_sign = sign(gsea$NES)
  gset = gsea %>% dplyr::select(inPathway) %>% tidyr::separate_rows(inPathway)
  gs = intersect(paste(gset$inPathway),names(rnk))
  es.data <- plotEnrichment(pathway=gs, stats=rnk, gseaParam = 1)
  x <- es.data$data$x
  y <- es.data$data$y
  xranks = sort(unname(as.vector(na.omit(match(gs, names(rnk))))))
  
  gs_ranks = data.frame(xranks=xranks)
  gs_scores = data.frame(x = x, y = y)
  ys = round( y[ -c(1, length(y)) ], 10)
  if(ES_sign > 0) {
    keep=seq(2,length(ys), 2)
  } else {
    keep=seq(1,length(ys), 2)
  }    
  es_data=data.frame(rank=xranks,running_es=ys[keep])
  
  # set positive/negative params    
  
  if( ES_sign > 0 ){
    
    all_ranks = data.frame(Index=1:length(rnk), Rank=rnk, LE=ifelse(lep==TRUE,'LE','Outside'), order=ifelse(lep==TRUE,2,1))
    
    if(add_top==TRUE & ntop > 0){
      if(sum(lep) < ntop) {
        ntop = sum(lep)                
        warning(sprintf("Max number of leading edge genes available is %g", sum(lep)))
      }
      top = sort(rnk[lep], decreasing=TRUE)[1:ntop]; top = paste(names(top), collapse='\n')
      topx = 1; topy= ymin
      nx=which(all_ranks$LE=="LE")[1]; ny = 0-yrange/30
      v=1; h=0; hn=0
      colorMargin = color.margins()['up']
      angletop=90
    }
    
  } else if ( ES_sign < 0) {
    
    all_ranks = data.frame(Index=1:length(rnk), Rank=rnk, LE=ifelse(len==TRUE,'LE','Outside'), order=ifelse(len==TRUE,2,1))
    
    if(add_top==TRUE & ntop > 0){
      if(sum(len) < ntop) {
        ntop = sum(len)
        warning(sprintf("Max number of leading edge genes available is %g", sum(len)))
      }
      top = sort(rnk[len],decreasing=FALSE)[1:ntop]; top = paste((names(top)), collapse='\n') # rev(names(top)) if angle=90
      topx = length(rnk); topy= ymax
      nx=which(all_ranks$LE=="LE")[sum(len)]; ny=0+yrange/30
      v=1; h=0; hn=1
      colorMargin = color.margins()['dn']
      angletop=-90
    } 
  }
  
  # set miscelanous
  
  #.. zero arrow
  df_arrow <- data.frame(x1 = zero+0.5, x2 = zero+0.5, y1 = 0, y2 = 0+yrange/16)
  
  #.. keep only gene set ranks
  if(image_size == 'reduced') { all_ranks$Rank = ifelse(all_ranks$LE=='LE', all_ranks$Rank, NA) } 
  
  #.. color ranks    
  qua = quantile(abs(rnk), 0.95)
  newrnk = ifelse(rnk > qua, qua, rnk); newrnk[newrnk < -qua] = -qua
  all_ranks$Ranklimit = newrnk
  
  #.. running score line color
  if(ES_colo == "green") {
    line_colo = "green2"
  } else {
    line_colo = colorMargin }
  
  #.. top line type
  topL = ifelse(ES_sign==1 , max(y), min(y))
  which_topL = which(y==topL)
  df_topL = data.frame(x1 = x[which_topL], y1 = topL, x2 = x[which_topL], y2=0) 
  
  #.. base text size
  base = 12
  
  # generate rank subplot (barplot + heatbar)
  
  p <- ggplot( all_ranks, aes( x=Index, y=Rank ) ) +
    
    geom_bar(stat="identity", aes(fill=LE), width=10, order=order, color=NA, show.legend=FALSE) + 
    
    scale_fill_manual(values=c("Outside"="#D8D8D855","LE"=paste(colorMargin))) +  
    
    scale_y_continuous(limits=c( min(rnk), max(rnk) )) +        
    
    xlab("Rank") + ylab("Gene score") + 
    
    annotate(geom="text", x=topx, y=topy, label=top, angle=angletop, color=colorMargin, hjust=h, size=cex_top, vjust=v) + 
    
    theme_bw() + theme( text = element_text(size = base+3), axis.text = element_text(size = base)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    
    annotate(geom="text", x=nx, y=ny, label=sprintf("N=%g",sum(all_ranks$LE=="LE")), color=colorMargin, hjust=hn, size=3) +
    
    geom_segment(data=df_arrow, aes(x = x1, y = y2, xend = x2, yend = y1), arrow =arrow(length = unit(0.2, "cm"), angle=25, type='closed'), size=0.3, color='black', inherit.aes=FALSE) +            
    annotate(geom="text", x=zero+0.5, y=0+yrange/11, label=paste('zero crossed at', zero+1), color="black", size=3)
  
  if(image_size == 'reduced'){
    
    p = p +  geom_segment(aes(xend = 1, y = 0, x = length(rnk), yend = 0), col="black", size=0.1) +        
      annotate(geom='text',x=-Inf,y= -Inf, label="+", hjust=-0.4, vjust=-0.2, size=7) +
      annotate(geom='text',x=Inf, y=-Inf, label="_", hjust=1.7, vjust=-1, size=5.5, fontface='bold')
  }
  
  if(length(rnk) >= 1000){
    
    p = p + scale_x_continuous(position='bottom', limits=c(0,length(rnk)), labels = function(l) {trans=l/1000; paste0(trans, "K")})
    
  } else {
    
    p = p + scale_x_continuous(position='bottom', limits=c(0,length(rnk)))
  }
  
  pRank <- p +  theme(plot.margin = unit(c(l=0,r=0.03,t=0,b=0), "npc"))
  
  # generate running ES sublot
  
  q <- ggplot(gs_scores, aes(x=x, y=y)) + geom_line(color=line_colo) + 
    
    theme_bw() +  theme( text = element_text(size = base+3), axis.text = element_text(size = base)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    
    scale_y_continuous(expand=expand_scale(add = c(0.1, 0)) ) +
    
    #annotate(geom='text',-Inf, Inf, label="+", hjust=-0.1, vjust=0, size=7) +            
    #annotate(geom='text',Inf, Inf, label="\u2013", hjust=1.2, vjust=0, size=4, fontface='bold') +
    
    ylab("Enrichment Score (ES)") + xlab("Rank") +
    
    geom_hline(yintercept=0)
  
  if(image_size == 'reduced') {
    
    q = q + geom_rug(data=gs_ranks, aes(x=xranks), sides='b', show.legend=FALSE, length = unit(0.04, "npc"), inherit.aes=FALSE, size=0.3)
    
  } else {
    
    q = q + geom_rug(data=gs_ranks, aes(x=xranks), sides='b', show.legend=FALSE, length = unit(0.08, "npc"), inherit.aes=FALSE, size=0.3) +
      
      geom_rug(data=all_ranks, aes(x=Index, color=Ranklimit), sides='b', show.legend=FALSE, length = unit(0.05, "npc"), inherit.aes=FALSE) +         
      scale_color_gradient2(high=color.margins()['up'], low=color.margins()['dn'], mid=color.margins()['md'], midpoint=0, limits=c(-1,1)*qua )
  }
  
  if(line_top == 'horizontal'){
    
    q = q + geom_hline(yintercept=topL, linetype='dashed')
    
  } else if(line_top == "coordinate") {
    
    q = q + geom_segment(data=df_topL, aes(x = x1, y = y2, xend = x2, yend = y1), size=0.3, color=line_colo, linetype='dashed', inherit.aes=FALSE)
  }
  
  if(length(rnk) >= 1000){
    
    q = q + scale_x_continuous(position='bottom', limits=c(0,length(rnk)+1), labels = function(l) {trans=l/1000; paste0(trans, "K")})
    
  } else { 
    
    q = q + scale_x_continuous(position='bottom', limits=c(0,length(rnk)+1))
    
  }
  pRunes <- q + theme(plot.margin = unit(c(l=0,r=0.03,t=0,b=0), "npc")) +        
    annotate(geom='text',x=-Inf,y= -Inf, label="+", hjust=-0.4, vjust=-0.2, size=7) +
    annotate(geom='text',x=Inf, y=-Inf, label="_", hjust=1.7, vjust=-1, size=5.5, fontface='bold')
  
  # grobs + output table
  
  all_ranks$Rank = round(ranks$genescores[match(rownames(all_ranks), ranks$geneid)],10)
  gsea_out = tibble::rownames_to_column(all_ranks, 'geneid') %>% dplyr::filter(Index %in% gs_ranks$xranks) %>% dplyr::mutate("leadingEdge"=ifelse(LE=='Outside',FALSE,TRUE)) %>% dplyr::inner_join(es_data, by=c('Index'='rank'))  %>% dplyr::rename('geneScore'='Rank', 'geneRank'="Index", "runningES"="running_es", "gene"='geneid') %>% dplyr::mutate(contrast=gsea$contrast, pathway=gsea$pathway) %>% dplyr::select(contrast,pathway, geneRank, runningES, gene, leadingEdge, geneScore) 
  
  return(list( pRunes = pRunes, pRank = pRank, gsea_out = gsea_out ))
}

## prep sample metadata
prep.metadata <- function(meta, gsea, samples, Sample, Group, transformation, dropREF, reference, own_palette, label_palette, filterSamples) {
  
  meta <- dplyr::select(meta, c(as.name(Sample),as.name(Group))) %>% dplyr::rename("Sample"=Sample, "Group"=Group)
  meta <- meta[ na.omit(match(samples, meta$Sample)), ]
  groups = unlist(strsplit(gsea$contrast, "-"))
  if (filterSamples) {
    meta <- meta %>% dplyr::filter(Sample %in% samples)
  } else {    
    meta <- meta %>% dplyr::mutate(Group = gsub(" ", "_", Group)) %>% dplyr::filter(Group %in% groups)
  }
  if (reference) { reference_group = groups[2] }
  
  if (!reference) {
    
    if (grepl('reference', transformation)) {
      stop("\nERROR: Reference phenotype equels FALSE, but gene transformation with reference phenotype selected (Heatmap parameters)\n")
      
    } else {
      meta$Condition <- rep('Experiment', nrow(meta))
    }
    
  } else if (reference) {
    
    if (! grepl('reference', transformation)) {
      stop("\nERROR: Reference phenotype equels TRUE, but gene transformation with reference phenotype not selected (Heatmap parameters)\n")
      
    } else {
      
      if (reference_group %in% meta$Group) {            
        meta$Condition <- ifelse(meta$Group==reference_group, 'Reference', 'Experiment')
        
      } else if (! reference_group %in% meta$Group) { 
        stop(sprintf("\nERROR: Reference phenotype (%s) not found \nSUGGESTION: check if input datasets are correct that is sample metadata (Sample id and Group id), gene expression (column names), and (GSEA contrasts)\n", reference))
      }
    }
  }  
  
  color_label <- factor(meta$Group)
  if (label_palette == "Custom") {
    colors <- own_palette
    if( length(colors) < length(levels(color_label)) & dropREF == TRUE) {
      colors <- c(colors,rep('grey'),1)
    } else if(length(colors) != length(levels(color_label)) ){
      stop('ERROR: select number of colors in Custom palette equal to the number of phenotype labels')
    }
  } else {        
    n = ifelse(is.element(label_palette, c('Paired','Set3')), 12, 8)
    n_color <- length(levels(color_label))
    colors <- colors <- brewer.pal(n, label_palette)[1:n_color]
  }
  levels(color_label) <- colors
  meta$Group_color <- paste(color_label)
  
  return(meta)
}

## transfer gene expression by row with R (input r data.frame)
transform.gex <- function(df, transformation, meta) {
  
  if ( ! all(colnames(df) == meta$sample) ) { 
    
    stop("\nERROR: column names in gene expression dataset and sample metadata are not matched correctly: contact template maintainer at michaloa@mail.nih.gov")
  }
  
  if (transformation == 'median centering') {
    mat <- t(apply(df, 1, function(y) y-median(y, na.rm=TRUE)) )
    
  } else if ( transformation == 'mean centering') {
    mat <- t(apply(df, 1, function(y) y-mean(y, na.rm=TRUE)) )
    
  } else if (transformation == 'z-score') {
    mat <- t(apply(df, 1, function(y) (y-mean(y, na.rm=TRUE))/sd(y, na.rm=TRUE)) )
    
  } else if (transformation == 'reference median centering' & !is.null(meta)) {
    mat <- df - apply(df[, meta$Condition == 'Reference'], 1, median, na.rm=TRUE)
    
  } else if (transformation == 'reference mean centering' & !is.null(meta) ) {
    mat <- df - apply(df[, meta$Condition == 'Reference'], 1, mean, na.rm=TRUE)
    
  } else {
    mat <- df
  } 
  
  return(mat)
}

## clustering order
prep.clustOrder <- function(df, linkage, distance, way) {
  
  rowv=FALSE
  colv=FALSE
  
  if( (way == 'rows'| way == 'rows and columns') ){
    
    if (distance=='1-Spearman') {    
      rowv = hclust(as.dist(1 - cor(t(df), method='spearman')), method = linkage)
    } else if(distance == "1-Pearson") {
      rowv = hclust(as.dist(1 - cor(t(df), method='pearson')), method = linkage)
    } else if (distance == 'Euclidean'){
      rowv = hclust(dist(df), method = linkage)
    } else { 
      rowv = FALSE
    }}
  
  if( (way == 'columns'| way == 'rows and columns') ){
    
    if (distance=='1-Spearman') {    
      colv = hclust(as.dist(1 - cor(df, method='spearman')), method = linkage)
    } else if(distance == "1-Pearson") {
      colv = hclust(as.dist(1 - cor(df, method='pearson')), method = linkage)
    } else if (distance == 'Euclidean'){
      colv = hclust(dist(t(df)), method = linkage)
    } else { 
      colv = FALSE
    }}
  
  return(list(rowv=rowv,colv=colv))
}

# plot heatmap
plot.heatmap <- function(gex, meta, heat_colors=NULL, limit, rowv, colv, transformation, show_rownames, show_colnames, show_coldend, show_rowdend, row_size, col_size, heatmap_leg, sample_leg) {
  
  meta_label = meta$Group
  meta_color = unlist(lapply(split(meta$Group_color, meta_label), unique))
  ha = HeatmapAnnotation(Class = meta_label, col=list(Class=meta_color), annotation_height = unit(rep(0.3,1), "cm"), annotation_legend_param = list(title_gp = gpar(fontsize = 7), grid_width=unit( c(0.3), "cm"), labels_gp = gpar(fontsize = 8)))
  
  pal = color.margins()
  if(is.null(limit)) {
    limit = quantile(abs(as.matrix(gex)), 0.95)
  } else { limit = limit }
  palette_function <- circlize::colorRamp2( c(-limit, 0, limit), c(pal['dn'],pal['md'],pal['up']), space="LAB")
  if(transformation=='z-score') { heat_legend ='SD' } else { heat_legend = 'log2' }
  
  # set clustering matrix
  
  h1 = Heatmap(gex, col=palette_function, cluster_rows=rowv, cluster_columns=colv, name = heat_legend, column_title = "", column_title_gp = gpar(fontsize = 12), column_title_side = "top", show_column_names = show_colnames, column_names_gp = gpar(fontsize = col_size), row_names_gp = gpar(fontsize = row_size), show_row_names = show_rownames, column_dend_height = unit(1, "cm"), column_dend_reorder = T, row_dend_reorder=T, show_column_dend=show_coldend, show_row_dend=show_rowdend, top_annotation=ha,heatmap_legend_param = list(color_bar = "continuous", title_gp = gpar(fontsize = 7),labels_gp = gpar(fontsize = 6)))
  
  gb_heatmap = grid.grabExpr(draw(h1, heatmap_legend_side='right',  annotation_legend_side='right' , show_heatmap_legend = heatmap_leg, show_annotation_legend = sample_leg) )
  
  return(gb_heatmap)
}     

# color ranks
color.ranks <- function(val) {
  posVal=val[val>=0]
  up=sequential_hcl(length(posVal),h=0,c.=c(180,0),l=c(30,90),power=1.5,gamma=NULL,fixup=TRUE,alpha=1)[rev(rank(posVal))]
  downVal=val[val<0]
  down= sequential_hcl(length(downVal),h=260,c.=c(90,0),l=c(30,90),power=2,gamma=NULL,fixup=TRUE,alpha=1)[rank(downVal)]
  return(c(up,down))
}

# set marginal colors (default: red, blue, whitesmoke, no transparency)
color.margins <- function(dn_hew=260, up_hew=0, md_hew=94, dn_c=90, up_c=180, md_c=0, dn_l=30, up_l=30, md_l=97, a=1) {
  dn = sequential_hcl(1, h = dn_hew, c. = c(dn_c), l = c(dn_l), fixup = TRUE, alpha = a)
  up = sequential_hcl(1, h = up_hew, c. = c(up_c), l = c(up_l), fixup = TRUE, alpha = a)
  md = sequential_hcl(1, h = md_hew, c. = c(md_c), l = c(md_l), fixup = TRUE, alpha = a)
  return(c(dn=dn,md=md,up=up))
}

# set font size in LE heatmap
find.fontsize <- function(n, max_size, min_size, stepdown, max_n, min_n) {
  
  numbers = seq(min_n, max_n, by=25)
  font_sizes = seq(max_size, min_size, by=-stepdown)
  index = which.min(abs(numbers-n))
  if (index > length(font_sizes)) { index = length(font_sizes) }
  return(font_sizes[index])
}

# adjust old GSEA table (released v67 or lower)
# adjust old GSEA table (released v67 or lower)
adjust.v67 <- function(
    input, gS, tC, sp, db,
    required_columns = c(
      "contrast","geneScore","fdr_correction_mode","collection","pathway",
      "pval","padj","ES","NES","nMoreExtreme","size","leadingEdge",
      "size_leadingEdge","inPathway"
    )
) {
  
  missing_columns = setdiff(required_columns, colnames(input))
  
  if (length(missing_columns) > 0) {
    
    cat("WARNING: Output from outdated 'Preranked GSEA [CCBR]' detected.\n")
    cat(sprintf(
      "\nWARNING: Missing columns to be added (%s)\n\t 'inPathway' column will include all genes annotated to a gene set - run the latest released version of 'Preranked GSEA [CCBR]'\n\t to return only genes mapped in the dataset.\n",
      paste(missing_columns, collapse = ", ")
    ))
    
    if (is.null(gS)) {
      stop("\nERROR: 'Gene score' parameter not provided (Advanced Parameters)\n")
    } else {
      
      if (length(gS) > 1) {
        gS = gS[1]
        cat(sprintf(
          "\nWARNING: too many values for 'Gene score' provided; only the first one used, '%s'\n",
          gS
        ))
      }
      
      has.underscore = substring(gS, 1, 1) == "_"
      
      if (!has.underscore) {
        gS = paste0("_", gS)
        cat(sprintf(
          "\nWARNING: 'Gene score' (Advanced Parameters) should start with underscore to match column naming convention in DEG table used for GSEA run (e.g. '_tstat');\n\t column 'geneScore' is assigned the provided value with '_' added ('%s');\n\t this template will fail if this is not an adequate specification.\n",
          gS
        ))
      }
    }
    
    if (is.null(tC)) {
      tC = "NULL"
      stop("\nERROR: 'Tested contrast' parameter not provided (Advanced Parameters)\n")
    } else {
      
      if (length(tC) > 1) {
        tC = tC[1]
        cat(sprintf(
          "\nWARNING: too many values for 'Tested contrast' provided; only the first one used, '%s'\n",
          tC
        ))
      }
      
      is.contrast = any(grepl("-", tC))
      
      if (!is.contrast) {
        stop(sprintf(
          "\nERROR: 'Tested contrast' (Advanced Parameters) should be specified exactly the same way as when Preranked GSEA was run (e.g. treated-control);\n\t provided value is '%s'",
          tC
        ))
      }
    }
    
    input <- input %>%
      dplyr::mutate(
        contrast = tC,
        geneScore = gS,
        fdr_correction_mode = "over all collections",
        size_leadingEdge = sapply(strsplit(leadingEdge, ","), length)
      )
    
    input <- input[, match(required_columns[required_columns != "inPathway"], colnames(input))]
    
    db <- db %>%
      dplyr::filter(species == sp) %>%
      dplyr::filter(collection %in% unique(input$collection)) %>%
      dplyr::filter(gene_set_name %in% unique(input$pathway)) %>%
      dplyr::collect()
    
    db <- db %>%
      dplyr::group_by(collection, gene_set_name, species) %>%
      dplyr::summarize(inPathway = paste(unique(gene_symbol), collapse = ",")) %>%
      dplyr::ungroup() %>%
      data.frame()
    
    input <- input %>%
      dplyr::left_join(db, by = c("collection" = "collection", "pathway" = "gene_set_name"))
    
    input <- input[, na.omit(match(c(required_columns, "species"), colnames(input)))]
    
    return(input)
  } else {
    
    return(input)
  }
}
