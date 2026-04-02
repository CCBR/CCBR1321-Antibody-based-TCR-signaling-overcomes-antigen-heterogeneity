input <- read.csv("ccbr1321_RawCountFile_rsemgenes.csv")

# clean counts
CleanRawCounts <- function(Remove_untitled_col) {
  
  library(stringr)
  library(tidyr)
  library(dplyr)
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters:
  raw_counts_matrix=Remove_untitled_col
  Data_type='Bulk RNAseq'
  gene_id_column='gene_id'
  samples_to_rename = c("")
  
  #Advanced Parameters:
  cleanup_column_names=TRUE
  split_gene_name = TRUE
  aggregate_rows_with_duplicate_gene_names=TRUE
  gene_name_column_to_use_for_collapsing_duplicates=''
  
  ###################################################################################
  removeVersion <- function(ids){
    return(unlist(lapply(stringr::str_split(ids, "[.]"), "[[",1)))
  }  
  
  print(Data_type)
  
  ##################################     
  ##### Sample Name Check
  ################################## 
  
  ## duplicate col name
  if(sum(duplicated(colnames(raw_counts_matrix)))!=0){
    print("Duplicate column names are not allowed, the following columns were duplicated.\n")
    colnames(raw_counts_matrix)[duplicated(colnames(raw_counts_matrix))]
    stop("Duplicated columns")
  }
  
  ##################################     
  ##### Manually rename samples
  ################################## 
  
  if (!is.null(samples_to_rename)) {
    if (samples_to_rename != c("")) {
      for (x in samples_to_rename) {
        old <- strsplit(x, ": ?")[[1]][1]
        new <- strsplit(x, ": ?")[[1]][2]
        colnames(raw_counts_matrix)[colnames(raw_counts_matrix)%in%old]=new
      }
    }
  }
  
  
  
  ##################################     
  ##### Cleanup Columns
  ##################################    
  if(cleanup_column_names){
    cl_og=colnames(raw_counts_matrix)
    ## convert special charchers to _
    cl2 <- gsub('-| |\\:','_',colnames(raw_counts_matrix))
    if (length(cl2[(cl2)!=colnames(raw_counts_matrix)])>0) {
      print('Columns had special characters relpaced with _ ')
      # (colnames(raw_counts_matrix)[(colnames(raw_counts_matrix))!=cl2])
      # print(cl2[(cl2)!=colnames(raw_counts_matrix)])
      colnames(raw_counts_matrix) = cl2
    }
    
    ## if names begin with number add X
    cl2=sub("^(\\d)", "X\\1", colnames(raw_counts_matrix))
    if (length(cl2[(cl2)!=colnames(raw_counts_matrix)])>0) {
      
      print('Columns started with numbers and an X was added to colname :')
      # (colnames(raw_counts_matrix)[(colnames(raw_counts_matrix))!=cl2])
      # print(cl2[(cl2)!=colnames(raw_counts_matrix)])
      colnames(raw_counts_matrix) = cl2
    }
    #print("Original Colnames:")
    #print(cl_og[(cl_og)!=colnames(df)])
    #print("Modified Colnames:")
    #print(colnames(df)[colnames(df)!=(cl_og)]%>%as.data.frame)
    
    #print("Final Colnames:")   
    
  }else{
    
    ## invalid name format
    if(any(make.names(colnames(raw_counts_matrix))!=colnames(raw_counts_matrix))){
      print("Error: The following counts matrix column names are not valid:\n")
      print(colnames(raw_counts_matrix)[make.names(colnames(raw_counts_matrix))!=colnames(raw_counts_matrix)])
      print("Likely causes are columns starting with numbers or other special characters eg spaces.\n")
      # stop("Bad column names.")
    }
    ## Names Contain dashes
    if(sum(grepl("-",colnames(raw_counts_matrix)))!=0){
      print("The sample names cannot contain dashes.")
      print(colnames(raw_counts_matrix)[grepl("-",colnames(raw_counts_matrix))])
      # stop("No dashes allowed in column names")
    }
  }
  
  
  ##################################    
  ## Split Ensemble + Gene name
  ##################################
  ## First check if Feature ID column  can be split by ",|_-:"
  ## Then check if one column contains Ensemble (regex '^ENS[A-Z]+[0-9]+')
  ##   check if Ensemble ID has version info and remove version
  ##   If one column contains Ensemble ID Assume other column is Gene names
  ## If Column does not contain Ensmeble ID name split columns Gene_ID_1 and Gene_ID_2
  print("")
  
  if(split_gene_name==T){
    Ensembl_ID=  str_split_fixed(raw_counts_matrix[,gene_id_column],'_|-|:|\\|',n=2)%>%data.frame()
    EnsCol= apply(Ensembl_ID, c(1,2), function(x) grepl('^ENS[A-Z]+[0-9]+', x))
    
    
    if(""%in%Ensembl_ID[,1]|""%in%Ensembl_ID[,2]){
      print(paste0("Not able to identify multiple id's in ", gene_id_column ))
      # colnames(df)[colnames(df)%in%clm]=gene_col
      if(Data_type=='Bulk RNAseq') { 
        colnames(raw_counts_matrix)[colnames(raw_counts_matrix)%in%gene_id_column]='Gene'
      }else if(Data_type=='Proteomics'){
        colnames(raw_counts_matrix)[colnames(raw_counts_matrix)%in%gene_id_column]='FeatureID'
      }else { print('incorrect Data Type'); incorrect_Data_Type }
    }else{
      ## at least one column must have all ensemble ids found in EnsCol 
      if (nrow(EnsCol[EnsCol[,1]==T,])==nrow(Ensembl_ID)|nrow(EnsCol[EnsCol[,2]==T,])==nrow(Ensembl_ID)){
        if(Data_type=='Bulk RNAseq') { 
          colnames(Ensembl_ID)[colSums(EnsCol)!=nrow(Ensembl_ID)]='Gene'
        }else if(Data_type=='Proteomics'){
          colnames(Ensembl_ID)[colSums(EnsCol)!=nrow(Ensembl_ID)]='FeatureID'
        }
        ## check if Ensmble column has version information
        if(grepl('^ENS[A-Z]+[0-9]+\\.[0-9]+$', Ensembl_ID[,colSums(EnsCol)==nrow(Ensembl_ID)])%>%sum()==nrow(Ensembl_ID)){
          colnames(Ensembl_ID)[colSums(EnsCol)==nrow(Ensembl_ID)]='Ensembl_ID_version'
          Ensembl_ID$Ensembl_ID=removeVersion(Ensembl_ID$Ensembl_ID_version)
        }else{
          colnames(Ensembl_ID)[colSums(EnsCol)==nrow(Ensembl_ID)]='Ensembl_ID'
        }
      }else{
        colnames(Ensembl_ID)=c('Feature_id_1','Feature_id_2')
        print("Could not determine ID formats from split 'Feature ID' Column")
        
      }
      raw_counts_matrix <- cbind(Ensembl_ID,raw_counts_matrix[,!colnames(raw_counts_matrix)%in%gene_id_column])
    }         
  }else{
    if(Data_type=='Bulk RNAseq') { 
      colnames(raw_counts_matrix)[colnames(raw_counts_matrix)%in%gene_id_column]='Gene'
    }else if(Data_type=='Proteomics'){
      colnames(raw_counts_matrix)[colnames(raw_counts_matrix)%in%gene_id_column]='FeatureID'
    }else { print('incorrect Data Type'); incorrect_Data_Type }
  }
  
  ##################################
  ## If duplicate gene aggregate information to single row
  ##################################   
  ## If user uses "Feature ID" column then switch to empty for appropriate behavor based on other parameters
  if(gene_name_column_to_use_for_collapsing_duplicates==gene_id_column){
    gene_name_column_to_use_for_collapsing_duplicates=""
  }
  
  if(gene_name_column_to_use_for_collapsing_duplicates==""&
     ('Feature_id_1'%in%colnames(raw_counts_matrix))==F){
    if(Data_type=='Bulk RNAseq') { 
      gene_name_column_to_use_for_collapsing_duplicates='Gene'
    }else if(Data_type=='Proteomics'){
      gene_name_column_to_use_for_collapsing_duplicates='FeatureID'
    }
  }  
  
  #geneids<-df[,gene_col]
  nums <- unlist(lapply(raw_counts_matrix, is.numeric)) 
  nums = names(nums[nums])
  print('')
  print('Columns that can be used to aggregate gene information' )
  print(raw_counts_matrix[,!names(raw_counts_matrix) %in% nums,drop=F]%>%colnames())
  
  print('')
  
  
  if(gene_name_column_to_use_for_collapsing_duplicates==""){
    
    if(split_gene_name==F){     
      ## If no additional Column name given for Aggregation then display Feature ID duplicates
      print(paste0("genes with duplicate IDs in ",gene_id_column,":")) 
      
      ## Print original Column name for user Reference then use new Column name to subset table
      if(Data_type=='Bulk RNAseq') { 
        gene_id_column='Gene'
      }else if(Data_type=='Proteomics'){
        gene_id_column='FeatureID'
      }
      raw_counts_matrix[duplicated(raw_counts_matrix[,gene_id_column]),gene_id_column]%>%unique()%>%as.character()%>%write( stdout())
      
    }else if(split_gene_name==T&grepl('Feature_id_1',colnames(raw_counts_matrix))==F){  
      if(Data_type=='Bulk RNAseq') { 
        gene_id_column='Gene'
      }else if(Data_type=='Proteomics'){
        gene_id_column='FeatureID'
      }
      print(paste0("genes with duplicate IDs in ",gene_id_column,":"))
      
      raw_counts_matrix[duplicated(raw_counts_matrix[,gene_name_column_to_use_for_collapsing_duplicates]),gene_name_column_to_use_for_collapsing_duplicates]%>%unique()%>%as.character()%>%write( stdout())
      
      
    }else if(split_gene_name==T&grepl('Feature_id_1',colnames(raw_counts_matrix))==T){  
      print(paste0("genes with duplicate IDs in ",'Feature_id_1',":"))
      
      raw_counts_matrix[duplicated(raw_counts_matrix[,'Feature_id_1']),'Feature_id_1']%>%unique()%>%as.character()%>%write( stdout())
      
      print(paste0("genes with duplicate IDs in ",'Feature_id_2',":"))
      
      raw_counts_matrix[duplicated(raw_counts_matrix[,'Feature_id_2']),'Feature_id_2']%>%unique()%>%as.character()%>%write( stdout())
      
    }
  }
  
  if(aggregate_rows_with_duplicate_gene_names == TRUE){
    
    print("Aggregating the counts for the same ID in different chromosome locations.")
    print("Column used to Aggregate duplicate IDs: ")
    print(gene_name_column_to_use_for_collapsing_duplicates)
    print("Number of rows before Collapse: ")
    print(nrow(raw_counts_matrix))
    
    if(sum(duplicated(raw_counts_matrix[,gene_name_column_to_use_for_collapsing_duplicates]))!=0){
      print("")
      print("Duplicate IDs: ")
      print(raw_counts_matrix[duplicated(raw_counts_matrix[,gene_name_column_to_use_for_collapsing_duplicates]),gene_name_column_to_use_for_collapsing_duplicates]%>%as.character%>%unique)
      
      dfagg=raw_counts_matrix[,c(gene_name_column_to_use_for_collapsing_duplicates,nums)]%>%group_by_at(gene_name_column_to_use_for_collapsing_duplicates)%>%summarise_all(sum)
      
      if (ncol(raw_counts_matrix[,!names(raw_counts_matrix) %in% nums, drop = FALSE])>1) {
        ## collapse non-numeric columns
        dfagg2=raw_counts_matrix[,!names(raw_counts_matrix) %in% nums]%>%group_by_at(gene_name_column_to_use_for_collapsing_duplicates)%>%summarise_all(paste,collapse=',')
        
        dfagg=merge(dfagg2,dfagg,by=eval(gene_name_column_to_use_for_collapsing_duplicates),sort = F)%>%as.data.frame()
      }
      dfout=dfagg
      print("Number of rows after Collapse: ")
      print(nrow(dfout))
    }else{
      print(paste0("no duplicated IDs in ",gene_name_column_to_use_for_collapsing_duplicates))
      dfout=raw_counts_matrix
    }
  }else{
    if(gene_name_column_to_use_for_collapsing_duplicates!=""){
      print("")
      print(paste0("Duplicate IDs in ",gene_name_column_to_use_for_collapsing_duplicates," Column:"))
      print(raw_counts_matrix[duplicated(raw_counts_matrix[,gene_name_column_to_use_for_collapsing_duplicates]),gene_name_column_to_use_for_collapsing_duplicates]%>%as.character%>%unique)
    }
    
    print("")
    print(paste0("If you desire to Aggregate row feature information select appropriate Column to use for collapsing duplicates"))
    
    dfout=raw_counts_matrix}
  
  return(dfout)
}   
cleaned_counts = CleanRawCounts(input)

colnames(cleaned_counts)[colnames(cleaned_counts) == "X8_hYP7_305_T_resent"] <- 
  "X8_Y_H_305_T_resent"

# filter counts
meta <- read.csv("Ccbr1321 metadata.csv")

Filtered_Counts <- function(CleanRawCounts, Ccbr1321_metadata) {
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(limma)
  library(amap)
  library(colorspace)
  library(dendsort)
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(gplots)
  library(gridExtra)
  library(gridGraphics)
  library(lattice)
  library(magrittr)
  library(plotly)
  library(RColorBrewer)
  library(RCurl)
  library(reshape2)
  library(stringr)
  library(tidyverse)
  library(tibble)
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters:
  counts_matrix = CleanRawCounts
  sample_metadata <- Ccbr1321_metadata   
  gene_names_column <- "Gene"
  
  columns_to_include <- c("X1_CD19_314_T_resent",
                          "X2_Y_Y_319_T_resent",
                          "X3_Y_Y_304_T_resent",
                          "X4_Y_Y_729_T_resent",
                          "X5_Y_Y_727_T_resent",
                          "X6_CD19_738_T_resent","X7_CD19_271_T_resent",
                          "X8_Y_H_305_T_resent","X9_CD19_324_T_resent",
                          "X10_hYP7_731_T_resent","X11_Y_H_314_T_resent",
                          "X12_hYP7_740_T_resent","X13_Y_H_737_T_resent",
                          "X14_hYP7_732_T_resent","X15_Y_H_730_T_resent",
                          "X16_hYP7_303_T_resent")
  
  sample_names_column <- "Sample"
  groups_column <- "Group"
  labels_column <- "Sample"
  
  #Filtering Parameters:
  outlier_samples_to_remove <- c()
  use_cpm_counts_to_filter <- TRUE
  Minimum_Count_Value_to_be_Considered_Nonzero <- 1
  Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total <- 1
  Use_Group_Based_Filtering <- TRUE
  Minimum_Number_of_Samples_with_Nonzero_Counts_in_a_Group <- 3
  
  #PCA Parameters:
  principal_component_on_x_axis<-1 
  principal_component_on_y_axis<-2 
  legend_position_for_PCA <- "top"
  point_size_for_pca<-1
  add_labels_to_PCA <- TRUE
  label_font_size <- 3
  label_offset_y_ <- 2
  label_offset_x_ <- 2
  samples_to_rename_manually <- c("")
  
  #Histogram Parameters:
  color_histogram_by_group <- FALSE     
  set_min_max_for_x_axis_for_histogram <- FALSE
  minimum_for_x_axis_for_histogram <- -1
  maximum_for_x_axis_for_histogram <- 1
  legend_position_for_histogram <- 'none'
  legend_font_size_for_histogram <- 10
  number_of_histogram_legend_columns <- 6
  
  #Visualization Parameters:
  colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  number_of_image_rows <- 2
  interactive_plots <- FALSE
  
  #TCGA:
  plot_correlation_matrix_heatmap <- F
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  getourrandomcolors<-function(k){
    seed=10
    n <- 2e3
    ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
    ourColorSpace <- as(ourColorSpace, "LAB")
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    return( unname(hex(LAB(km$centers))))
  }
  
  make_heatmap <- function(counts_matrix, metadata,colorval) {
    mat <- as.matrix(counts_matrix) 
    tcounts=t(mat)
    tcounts=merge(metadata,tcounts,by.x=sample_names_column,by.y='row.names')
    rownames(tcounts)=tcounts[,labels_column]
    tcounts=tcounts[,!colnames(tcounts)%in%colnames(metadata)]
    d=Dist(tcounts,method="correlation",diag=TRUE)
    dend = rev(dendsort(as.dendrogram(hclust( d,method="average"))))
    m=as.matrix(d)
    sample_metadata <- metadata
    rownames(sample_metadata) = sample_metadata[[labels_column]]
    idx = as.factor(sample_metadata[rownames(m),groups_column])
    col = colorval
    cols <- col[idx]
    new.palette=colorRampPalette(c("blue","green","yellow"),space="rgb")
    
    mk<-function(){
      if(length(colnames(m))>20){
        par(mar=c(0,0,0,0))
        heatmap.2(m,
                  labRow = NA, 
                  labCol = NA,
                  col=new.palette(20),
                  trace="none",
                  colRow = col[idx], 
                  colCol = col[idx],
                  rowDendrogram=dend,
                  colDendrogram=dend,
                  RowSideColors = col[idx],
                  ColSideColors = col[idx],
                  dendrogram = "row",
                  cexRow=3,
                  cexCol=3,
                  margins=c(0,0),   
                  lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                  lhei=c(.2,4,2), 
                  lwid=c(1, .2,4 ), 
                  key.par=list(mgp=c(1.75, .5, 0), 
                               mar=c(7, 2, 3.5, 0), 
                               cex.axis=.1, 
                               cex.lab=3, 
                               cex.main=1, 
                               cex.sub=1),
                  key.xlab = "Correlation",
                  key.ylab="Count",
                  key.title=" ")       
      } else {
        heatmap.2(m,col=new.palette(20),
                  trace="none",
                  colRow = col[idx], 
                  colCol = col[idx],
                  rowDendrogram=dend,
                  colDendrogram=dend,
                  RowSideColors = col[idx],
                  ColSideColors = col[idx],
                  dendrogram = "row",
                  cexRow=3,cexCol=3,margins=c(4,1),  
                  lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                  lhei=c( .2,4,2), 
                  lwid=c(1, .2,4),
                  key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),
                  key.xlab = "Correlation",
                  key.ylab="Count",
                  key.title=" ")
      }
    }
    
    tg<-mk()
    grid.echo(mk)
    gh1<-grid.grab()
    mklegend<-function(){
      plot.new()
      legend(x="top", legend=levels(idx), col=col[as.factor(levels(idx))],pch=15,x.intersp=3,bty ="n",cex=2)
    }
    grid.echo(mklegend )
    gh2<-grid.grab()
    lay <- c(1,3)
    grid.newpage()
    grid.arrange(gh1,gh2,nrow=1,widths=c(unit(1000, "bigpts"),unit(300, "bigpts")))
    gh<-grid.grab()
    return(gh)
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  samples_to_include=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]]
  anno_col=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]==F]
  
  samples_to_include <- samples_to_include[! samples_to_include %in% outlier_samples_to_remove]
  samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
  samples_to_include <- samples_to_include[samples_to_include != "Gene"]
  samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
  samples_to_include <- samples_to_include[samples_to_include %in% sample_metadata[[sample_names_column]]]
  
  ##create unique rownames to correctly add back Annocolumns at end of template
  counts_matrix[,gene_names_column]=paste0(counts_matrix[,gene_names_column],'_',1:nrow(counts_matrix))
  
  anno_col=c(anno_col,gene_names_column)%>%unique
  anno_tbl=counts_matrix[,anno_col,drop=F]%>%as.data.frame
  
  df <- counts_matrix[,c(gene_names_column,samples_to_include)]
  gene_names <- NULL
  gene_names$GeneID <- counts_matrix[,gene_names_column]
  
  #print(colnames(df.final))
  
  ### This code block does input data validation
  
  sample_metadata <- sample_metadata[match(colnames(df),sample_metadata[[sample_names_column]]),] #First match sample metadata to counts matrix
  sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
  sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] #Remove empty columns
  rownames(sample_metadata) <- sample_metadata[[sample_names_column]]
  
  
  ### Remove specal characters from Metadata Column. Replace with _
  sample_metadata[,groups_column]=gsub('-| |!|\\*|\\.',"_",sample_metadata[,groups_column])
  
  
  #### remove low count genes ########
  
  df <- df[complete.cases(df),]
  ## duplicate Rows should be removed in Clean_Raw_Counts template
  #df %>% dplyr::group_by(.data[[gene_names_column]]) %>% summarise_all(sum) %>% as.data.frame() -> df
  print(paste0("Number of features before filtering: ", nrow(df)))
  
  ## USE CPM Transformation
  if (use_cpm_counts_to_filter == TRUE){
    trans.df=df
    trans.df[, -1]=edgeR::cpm(as.matrix(df[, -1]))
    counts_label="Filtered Counts (CPM)"
  } else {
    trans.df=df
    counts_label="Filtered Counts"
    
  }
  
  if (Use_Group_Based_Filtering == TRUE) {
    rownames(trans.df) <- trans.df[,gene_names_column]
    trans.df[,gene_names_column] <- NULL
    
    counts <- trans.df >  Minimum_Count_Value_to_be_Considered_Nonzero # boolean matrix
    
    tcounts <- as.data.frame(t(counts))
    colnum <- dim(counts)[1] # number of genes
    tcounts <- merge(sample_metadata[groups_column], tcounts, by="row.names")
    tcounts$Row.names <- NULL
    melted <- melt(tcounts, id.vars=groups_column)
    tcounts.tot <- dplyr::summarise(dplyr::group_by_at(melted, c(groups_column, "variable")), sum=sum(value))
    tcounts.tot %>% tidyr::spread(variable, sum) -> tcounts.group
    colSums(tcounts.group[(1:colnum+1)]>=Minimum_Number_of_Samples_with_Nonzero_Counts_in_a_Group) >= 1 -> tcounts.keep 
    df.filt <- trans.df[tcounts.keep, ]
    df.filt %>% rownames_to_column(gene_names_column) -> df.filt
  } else {
    
    trans.df$isexpr1 <- rowSums(as.matrix(trans.df[, -1]) > Minimum_Count_Value_to_be_Considered_Nonzero) >= Minimum_Number_of_Samples_with_Nonzero_Counts_in_Total
    
    df.filt <- as.data.frame(trans.df[trans.df$isexpr1, ])
  }
  
  #colnames(df.filt)[colnames(df.filt)==gene_names_column] <- "Gene"
  print(paste0("Number of features after filtering: ", nrow(df.filt)))
  
  ######## Start PCA ###############
  
  edf <- log((as.matrix(df.filt[,samples_to_include]+0.5)))
  rownames(edf) <- df.filt[,1]
  tedf <- t(edf)
  tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  tedf <- tedf[, apply(tedf, 2, var) != 0]
  pca <- prcomp(tedf, scale.=T)
  
  pcx <- paste0("PC",principal_component_on_x_axis)
  pcy <- paste0("PC",principal_component_on_y_axis)
  pca.df <- as.data.frame(pca$x) %>% dplyr::select(.data[[pcx]], .data[[pcy]])
  pca.df$group <- sample_metadata[[groups_column]]
  pca.df$sample <- sample_metadata[[labels_column]]
  perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
  perc.var <- formatC(perc.var,format = "g",digits=4)
  pc.x.lab <- paste0(pcx," ", perc.var[principal_component_on_x_axis],"%")
  pc.y.lab <- paste0(pcy," ", perc.var[principal_component_on_y_axis],"%")
  labelpos <- pca.df
  labelpos$mean_y <- pca.df[[pcy]]+label_offset_y_
  labelpos$mean_x <- pca.df[[pcx]]+label_offset_x_
  pca.df$xdata <- pca.df[[pcx]]
  pca.df$ydata <- pca.df[[pcy]]
  
  # Manual changes to sample names
  replacements = samples_to_rename_manually
  
  if (!is.null(samples_to_rename_manually)) {
    if (replacements != c("")) {
      for (x in replacements) {
        old <- strsplit(x, ": ?")[[1]][1]
        new <- strsplit(x, ": ?")[[1]][2]
        pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
      }
    }
  }
  
  colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
  names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  if(length(colors_for_plots) == 0){
    colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  }
  colorval <- colorlist[colors_for_plots]
  colorval <- unname(colorval) #remove names which affect ggplot
  
  if (length(unique(sample_metadata[[groups_column]])) > length(colorval)) {
    ## Original color-picking code.
    k=length(unique(sample_metadata[[groups_column]]))-length(colorval)
    more_cols<- getourrandomcolors(k) 
    colorval <- c(colorval , more_cols)
  }
  
  # if (add_labels_to_PCA == TRUE){
  #   g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
  #     theme_bw() +
  #     theme(legend.title=element_blank()) +
  #     theme(legend.position=legend_position_for_PCA) +
  #     geom_point(aes(color=group), size=point_size_for_pca) +
  #     geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
  #                                  label=sample, color=group, vjust="inward", hjust="inward"), size=label_font_size, show.legend=FALSE) +
  #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #           panel.background = element_blank()) +
  #     scale_colour_manual(values = colorval) +
  #     xlab(pc.x.lab) + ylab(pc.y.lab)
  # } else {
  #   g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
  #     theme_bw() +
  #     theme(legend.title=element_blank()) +
  #     theme(legend.position=legend_position_for_PCA) +
  #     geom_point(aes(color=group,text=sample), size=point_size_for_pca) +
  #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #           panel.background = element_blank()) +
  #     scale_colour_manual(values = colorval) +
  #     xlab(pc.x.lab) + ylab(pc.y.lab)    
  # }
  
  par(mfrow = c(2,1))
  
  df.m <- melt(edf,id.vars=c(gene_names_column))
  df.m = dplyr::rename(df.m,sample=Var2)
  
  # if(set_min_max_for_x_axis_for_histogram == TRUE){
  #   xmin = minimum_for_x_axis_for_histogram
  #   xmax = maximum_for_x_axis_for_histogram
  # } else {
  #   xmin = min(df.m$value)
  #   xmax = max(df.m$value)
  # }
  # 
  # if(color_histogram_by_group == TRUE){
  #   df.m %>% mutate(colgroup = sample_metadata[sample,groups_column]) -> df.m
  #   df.m = df.m[complete.cases(df.m[, "colgroup"]),]
  #   df.m$colgroup = gsub("\\s","_",df.m$colgroup)
  #   df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
  #   #print(unique(df.m$sample))
  #   
  #   # plot Density 
  #   g2 = ggplot(df.m, aes(x=value, group=sample)) + 
  #     geom_density(aes(colour = colgroup)) +
  #     xlab(counts_label) + ylab("Density") +
  #     theme_bw() +
  #     theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) + 
  #     ggtitle("Frequency Histogram") +
  #     xlim(xmin,xmax) +
  #     #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
  #     scale_colour_manual(values=colorval)
  #   guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
  # } else {
  #   
  #   df.m$sample = sample_metadata[df.m$sample,labels_column]
  #   n=length(unique(df.m$sample))
  #   cols<- getourrandomcolors(n) 
  #   
  #   g2 = ggplot(df.m, aes(x=value, group=sample)) + 
  #     geom_density(aes(colour = sample )) +
  #     xlab(counts_label) + ylab("Density") +
  #     theme_bw() +
  #     theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) +  
  #     ggtitle("Frequency Histogram") +
  #     xlim(xmin,xmax) +
  #     #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
  #     scale_colour_manual(values=cols) +
  #     guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
  # }
  
  #dev.off()
  
  imageWidth = 3000
  imageHeight = 1500*2
  dpi = 300
  
  # if(plot_correlation_matrix_heatmap == TRUE){
  #   if(interactive_plots == TRUE){
  #     p1=(g)%>%ggplotly(tooltip = c("sample","group"))
  #     p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = c("sample"))
  #     fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
  #     fig=fig %>% layout(title = 'Interactive PCA and Histogram')
  #     print(fig)
  #   } else {
  #     require(gridExtra)
  #     gh<-make_heatmap(df.filt[,samples_to_include],sample_metadata,colorval)
  #     grid.arrange(g,g2,gh, nrow=number_of_image_rows)
  #     #dev.off()
  #   }  
  # } else {
  #   if(interactive_plots == TRUE){
  #     p1=(g)%>%ggplotly(tooltip = c("sample","group"))
  #     p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = "sample" )
  #     fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
  #     fig=fig %>% layout(title = 'Interactive PCA and Histogram')
  #     print(fig)
  #   } else {
  #     grid.arrange(g,g2, nrow=number_of_image_rows)
  #     #dev.off()
  #   }
  # }
  
  df %>% filter(.data[[gene_names_column]] %in% df.filt[,gene_names_column]) -> df.final
  # colnames(df.final)[colnames(df.final)==gene_names_column] <- "Gene"
  
  # print('')
  # print('Sample Columns')
  # print(colnames(df.final[,!colnames(df.final)%in%gene_names_column]))
  # print('Annotation Columns')
  # print(colnames(anno_tbl))
  
  df.final=merge(anno_tbl,df.final,by=gene_names_column,all.y=T)
  df.final[,gene_names_column]=gsub('_[0-9]+$',"",df.final[,gene_names_column])
  
  return(df.final)
}
filtered_counts <- Filtered_Counts(CleanRawCounts = cleaned_counts, 
                                   Ccbr1321_metadata = meta)

# Normalization [CCBR] (afc2524c-9bae-4873-98c1-e06ef8f4632b): v260
NormalizedCounts <- function(Filtered_Counts, Ccbr1321_metadata) {
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(limma)
  library(tidyverse)
  library(edgeR)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(RColorBrewer)
  library(colorspace)
  library(stringr)
  library(RCurl)
  library(reshape2)
  library(gridExtra)
  library(amap)
  library(lattice)
  library(gplots)
  #library(gridGraphics)
  library(dendsort)
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters:
  counts_matrix <- Filtered_Counts
  sample_metadata <- Ccbr1321_metadata
  gene_names_column <- "Gene"
  columns_to_include = c("X1_CD19_314_T_resent","X2_Y_Y_319_T_resent",
                         "X3_Y_Y_304_T_resent","X4_Y_Y_729_T_resent",
                         "X5_Y_Y_727_T_resent","X6_CD19_738_T_resent",
                         "X7_CD19_271_T_resent","X8_Y_H_305_T_resent",
                         "X9_CD19_324_T_resent","X10_hYP7_731_T_resent",
                         "X11_Y_H_314_T_resent","X12_hYP7_740_T_resent",
                         "X13_Y_H_737_T_resent","X14_hYP7_732_T_resent",
                         "X15_Y_H_730_T_resent","X16_hYP7_303_T_resent")
  sample_names_column <- "Sample"
  groups_column <- "Group"
  labels_column <- "Sample"
  
  #Normalization Parameters:
  input_in_log_counts <- FALSE
  normalization_method <- "quantile"
  
  #PCA parameters:
  samples_to_rename_manually_on_pca <- c("")
  add_labels_to_pca <- TRUE
  principal_component_on_x_axis_for_pca <- 1
  principal_component_on_y_axis_for_pca <- 2
  legend_position_for_pca <- "top"
  label_offset_x_for_pca <- 2
  label_offset_y_for_pca <- 2
  label_font_size_for_pca <- 3
  point_size_for_pca <- 2
  
  #Histogram parameters:
  color_histogram_by_group <- FALSE
  set_min_max_for_x_axis_for_histogram <- FALSE 
  minimum_for_x_axis_for_histogram <- -1
  maximum_for_x_axis_for_histogram <- 1
  legend_font_size_for_histogram <- 10
  legend_position_for_histogram <- "none"
  number_of_histogram_legend_columns <- 6
  
  #Visualization Parameters:
  number_of_image_rows <- 2
  colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  make_plots_interactive <- FALSE
  plot_correlation_matrix_heatmap <- TRUE
  
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  getourrandomcolors<-function(k){
    seed=10
    n <- 2e3
    ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
    ourColorSpace <- as(ourColorSpace, "LAB")
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    return( unname(hex(LAB(km$centers))))
  }
  
  make_heatmap <- function(counts_matrix, metadata,colorval) {
    mat <- as.matrix(counts_matrix) 
    tcounts=t(mat)
    tcounts=merge(metadata,tcounts,by.x=sample_names_column,by.y='row.names')
    rownames(tcounts)=tcounts[,labels_column]
    tcounts=tcounts[,!colnames(tcounts)%in%colnames(metadata)]
    d=Dist(tcounts,method="correlation",diag=TRUE)
    dend = rev(dendsort(as.dendrogram(hclust( d,method="average"))))
    m=as.matrix(d)
    sample_metadata <- metadata
    rownames(sample_metadata) = sample_metadata[[labels_column]]
    idx = as.factor(sample_metadata[rownames(m),groups_column])
    col = colorval
    cols <- col[idx]
    new.palette=colorRampPalette(c("blue","green","yellow"),space="rgb")
    
    mk<-function(){
      if(length(colnames(m))>20){
        par(mar=c(0,0,0,0))
        heatmap.2(m,
                  labRow = NA, 
                  labCol = NA,
                  col=new.palette(20),
                  trace="none",
                  colRow = col[idx], 
                  colCol = col[idx],
                  rowDendrogram=dend,
                  colDendrogram=dend,
                  RowSideColors = col[idx],
                  ColSideColors = col[idx],
                  dendrogram = "row",
                  cexRow=3,
                  cexCol=3,
                  margins=c(0,0),   
                  lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                  lhei=c(.2,4,2), 
                  lwid=c(1, .2,4 ), 
                  key.par=list(mgp=c(1.75, .5, 0), 
                               mar=c(7, 2, 3.5, 0), 
                               cex.axis=.1, 
                               cex.lab=3, 
                               cex.main=1, 
                               cex.sub=1),
                  key.xlab = "Correlation",
                  key.ylab="Count",
                  key.title=" ")       
      } else {
        heatmap.2(m,col=new.palette(20),
                  trace="none",
                  colRow = col[idx], 
                  colCol = col[idx],
                  rowDendrogram=dend,
                  colDendrogram=dend,
                  RowSideColors = col[idx],
                  ColSideColors = col[idx],
                  dendrogram = "row",
                  cexRow=3,cexCol=3,margins=c(4,1),  
                  lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                  lhei=c( .2,4,2), 
                  lwid=c(1, .2,4),
                  key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),
                  key.xlab = "Correlation",
                  key.ylab="Count",
                  key.title=" ")
      }
    }
    
    tg<-mk()
    grid.echo(mk)
    gh1<-grid.grab()
    mklegend<-function(){
      plot.new()
      legend(x="top", legend=levels(idx), col=col[as.factor(levels(idx))],pch=15,x.intersp=3,bty ="n",cex=2)
    }
    grid.echo(mklegend )
    gh2<-grid.grab()
    lay <- c(1,3)
    grid.newpage()
    grid.arrange(gh1,gh2,nrow=1,widths=c(unit(1000, "bigpts"),unit(300, "bigpts")))
    gh<-grid.grab()
    return(gh)
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  samples_to_include=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]]
  anno_col=columns_to_include[columns_to_include%in%sample_metadata[,sample_names_column,drop=T]==F]
  
  samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
  samples_to_include <- samples_to_include[samples_to_include != "Gene"]
  samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
  
  ##create unique rownames to correctly add back Annocolumns at end of template
  counts_matrix[,gene_names_column]=paste0(counts_matrix[,gene_names_column],'_',1:nrow(counts_matrix))
  
  anno_col=c(anno_col,gene_names_column)%>%unique
  anno_tbl=counts_matrix[,anno_col,drop=F]%>%as.data.frame
  
  
  df.filt <- counts_matrix[,samples_to_include]
  gene_names <- NULL
  gene_names$GeneID <- counts_matrix[,1]
  
  
  
  sample_metadata <- sample_metadata[match(colnames(df.filt),sample_metadata[[sample_names_column]]),] #First match sample metadata to counts matrix
  sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
  sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] #Remove empty columns
  rownames(sample_metadata) <- sample_metadata[[sample_names_column]]
  
  df.filt <- df.filt[,match(sample_metadata[[sample_names_column]],colnames(df.filt))] #Match counts matrix columns to sample metadata
  
  #If input is in log space, linearize
  if(input_in_log_counts == TRUE){
    x <- DGEList(counts=2^df.filt, genes=gene_names)
  } else {
    x <- DGEList(counts=df.filt, genes=gene_names)     
  }
  
  v <- voom(x,normalize=normalization_method)
  rownames(v$E) <- v$genes$GeneID
  as.data.frame(v$E) %>% rownames_to_column(gene_names_column) -> df.voom
  print(paste0("Total number of features included: ", nrow(df.voom)))
  
  #Start PCA Plot:
  
  edf <- v$E
  tedf <- t(edf)
  tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  tedf <- tedf[, apply(tedf, 2, var) != 0]
  pca <- prcomp(tedf, scale.=T)
  
  pcx <- paste0("PC",principal_component_on_x_axis_for_pca)
  pcy <- paste0("PC",principal_component_on_y_axis_for_pca)
  pca.df <- as.data.frame(pca$x) %>% dplyr::select(.data[[pcx]], .data[[pcy]])
  pca.df$group <- sample_metadata[[groups_column]]
  pca.df$sample <- sample_metadata[[labels_column]]
  perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
  perc.var <- formatC(perc.var,format = "g",digits=4)
  pc.x.lab <- paste0(pcx," ", perc.var[principal_component_on_x_axis_for_pca],"%")
  pc.y.lab <- paste0(pcy," ", perc.var[principal_component_on_y_axis_for_pca],"%")
  labelpos <- pca.df
  labelpos$mean_y <- pca.df[[pcy]]+label_offset_y_for_pca
  labelpos$mean_x <- pca.df[[pcx]]+label_offset_x_for_pca
  pca.df$xdata <- pca.df[[pcx]]
  pca.df$ydata <- pca.df[[pcy]]
  
  # Manual changes to sample names
  replacements = samples_to_rename_manually_on_pca
  
  if (!is.null(replacements)) {
    if (replacements != c("")) {
      for (x in replacements) {
        old <- strsplit(x, ": ?")[[1]][1]
        new <- strsplit(x, ": ?")[[1]][2]
        pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
      }
    }
  }
  
  colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
  names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  if(length(colors_for_plots) == 0){
    colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  }
  colorval <- colorlist[colors_for_plots]
  colorval <- unname(colorval) #remove names which affect ggplot
  
  if (length(unique(sample_metadata[[groups_column]])) > length(colorval)) {
    ## Original color-picking code.
    k=length(unique(sample_metadata[[groups_column]]))-length(colorval)
    more_cols<- getourrandomcolors(k) 
    colorval <- c(colorval , more_cols)
  }
  
  if (add_labels_to_pca == TRUE){
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position=legend_position_for_pca) +
      geom_point(aes(color=group), size=point_size_for_pca) +
      geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
                                   label=sample, color=group, vjust="inward", hjust="inward"), size=label_font_size_for_pca, show.legend=FALSE) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      scale_colour_manual(values = colorval) +
      xlab(pc.x.lab) + ylab(pc.y.lab)
  } else {
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      theme(legend.position=legend_position_for_pca) +
      geom_point(aes(color=group,text=sample), size=point_size_for_pca) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      scale_colour_manual(values = colorval) +
      xlab(pc.x.lab) + ylab(pc.y.lab)    
  }
  
  par(mfrow = c(2,1))
  
  #Histogram Plot:
  
  df.m <- melt(edf,id.vars=c(gene_names_column))
  df.m = dplyr::rename(df.m,sample=Var2)
  
  if(set_min_max_for_x_axis_for_histogram == TRUE){
    xmin = minimum_for_x_axis_for_histogram
    xmax = maximum_for_x_axis_for_histogram
  } else {
    xmin = min(df.m$value)
    xmax = max(df.m$value)
  }
  
  if(color_histogram_by_group == TRUE){
    df.m %>% mutate(colgroup = sample_metadata[sample,groups_column]) -> df.m
    df.m = df.m[complete.cases(df.m[, "colgroup"]),]
    df.m$colgroup = gsub("\\s","_",df.m$colgroup)
    df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
    #print(unique(df.m$sample))
    
    # plot Density 
    g2 = ggplot(df.m, aes(x=value, group=sample)) + 
      geom_density(aes(colour = colgroup)) +
      xlab("Filtered Counts") + ylab("Density") +
      theme_bw() +
      theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) + 
      ggtitle("Frequency Histogram") +
      xlim(xmin,xmax) +
      #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
      scale_colour_manual(values=colorval)
  } else {
    
    df.m$sample = sample_metadata[df.m$sample,labels_column]
    n=length(unique(df.m$sample))
    cols<- getourrandomcolors(n) 
    
    g2 = ggplot(df.m, aes(x=value, group=sample)) + 
      geom_density(aes(colour = sample)) +
      xlab("Filtered Counts") + ylab("Density") +
      theme_bw() +
      theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) +  
      ggtitle("Frequency Histogram") +
      xlim(xmin,xmax) +
      #scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
      scale_colour_manual(values=cols)#+
    guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
  }
  
  # dev.off()
  
  imageWidth = 3000
  imageHeight = 1500*2
  dpi = 300
  
  # if(plot_correlation_matrix_heatmap == TRUE){
  #   if(make_plots_interactive == TRUE){
  #     p1=(g)%>%ggplotly(tooltip = c("sample","group"))
  #     p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = c("sample"))
  #     fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
  #     fig=fig %>% layout(title = 'Interactive PCA and Histogram')
  #     print(fig)
  #   } else {
  #     require(gridExtra)
  #     gh<-make_heatmap(df.filt,sample_metadata,colorval)
  #     grid.arrange(g,g2,gh, nrow=number_of_image_rows)
  #     # dev.off()
  #   }  
  # } else {
  #   if(make_plots_interactive == TRUE){
  #     p1=(g)%>%ggplotly(tooltip = c("sample","group"))
  #     p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = "sample" )
  #     fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
  #     fig=fig %>% layout(title = 'Interactive PCA and Histogram')
  #     print(fig)
  #   } else {
  #     grid.arrange(g,g2, nrow=number_of_image_rows)
  #     # dev.off()
  #   }
  # }    
  
  print("Sample columns")
  print(colnames(df.voom)[!colnames(df.voom)%in%gene_names_column])
  print("Feature Columns")
  print(colnames(anno_tbl))
  
  df.voom=merge(anno_tbl,df.voom,by=gene_names_column,all.y=T)
  df.voom[,gene_names_column]=gsub('_[0-9]+$',"",df.voom[,gene_names_column])
  
  return(df.voom)
}
normalized_counts_input <- NormalizedCounts(Filtered_Counts = filtered_counts, Ccbr1321_metadata = meta)

# Gene Boxplot with Statistics [CCBR] [Beta] (b76846e8-1bfe-4027-8088-11b1b8d3b0cc): v35
Tcm_box <- function( NormalizedCounts, Ccbr1321_metadata) {
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  current_packages <- c("ggplot2","tidyr","dplyr","tibble","magrittr",
                        "reshape2","ggbeeswarm","RColorBrewer","stringr",
                        "l2p","l2psupp")
  lapply(current_packages, function(x) library(x, character.only = T, quietly = T))    
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  # Basic parameters
  normalized_counts = NormalizedCounts
  sample_metadata = Ccbr1321_metadata
  gene_column = "Gene"
  sample_column = "Sample"
  #genes = c("MAP2K6","UHRF1","PPP2R5D","EIF2B3")
  genes = c("GLUT1","SLC2A1","HK2","LDHA")
  run_anova = TRUE
  category_column = "Group" 
  categories = c("Y_Y","Y_H","hYP7_CAR_T","CD19")
  category_labels = c("Y_Y","Y_H","hYP7_CAR_T","CD19")
  pval_column_suffix = "_adjpval"
  
  # Visualization parameters
  draw_jitterplot = TRUE
  plot_type = "Box plot"
  add_annotations = TRUE
  x_axis_title = "none"
  x_axis_title_size = 20
  y_axis_title = "none"
  y_axis_title_size = 20
  y_axis_text_size = 20
  colors_to_use = c("Deep Red","Vivid Blue","Green","Purple","Bright Orange","Yellow","Brown","Pink","Grey","Burnt Orange","Teal Green","Soft Purple","Hot Pink","Leaf Green","Mustard Yellow")
  shapes_to_use = c("filled circle")
  significance_colors_to_use = c("Black")
  plot_width = 0.3
  dot_size = 2
  title = "auto"
  title_position = "3"
  title_size = 20
  title_face = "italic"
  legend_text_size = 20
  legend_text_face = "plain"
  legend_position = "right"
  legend_title = "auto"
  legend_title_size = 20
  margin_left = 6
  margin_right = 6
  margin_top = 1
  margin_bottom = 1
  
  # Jitter Plot parameters
  jitter_width = 0.10
  jitter_height = 0.1
  
  # Beeswarm parameters
  dodge_width = 0.75
  beeswarm_spread = 3
  beeswarm_method = "density"
  
  # Advanced parameters
  sum_duplicates = TRUE
  
  ##---------------------------------- ##
  ## Parameter Misspecification Errors ##
  ## --------------------------------- ##
  
  if(length(colors_to_use) != 1 && length(categories) > length(colors_to_use)){
    stop("ERROR: Need more colors to accomodate more categories. Other choice is to choose 1 color to use for all categories")
  }
  
  if(length(shapes_to_use) != 1 && length(categories) > length(shapes_to_use)){
    stop("ERROR: Need more shapes to accomodate more categories. Other choice is to choose 1 shape to use for all categories")
  }
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  asterisks <- function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("")
  }
  
  anncolor <- function(p) {
    if (p < 0.001) return(significance_colors[1])
    if (p < 0.01) return(significance_colors[2])
    if (p < 0.05) return(significance_colors[3])
    return("")
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  genes <- updategenes(genes)
  
  # Extract old and new gene names
  old_gene_names <- genes$oldname
  new_gene_names <- genes$newname
  
  # Compare old and new gene names, and print the new ones if they differ
  for (i in seq_along(new_gene_names)) {
    if (new_gene_names[i] != old_gene_names[i]) {
      cat("Old name:", old_gene_names[i], "-> New name:", new_gene_names[i], "\n")
    }
  }
  
  genes <- genes$newname
  newgenes <-  updategenes(normalized_counts[[gene_column]])
  normalized_counts[[gene_column]] <- newgenes$newname
  
  ##Clean up the normalized counts data:
  numeric_columns <- names(normalized_counts)[sapply(normalized_counts, is.numeric)]
  normcounts_filt <- normalized_counts %>% filter(.data[[gene_column]] != "---") %>% select(c(gene_column,numeric_columns))
  if(sum_duplicates == TRUE){
    normcounts_filt <- normcounts_filt %>%
      group_by(.data[[gene_column]]) %>%
      summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
      ungroup()
  } else {
    normcounts_max <- rowSums(as.matrix(normcounts_filt[,-1]))
    normcounts_filt$max <- normcounts_max 
    normcounts_filt <- normcounts_filt %>% group_by(.data[[gene_column]]) %>% filter(max == max(max))
    normcounts_filt$max <- NULL
  }
  normcounts_filt <- as.data.frame(normcounts_filt)
  anova.list <- list()
  
  sample_metadata <- sample_metadata %>% filter(.data[[category_column]] %in% categories)
  sample_metadata <- sample_metadata %>% filter(.data[[sample_column]] %in% colnames(normcounts_filt))
  
  categ_filt <- unique(sample_metadata[[category_column]])
  if(length(categ_filt) < 2){
    stop(sprintf("error: Need more than 1 category to compare: %s", categ_filt))    
  }
  cnames <- colnames(normcounts_filt)[colnames(normcounts_filt) %in% sample_metadata[[sample_column]]]
  normcounts_filt <- normcounts_filt %>% 
    select(c(gene_column,cnames)) %>%
    filter(.data[[gene_column]] %in% genes)
  
  #Set up colors and shapes:
  colorlist <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999","#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#F0027F", "#8DD3C7","#000000")
  names(colorlist) <- c("Deep Red","Vivid Blue","Green","Purple","Bright Orange","Yellow","Brown","Pink","Grey","Burnt Orange","Teal Green","Soft Purple","Hot Pink","Leaf Green","Mustard Yellow","Bronze","Dark Gray","Bright Magenta","Light Aqua","Black")
  
  shapelist <- seq(0,25,1)
  names(shapelist) <- c("square","circle","triangle point up","plus","cross","diamond","triangle point down","square cross","star","diamond plus","circle plus","triangles up and down","square plus","circle cross","square and triangle down","filled square","filled circle","filled triangle point-up","filled diamond","solid circle","bullet (smaller circle)","filled circle blue","filled square blue","filled diamond blue","filled triangle point-up blue","filled triangle point down blue")
  
  if(length(colors_to_use) == 1){
    colors <- rep(colorlist[colors_to_use], length(categ_filt))
  } else {
    colors <- colorlist[colors_to_use]
  }
  
  if(length(shapes_to_use) == 1){
    shapes <- rep(shapelist[shapes_to_use], length(categ_filt))
  } else {
    shapes <- shapelist[shapes_to_use]
  }
  
  if(length(significance_colors_to_use) == 1){
    significance_colors <- rep(colorlist[significance_colors_to_use], 3)
  } else {
    significance_colors <- colorlist[significance_colors_to_use]
  }
  
  #Draw Plots:
  set.seed(42)
  for(g in genes){
    tryCatch({
      if(g %in% normcounts_filt[[gene_column]]){
        counts.filt <- normcounts_filt %>% filter(.data[[gene_column]] == g) %>%  select(-1) %>% as.numeric()  
        names(counts.filt) <- colnames(normcounts_filt[,-1])
        counts.filt <- counts.filt[sample_metadata[[sample_column]]]
        rownames(sample_metadata) = sample_metadata[[sample_column]]
        group = sample_metadata[[category_column]][match(names(counts.filt),rownames(sample_metadata))]  # Add Group here
        df.m <- data.frame(value = counts.filt, group = group)
        #df.m$group <- factor(gsub("_"," ",df.m$group), levels = gsub("_"," ",categories))
        df.m$group <- factor(df.m$group, levels = categories)
        
        if(run_anova == TRUE){
          anova <- aov(df.m$value ~ df.m$group)
          anova.df <- as.data.frame(TukeyHSD(anova,ordered=TRUE)$`df.m$group`)
          anova.df %>% rownames_to_column("GroupComp") -> anova.df
          colnames(anova.df) <- str_replace(colnames(anova.df), " ", "_")
          anova.df %>% arrange(p_adj) -> anova.df
          stats.tab <- anova.df %>% mutate(diff = case_when(diff < 1 ~ -1/diff , TRUE ~ diff))
          anova.gene.df <- cbind(gene = g, anova.df)
          anova.list[[g]] <- anova.gene.df
        } else {
          pvalcols <- colnames(normalized_counts)[grepl(pval_column_suffix,colnames(normalized_counts))]
          fchangecols <- colnames(normalized_counts)[grepl("_FC",colnames(normalized_counts))]
          pvals <- normalized_counts %>% filter(.data[[gene_column]] == g) %>%  select(pvalcols) %>% as.numeric()
          pvalcols <- gsub(pval_column_suffix, "", pvalcols)
          FC <- normalized_counts %>% filter(.data[[gene_column]] == g) %>%  select(fchangecols) %>% as.numeric()
          fchangecols <- str_remove(fchangecols, "_FC")
          stats.tab <- data.frame("GroupComp"=pvalcols,
                                  "diff" = FC,
                                  "p_adj"=pvals)
          stats.tab <- stats.tab %>% filter(str_count(GroupComp, "-") <= 2)
          stats.tab <- stats.tab %>% arrange(p_adj)
        }
        #Set up colors:
        names(colors) <- levels(df.m$group)
        df.m$colors <- factor(colors[df.m$group])
        
        #Set up shapes:
        names(shapes) <- levels(df.m$group)
        df.m$shapes <- factor(shapes[df.m$group])
        
        #Set up category labels:
        names(category_labels) <- levels(df.m$group)
        df.m$category_labels <- factor(category_labels[df.m$group])
        
        gp <- ggplot(df.m, aes(x = group, 
                               y = value, 
                               shape = shapes, 
                               color = group))  # Use group for color
        
        if (draw_jitterplot == TRUE) {
          gp <- gp + geom_jitter(position = position_jitter(width = jitter_width,
                                                            height = jitter_height), 
                                 size = dot_size, 
                                 show.legend = c(shape = FALSE, 
                                                 color = TRUE)) +
            theme_bw() +
            scale_color_manual(values = colors, 
                               labels = category_labels,
                               name = NULL)+
            guides(shape = "none") 
        } else {
          gp <- gp + geom_beeswarm(cex = beeswarm_spread, 
                                   size = dot_size, 
                                   dodge.width = dodge_width, 
                                   priority = beeswarm_method,
                                   show.legend = c(shape = FALSE, 
                                                   color = TRUE)) +
            theme_bw() +
            scale_color_manual(values = colors, 
                               labels = category_labels,
                               name = NULL) +
            guides(shape = "none")
        }
        
        # Update legend and labels
        if(y_axis_title == "auto"){
          y_axis_title_value <- g
        } else {
          y_axis_title_value <- ""
        }
        
        if(x_axis_title == "auto"){
          x_axis_title_value <- category_column
        } else {
          x_axis_title_value <- ""
        }
        
        if(title == "auto"){
          title_value <- g
        } else {
          title_value <- ""
        }
        
        gp <- gp + labs(title = title_value, 
                        x = x_axis_title_value, 
                        y = y_axis_title_value) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.title.y = element_text(size = y_axis_title_size, vjust = 5),
            axis.text.y = element_text(size = y_axis_text_size), 
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5,
                                      vjust = title_position,
                                      face = title_face, 
                                      size = title_size), 
            legend.key.size = unit(2, "line"),
            legend.text = element_text(size = legend_text_size,
                                       face = legend_text_face),
            legend.position = legend_position,
            legend.title = element_blank()) +
          scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
          expand_limits(y = c(min(df.m$value - 1, na.rm = TRUE), max(df.m$value + 1, na.rm = TRUE))) 
        
        # Add boxplot if enabled
        if (plot_type == "Box plot") {      
          gp <- gp + geom_boxplot(outlier.shape = NA, 
                                  alpha = 0.3, 
                                  width = plot_width,
                                  position = position_dodge(width = 0.75))
          base_y_position <- max(df.m$value) + 0.2 * max(df.m$value)
          spacing <- 0.2 * max(df.m$value) # Adjust as needed for spacing
        } else if (plot_type == "Violin plot") { 
          gp <- gp + geom_violin(trim = FALSE, 
                                 alpha = 0.3, 
                                 width = plot_width,
                                 position = position_dodge(width = 0.75),
                                 scale = "width")
          # Build the ggplot object to access the computed data
          gp_build <- ggplot_build(gp)
          
          # Extract the data used for the violin plot
          violin_data <- gp_build$data[[3]] 
          
          # Find the maximum y-value in the computed violin data
          base_y_position <- max(violin_data$y, na.rm = TRUE) 
          spacing <- 0.2 * base_y_position # Adjust as needed for spacing
          base_y_position <- base_y_position + 0.2 * base_y_position    
        }
        
        #Set the margins:
        gp <- gp + theme(plot.margin = margin(t = margin_top, 
                                              r = margin_right, 
                                              b = margin_bottom, 
                                              l = margin_left, 
                                              unit = "cm"))
        
        if(add_annotations == TRUE){
          # Base y_position for annotations
          
          for(i in 1:nrow(stats.tab)){
            if(stats.tab$p_adj[i] < 0.05){
              current_y_position <- base_y_position + (i-1) * spacing
              grp <- unlist(strsplit(stats.tab$GroupComp[i], "-"))
              if(run_anova == FALSE){ #When sample metadata column needs renaming to match the DEG column names
                grp <- names(category_labels)[match(grp, names(category_labels))]
              }
              xpos1 <- which(levels(df.m$group) == grp[1])
              xpos2 <- which(levels(df.m$group) == grp[2])
              xpos_mid <- mean(c(xpos1, xpos2))
              group <- factor(grp[1], levels=categories)
              
              current_y_position <- base_y_position + (i-1) * spacing
              anno_data <- data.frame(xstart = xpos1, 
                                      xend = xpos2, 
                                      y = current_y_position, 
                                      label = asterisks(stats.tab$p_adj[i]),
                                      colour = anncolor(stats.tab$p_adj[i]),
                                      group = group)
              
              add_gp_annotations <- function(gplot){ 
                plot_build <- ggplot_build(gplot)
                
                # Extract the y-axis limits (ylim)
                ylim_values <- plot_build$layout$panel_params[[1]]$y.range
                y_text_pos <- ylim_values[2]*0.05 + anno_data$y
                
                gplot <- gplot + geom_segment(data = anno_data, 
                                              aes(x = xstart, 
                                                  xend = xend, 
                                                  y = y, 
                                                  yend = y),  
                                              linewidth = 1, 
                                              inherit.aes = FALSE,
                                              show.legend = FALSE) +
                  geom_text(data = anno_data, 
                            aes(x = (xstart + xend)/2, 
                                y = y_text_pos, 
                                label = label), 
                            size = 10,
                            inherit.aes = FALSE,
                            show.legend = FALSE) + 
                  guides(fill = guide_legend(reverse = FALSE))
                return(gplot)
              }
              gp <- add_gp_annotations(gp)
            }
          }
        }
        
        print(gp)
      } else {
        cat(sprintf("Gene %s is not in dataset \n", g))
        next
      }}, error = function(msg){
        print(paste0(g," has zero expression across samples and cannot be analyzed"))
      })    
    
    ggsave(g,gp)
  }
  
  combined_matrix <- do.call(rbind, anova.list)
  print(combined_matrix)
  
  results <- list(normcounts_filt, combined_matrix)
  
  return(results)
}
boxplot <- Tcm_box(NormalizedCounts = normalized_counts_input, Ccbr1321_metadata = meta)

anova_boxplot_results <- boxplot
rownames(anova_boxplot_results) <- NULL

write.csv(anova_boxplot_results, "pptx3_anova_boxplot_results.csv")

# DEG
DEGAnalysis <- function(Filtered_Counts, Ccbr1321_metadata) {
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(limma)
  library(tidyverse)
  library(edgeR)
  library(stringr)
  library(grid)
  library(gridExtra)
  
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters:
  counts_matrix <- Filtered_Counts 
  sample_metadata <- Ccbr1321_metadata
  gene_names_column="Gene"
  sample_name_column<-"Sample"
  columns_to_include = c("X1_CD19_314_T_resent",
                         "X2_Y_Y_319_T_resent",
                         "X3_Y_Y_304_T_resent",
                         "X4_Y_Y_729_T_resent",
                         "X5_Y_Y_727_T_resent",
                         "X6_CD19_738_T_resent","X7_CD19_271_T_resent",
                         "X8_Y_H_305_T_resent","X9_CD19_324_T_resent",
                         "X10_hYP7_731_T_resent","X11_Y_H_314_T_resent",
                         "X12_hYP7_740_T_resent","X13_Y_H_737_T_resent",
                         "X14_hYP7_732_T_resent","X15_Y_H_730_T_resent",
                         "X16_hYP7_303_T_resent")
  contrast_variable_column<-c("Group")
  contrasts<-c("hYP7_CAR_T-CD19","Y_Y-hYP7_CAR_T","Y_Y-CD19",
               "Y_H-CD19","Y_Y-Y_H")
  covariates_columns=c("Group")
  
  #Advanced Parameters:
  input_in_log_counts <- FALSE
  return_mean_and_sd<-FALSE
  return_normalized_counts<-TRUE
  normalization_method<-"quantile"
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  if(make.names(colnames(counts_matrix))!=colnames(counts_matrix)){
    print("Error: The following counts matrix column names are not valid:\n")
    print(colnames(counts_matrix)[make.names(colnames(counts_matrix))!=colnames(counts_matrix)])
    
    print("Likely causes are columns starting with numbers or other special characters eg spaces.")
    stop("Bad column names.")
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  samples_to_include=columns_to_include[columns_to_include%in%sample_metadata[,sample_name_column,drop=T]]
  anno_col=columns_to_include[columns_to_include%in%sample_metadata[,sample_name_column,drop=T]==F]
  
  samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
  samples_to_include <- samples_to_include[samples_to_include != "Gene"]
  samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
  
  
  ##create unique rownames to correctly add back Annocolumns at end of template
  counts_matrix[,gene_names_column]=paste0(counts_matrix[,gene_names_column],'_',1:nrow(counts_matrix))
  
  anno_col=c(anno_col,gene_names_column)%>%unique
  anno_tbl=counts_matrix[,anno_col,drop=F]%>%as.data.frame
  
  df.m <- counts_matrix[,c(gene_names_column,samples_to_include)]
  gene_names <- NULL
  gene_names$GeneID <- counts_matrix[,gene_names_column]
  
  ### This code block does input data validation
  sample_metadata <- sample_metadata[match(colnames(df.m),sample_metadata[,sample_name_column]),]
  sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ]
  df.m <- df.m[,match(sample_metadata[,sample_name_column],colnames(df.m))]
  
  #Create DGEList object from counts
  if(input_in_log_counts == TRUE){
    x <- DGEList(counts=2^df.m, genes=gene_names)
  } else {
    x <- DGEList(counts=df.m, genes=gene_names) 
  }
  
  #Put covariates in order 
  covariates_columns=covariates_columns[order(covariates_columns!=contrast_variable_column)]
  
  for(ocv in covariates_columns){
    sample_metadata[,ocv]=gsub(" ","_",sample_metadata[,ocv])
  }
  
  contrasts=gsub(" ","_",contrasts)
  cov <- covariates_columns[!covariates_columns %in% contrast_variable_column]
  
  #Combine columns if 2-factor analysis
  if(length(contrast_variable_column)>1){
    sample_metadata %>% dplyr::mutate(contmerge = paste0(.data[[contrast_variable_column[1]]],".",.data[[contrast_variable_column[2]]])) -> sample_metadata
  } else {
    sample_metadata %>% dplyr::mutate(contmerge = .data[[contrast_variable_column]]) -> sample_metadata
  }
  
  contrast_var <- factor(sample_metadata$contmerge)
  
  if(length(cov) >0){
    dm.formula <- as.formula(paste("~0 +", paste("contmerge", paste(cov, sep="+", collapse="+"),sep="+")))
    design=model.matrix(dm.formula, sample_metadata)
    colnames(design) <- gsub("contmerge","",colnames(design))
  } else {
    dm.formula <- as.formula(~0 + contmerge)
    design=model.matrix(dm.formula, sample_metadata)
    colnames(design) <- levels(contrast_var)
  }
  
  #colnames(design) <- str_replace_all(colnames(design), contrast_variable_column, "")
  
  if (normalization_method %in% c("TMM","TMMwzp","RLE","upperquartile")){
    x <- calcNormFactors(x, method = normalization_method) 
    rownames(x) <- x$genes$GeneID
    v <- voom(x,design=design,normalize="none")
  } else {
    v <- voom(x,design=design,normalize=normalization_method,save.plot = TRUE)
  }
  
  rownames(v$E) <- v$genes$GeneID
  as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
  fit <- lmFit(v, design)
  cm <- makeContrasts(contrasts = contrasts, levels=design)
  
  #Print Mean-variance Plot
  sx <- v$voom.xy$x
  sy <- v$voom.xy$y
  xyplot <- as.data.frame(cbind(sx,sy))
  voomline <- as.data.frame(cbind(x=v$voom.line$x,y=v$voom.line$y))
  
  g <- ggplot() +
    geom_point(data=xyplot, aes(x=sx,y=sy),size=1) +
    theme_bw() +
    geom_smooth(data=voomline, aes(x=x,y=y),color = "red") +
    ggtitle("voom: Mean-variance trend") +
    xlab(v$voom.xy$xlab) + ylab(v$voom.xy$ylab) + 
    theme(axis.title=element_text(size=12),
          plot.title = element_text(size = 14, face = "bold",hjust = 0.5))
  
  #Print out sample numbers:
  #
  sampsize <- colSums(design)
  titleval <- "Please note Sample size:"
  titletext <- paste(names(sampsize), sampsize, sep = "=", collapse = " \n ") 
  titleall <- paste(titleval,"\n",titletext,"\n\n\n")
  
  contrast <- colnames(cm)
  connames <- strsplit(contrast,"-")
  connames <- lapply(connames,function(x) {gsub("\\(","",gsub("\\)","",x))})
  contrastsize <- lapply(connames,function(x) sampsize[unlist(x)])
  footnotetext <- paste(contrast, contrastsize, sep = " : ", collapse = "\n") 
  footnotetext <- paste("\n\n\nContrasts:\n",footnotetext)
  
  textall <- textGrob(paste0(titleall, footnotetext),gp=gpar(fontsize=10))
  
  #Run Contrasts
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  logFC = fit2$coefficients
  colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
  tstat = fit2$t
  colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
  FC = 2^fit2$coefficients
  FC = ifelse(FC<1,-1/FC,FC)
  colnames(FC)=paste(colnames(FC),"FC",sep="_")
  pvalall=fit2$p.value
  colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
  pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
  colnames(pvaladjall)=paste(colnames(fit2$coefficients),"adjpval",sep="_")
  
  
  if(return_mean_and_sd == TRUE){
    tve <- t(v$E)        
    mean.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% dplyr::mutate(group=sample_metadata[sample_metadata[,sample_name_column]==Sample,contrast_variable_column]) %>% group_by(group) %>% summarise_all(funs(mean)) %>% as.data.frame()
    mean.df[,-c(1,2)] %>% as.matrix() %>% t() -> mean
    colnames(mean) <- mean.df[,1]
    colnames(mean)=paste(colnames(mean),"mean", sep="_")
    colnames(mean) = gsub("\\.", "_", colnames(mean))
    
    sd.df <- as.data.frame(tve) %>% rownames_to_column("Sample") %>% dplyr::mutate(group=sample_metadata[sample_metadata[,sample_name_column]==Sample,contrast_variable_column]) %>% group_by(group) %>% summarise_all(funs(sd)) %>% as.data.frame()
    sd.df[,-c(1,2)] %>% as.matrix() %>% t() -> sd
    colnames(sd) <- sd.df[,1]
    colnames(sd)=paste(colnames(sd), "sd",sep="_")
    colnames(sd) = gsub("\\.", "_", colnames(sd))
    finalres=as.data.frame(cbind(mean, sd,  FC, logFC, tstat, pvalall, pvaladjall)) 
  } else {
    finalres=as.data.frame(cbind(FC, logFC, tstat, pvalall, pvaladjall))
  }
  
  if(return_normalized_counts == TRUE){
    finalres = as.data.frame(cbind(finalres, v$E))
  }
  
  finalres %>% rownames_to_column("Gene") -> finalres
  print(paste0("Total number of genes included: ", nrow(finalres)))
  
  getgenelists <- function(FClimit,pvallimit,pval){
    upreggenes <- list()
    downreggenes <- list()
    for(i in 1:length(contrasts)){
      if(pval == "pval"){
        finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvalall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
        finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvalall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]]        
      } else {
        finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] > FClimit & .data[[colnames(pvaladjall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> upreggenes[[i]] 
        finalres %>% dplyr::filter(.data[[colnames(FC)[i]]] < -FClimit & .data[[colnames(pvaladjall)[i]]] < pvallimit) %>% pull(Gene) %>% length() -> downreggenes[[i]] 
      }
    }
    names(upreggenes) <- contrasts
    names(downreggenes) <- contrasts
    allreggenes <- rbind(unlist(upreggenes),unlist(downreggenes))
    rownames(allreggenes) <- c(paste0("upreg>",FClimit, ", ",pval,"<",pvallimit),paste0("downreg<-",FClimit, ", ",pval,"<",pvallimit))
    return(allreggenes)
  }
  
  FCpval1 <- getgenelists(FClimit = 1.2, pvallimit = 0.05,"pval")
  FCpval2 <- getgenelists(FClimit = 1.2, pvallimit = 0.01,"pval")
  FCadjpval1 <- getgenelists(FClimit = 1.2, pvallimit = 0.05,"adjpval")
  FCadjpval2 <- getgenelists(FClimit = 1.2, pvallimit = 0.01,"adjpval")
  
  wraplines <- function(y){
    j = unlist(strsplit(y,"-"))
    k = strwrap(j, width = 10)
    l = paste(k,collapse="\n-")
    return(l)
  }
  
  pvaltab <- rbind(FCpval1,FCpval2,FCadjpval1,FCadjpval2)
  colnames(pvaltab) <- sapply(colnames(pvaltab), function(x) wraplines(x))
  table2 <- tableGrob(pvaltab, theme=ttheme_default(base_size = 10))
  table2$layout$clip <- "off"
  
  layout <- rbind(c(1,2),
                  c(1,2),
                  c(3,3))
  
  
  #Printing all together (tables and plot)
  grid.newpage()
  grid.arrange(textall, g, table2, layout_matrix=layout)
  
  #Printing in brand new multiviz
  grid.newpage()
  print(g)
  grid.newpage()
  grid.draw(textall)
  grid.newpage()
  grid.draw(table2)
  
  ### add back Anno columns and Remove row number from Feature Column
  colnames(finalres)[colnames(finalres)%in%"Gene"]=gene_names_column
  
  finalres=merge(anno_tbl,finalres,by=gene_names_column,all.y=T)
  finalres[,gene_names_column]=gsub('_[0-9]+$',"",finalres[,gene_names_column])
  
  call_me_alias<-colnames(finalres)
  colnames(finalres)<-gsub("\\(|\\)","", call_me_alias)
  df.final<-finalres
  
  return(df.final) 
}

deg_res <- DEGAnalysis(Filtered_Counts = filtered_counts, Ccbr1321_metadata = meta)

write.csv(deg_res, "pptx1_deg_results.csv")
# L2P Analysis for Multiple Comparisons [CCBR] (faa87d26-5559-4852-ac45-4c576502c58c): v83
L2PMultiple <- function(DEGAnalysis) {
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  #libs <- c("l2p","l2psupp","dplyr","tidyr","magrittr","ggplot2","stringr","RCurl","plyr")
  libs <- c("l2p","l2psupp","dplyr","tidyr","magrittr","ggplot2","stringr","RCurl","tidyverse")
  lapply(libs, function(x) library(x, character.only = T, quietly = T))
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Primary Inputs:
  deg_table <- DEGAnalysis
  
  #Basic Parameters:
  gene_names_column <- "Gene"
  species = tolower("Human")
  collections_to_include = c("H","GO")
  select_by_rank <- F
  top_pathways <- 30
  number_of_significant_events <- 1
  
  #Genelist selected by t-statistic rank parameters:
  t_statistic_columns <- c("hYP7_CAR_T-CD19_tstat",
                           "Y_Y-hYP7_CAR_T_tstat",
                            "Y_Y-CD19_tstat",
                            "Y_H-CD19_tstat",
                            "Y_Y-Y_H_tstat") 
  select_top_percentage_of_genes <- TRUE
  select_top_genes <- 500
  
  #Genelist selected by fold-change and pval parameters:    
  significance_columns <- c("hYP7_CAR_T-CD19_pval","Y_Y-hYP7_CAR_T_pval",
                            "Y_Y-CD19_pval","Y_H-CD19_pval","Y_Y-Y_H_pval") 
  significance_threshold <- 0.05
  fold_change_columns <- c("hYP7_CAR_T-CD19_FC","Y_Y-hYP7_CAR_T_FC",
                           "Y_Y-CD19_FC","Y_H-CD19_FC","Y_Y-Y_H_FC")
  fold_change_threshold <- 1.2
  
  #Visual Parameters:
  plot_bubble_size <- "pval"
  plot_bubble_color <- "enrichment_score"
  pathway_axis_label_max_length <- 70
  pathway_axis_label_font_size <- 5
  
  #Advanced Parameters:
  use_built_in_gene_universe = FALSE
  minimum_pathway_hit_count = 15 
  pathway_size_limit <- 500
  p_value_limit <- 0.05
  use_fdr_for_significance <- FALSE
  
  # duplicated pathways
  pathways_to_remove <- c("Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
                          "organic acid catabolic process",
                          # "positive regulation of nervous system development",
                          # "positive regulation of neurogenesis",
                          "small molecule catabolic process",
                          "sterol metabolic process",
                          "cellular amine metabolic process",
                          "amine metabolic process")
  rename_groups <- c()
  
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  if(select_by_rank == TRUE && is.null(t_statistic_columns)){
    stop("ERROR: Choose t-statistics columns since selecting by rank")
  }
  
  if((select_by_rank == FALSE && is.null(fold_change_columns) == TRUE) | (select_by_rank == FALSE && is.null(significance_columns)))  {
    stop("ERROR: Choose fold change and p-value columns and ensure both are in the same order")   
  }
  
  if(select_by_rank == FALSE && !is.null(fold_change_threshold) && !is.null(significance_columns)) {
    FCgroups <- gsub("_FC|_logFC|avg_logFC_|avg_log2FC_","",fold_change_columns)
    pvalgroups <- gsub("_pval|_adjpval|p_val_|p_val_adj_","",significance_columns)
    if(!identical(FCgroups,pvalgroups)){
      stop("ERROR: Make sure fold change and pval columns are in the same group order")
    }
  }
  
  if(sum(grepl("logFC|avg_log2FC",fold_change_columns) | sum(grepl("_FC",fold_change_columns))) != length(fold_change_columns)){
    stop("ERROR: Make sure fold change columns are consistently either log fold change or fold change")  
  }
  
  if(sum(grepl("_pval",significance_columns) | sum(grepl("_adjpval",significance_columns)) | sum(grepl("p_val",significance_columns))) != length(significance_columns)){
    stop("ERROR: Make sure significance columns are consistently either pval or adjpval")  
  }
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  return_org_genes <- function(l2pout){
    l2pgenes <- as.list(l2pout$genes)
    l2pgenes <- lapply(l2pgenes, function(x) unlist(strsplit(x," ")))
    l2pgenesnew <- lapply(l2pgenes,function(a) o2o(a,"human",species))
    l2pout$orig_genes <- l2pgenesnew
    l2pout$orig_genes <- sapply(l2pout$orig_genes, paste, collapse=' ')
    l2pout %>% arrange(pval) -> l2pout
    return(l2pout)
  }
  
  return_orig_genes <- function(l2pout, new_gene_names){
    l2pgenes <- as.list(l2pout$genesinpathway)
    l2pgenes <- lapply(l2pgenes, function(x) unlist(strsplit(x," ")))
    l2pgenesorig <- lapply(l2pgenes,function(a) names(new_gene_names[a]))
    l2pout$orig_genes <- l2pgenesorig
    l2pout$orig_genes <- sapply(l2pout$orig_genes, paste, collapse=' ')
    l2pout %>% arrange(pval) -> l2pout
    return(l2pout)
  }
  
  run_l2p <- function(genes_to_include) {   
    genes_to_include <- as.vector(unique(unlist(genes_to_include)))
    
    if(species == "human"){
      new_gene_names <- sapply(genes_to_include, function(x) updategenes(x, trust=1))
      
      #Print out genes in genelist that are updated:
      updated_genes_idx <- sapply(seq_along(new_gene_names), function(i) names(new_gene_names)[i] != new_gene_names[[i]])
      updated_genes <- new_gene_names[updated_genes_idx]
      updated_genes_num <- sum(updated_genes_idx)
      
      print(paste("Number of updated genes:",updated_genes_num))
      
      cat("Original:Updated\n")
      print(sapply(seq_along(updated_genes), function(i) print(paste0(names(updated_genes)[i],":",updated_genes[i]))))
      
      #Set up genelist and gene universe using updated genes
      genes_to_include <- as.character(new_gene_names)
      gene_universe <- updategenes(gene_universe,trust=1)
      lastgene <- names(tail(new_gene_names,1))
      
      # genes_to_include <- updategenes(genes_to_include,trust=1) 
      # genes_universe <- updategenes(genes_universe,trust=1)
    }
    
    if (species != "human") {
      #Get homologs for genelist:
      orth_genes <- sapply(genes_to_include, function(x) o2o(x, species,"human")[1])
      no_orth <- names(orth_genes[unlist(lapply(orth_genes, function(x) is.na(x)))])
      orth <- orth_genes[unlist(lapply(orth_genes, function(x) !is.na(x)))]
      
      #Print out numbers of genes with homologs in genelist:
      num_no_orth <- length(no_orth)
      perc_num_no_orth <- formatC((length(no_orth)/length(genes_to_include))*100,digits=2,format="f")
      num_orth <- length(orth)
      perc_num_orth <- formatC((length(orth)/length(genes_to_include))*100,digits=2,format="f")
      
      cat(paste("\n\nNumber of genes in the genelist without a homologue:",num_no_orth,", Percentage:",perc_num_no_orth,"%\n"))
      print(no_orth)
      
      cat(paste("\n\nNumber of genes in the genelist with a homologue:",num_orth,",Percentage:",perc_num_orth,"%\n"))
      cat("\nGene:Homolog\n")
      print(sapply(seq_along(orth), function(i) paste0(names(orth)[i],":",orth[i])))
      
      #Set up genelist using gene homologs
      lastgene <- names(tail(orth,1))
      genes_to_include <- as.character(unlist(orth))
      
      #Get homologs for gene universe:
      orth_gene_universe <- sapply(gene_universe, function(x) o2o(x, species,"human")[1])
      no_orth_gu <- names(orth_gene_universe[unlist(lapply(orth_gene_universe, function(x) is.na(x)))])
      orth_gu <- orth_gene_universe[unlist(lapply(orth_gene_universe, function(x) !is.na(x)))]
      
      gene_universe <- as.character(unlist(orth_gu))
      
      #Print out numbers of genes with homologs in genelist:
      num_no_orth <- length(no_orth_gu)
      num_orth <- length(orth_gu)
      cat(paste("\n\nNumber of genes in the gene universe without a homologue:",num_no_orth,"\n"))
      #print(no_orth_gu)
      cat(paste("Number of genes in the gene universe with a homologue:",num_orth,"\n"))
      #print(gene_universe)
      
      #genes_to_include <- o2o(genes_to_include,species,"human")
      #genes_universe <- o2o(unique(genes_universe),species,"human")
    }
    
    if (use_built_in_gene_universe == TRUE) {
      x <- l2p(genes_to_include, categories=collections_to_include)
      print("Using built-in gene universe.")
    } else {
      x <- l2p(genes_to_include, categories=collections_to_include, universe=gene_universe)
    }
    
    if (species != "human"){
      x <- as.data.frame(return_org_genes(x))
    } else {
      x <- as.data.frame(return_orig_genes(x, new_gene_names))
    }
    
    return(x)
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  genelists <- list()
  lists <- list()
  if(select_by_rank == TRUE){
    compnum <- length(t_statistic_columns)
    groups <- gsub("_tstat","",t_statistic_columns)
    for (i in 1:compnum){
      deg_table %>% dplyr::select(.data[[gene_names_column]],.data[[t_statistic_columns[i]]]) -> genesmat
      if(select_top_percentage_of_genes == TRUE){
        numselect <- ceiling(0.1*dim(deg_table)[1])
      } else {
        numselect <- select_top_genes
      }
      genesmat %>% dplyr::filter(!is.na(.data[[t_statistic_columns[i]]])) %>% dplyr::arrange(desc(.data[[t_statistic_columns[i]]])) -> genesmat
      genesmat %>% head(numselect) %>% pull(.data[[gene_names_column]]) -> lists[[1]]
      genesmat %>% dplyr::filter(!is.na(.data[[t_statistic_columns[i]]])) %>% dplyr::arrange(.data[[t_statistic_columns[i]]]) -> genesmat
      genesmat %>% head(numselect) %>% pull(.data[[gene_names_column]]) -> lists[[2]]   
      genelists[[i]] <- list(lists[[1]],lists[[2]])
    } 
  } else {
    compnum <- length(fold_change_columns)
    if(sum(grepl("logFC|avg_log2FC",fold_change_columns)) == compnum){
      groups <- gsub("_logFC|avg_log2FC_","",fold_change_columns)
      for (i in 1:compnum){
        deg_table %>% dplyr::select(.data[[gene_names_column]],.data[[fold_change_columns[i]]],.data[[significance_columns[i]]]) -> genesmat
        logFC_threshold <- log2(fold_change_threshold) #Upregulated genes
        genesmat %>% dplyr::arrange(.data[[significance_columns[i]]]) %>% dplyr::filter(.data[[significance_columns[i]]] <= significance_threshold & .data[[fold_change_columns[i]]] >= logFC_threshold) %>% pull(gene_names_column) -> lists[[1]]
        logFC_threshold <- -1*log2(fold_change_threshold) #Downregulated genes
        genesmat %>% dplyr::arrange(.data[[significance_columns[i]]]) %>% dplyr::filter(.data[[significance_columns[i]]] <= significance_threshold & .data[[fold_change_columns[i]]] <= logFC_threshold) %>% pull(gene_names_column) -> lists[[2]]
        genelists[[i]] <- list(lists[[1]],lists[[2]])
      }
    } else {
      groups <- gsub("_FC","",fold_change_columns) 
      for (i in 1:compnum){
        deg_table %>% dplyr::select(.data[[gene_names_column]],.data[[fold_change_columns[i]]],.data[[significance_columns[i]]]) -> genesmat
        genesmat %>% dplyr::arrange(.data[[significance_columns[i]]]) %>% dplyr::filter(.data[[significance_columns[i]]] <= significance_threshold & .data[[fold_change_columns[i]]] >= fold_change_threshold) %>% pull(gene_names_column) -> lists[[1]]
        genesmat %>% dplyr::arrange(.data[[significance_columns[i]]]) %>% dplyr::filter(.data[[significance_columns[i]]] <= significance_threshold & .data[[fold_change_columns[i]]] <= -1*fold_change_threshold) %>% pull(gene_names_column) -> lists[[2]]
        genelists[[i]] <- list(lists[[1]],lists[[2]])
      } 
    }
  }
  
  names(genelists) <- groups 
  
  #Error messaging for low number of genes in genelist:
  genelengths <- lapply(genelists, function(x) lapply(x, length))
  genelistnums <- unlist(lapply(genelengths, function(x) lapply(x, function(x) x[1])))
  names(genelistnums) <- gsub("1$","_upregulated",names(genelistnums))
  names(genelistnums) <- gsub("2$","_downregulated",names(genelistnums))
  lowgenenums <- names(genelistnums)[genelistnums < 100]
  if(sum(genelistnums < 100) > 0){
    warning(sprintf("The size of these genelists are below a desired size threshold of 100: Try loosening the criteria for p-values up to 0.15 and/or fold change to 1.2"))
    print(lowgenenums) 
  }
  
  reg.table <- list()
  for(i in 1:length(genelists)){
    genesize <- lapply(genelists[[i]],function(x) length(x))
    lastgene <- lapply(genelists[[i]],function(x) tail(x,1))
    if(is.null(fold_change_columns[i])){
      coln <- colnames(deg_table)
      #colidx <- coln[grep(groups[i],coln)]
      #FCcol <- colidx[grep("_FC",colidx)]
      FCcol <- paste0(groups[i],"_FC")
    } else {
      FCcol <- fold_change_columns[i]            
    }
    if(is.null(significance_columns[i])){
      coln <- colnames(deg_table)
      #colidx <- coln[grep(groups[i],coln)]
      #pvalcol <- colidx[grep("_pval",colidx)]
      pvalcol <- paste0(groups[i],"_pval")
    } else {
      pvalcol <- significance_columns[i]         
    }
    
    deg_table %>% dplyr::select(.data[[gene_names_column]],.data[[FCcol]],.data[[pvalcol]]) %>% mutate(Group = names(genelists)[i]) %>% mutate(Direction = "Upregulated") %>% mutate(Listsize = genesize[[1]]) %>% dplyr::filter(.data[[gene_names_column]] == lastgene[1]) -> upreggene
    deg_table %>% dplyr::select(.data[[gene_names_column]],.data[[FCcol]],.data[[pvalcol]]) %>% mutate(Group = names(genelists)[i]) %>% mutate(Direction = "Downregulated") %>% mutate(Listsize = genesize[[2]]) %>% dplyr::filter(.data[[gene_names_column]] == lastgene[2]) -> downreggene
    reg.table[[i]] <- rbind(upreggene,downreggene)
    colnames(reg.table[[i]]) <- c("Last_Gene_in_list","FC","pval or adjpval","Group","Direction","Listsize")
  }
  
  combined_table <- bind_rows(reg.table)
  cat("Check for significance of selected genelists:\n\n")
  print(combined_table)
  gene_universe = as.vector(unique(unlist(deg_table[gene_names_column])))
  
  l2presults <- list()
  for (i in 1:length(genelists)){
    l2presults[[i]] <- lapply(genelists[[i]], function(x) {run_l2p(x)})
  }
  
  #Error messaging for zero significant pathway results:
  names(l2presults) <- groups
  l2plength <- lapply(l2presults, function(x) lapply(x, dim))
  l2poutpaths <- unlist(lapply(l2plength, function(x) lapply(x, function(x) x[1])))
  names(l2poutpaths) <- gsub("1$","_upregulated",names(l2poutpaths))
  names(l2poutpaths) <- gsub("2$","_downregulated",names(l2poutpaths))
  nopaths <- names(l2poutpaths)[l2poutpaths== 0]
  
  if(sum(l2poutpaths == 0) > 0){
    stop(sprintf("ERROR: At least one of the l2p results shows no significant pathways, probably because of too few genes in genelist: %s", nopaths)) 
  }
  
  colname = unlist(names(genelists))
  pathlist <- list()
  
  for (i in 1:length(l2presults)){
    paths <- list()
    if(use_fdr_for_significance == TRUE){
      paths <- lapply(l2presults[[i]], function(x) {x %>%
          dplyr::filter(number_hits > minimum_pathway_hit_count) %>%
          dplyr::filter(fdr < p_value_limit) %>%
          head(top_pathways) %>% dplyr::select(pathway_name)})    
    } else {
      paths <- lapply(l2presults[[i]], function(x) {x %>%
          dplyr::filter(number_hits > minimum_pathway_hit_count) %>%
          dplyr::filter(pval < p_value_limit) %>%
          head(top_pathways) %>% dplyr::select(pathway_name)})
    }
    pathlist[[i]] <- unlist(lapply(paths,function(x) {unlist(x, use.names=FALSE)}))
  }
  
  path.all <-  data.frame(pathwayname = unlist(pathlist))
  
  path.all %>% group_by(pathwayname) %>% 
    tally() %>% 
    arrange(dplyr::desc(n)) %>% 
    dplyr::filter(n>=number_of_significant_events) %>% 
    dplyr::pull(pathwayname) -> path.select
  
  pathmerge <- list()
  for (i in 1:length(l2presults)){
    pathselect <- list()
    pathselect <- lapply(l2presults[[i]], function(x) {
      dplyr::filter(x, pathway_name %in% path.select) %>% 
        mutate(total = number_hits + number_misses) %>%  
        dplyr::distinct(pathway_name, .keep_all = TRUE) %>% 
        dplyr::filter(total < pathway_size_limit) %>%
        dplyr::mutate(percent_gene_hits_per_pathway = percent_gene_hits_per_pathway) %>%
        dplyr::select(pathway_name,pathway_id,category,enrichment_score,number_hits,total,percent_gene_hits_per_pathway,pval,fdr,genesinpathway,orig_genes) 
    })                              
    
    pathselect.merge <- merge(pathselect[[1]],pathselect[[2]],by="pathway_name",all=TRUE) %>% 
      dplyr::mutate_if(is.numeric, tidyr::replace_na, 0) %>%
      dplyr::mutate(net_enrichment_score = enrichment_score.x-enrichment_score.y) %>% 
      dplyr::mutate(net_number_hits = number_hits.x - number_hits.y) %>%
      dplyr::mutate(enrichment_score = case_when(net_enrichment_score>0 ~enrichment_score.x,
                                                 net_enrichment_score < 0 ~ -1*enrichment_score.y,TRUE ~ 0)) %>% 
      dplyr::mutate(categ = case_when(net_enrichment_score > 0 ~category.x,
                                      net_enrichment_score < 0 ~category.y, TRUE ~ "NA")) %>%
      dplyr::mutate(number_hits = case_when(net_enrichment_score > 0 ~number_hits.x,
                                            net_enrichment_score < 0 ~ number_hits.y, TRUE ~ as.integer(0))) %>%
      dplyr::mutate(percent_gene_hits_per_pathway = case_when(net_enrichment_score > 0 ~percent_gene_hits_per_pathway.x,
                                                              net_enrichment_score < 0 ~ -1*percent_gene_hits_per_pathway.y, TRUE ~ 0)) %>%
      dplyr::mutate(pval = case_when(net_enrichment_score > 0 ~pval.x,
                                     net_enrichment_score < 0 ~ pval.y, 
                                     TRUE ~ 0)) %>%
      dplyr::mutate(fdr = case_when(net_enrichment_score > 0 ~fdr.x,
                                    net_enrichment_score < 0 ~ fdr.y,
                                    TRUE ~ 0)) %>%
      dplyr::mutate(genes = case_when(net_enrichment_score > 0 ~genesinpathway.x,
                                      net_enrichment_score < 0 ~genesinpathway.y)) %>%
      dplyr::mutate(orig_genes = case_when(net_enrichment_score > 0 ~orig_genes.x,
                                           net_enrichment_score < 0 ~orig_genes.y)) %>%
      dplyr::mutate(pathway_id = case_when(net_enrichment_score > 0 ~pathway_id.x,
                                           net_enrichment_score < 0 ~pathway_id.y)) %>%
      dplyr::mutate(group=colname[i])
    pathmerge[[i]] <- pathselect.merge
  }
  
  pathall <- bind_rows(pathmerge) %>% select(pathway_name,pathway_id,categ,number_hits, percent_gene_hits_per_pathway,enrichment_score,pval,fdr,net_number_hits,net_enrichment_score,genes,orig_genes,group) %>% arrange(pval)
  
  if(select_by_rank == TRUE){
    grouplevel <- gsub("_tstat","",t_statistic_columns)
  } else {
    grouplevel <- gsub("_FC|_logFC|avg_log2FC_","",fold_change_columns)
  }
  pathall$group <- factor(pathall$group,levels=grouplevel)
  pathall %>% dplyr::filter(!pathway_name %in% pathways_to_remove) -> pathall
  if(length(rename_groups)>0){
    names(rename_groups) <- grouplevel
    pathall$group <- rename_groups[pathall$group]
    pathall$group <- factor(pathall$group,levels=rename_groups)
  }
  
  pathall2 <- pathall
  pathall2$pathway_name2 <- str_to_upper(pathall2$pathway_name)
  pathall2$pathway_name2 <- gsub("_"," ",pathall2$pathway_name2)
  pathall2 %>% dplyr::mutate(pathway_name2 = stringr::str_wrap(pathway_name2,pathway_axis_label_max_length)) -> pathall2
  maxabscore <- max(abs(range(pathall2[[plot_bubble_color]])))
  maxscore = maxabscore
  minscore = -1*maxabscore
  
  col1 <- sym(plot_bubble_size)
  col2 <- sym(plot_bubble_color)
  if(plot_bubble_size %in% c("pval","fdr")){
    g <- ggplot(pathall2, aes(x = group, y = reorder(pathway_name2,enrichment_score), size = -log10(!!col1), colour = !!col2)) + 
      geom_point() +
      scale_size_continuous(range = c(1,50)) +
      theme_bw() + 
      ylab("Pathways") +
      xlab("Treatment") +
      scale_colour_gradient2(limits=c(minscore, maxscore),midpoint = 0,low="blue",mid="darkgrey", high="tomato",oob = scales::squish) +
      #scale_colour_gradient(low="darkblue", high="tomato") +
      scale_size(range = c(0, 10)) +
      scale_y_discrete(expand = c(0.05, 0.05)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),axis.text.y = element_text(colour="black", size = 14)) + 
      guides(
        size = guide_legend(override.aes = list(color = "black", fill = "white",
                                                shape = 1))
      )
  } else {  
    g <- ggplot(pathall2, aes(x = group, y = reorder(pathway_name2,enrichment_score), size = !!col1, colour = !!col2)) + 
      geom_point() +
      scale_size_continuous(range = c(1,50)) +
      theme_bw() + 
      ylab("Pathways") +
      xlab("Treatment") +
      scale_colour_gradient2(limits=c(minscore, maxscore),midpoint = 0,low="darkblue",mid="grey", high="tomato",oob = scales::squish) +
      scale_size(range = c(0, 10)) +
      scale_y_discrete(expand = c(0.05, 0.05)) +
      guides(
        size = guide_legend(override.aes = list(color = "black", fill = "white",
                                                shape = 1))
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_text(colour="black", size = pathway_axis_label_font_size))
  } 
  
  print(g)
  ggsave("ccbr1321_dotplots.pdf",g, width = 12, height = 16, dpi = 500)
  return(pathall)
}

l2p_res = L2PMultiple(DEGAnalysis = deg_res)

write.csv(l2p_res, "l2p_res.csv")

# Expression Heatmap [CCBR] (89a32987-10d9-4233-91c3-e9adf3dcc517): v550
heatmap_1 <- function(Filtered_Counts,Ccbr1321_metadata) {
  ## This function uses pheatmap to draw a heatmap, scaling first by rows
  ## (with samples in columns and genes in rows)
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(colorspace)
  library(dendsort)
  library(ComplexHeatmap)
  library(dendextend)
  library(tibble)
  library(stringr)
  library(RColorBrewer)
  library(dplyr)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(gridGraphics)
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters:
  counts_matrix <- Filtered_Counts
  sample_metadata <- Ccbr1321_metadata
  gene_column_name <- "Gene"
  group_columns <- c("Group")
  sample_name_column <- "Sample"
  samples_to_include = c("X1_CD19_314_T_resent","X2_Y_Y_319_T_resent",
                         "X3_Y_Y_304_T_resent","X4_Y_Y_729_T_resent",
                         "X5_Y_Y_727_T_resent","X6_CD19_738_T_resent",
                         "X7_CD19_271_T_resent","X8_Y_H_305_T_resent",
                         "X9_CD19_324_T_resent","X10_hYP7_731_T_resent",
                         "X11_Y_H_314_T_resent","X12_hYP7_740_T_resent",
                         "X13_Y_H_737_T_resent","X14_hYP7_732_T_resent",
                         "X15_Y_H_730_T_resent","X16_hYP7_303_T_resent")
  include_all_genes <- FALSE
  filter_top_genes_by_variance = TRUE
  top_genes_by_variance_to_include <- 0
  specific_genes_to_include_in_heatmap = c(
    "TMEM43", "AGPAT4", "ZBTB43", "AC020916.1", "NPIPB12", "CD69", "PDCD1", "RGPD1",
    "TNF", "DUSP10", "IFNG", "LEF1", "CDC7", "MKI67", "EOMES", "SLAMF7", "PECAM1",
    "GZMK", "LAG3", "ITGAX", "ATP8B4", "BCL2L11", "LRMP", "TIGIT", "CSF1",
    "HLA-DQA1", "ADORA2A", "FASLG", "IL2RA", "JAKMIP1", "M6PR", "BAB37", "FAS",
    "ICOS", "CD4", "HAVCR2", "GZMB", "CD8A", "CCL5", "GZMA", "GZMB", "CD8B",
    "CCNB2", "TOP2A", "MTHFD1", "SM3BP1", "HIMS3", "HJURP", "RPA1", "LRRC20",
    "XRCC5", "CCNB1", "KLF6", "FBXO10", "TOX4", "SLC12A7", "PLEKHG3", "CYP1B1",
    "CORO1B", "CLU", "DHCR24", "H1F0", "KLF1C", "MYC", "TRIB1", "IDH1", "APLP2",
    "FASN", "NADK2", "ATP1B1", "HMGCS1", "GLUD1", "SCD"
  )
  
  specific_genes_to_include_in_heatmap <- stringr::str_flatten(gsub(","," ",specific_genes_to_include_in_heatmap), collapse = " ")
  
  
  #Visual Parameters:
  heatmap_color_scheme <- "Bu Wt Rd"
  autoscale_heatmap_color <- TRUE
  set_min_heatmap_color <- -2
  set_max_heatmap_color <- 2
  group_colors <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  assign_group_colors <- FALSE
  assign_color_to_sample_groups <- c()
  legend_font_size <- 10 
  display_gene_names <- TRUE
  gene_name_font_size <- 6
  display_sample_names <- TRUE
  sample_name_font_size <- 8
  display_dendrograms <- TRUE
  reorder_dendrogram <- FALSE
  reorder_dendrogram_order <- c()
  manually_rename_samples <- FALSE
  samples_to_rename <- c("")
  display_numbers <- FALSE
  aspect_ratio <- "Auto"
  
  #Advanced Parameters
  distance_metric <- "correlation"
  clustering_method <- "average"
  center_and_rescale_expression <- TRUE
  cluster_genes <- TRUE
  cluster_samples <- TRUE
  arrange_sample_columns <- FALSE
  order_by_gene_expression <- FALSE
  gene_to_order_columns <- " "
  gene_expression_order <- "low_to_high"
  return_z_scores <- TRUE
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  if(include_all_genes == TRUE && filter_top_genes_by_variance == TRUE){
    stop("ERROR: Choose only one of 'Include all genes' or 'Filter top genes by variance' as TRUE")
  }
  
  if((cluster_samples == TRUE && arrange_sample_columns == TRUE) | (arrange_sample_columns == TRUE && order_by_gene_expression == TRUE) | 
     (arrange_sample_columns == TRUE && cluster_samples == TRUE) | (cluster_samples == FALSE && arrange_sample_columns == FALSE && order_by_gene_expression == FALSE)) {
    stop("ERROR: Choose only one of 'Cluster Samples', 'Arrange sample columns', or 'order by gene expression' as TRUE")   
  }
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  getourrandomcolors<-function(k){
    seed=10
    n <- 2e3
    ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
    ourColorSpace <- as(ourColorSpace, "LAB")
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    return( unname(hex(LAB(km$centers))))
  }
  
  ## Begin pal() color palette function∂:
  pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
    if (n < 1L) {
      return(character(0L))
    }
    h <- rep(h, length.out = 2L)
    c <- c[1L]
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, -1, length = n)
    rval <- hex(
      polarLUV(
        L = l[2L] - diff(l) * abs(rval)^power[2L], 
        C = c * abs(rval)^power[1L],
        H = ifelse(rval > 0, h[1L], h[2L])
      ),
      fixup=fixup, ...
    )
    if (!missing(alpha)) {
      alpha <- pmax(pmin(alpha, 1), 0)
      alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                      width = 2L, upper.case = TRUE)
      rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  } 
  # End pal() color palette function:
  
  ## Begin doheatmap() function:
  doheatmap <- function(dat, clus, clus2, ht, rn, cn, col, dispnum) {
    #require(pheatmap)
    #require(dendsort)
    col.pal <- np[[col]]
    if (FALSE) {
      col.pal = rev(col.pal)
    }
    # Define metrics for clustering
    drows1 <- distance_metric
    dcols1 <- distance_metric
    minx = min(dat)
    maxx = max(dat)
    if (autoscale_heatmap_color) {
      breaks = seq(minx, maxx, length=100)
      legbreaks = seq(minx, maxx, length=5)
    } else {
      breaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=100)
      legbreaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=5)
    }
    breaks = sapply(breaks, signif, 4)
    legbreaks = sapply(legbreaks, signif, 4)
    # Run cluster method using 
    hc = hclust(dist(t(dat)), method=clustering_method)
    hcrow = hclust(dist(dat), method=clustering_method)
    if (FALSE) {
      sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
    } else {
      sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
    }
    if (clus) {
      colclus <- sort_hclust(hc)
    } else {
      colclus = FALSE
    }
    if (clus2) {
      rowclus <- sort_hclust(hcrow)
    } else {
      rowclus = FALSE
    }
    if (display_dendrograms) {
      treeheight <- 25
    } else {
      treeheight <- 0
    }
    
    hm.parameters <- list(
      dat, 
      color=col.pal,
      legend_breaks=legbreaks,
      legend=TRUE,
      scale="none",
      treeheight_col=treeheight,
      treeheight_row=treeheight,
      kmeans_k=NA,
      breaks=breaks,
      display_numbers=dispnum,
      number_color = "black",
      fontsize_number = 8,
      height=80,
      cellwidth = NA, 
      cellheight = NA, 
      fontsize= legend_font_size,   
      fontsize_row=gene_name_font_size,
      fontsize_col=sample_name_font_size,
      show_rownames=rn, 
      show_colnames=cn,
      cluster_rows=rowclus, 
      cluster_cols=clus,
      clustering_distance_rows=drows1, 
      clustering_distance_cols=dcols1,
      annotation_col = annotation_col,
      annotation_colors = annot_col,
      labels_col = labels_col
    )
    mat = t(dat)
    callback = function(hc, mat) {
      dend = rev(dendsort(as.dendrogram(hc)))
      if(reorder_dendrogram == TRUE) {
        dend %>% dendextend::rotate(reorder_dendrogram_order) -> dend
      } else {
        dend %>% dendextend::rotate(c(1:nobs(dend))) 
      }
      as.hclust(dend)
    }
    
    ## Make Heatmap
    phm <- do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
    
  }
  # End doheatmap() function.
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  ## Build different color spectra options for heatmap:
  np0 = pal(100) 
  np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1) # Blue to Red
  np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)) # Red to Vanilla
  np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) # Violet to Pink
  np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)) #Red to yellow to blue
  np5 = colorRampPalette(c("steelblue","white", "red"))(100)  # Steelblue to White to Red
  
  ## Gather list of color spectra and give them names for the GUI to show.
  np = list(np0, np1, np2, np3, np4, np5)
  names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")
  
  ## Parse input counts matrix. Subset by samples.
  df1 <- counts_matrix
  # Swap out Gene Name column name, if it's not 'Gene'.
  if(gene_column_name != "Gene"){
    # Drop original Gene column
    df1 = df1[,!(colnames(df1)%in% c("Gene")) ]
    # Rename column to Gene
    colnames(df1)[which(colnames(df1) == gene_column_name)] <- 'Gene'
  }
  # Get sample columns
  samples_to_include <- samples_to_include[samples_to_include != gene_column_name]
  samples_to_include <- samples_to_include[samples_to_include != "Gene"]
  samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
  
  # Build new counts matrix containing only sample subset chosen by user.
  df1 <- df1[,append("Gene", samples_to_include)]
  df.orig = df1
  df.orig %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean)) -> df
  df.mat = df[ , (colnames(df) != "Gene" )] %>% as.data.frame
  df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
  row.names(df.mat) <- df$Gene
  rownames(df.mat) <- str_wrap(rownames(df.mat),30) #for really long geneset names
  df.mat <- as.data.frame(df.mat)
  
  ## Subset counts matrix by genes.
  # Toggle to include all genes in counts matrix (in addition to any user-submitted gene list).
  if (include_all_genes == FALSE) {
    # Add user-submitted gene list (optional).
    genes_to_include_parsed = c()
    genes_to_include_parsed = strsplit(specific_genes_to_include_in_heatmap, " ")[[1]]
    df.mat[genes_to_include_parsed,] -> df.final.extra.genes
    if(filter_top_genes_by_variance == TRUE) {
      # Want to filter all genes by variance.
      df.final = as.matrix(df.mat)
      var <- matrixStats::rowVars(df.final)
      df <- as.data.frame(df.final)
      rownames(df) <- rownames(df.final)
      df.final <- df
      df.final$var <- var
      df.final %>% rownames_to_column("Gene") -> df.final 
      df.final %>% dplyr::arrange(desc(var)) -> df.final
      df.final.extra.genes = dplyr::filter(df.final, Gene %in% genes_to_include_parsed)
      df.final = df.final[1:top_genes_by_variance_to_include,]
      df.final = df.final[complete.cases(df.final),]
      # Rbind user gene list to variance-filtered gene list and deduplicate.
      df.final <- rbind(df.final, df.final.extra.genes)
      df.final <- df.final[!duplicated(df.final),] 
      rownames(df.final) <- df.final$Gene
      df.final$Gene <- NULL
      df.final$var <- NULL
    } else {
      # Want to use ONLY user-provided gene list.
      df.final <- df.final.extra.genes
      df.final <- df.final[!duplicated(df.final),]
      # Order genes in heatmap by user-submitted order of gene names.
      df.final <- df.final[genes_to_include_parsed,]
      #df.final$Gene <- NULL
    }
  } else {
    df.final <- df.mat
    df.final$Gene <- NULL
  }
  
  ## Optionally apply centering and rescaling (default TRUE).
  if (center_and_rescale_expression == TRUE) {
    tmean.scale = t(scale(t(df.final)))
    tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
    tmean.scale = na.omit(tmean.scale)
  } else {
    tmean.scale = df.final
  }
  
  if(order_by_gene_expression == TRUE){
    gene_to_order_columns <- gsub(" ","",gene_to_order_columns)
    if(gene_expression_order == "low_to_high"){
      tmean.scale <- tmean.scale[,order(tmean.scale[gene_to_order_columns,])] #order from low to high 
    } else{
      tmean.scale <- tmean.scale[,order(-tmean.scale[gene_to_order_columns,])] #order from high to low  
    }
  }
  
  df.final <- as.data.frame(tmean.scale)
  
  ## Parse input sample metadata and add annotation tracks to top of heatmap.
  annot <- sample_metadata
  # Filter to only samples user requests.
  annot %>% dplyr::filter(.data[[sample_name_column]] %in% samples_to_include) -> annot
  # Arrange sample options.
  if(arrange_sample_columns) {
    annot %>% dplyr::arrange_(.dots=group_columns) -> annot
    df.final <- df.final[,match(annot[[sample_name_column]],colnames(df.final))] 
  }
  # Build subsetted sample metadata table to use for figure.
  
  colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
  names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  group_colors <- colorlist[group_colors]
  
  annot %>% dplyr::select(group_columns) -> annotation_col    
  annotation_col = as.data.frame(unclass(annotation_col))
  annotation_col[] <- lapply(annotation_col,factor)
  x <- length(unlist(lapply(annotation_col,levels)))
  if(x>length(group_colors)){
    k=x-length(group_colors)
    more_cols<- getourrandomcolors(k) 
    group_colors <- c(group_colors, more_cols)
  }
  rownames(annotation_col) <- annot[[sample_name_column]]
  annot_col = list()
  b=1
  i=1
  while (i <= length(group_columns)){
    nam <- group_columns[i]
    grp <- as.factor(annotation_col[,i])
    c <- b+length(levels(grp))-1
    col = group_colors[b:c]
    names(col) <- levels(grp)
    assign(nam,col)
    annot_col = append(annot_col,mget(nam))
    b = c+1
    i=i+1
  }
  
  if(assign_group_colors == TRUE){
    colassign <- assign_color_to_sample_groups
    groupname <- c()
    groupcol <- c() 
    for (i in 1:length(colassign)) {
      groupname[i] <- strsplit(colassign[i], ": ?")[[1]][1]
      groupcol[i] <- strsplit(colassign[i], ": ?")[[1]][2]
    }
    annot_col[[1]][groupname] <- groupcol
  }
  
  ## Setting labels_col for pheatmap column labels.
  if (manually_rename_samples == TRUE) {
    # Use user-provided names to rename samples.
    replacements = samples_to_rename
    old <- c()
    new <- c()
    labels_col <- colnames(df.final)
    for (i in 1:length(replacements)) {
      old <- strsplit(replacements[i], ": ?")[[1]][1]
      new <- strsplit(replacements[i], ": ?")[[1]][2]
      old=gsub("^[[:space:]]+|[[:space:]]+$","",old)
      new=gsub("^[[:space:]]+|[[:space:]]+$","",new)
      labels_col[labels_col==old]=new           
    }
  } else {
    ## Use original column names for samples.
    labels_col <- colnames(df.final)
  }
  
  ## Print number of genes to log.
  print(paste0("The total number of genes in heatmap: ", nrow(df.final)))
  
  # manually remove ALB, don't know why it is showing up
  df.final = df.final[rownames(df.final) != 'ALB',]
  
  ## Make the final heatmap.
  p <- doheatmap(dat=df.final, clus=cluster_samples, clus2=cluster_genes, ht=50, rn=display_gene_names, cn=display_sample_names, col=heatmap_color_scheme, dispnum=display_numbers)
  p@matrix_color_mapping@name <- " "
  p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
  p@column_title_param$gp$fontsize <- 10
  print(p)
  
  ## If user sets toggle to TRUE, return Z-scores.
  ## Else return input counts matrix by default (toggle FALSE).
  ## Returned matrix includes only genes & samples used in heatmap.
  if(return_z_scores){
    df.new <- data.frame(tmean.scale) # Convert to Z-scores.
    df.new %>% rownames_to_column("Gene") -> df.new
    return(df.new)
  } else {
    df.final %>% rownames_to_column("Gene") -> df.new
    return(df.new)
  }
}
heatmap_2 <- function(Filtered_Counts,Ccbr1321_metadata) {
  ## This function uses pheatmap to draw a heatmap, scaling first by rows
  ## (with samples in columns and genes in rows)
  
  ## --------- ##
  ## Libraries ##
  ## --------- ##
  
  library(colorspace)
  library(dendsort)
  library(ComplexHeatmap)
  library(dendextend)
  library(tibble)
  library(stringr)
  library(RColorBrewer)
  library(dplyr)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(gridGraphics)
  
  ## -------------------------------- ##
  ## User-Defined Template Parameters ##
  ## -------------------------------- ##
  
  #Basic Parameters:
  counts_matrix <- Filtered_Counts
  sample_metadata <- Ccbr1321_metadata
  gene_column_name <- "Gene"
  group_columns <- c("Group")
  sample_name_column <- "Sample"
  samples_to_include = c("X1_CD19_314_T_resent","X2_Y_Y_319_T_resent",
                         "X3_Y_Y_304_T_resent","X4_Y_Y_729_T_resent",
                         "X5_Y_Y_727_T_resent","X6_CD19_738_T_resent",
                         "X7_CD19_271_T_resent","X8_Y_H_305_T_resent",
                         "X9_CD19_324_T_resent","X10_hYP7_731_T_resent",
                         "X11_Y_H_314_T_resent","X12_hYP7_740_T_resent",
                         "X13_Y_H_737_T_resent","X14_hYP7_732_T_resent",
                         "X15_Y_H_730_T_resent","X16_hYP7_303_T_resent")
  include_all_genes <- FALSE
  filter_top_genes_by_variance = TRUE
  top_genes_by_variance_to_include <- 0
  specific_genes_to_include_in_heatmap = genes <- c(
    "PLEKHF2", "TBL1X", "FBXO10", "SS18L1", "CSPP1", "FLCN", "NATD1", "LAMTOR3",
    "ABCD4", "UBE2W", "CRAMP1", "PPP2R5B", "ARRDC3", "SLC3A2", "AC012321.1",
    "GK", "GPAT3", "AC104825.1", "PPP1R15A", "AC016831.1", "RAPGEF2", "FOS",
    "SMARCD3", "MARCH3", "AC093525.6", "UAP1L1", "NEAT1", "AC022144.1", "HOMER1",
    "KLF6", "WNT5B", "VCAM1", "DPH6", "CHAF1B", "XRCC3", "POLD1", "LRRC20",
    "HPF1", "EIF2B3", "TMA16", "CDPF1", "FXN", "USP5", "MR1I", "MAP2K6", "HMGXB4",
    "HNRNPD", "ACAA2", "ATP5F1A", "CWC27", "DNASE2", "NSMCE2", "PPP2R5D",
    "ZCRB1", "SLC25A5", "LINS2", "C11orf49", "RPA1", "NASP", "UHRF1", "MCM4",
    "MTHFD1", "DHFR", "LIG1", "MCM5", "RMI2", "GINS2", "FEN1", "MCM7"
  )
  
  specific_genes_to_include_in_heatmap <- stringr::str_flatten(gsub(","," ",specific_genes_to_include_in_heatmap), collapse = " ")
  
  
  #Visual Parameters:
  heatmap_color_scheme <- "Bu Wt Rd"
  autoscale_heatmap_color <- TRUE
  set_min_heatmap_color <- -2
  set_max_heatmap_color <- 2
  group_colors <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  assign_group_colors <- FALSE
  assign_color_to_sample_groups <- c()
  legend_font_size <- 10 
  display_gene_names <- TRUE
  gene_name_font_size <- 6
  display_sample_names <- TRUE
  sample_name_font_size <- 8
  display_dendrograms <- TRUE
  reorder_dendrogram <- FALSE
  reorder_dendrogram_order <- c()
  manually_rename_samples <- FALSE
  samples_to_rename <- c("")
  display_numbers <- FALSE
  aspect_ratio <- "Auto"
  
  #Advanced Parameters
  distance_metric <- "correlation"
  clustering_method <- "average"
  center_and_rescale_expression <- TRUE
  cluster_genes <- TRUE
  cluster_samples <- TRUE
  arrange_sample_columns <- FALSE
  order_by_gene_expression <- FALSE
  gene_to_order_columns <- " "
  gene_expression_order <- "low_to_high"
  return_z_scores <- TRUE
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  if(include_all_genes == TRUE && filter_top_genes_by_variance == TRUE){
    stop("ERROR: Choose only one of 'Include all genes' or 'Filter top genes by variance' as TRUE")
  }
  
  if((cluster_samples == TRUE && arrange_sample_columns == TRUE) | (arrange_sample_columns == TRUE && order_by_gene_expression == TRUE) | 
     (arrange_sample_columns == TRUE && cluster_samples == TRUE) | (cluster_samples == FALSE && arrange_sample_columns == FALSE && order_by_gene_expression == FALSE)) {
    stop("ERROR: Choose only one of 'Cluster Samples', 'Arrange sample columns', or 'order by gene expression' as TRUE")   
  }
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  getourrandomcolors<-function(k){
    seed=10
    n <- 2e3
    ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
    ourColorSpace <- as(ourColorSpace, "LAB")
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    return( unname(hex(LAB(km$centers))))
  }
  
  ## Begin pal() color palette function∂:
  pal = function (n, h=c(237, 43), c=100, l=c(70, 90), power=1, fixup=TRUE, gamma=NULL, alpha=1, ...) {
    if (n < 1L) {
      return(character(0L))
    }
    h <- rep(h, length.out = 2L)
    c <- c[1L]
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, -1, length = n)
    rval <- hex(
      polarLUV(
        L = l[2L] - diff(l) * abs(rval)^power[2L], 
        C = c * abs(rval)^power[1L],
        H = ifelse(rval > 0, h[1L], h[2L])
      ),
      fixup=fixup, ...
    )
    if (!missing(alpha)) {
      alpha <- pmax(pmin(alpha, 1), 0)
      alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                      width = 2L, upper.case = TRUE)
      rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  } 
  # End pal() color palette function:
  
  ## Begin doheatmap() function:
  doheatmap <- function(dat, clus, clus2, ht, rn, cn, col, dispnum) {
    #require(pheatmap)
    #require(dendsort)
    col.pal <- np[[col]]
    if (FALSE) {
      col.pal = rev(col.pal)
    }
    # Define metrics for clustering
    drows1 <- distance_metric
    dcols1 <- distance_metric
    minx = min(dat)
    maxx = max(dat)
    if (autoscale_heatmap_color) {
      breaks = seq(minx, maxx, length=100)
      legbreaks = seq(minx, maxx, length=5)
    } else {
      breaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=100)
      legbreaks = seq(set_min_heatmap_color, set_max_heatmap_color, length=5)
    }
    breaks = sapply(breaks, signif, 4)
    legbreaks = sapply(legbreaks, signif, 4)
    # Run cluster method using 
    hc = hclust(dist(t(dat)), method=clustering_method)
    hcrow = hclust(dist(dat), method=clustering_method)
    if (FALSE) {
      sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
    } else {
      sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
    }
    if (clus) {
      colclus <- sort_hclust(hc)
    } else {
      colclus = FALSE
    }
    if (clus2) {
      rowclus <- sort_hclust(hcrow)
    } else {
      rowclus = FALSE
    }
    if (display_dendrograms) {
      treeheight <- 25
    } else {
      treeheight <- 0
    }
    
    hm.parameters <- list(
      dat, 
      color=col.pal,
      legend_breaks=legbreaks,
      legend=TRUE,
      scale="none",
      treeheight_col=treeheight,
      treeheight_row=treeheight,
      kmeans_k=NA,
      breaks=breaks,
      display_numbers=dispnum,
      number_color = "black",
      fontsize_number = 8,
      height=80,
      cellwidth = NA, 
      cellheight = NA, 
      fontsize= legend_font_size,   
      fontsize_row=gene_name_font_size,
      fontsize_col=sample_name_font_size,
      show_rownames=rn, 
      show_colnames=cn,
      cluster_rows=rowclus, 
      cluster_cols=clus,
      clustering_distance_rows=drows1, 
      clustering_distance_cols=dcols1,
      annotation_col = annotation_col,
      annotation_colors = annot_col,
      labels_col = labels_col
    )
    mat = t(dat)
    callback = function(hc, mat) {
      dend = rev(dendsort(as.dendrogram(hc)))
      if(reorder_dendrogram == TRUE) {
        dend %>% dendextend::rotate(reorder_dendrogram_order) -> dend
      } else {
        dend %>% dendextend::rotate(c(1:nobs(dend))) 
      }
      as.hclust(dend)
    }
    
    ## Make Heatmap
    phm <- do.call("pheatmap", c(hm.parameters, list(clustering_callback=callback)))
    
  }
  # End doheatmap() function.
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  ## Build different color spectra options for heatmap:
  np0 = pal(100) 
  np1 = diverge_hcl(100, c=100, l=c(30, 80), power=1) # Blue to Red
  np2 = heat_hcl(100, c=c(80, 30), l=c(30, 90), power=c(1/5, 2)) # Red to Vanilla
  np3 = rev(heat_hcl(100, h=c(0, -100), c=c(40, 80), l=c(75, 40), power=1)) # Violet to Pink
  np4 = rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(100)) #Red to yellow to blue
  np5 = colorRampPalette(c("steelblue","white", "red"))(100)  # Steelblue to White to Red
  
  ## Gather list of color spectra and give them names for the GUI to show.
  np = list(np0, np1, np2, np3, np4, np5)
  names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Bu Yl Rd","Bu Wt Rd")
  
  ## Parse input counts matrix. Subset by samples.
  df1 <- counts_matrix
  # Swap out Gene Name column name, if it's not 'Gene'.
  if(gene_column_name != "Gene"){
    # Drop original Gene column
    df1 = df1[,!(colnames(df1)%in% c("Gene")) ]
    # Rename column to Gene
    colnames(df1)[which(colnames(df1) == gene_column_name)] <- 'Gene'
  }
  # Get sample columns
  samples_to_include <- samples_to_include[samples_to_include != gene_column_name]
  samples_to_include <- samples_to_include[samples_to_include != "Gene"]
  samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
  
  # Build new counts matrix containing only sample subset chosen by user.
  df1 <- df1[,append("Gene", samples_to_include)]
  df.orig = df1
  df.orig %>% dplyr::group_by(Gene) %>% summarise_all(funs(mean)) -> df
  df.mat = df[ , (colnames(df) != "Gene" )] %>% as.data.frame
  df %>% dplyr::mutate(Gene = stringr::str_replace_all(Gene, "_", " ")) -> df
  row.names(df.mat) <- df$Gene
  rownames(df.mat) <- str_wrap(rownames(df.mat),30) #for really long geneset names
  df.mat <- as.data.frame(df.mat)
  
  ## Subset counts matrix by genes.
  # Toggle to include all genes in counts matrix (in addition to any user-submitted gene list).
  if (include_all_genes == FALSE) {
    # Add user-submitted gene list (optional).
    genes_to_include_parsed = c()
    genes_to_include_parsed = strsplit(specific_genes_to_include_in_heatmap, " ")[[1]]
    df.mat[genes_to_include_parsed,] -> df.final.extra.genes
    if(filter_top_genes_by_variance == TRUE) {
      # Want to filter all genes by variance.
      df.final = as.matrix(df.mat)
      var <- matrixStats::rowVars(df.final)
      df <- as.data.frame(df.final)
      rownames(df) <- rownames(df.final)
      df.final <- df
      df.final$var <- var
      df.final %>% rownames_to_column("Gene") -> df.final 
      df.final %>% dplyr::arrange(desc(var)) -> df.final
      df.final.extra.genes = dplyr::filter(df.final, Gene %in% genes_to_include_parsed)
      df.final = df.final[1:top_genes_by_variance_to_include,]
      df.final = df.final[complete.cases(df.final),]
      # Rbind user gene list to variance-filtered gene list and deduplicate.
      df.final <- rbind(df.final, df.final.extra.genes)
      df.final <- df.final[!duplicated(df.final),] 
      rownames(df.final) <- df.final$Gene
      df.final$Gene <- NULL
      df.final$var <- NULL
    } else {
      # Want to use ONLY user-provided gene list.
      df.final <- df.final.extra.genes
      df.final <- df.final[!duplicated(df.final),]
      # Order genes in heatmap by user-submitted order of gene names.
      df.final <- df.final[genes_to_include_parsed,]
      #df.final$Gene <- NULL
    }
  } else {
    df.final <- df.mat
    df.final$Gene <- NULL
  }
  
  ## Optionally apply centering and rescaling (default TRUE).
  if (center_and_rescale_expression == TRUE) {
    tmean.scale = t(scale(t(df.final)))
    tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
    tmean.scale = na.omit(tmean.scale)
  } else {
    tmean.scale = df.final
  }
  
  if(order_by_gene_expression == TRUE){
    gene_to_order_columns <- gsub(" ","",gene_to_order_columns)
    if(gene_expression_order == "low_to_high"){
      tmean.scale <- tmean.scale[,order(tmean.scale[gene_to_order_columns,])] #order from low to high 
    } else{
      tmean.scale <- tmean.scale[,order(-tmean.scale[gene_to_order_columns,])] #order from high to low  
    }
  }
  
  df.final <- as.data.frame(tmean.scale)
  
  ## Parse input sample metadata and add annotation tracks to top of heatmap.
  annot <- sample_metadata
  # Filter to only samples user requests.
  annot %>% dplyr::filter(.data[[sample_name_column]] %in% samples_to_include) -> annot
  # Arrange sample options.
  if(arrange_sample_columns) {
    annot %>% dplyr::arrange_(.dots=group_columns) -> annot
    df.final <- df.final[,match(annot[[sample_name_column]],colnames(df.final))] 
  }
  # Build subsetted sample metadata table to use for figure.
  
  colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
  names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
  group_colors <- colorlist[group_colors]
  
  annot %>% dplyr::select(group_columns) -> annotation_col    
  annotation_col = as.data.frame(unclass(annotation_col))
  annotation_col[] <- lapply(annotation_col,factor)
  x <- length(unlist(lapply(annotation_col,levels)))
  if(x>length(group_colors)){
    k=x-length(group_colors)
    more_cols<- getourrandomcolors(k) 
    group_colors <- c(group_colors, more_cols)
  }
  rownames(annotation_col) <- annot[[sample_name_column]]
  annot_col = list()
  b=1
  i=1
  while (i <= length(group_columns)){
    nam <- group_columns[i]
    grp <- as.factor(annotation_col[,i])
    c <- b+length(levels(grp))-1
    col = group_colors[b:c]
    names(col) <- levels(grp)
    assign(nam,col)
    annot_col = append(annot_col,mget(nam))
    b = c+1
    i=i+1
  }
  
  if(assign_group_colors == TRUE){
    colassign <- assign_color_to_sample_groups
    groupname <- c()
    groupcol <- c() 
    for (i in 1:length(colassign)) {
      groupname[i] <- strsplit(colassign[i], ": ?")[[1]][1]
      groupcol[i] <- strsplit(colassign[i], ": ?")[[1]][2]
    }
    annot_col[[1]][groupname] <- groupcol
  }
  
  ## Setting labels_col for pheatmap column labels.
  if (manually_rename_samples == TRUE) {
    # Use user-provided names to rename samples.
    replacements = samples_to_rename
    old <- c()
    new <- c()
    labels_col <- colnames(df.final)
    for (i in 1:length(replacements)) {
      old <- strsplit(replacements[i], ": ?")[[1]][1]
      new <- strsplit(replacements[i], ": ?")[[1]][2]
      old=gsub("^[[:space:]]+|[[:space:]]+$","",old)
      new=gsub("^[[:space:]]+|[[:space:]]+$","",new)
      labels_col[labels_col==old]=new           
    }
  } else {
    ## Use original column names for samples.
    labels_col <- colnames(df.final)
  }
  
  ## Print number of genes to log.
  print(paste0("The total number of genes in heatmap: ", nrow(df.final)))
  
  # manually remove ALB, don't know why it is showing up
  df.final = df.final[rownames(df.final) != 'ALB',]
  
  ## Make the final heatmap.
  p <- doheatmap(dat=df.final, clus=cluster_samples, clus2=cluster_genes, ht=50, rn=display_gene_names, cn=display_sample_names, col=heatmap_color_scheme, dispnum=display_numbers)
  p@matrix_color_mapping@name <- " "
  p@matrix_legend_param$at <- as.numeric(formatC(p@matrix_legend_param$at, 2))
  p@column_title_param$gp$fontsize <- 10
  print(p)
  
  ## If user sets toggle to TRUE, return Z-scores.
  ## Else return input counts matrix by default (toggle FALSE).
  ## Returned matrix includes only genes & samples used in heatmap.
  if(return_z_scores){
    df.new <- data.frame(tmean.scale) # Convert to Z-scores.
    df.new %>% rownames_to_column("Gene") -> df.new
    return(df.new)
  } else {
    df.final %>% rownames_to_column("Gene") -> df.new
    return(df.new)
  }
}

heatmap_pptx1_1 <- heatmap_1(Filtered_Counts = filtered_counts, 
                             Ccbr1321_metadata = meta)

write.csv(heatmap_pptx1_1, "heatmap_pptx1_1.csv")

heatmap_pptx1_2 <- heatmap_2(Filtered_Counts = filtered_counts, 
                             Ccbr1321_metadata = meta)

write.csv(heatmap_pptx1_2, "heatmap_pptx1_2.csv")


