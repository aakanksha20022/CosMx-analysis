----#FUNCTIONS------
#Function to process the expression matrix
process_em= function(em_path)
{
  #Reading in the em matrix
  em_mtx= read.csv(em_path, row.names=1, header=T)
  
  #cleaning column names
  colnames(em_mtx) <- gsub("\\.", "_", colnames(em_mtx))
  
  #Mapping gene symbols to ensembl IDs
  ensembl_ids <- mapIds(org.Hs.eg.db, keys = rownames(em_mtx), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  na_genes <- which(is.na(ensembl_ids))
  if (length(na_genes) > 0) {
    cat("Unmapped genes:", rownames(em_mtx)[na_genes], "\n")
  }
  
  # Ensuring unique Ensembl IDs by using make.names
  ensembl_ids_unique <- make.names(ensembl_ids, unique = TRUE)
  # Updating rownames of epithelial_mtx with unique Ensembl IDs
  rownames(em_mtx)= ensembl_ids_unique
  
  # makes sure the batch and group information is a factor
  batch = cut(colSums(em_mtx),2)
  
  # Correct the EM using CombatSeq
  sample_group <- col_data$sample_group
  counts_corrected = ComBat_seq(as.matrix(em_mtx), batch=batch, group=sample_group)
  return (counts_corrected)
}

#Function to perform de
do_de= function(counts_corrected, col_data)
{
  # Create DESeq2 dataset object
  dds <- DESeqDataSetFromMatrix(countData = counts_corrected,
                                colData = col_data,
                                design = ~ sample_group + patient_id)
  
  # Perform pre-filtering for gene counts
  keep_gene=rowSums(counts(dds) >= 10) >= 3
  dds=dds[keep_gene,]
  
  # Run DESeq2 pipeline
  dds <- DESeq(dds)
  return (dds)
}
  
# Extract results
get_res= function(dds)
{
  res <- results(dds, c('sample_group' , 'positive', 'negative'))
  
  # Sort results by adjusted p-value (FDR)
  res <- res[order(res$padj), ]
  # Return as a list
  return(res)
}

#function to format de table
format_de= function(res)
{
  res= res[, !colnames(res) %in% "baseMean"]
  res= res[, !colnames(res) %in% "lfcSE"]
  res= res[, !colnames(res) %in% "stat"]
  res$ID= rownames(res)
  res <- res[, c(ncol(res), 1:(ncol(res)-1))]
  rownames(res)= NULL
  names(res)= c('ID', 'Log2Fold', 'P', 'P.Adj')
  return(res)
}

#Function to format the em table
format_norm_em= function(dds)
{
  resdata = round(as.data.frame(counts(dds, normalized=TRUE)),2)
  resdata_IDs = resdata
  resdata_IDs$IDs = row.names(resdata)
  resdata_IDs = resdata_IDs[,ncol(resdata_IDs)]
  resdata =  cbind(resdata_IDs,resdata)
  colnames(resdata)[1] = "ID"
  return (resdata)
}