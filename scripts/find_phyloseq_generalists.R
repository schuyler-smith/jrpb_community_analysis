find_generalists <- function(phyloseq_obj, frequency = 0.2, treatments = NULL, drop_samples = FALSE){
  require(phyloseq)
  require(data.table)
  if(!(is.null(treatments))){
    treatment <- colnames(sample_data(phyloseq_obj)[,treatments])
    Treatment_Groups <- setDT(as(sample_data(phyloseq_obj)[,colnames(sample_data(phyloseq_obj)) %in% treatment], "data.frame"))
    # Treatment_Groups[, Treatment_Group := .GRP, by = Treatment_Groups]
    eval(parse(text=paste0("Treatment_Groups[, Treatment_Group := paste(", paste(treatment, collapse = ", "), ", sep = '-'), by = Treatment_Groups]")))
    sample_data(phyloseq_obj) <- cbind(sample_data(phyloseq_obj), data.frame(Treatment_Groups[,"Treatment_Group"]))
    Treatment_Groups <- unique(Treatment_Groups[,"Treatment_Group"])
    phyloseq_obj <- phyloseq(otu_table(phyloseq_obj), tax_table(phyloseq_obj), sample_data(phyloseq_obj))
    phyloseq_obj <- do.call(merge_phyloseq, 
                            apply(array(Treatment_Groups), 1, FUN = function(i){
                              sub_phy <- eval(parse(text=paste0("subset_samples(phyloseq_obj, Treatment_Group == '",i,"')")))
                              cutoff <- floor(ncol(otu_table(sub_phy)) * frequency)
                              sub_phy <- filter_taxa(sub_phy, function(x){sum(x != 0) >= cutoff}, TRUE)
                              return(sub_phy)})
    )
  }else {
    cutoff <- floor(ncol(otu_table(phyloseq_obj)) * frequency)
    phyloseq_obj <- filter_taxa(phyloseq_obj, function(x){sum(x != 0) >= cutoff}, TRUE)
  }
  if(drop_samples == TRUE){
    phyloseq_obj <- subset_samples(phyloseq_obj, sample_sums(phyloseq_obj) > 0)
  }
  return(phyloseq_obj)
}
