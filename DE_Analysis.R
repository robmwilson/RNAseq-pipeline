library("edgeR")
library("AnnotationDbi")
library("org.Hs.eg.db")

load_counts <- function(counts_path, info_path){
  
  # Load an edgeR DGEList from two csvs containing transcript counts and sample info
  samples = read.csv(info_path, row.names=1)
  samples$Treatment = relevel(samples$Treatment, ref="DMSO")
  counts = read.csv(counts_path, row.names="Symbol")
  genes = data.frame(Description=mapIds(org.Hs.eg.db,
                                        keys=row.names(counts), 
                                        keytype='SYMBOL', 
                                        column='GENENAME'))
  return(DGEList(counts, samples=samples, genes=genes, remove.zeros=TRUE))
}


make_GLM <- function(data, count_threshold=5){

  # Throw out any genes that don't meet threshold counts per million in all samples
  cpmlimit = count_threshold*min(data$samples$lib.size)*count_threshold/1000000
  sprintf("CPM Limit is assigned as %f", cpmlimit)
  data = data[rowSums(cpm(data)<cpmlimit) == 0, , keep.lib.sizes=FALSE]

  # Make a binary design matrix for GLM on the basis of treatment only. 
  design = model.matrix(~data$samples$Treatment, data=data$samples)
  colnames(design) = levels(data$samples$Treatment)
  
  # Calculate normalization factors by TMM, so that samples are scaled to a normalized number of reads
  data = calcNormFactors(data)
  
  # Estimate gene dispersion according to the design matrix.
  data = estimateDisp(data, design)
  
  # Fit a general linear model, i.e. a system of equations that predict how gene expression relates to treatment
  fit = glmFit(data, design)
  
  return(fit)
}
