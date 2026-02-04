#' @title  Construct a Multi Assay Experiment for ELMER analysis
#' @description 
#' This function will receive a gene expression and DNA methylation data objects 
#' and create a Multi Assay Experiment.
#' @param met A Summarized Experiment with one assay containing beta-values, 
#' a matrix or path of rda file only containing the data.
#' @param exp A Summarized Experiment with one assay, or
#' a matrix or path of rda file only containing the data. Rownames should be 
#' either Ensembl gene id (ensembl_gene_id) or gene symbol (external_gene_name)
#' @param genome Which is the default genome to make gene information. Options hg19 and hg38
#' @param colData A DataFrame or data.frame of the phenotype data for all participants. Must have column primary (sample ID).
#' @param sampleMap  A DataFrame or data.frame of the matching samples and colnames
#'  of the gene expression and DNA methylation matrix. This should be used if your matrix
#'  have different columns names. 
#'  This object must have following columns: 
#'  assay ("DNA methylation" and "Gene expression"), primary (sample ID) and colname (names of the columns of the matrix).
#' @param linearize.exp Take log2(exp + 1) in order to linearize relation between methylation and expression  
#' @param met.platform DNA methylation platform "450K" or "EPIC"
#' @param TCGA A logical. FALSE indicate data is not from TCGA (FALSE is default). 
#' TRUE indicates data is from TCGA and sample section will automatically filled in.
#' @param filter.probes A GRanges object contains the coordinate of probes which locate 
#'  within promoter regions or distal feature regions such as union enhancer from REMC and FANTOM5.
#'  See \code{\link{get.feature.probe}} function.
#' @param filter.genes List of genes ensemble ids to filter from object  
#' @param save If TRUE, MAE object will be saved into a file named as the argument save.file if this was set, otherwise as mae_genome_met.platform.rda.
#' @param save.filename Name of the rda file to save the object (must end in .rda)
#' @param met.na.cut Define the percentage of NA that the line should have to remove the probes 
#' for humanmethylation platforms.  
#' @return A MultiAssayExperiment object
#' @export 
#' @importFrom MultiAssayExperiment MultiAssayExperiment 
#' @importFrom SummarizedExperiment SummarizedExperiment makeSummarizedExperimentFromDataFrame assay assay<-
#' @examples
#' # NON TCGA example: matrices has different column names
#' gene.exp <- S4Vectors::DataFrame(
#'    sample1.exp = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
#'    sample2.exp = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3)
#' )
#'                   
#' dna.met <- S4Vectors::DataFrame(
#'    sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
#'    sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9)
#' )
#' sample.info <- S4Vectors::DataFrame(
#'   primary =  c("sample1","sample2"), 
#'   sample.type = c("Normal", "Tumor")
#'  )
#'  
#' sampleMap <- S4Vectors::DataFrame(
#'     assay = c("Gene expression","DNA methylation","Gene expression","DNA methylation"),
#'     primary = c("sample1","sample1","sample2","sample2"), 
#'     colname = c("sample1.exp","sample1.met","sample2.exp","sample2.met")
#' )
#' 
#' mae <- createMAE(
#'     exp = gene.exp, 
#'     met = dna.met, 
#'     sampleMap = sampleMap, 
#'     met.platform ="450K",
#'     colData = sample.info, 
#'     genome = "hg38"
#' ) 
#'                  
#' # You can also use sample Mapping and Sample information tables from a tsv file
#' # You can use the createTSVTemplates function to create the tsv files
#' readr::write_tsv(as.data.frame(sampleMap), path = "sampleMap.tsv")
#' readr::write_tsv(as.data.frame(sample.info), path = "sample.info.tsv")
#' 
#' mae <- createMAE(
#'    exp = gene.exp, 
#'    met = dna.met, 
#'    sampleMap = "sampleMap.tsv", 
#'    met.platform ="450K",
#'    colData = "sample.info.tsv", 
#'    genome = "hg38"
#' ) 
#'                  
#' # NON TCGA example: matrices has same column names
#' gene.exp <- S4Vectors::DataFrame(sample1 = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
#'                   sample2 = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3))
#' dna.met <- S4Vectors::DataFrame(sample1 = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2=  c("cg14324200"=0.3,"cg23867494"=0.9))
#' sample.info <- S4Vectors::DataFrame(primary =  c("sample1","sample2"), 
#'                                     sample.type = c("Normal", "Tumor"))
#' sampleMap <- S4Vectors::DataFrame(
#'      assay = c("Gene expression","DNA methylation","Gene expression","DNA methylation"),
#'     primary = c("sample1","sample1","sample2","sample2"), 
#'     colname = c("sample1","sample1","sample2","sample2")
#' )
#' mae <- createMAE(
#'   exp = gene.exp, 
#'   met = dna.met, 
#'   sampleMap = sampleMap, 
#'   met.platform ="450K",
#'   colData = sample.info, 
#'   genome = "hg38"
#' ) 
#' 
#' \dontrun{
#'    # TCGA example using TCGAbiolinks
#'    # Testing creating MultyAssayExperiment object
#'    # Load library
#'    library(TCGAbiolinks)
#'    library(SummarizedExperiment)
#'    
#'    samples <- c(
#'       "TCGA-BA-4074", "TCGA-BA-4075", "TCGA-BA-4077", "TCGA-BA-5149",
#'        "TCGA-UF-A7JK", "TCGA-UF-A7JS", "TCGA-UF-A7JT", "TCGA-UF-A7JV"
#'    )
#'    
#'    #1) Get gene expression matrix
#'    query.exp <- GDCquery(
#'        project = "TCGA-HNSC", 
#'        data.category = "Transcriptome Profiling", 
#'        data.type = "Gene Expression Quantification", 
#'        workflow.type = "STAR - Counts",
#'        barcode = samples
#'    )
#'    
#'    GDCdownload(query.exp)
#'    exp.hg38 <- GDCprepare(query = query.exp)
#'    
#'    # DNA Methylation
#'    query.met <- GDCquery(
#'       project = "TCGA-HNSC",
#'       data.category = "DNA Methylation", 
#'       data.type = "Methylation Beta Value",
#'       barcode = samples,
#'       platform = "Illumina Human Methylation 450"
#'    )
#'    
#'    GDCdownload(query.met)
#'    met <- GDCprepare(query = query.met)
#'    
#'    distal.enhancer <- get.feature.probe(
#'       genome = "hg19",
#'       met.platform = "450k"
#'    )
#'    
#'    
#'    mae.hg38 <- createMAE(
#'      exp = exp.hg38, met = met, 
#'      TCGA = TRUE, genome = "hg38",  
#'      filter.probes = distal.enhancer
#'    )
#'    values(getExp(mae.hg38))
#'    
#'    # Consisering it is TCGA and not SE
#'    mae.hg19.test <- createMAE(
#'      exp = assay(exp.hg19), met =  assay(met), 
#'      TCGA = TRUE, genome = "hg19",  
#'      filter.probes = distal.enhancer
#'    )
#'    
#'    mae.hg38 <- createMAE(
#'      exp = assay(exp.hg38), 
#'      met = assay(met), 
#'      TCGA = TRUE, 
#'      genome = "hg38",  
#'      filter.probes = distal.enhancer
#'    )
#'    values(getExp(mae.hg38))
#'    
#'    # Consisering it is not TCGA and SE
#'    # DNA methylation and gene expression Objects should have same sample names in columns
#'    not.tcga.exp <- exp.hg19 
#'    colnames(not.tcga.exp) <- substr(colnames(not.tcga.exp),1,15)
#'    not.tcga.met <- met 
#'    colnames(not.tcga.met) <- substr(colnames(not.tcga.met),1,15)
#'    
#'    phenotype.data <- data.frame(
#'       row.names = colnames(not.tcga.exp),
#'       primary = colnames(not.tcga.exp),
#'       samples = colnames(not.tcga.exp), 
#'       group = c(rep("group1",4),rep("group2",4))
#'    )
#'    
#'    distal.enhancer <- get.feature.probe(genome = "hg19",met.platform = "450k")
#'    mae.hg19 <- createMAE(
#'       exp = not.tcga.exp, 
#'       met =  not.tcga.met, 
#'       TCGA = FALSE, 
#'       filter.probes = distal.enhancer,
#'       genome = "hg19", 
#'       colData = phenotype.data
#'    )
#' }
#' createMAE
createMAE <- function (
  exp, 
  met, 
  colData, 
  sampleMap,
  linearize.exp = FALSE,
  filter.probes = NULL,
  met.na.cut = 0.2,
  filter.genes = NULL,
  met.platform = "450K",
  genome = NULL,
  save = TRUE,
  save.filename,
  TCGA = FALSE
) {
  
  if(missing(genome)) stop("Please specify the genome (hg38, hg19)")
  
  
  # Check if input are path to rda files
  if(is.character(exp)) exp <- get(load(exp))
  if(is.character(met)) met <- get(load(met))
  
  if(is(exp,"RangedSummarizedExperiment")){
    
    
    if(length(SummarizedExperiment::assays(exp)) > 1) {
      if("tpm_unstrand" %in% names(SummarizedExperiment::assays(exp))){
        message("Multiple assay found. Selecting tpm_unstrand")
        exp <- SummarizedExperiment::assays(exp)$tpm_unstrand
      } else {
        stop("RNA-seq object has more than one assay. Please make sure it is only")
      }
    }
  }
  
  suppressMessages({
    
    if(!missing(colData)) { 
      if(is.character(colData)) { 
        colData <- as.data.frame(read_tsv(colData))
      }
      if (!"primary" %in% colnames(colData)) stop("No primary column in colData input")
      rownames(colData) <- colData$primary
    }
    if(!missing(sampleMap)) { 
      if(is.character(sampleMap)) sampleMap <- read_tsv(sampleMap)
      if (!all(c("assay","colname","primary") %in% colnames(sampleMap))) 
        stop("All assay, primary and colname columns should be in sampleMap input")
    }
  })  
  
  # Expression data must have the ensembl_gene_id (Ensemble ID) and external_gene_name (Gene Symbol)
  required.cols <- c("external_gene_name", "ensembl_gene_id")
  # If my input is a data frame we will need to add metadata information for the ELMER analysis steps
  if(!is(exp,"RangedSummarizedExperiment")){
    exp <- makeSummarizedExperimentFromGeneMatrix(exp, genome)
  }
  # Add this here ?
  if(linearize.exp) assay(exp) <- log2(assay(exp) + 1)
  
  if(!is(met,"RangedSummarizedExperiment")){
    met <- makeSummarizedExperimentFromDNAMethylation(
      met = met, genome = genome, met.platform = met.platform
    )
  }
  met <- met[rowMeans(is.na(assay(met))) < met.na.cut, ]
  
  # Select the regions from DNA methylation that overlaps enhancer.
  if(!is.null(filter.probes)){
    if(is.character(filter.probes)){
      filter.probes <- get(load(filter.probes))
    }
  } 
  if(!is.null(filter.probes) & !is.null(met)){
    met <- met[rownames(met) %in% names(filter.probes),]
  }
  if(!is.null(filter.genes) & !is.null(exp)){
    exp <- exp[rownames(exp) %in% names(filter.genes),]
  } 
  
  # We will need to check if the fields that we need exists.
  # Otherwise we will need to create them
  if(!is(exp,"RangedSummarizedExperiment")){
    required.cols <- required.cols[!required.cols %in% colnames(values(exp))]
    if(length(required.cols) > 0) {
      gene.info <- get.GRCh(genome)
      colnames(gene.info)[grep("external_gene", colnames(gene.info))] <- "external_gene_name"
      if(all(grepl("ENSG",rownames(exp)))) {
        extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$ensembl_gene_id),required.cols])
        colnames(extra) <- required.cols
        values(exp) <- cbind(values(exp),extra)
      } else {
        stop("Please the gene expression matrix should receive ENSEMBLE IDs")
      }
    }
  } 
  if(TCGA){
    message("Checking if samples have both DNA methylation and Gene expression and if they are in the same order...")
    # If it is not TCGA we will assure the sample has both DNA methylation and gene expression
    ID <- intersect(substr(colnames(met),1,16), substr(colnames(exp),1,16))
    
    # Get only samples with both DNA methylation and Gene expression
    met <- met[,match(ID,substr(colnames(met),1,16))]
    exp <- exp[,match(ID,substr(colnames(exp),1,16))]
    stopifnot(all(substr(colnames(exp),1,16) == substr(colnames(met),1,16)))
    stopifnot(ncol(exp) == ncol(met))
    
    # Get clinical information
    if(missing(colData)) {
      colData <- TCGAbiolinks::colDataPrepare(colnames(met))
      # This will keep the same strategy the old ELMER version used:
      # Every type of tumor samples (starts with T) will be set to tumor and
      # every type of normal samples   (starts with N) will be set to normal 
      # See : https://github.com/lijingya/ELMER/blob/3e050462aa41c8f542530ccddc8fa607207faf88/R/Small.R#L8-L48
      colData$TN <- NA
      colData[grep("^N",colData$shortLetterCode),"TN"] <- "Normal" 
      colData[grep("^T",colData$shortLetterCode),"TN"] <- "Tumor" 
      
      colData$barcode <- NULL
      colData <- colData[!duplicated(colData),]      
      rownames(colData) <- colData$sample
    } 
    if(missing(sampleMap)) {
      sampleMap <- DataFrame(assay = c(rep("DNA methylation", length(colnames(met))), rep("Gene expression", length(colnames(exp)))),
                             primary = substr(c(colnames(met),colnames(exp)),1,16),
                             colname = c(colnames(met),colnames(exp)))
    }
    
    message("Creating MultiAssayExperiment")
    mae <- MultiAssayExperiment(
      experiments = list("DNA methylation" = met,"Gene expression" = exp),
      colData = colData,   
      sampleMap = sampleMap,
      metadata = list(TCGA= TRUE, genome = genome, met.platform = met.platform )
    )
  } else {
    
    if(missing(colData)){
      message <- paste(
        "Please set colData argument. A data frame with samples", 
        "information. All rownames should be colnames of DNA",
        "methylation and gene expression. An example is showed",
        "in MultiAssayExperiment documentation",
        "(access it with ?MultiAssayExperiment)"
      )
      stop(message)
    }
    
    if(missing(sampleMap)){
      # Check that we have the same number of samples
      message("Removing samples not found in both DNA methylation and gene expression (we are considering the names of the gene expression and DNA methylation columns to be the same) ")
      ID <- intersect(colnames(met), colnames(exp))
      met <- met[,match(ID,colnames(met))]
      exp <- exp[,match(ID,colnames(exp))]
      
      if(!all(colnames(exp) == colnames(met))) 
        stop("Error DNA methylation matrix and gene expression matrix are not in the same order")
      
      colData <- colData[match(ID,rownames(colData)),,drop = FALSE]
      sampleMap <- DataFrame(
        assay = c(rep("DNA methylation", length(colnames(met))), 
                  rep("Gene expression", length(colnames(exp)))),
        primary = c(colnames(met),colnames(exp)),
        colname = c(colnames(met),colnames(exp))
      )
      mae <- MultiAssayExperiment(
        experiments = list("DNA methylation" = met,"Gene expression" = exp),
        colData = colData,
        sampleMap = sampleMap,
        metadata = list(TCGA = FALSE, genome = genome, met.platform = met.platform )
      )
    } else {
      # Check that we have the same number of samples
      if(!all(c("primary","colname") %in% colnames(sampleMap))) 
        stop("sampleMap should have the following columns: primary (sample ID) and colname(DNA methylation and gene expression sample [same as the colnames of the matrix])")
      #if(!any(rownames(colData) %in% sampleMap$primary))
      #  stop("colData row names should be mapped to sampleMap primary column ")
      # Find which samples are DNA methylation and gene expression
      sampleMap.met <- sampleMap[sampleMap$assay %in% "DNA methylation",,drop = FALSE]
      sampleMap.exp <- sampleMap[sampleMap$assay %in% "Gene expression",,drop = FALSE]
      
      # Which ones have both DNA methylation and gene expression ?
      commun.samples <- intersect(sampleMap.met$primary,sampleMap.exp$primary)
      
      # Remove the one that does not have both data
      sampleMap.met <- sampleMap.met[match(sampleMap.met$primary,commun.samples),,drop = FALSE]
      sampleMap.exp <- sampleMap.exp[match(sampleMap.exp$primary,commun.samples),,drop = FALSE]
      
      # Ordering samples to be matched
      met <- met[,sampleMap.met$colname,drop = FALSE]
      exp <- exp[,sampleMap.exp$colname,drop = FALSE]
      
      if(!all(sampleMap.met$primary == sampleMap.exp$primary)) 
        stop("Error DNA methylation matrix and gene expression matrix are not in the same order")
      
      colData <- colData[match(commun.samples,colData$primary),,drop = FALSE]
      sampleMap <- DataFrame(
        assay = c(rep("DNA methylation", length(colnames(met))), 
                  rep("Gene expression", length(colnames(exp)))),
        primary = commun.samples,
        colname = c(colnames(met),colnames(exp))
      )
      mae <- MultiAssayExperiment(
        experiments=list("DNA methylation" = met,
                         "Gene expression" = exp),
        colData = colData,
        sampleMap = sampleMap,
        metadata = list(TCGA=FALSE, genome = genome, met.platform = met.platform )
      )
    }
  }
  if(save) {
    if(missing(save.filename)) save.filename <- paste0("mae_",genome,"_",met.platform,".rda")
    save(mae, file = save.filename,compress = "xz")
    message("MAE saved as: ", save.filename)
  }
  return(mae)
}

makeSummarizedExperimentFromGeneMatrix <- function(exp, genome = genome){
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from gene expression input")
  gene.info <- get.GRCh(genome)
  gene.info$chromosome_name <- paste0("chr",gene.info$chromosome_name)
  colnames(gene.info)[grep("external_gene", colnames(gene.info))] <- "external_gene_name"
  gene.info$strand[gene.info$strand == 1] <- "+"
  gene.info$strand[gene.info$strand == -1] <- "-"
  exp <- as.data.frame(exp)
  required.cols <- c("external_gene_name", "ensembl_gene_id")
  
  if(all(grepl("ENSG",rownames(exp)))) {
    exp$ensembl_gene_id <- gsub("\\.[0-9]*$","",rownames(exp))
    aux <- merge(exp, gene.info, by = "ensembl_gene_id", sort = FALSE)
    aux <- aux[!duplicated(aux$ensembl_gene_id),]
    rownames(aux) <- aux$ensembl_gene_id
    aux[,grep("entrezgene",colnames(aux))] <- NULL
    exp <- makeSummarizedExperimentFromDataFrame(
      aux[,!grepl("external_gene_name|ensembl_gene_id|entrezgene",colnames(aux))],    
      start.field = "start_position",
      end.field = c("end_position")
    )
    extra <- as.data.frame(gene.info[match(rownames(exp),gene.info$ensembl_gene_id),required.cols])
    colnames(extra) <- required.cols
    values(exp) <- cbind(values(exp),extra)
  } else {
    stop("Please the gene expression matrix should receive ENSEMBLE IDs (ENSG)")
  }
  return(exp)
}

#' @importFrom downloader download
#' @importFrom S4Vectors DataFrame
makeSummarizedExperimentFromDNAMethylation <- function(met, genome, met.platform) {
  message("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-")
  message("Creating a SummarizedExperiment from DNA methylation input")
  
  # Instead of looking on the size, it is better to set it as a argument as the annotation is different
  annotation <-   getInfiniumAnnotation(plat = met.platform, genome = genome)
  
  rowRanges <- annotation[names(annotation) %in% rownames(met),,drop=FALSE]
  
  # Remove masked probes, besed on the annotation
  rowRanges <- rowRanges[!rowRanges$MASK_general]
  
  colData <-  DataFrame(samples = colnames(met))
  met <- met[rownames(met) %in% names(rowRanges),,drop = FALSE]
  met <- met[names(rowRanges),,drop = FALSE]
  assay <- data.matrix(met)
  met <- SummarizedExperiment(
    assays = assay,
    rowRanges = rowRanges,
    colData = colData
  )
  return(met)
}


# getInfiniumAnnotation <- function(plat = "450K", genome = "hg38"){
#   if(tolower(genome) == "hg19" & toupper(plat) == "450K" ) return(sesameData::sesameDataGet("HM450.hg19.manifest"))
#   if(tolower(genome) == "hg19" & toupper(plat) == "EPIC" ) return(sesameData::sesameDataGet("EPIC.hg19.manifest"))
#   if(tolower(genome) == "hg38" & toupper(plat) == "450K" ) return(sesameData::sesameDataGet("HM450.hg38.manifest"))
#   if(tolower(genome) == "hg38" & toupper(plat) == "EPIC" ) return(sesameData::sesameDataGet("EPIC.hg19.manifest"))
# }

#' @title Get DNA methylation array metadata from SesameData
#' @noRd
#' @importFrom stringr str_c
getInfiniumAnnotation <- function(
  genome = c("hg38","hg19"),
  plat = c("450K","EPIC")
){
  genome <- match.arg(genome)
  arrayType <- match.arg(plat)
  
  check_package("sesameData")
  check_package("ExperimentHub")
  check_package("AnnotationHub")
  
  message("Accessing DNAm annotation from sesame package for: ", genome, " - ",arrayType)
  manifest <-  str_c(
    ifelse(arrayType == "450K","HM450","EPIC"),
    ".",
    genome,
    ".manifest"
  )
  ehub <- ExperimentHub::ExperimentHub()
  query <- AnnotationHub::query(ehub, c("sesameData",manifest))
  query <- query[query$title == manifest,]
  ah_id <- query$ah_id[query$rdatadateadded == max(as.Date(query$rdatadateadded))]
  ehub[[ah_id]]
}

#' @title register cores
#' @param package Package name
#' @noRd
check_package <- function(package){
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(package, " package is needed for this function to work. Please install it.",
         call. = FALSE)
  }
}

getdata <- function(...)
{
  e <- new.env()
  name <- data(..., package = "ELMER.data",envir = e)[1]
  e[[ls(envir = e)[1]]]
}

#' Create examples files for Sample mapping and information used in createMAE function
#' @description 
#' This function will receive the DNA methylation and gene expression matrix and will create
#' some examples of table for the argument colData and sampleMap used in ceeateMae function.
#' @param met DNA methylation matrix or Summarized Experiment
#' @param exp Gene expression matrix or Summarized Experiment
#' @examples 
#' gene.exp <- S4Vectors::DataFrame(sample1.exp = c("ENSG00000141510"=2.3,"ENSG00000171862"=5.4),
#'                   sample2.exp = c("ENSG00000141510"=1.6,"ENSG00000171862"=2.3))
#' dna.met <- S4Vectors::DataFrame(sample1.met = c("cg14324200"=0.5,"cg23867494"=0.1),
#'                        sample2.met =  c("cg14324200"=0.3,"cg23867494"=0.9))
#' createTSVTemplates(met = dna.met, exp = gene.exp)                       
#' @importFrom readr write_tsv
#' @export
createTSVTemplates <- function(met, exp) {
  assay <- c(rep("DNA methylation", ncol(met)),
             rep("Gene expression", ncol(exp)))
  primary <- rep("SampleX", ncol(met) + ncol(exp))
  colname <- c(colnames(met),colnames(exp))
  sampleMap <- data.frame(assay,primary,colname)
  message("===== Sample mapping example file ======")
  message("Saving example file as elmer_example_sample_mapping.tsv.")
  message("Please, fill primary column correctly")
  write_tsv(sampleMap,path = "elmer_example_sample_mapping.tsv")
  
  colData <- data.frame(primary = paste0("sample",1:ncol(met)), group = rep("To be filled",ncol(met)))
  message("===== Sample information example file ======")
  message("Saving example file as elmer_example_sample_metadata.tsv.")
  message("Please, fill primary column correctly, also you can add new columns as the example group column.")
  write_tsv(colData,path = "elmer_example_sample_metadata.tsv")
}

# splitmatix 
# @param x A matrix 
# @param by A character specify if split the matix by row or column.
# @return A list each of which is the value of each row/column in the matrix.
splitmatrix <- function(x,by="row") {
  if(by %in% "row"){
    out <- split(x, rownames(x))
  }else if (by %in% "col"){
    out <- split(x, colnames(x))
  }
  return(out)
}


#' lable linear regression formula 
#' @param df A data.frame object contains two variables: dependent 
#' variable (Dep) and explanation variable (Exp).
#' @param Dep A character specify dependent variable. The first column 
#' will be dependent variable as default.
#' @param Exp A character specify explanation variable. The second column 
#' will be explanation variable as default.
#' @return A linear regression formula
#' @importFrom stats coef lm
lm_eqn = function(df,Dep,Exp){
  if(missing(Dep)) Dep <- colnames(df)[1]
  if(missing(Exp)) Exp <- colnames(df)[2]
  m = lm(df[,Dep] ~ df[,Exp]);
  eq <- substitute(italic(y) == a + (b) %.% italic(x)*"\n"~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}



#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @description 
#' getTSS to fetch GENCODE gene annotation (transcripts level) from Bioconductor package biomaRt
#' If upstream and downstream are specified in TSS list, promoter regions of GENCODE gene will be generated.
#' @param TSS A list. Contains upstream and downstream like TSS=list(upstream, downstream).
#'  When upstream and downstream is specified, coordinates of promoter regions with gene annotation will be generated.
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @return GENCODE gene annotation if TSS is not specified. Coordinates of GENCODE gene promoter regions if TSS is specified.
#' @examples
#' # get GENCODE gene annotation (transcripts level)
#' \dontrun{
#'     getTSS <- getTSS()
#'     getTSS <- getTSS(genome.build = "hg38", TSS=list(upstream=1000, downstream=1000))
#' }
#' @export
#' @author Lijing Yao (maintainer: lijingya@usc.edu)
#' @import Seqinfo
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom biomaRt useEnsembl
getTSS <- function(genome = "hg38",
                   TSS = list(upstream = NULL, downstream = NULL)){
  
  if (tolower(genome) == "hg38") {
    tss <- getdata("Human_genes__GRCh38_p12__tss")
  } else {
    tss <- getdata("Human_genes__GRCh37_p13__tss")
  }
  
  # tries <- 0L
  # msg <- character()
  # while (tries < 3L) {
  #   tss <- tryCatch({
  #     host <- ifelse(genome == "hg19",  "grch37.ensembl.org","www.ensembl.org")
  #     message("Accessing ", host, " to get TSS information")
  #     
  #     ensembl <- tryCatch({
  #       useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", host =  host)
  #     },  error = function(e) {
  #       message(e)
  #       for(mirror in c("asia","useast","uswest")){
  #         x <- useEnsembl("ensembl",
  #                         dataset = "hsapiens_gene_ensembl",
  #                         mirror = mirror,
  #                         host =  host)
  #         if(class(x) == "Mart") {
  #           return(x)
  #         }
  #       }
  #       return(NULL)
  #     })
  #     
  #     if(is.null(host)) {
  #       message("Problems accessing ensembl database")
  #       return(NULL)
  #     }
  #     attributes <- c("chromosome_name",
  #                     "start_position",
  #                     "end_position", "strand",
  #                     "ensembl_gene_id", 
  #                     "transcription_start_site",
  #                     "transcript_start",
  #                     "ensembl_transcript_id",
  #                     "transcript_end",
  #                     "external_gene_name")
  #     chrom <- c(1:22, "X", "Y","M","*")
  #     db.datasets <- listDatasets(ensembl)
  #     description <- db.datasets[db.datasets$dataset=="hsapiens_gene_ensembl",]$description
  #     message(paste0("Downloading transcripts information from ", ensembl@host, ". Using: ", description))
  #     
  #     filename <-  paste0(gsub("[[:punct:]]| ", "_",description),"_tss.rda")
  #     if(!file.exists(filename)) {
  #       tss <- getBM(attributes = attributes, filters = c("chromosome_name"), values = list(chrom), mart = ensembl)
  #       tss <- tss[!duplicated(tss$ensembl_transcript_id),]
  #       save(tss, file = filename, compress = "xz")
  #     } else {
  #       message("Loading from disk")
  #       tss <- get(load(filename))  
  #     } 
  tss$chromosome_name <-  paste0("chr", tss$chromosome_name)
  tss$strand[tss$strand == 1] <- "+" 
  tss$strand[tss$strand == -1] <- "-" 
  tss <- makeGRangesFromDataFrame(tss,
                                  strand.field = "strand",
                                  start.field = "transcript_start", 
                                  end.field = "transcript_end", 
                                  keep.extra.columns = TRUE)
  
  if (!is.null(TSS$upstream) & !is.null(TSS$downstream)) 
    tss <- promoters(tss, upstream = TSS$upstream, downstream = TSS$downstream)
  #    tss
  #  }, error = function(e) {
  #    msg <<- conditionMessage(e)
  #    tries <<- tries + 1L
  #  })
  #  if(!is.null(tss)) break
  #}
  #if (tries == 3L) stop("failed to get URL after 3 tries:", "\n  error: ", msg)
  
  return(tss)
}

#' @title  Get human TF list from the UNiprot database
#' @description This function gets the last version of human TF list from the UNiprot database
#' @importFrom readr read_tsv
#' @return A data frame with the ensemble gene id.
getTF <- function(){
  print.header("Downloading TF list from Lambert, Samuel A., et al. The human transcription factors. Cell 172.4 (2018): 650-665.","subsection")
  # human.TF <- readr::read_table("http://humantfs.ccbr.utoronto.ca/download/v_1.01/TFs_Ensembl_v_1.01.txt",col_names = F)
  # colnames(human.TF) <- "ensembl_gene_id"
  human.TF <- getdata("human.TF")
  return(human.TF)
}

#' @importFrom biomaRt getBM useMart listDatasets
get.GRCh <- function(genome = "hg19", genes = NULL, as.granges = FALSE) {
  
  if (tolower(genome) == "hg38") {
    gene.location <- getdata("Human_genes__GRCh38_p12")
  } else {
    gene.location <- getdata("Human_genes__GRCh37_p13")
  }
  
  if (!is.null(genes))  
    gene.location <- gene.location[match(genes,gene.location$entrezgene),]
  
  if (as.granges) {
    gene.location$strand[gene.location$strand == 1] <- "+"
    gene.location$strand[gene.location$strand == -1] <- "-"
    gene.location$chromosome_name <- paste0("chr", gene.location$chromosome_name)
    gene.location <- makeGRangesFromDataFrame(gene.location, 
                                              seqnames.field = "chromosome_name", 
                                              start.field = "start_position", 
                                              end.field = "end_position", 
                                              keep.extra.columns = TRUE)
  }
  return(gene.location)
}



#' @title Get family of transcription factors
#' @description 
#' This will output a list each TF motif and TFs that binding the motis. Multiple TFs may
#' recognize a same motif such as TF family.  
#' The association between each motif famil and transcription factor was created using the 
#' (HOCOMOCO)[https://hocomoco14.autosome.org/human/mono?full=true] which TF structural families 
#' was created according to TFClass [@wingender2014tfclass]
#' This data is stored as a list whose elements 
#' are motifs and contents for each element are TFs which recognize the same motif that
#' is the name of the element. This data is used in function get.TFs in \pkg{ELMER} 
#' to identify the real regulator TF whose motif is enriched in a given set of probes 
#' and expression associate with average DNA methylation of these motif sites.
#' @importFrom rvest html_table  %>%
#' @importFrom xml2 read_html 
#' @export
#' @param classification Select if we will use Family classification or sub-family
#' @return A list of TFs and its family members
createMotifRelevantTfs <- function(classification = "family"){
  message("Retrieving TFClass ", classification," classification from ELMER.data.")
  if(classification == "family") motif.relevant.TFs <- getdata("TF.family")
  if(classification == "subfamily") motif.relevant.TFs <- getdata("TF.subfamily")
  return(motif.relevant.TFs)
}

#' @title Filtering probes
#' @description 
#' This function has some filters to the DNA methylation data
#' in each it selects probes to avoid correlations due to non-cancer 
#' contamination and for additional stringency.
#'  \itemize{
#' \item Filter 1: We usually call locus unmethylated when the methylation value < 0.3 and methylated when the methylation value > 0.3. 
#'       Therefore Meth_B is the percentage of methylation value > K. 
#'       Basically, this step will make sure we have at least a percentage of beta values lesser than K and n percentage of beta values greater K. 
#'       For example, if percentage is 5\%, the number of samples 100 and K = 0.3, 
#'       this filter will select probes that we have at least 5 (5\% of 100\%) samples have beta values > 0.3 and at least 5 samples have beta values < 0.3.
#'       This filter is importante as true promoters and enhancers usually have a pretty low value (of course purity can screw that up).  
#'       we often see lots of PMD probes across the genome with intermediate values like 0.4.  
#'       Choosing a value of 0.3 will certainly give some false negatives, but not compared to the number of false positives we thought we might get without this filter.
#' }
#' @references 
#' Yao, Lijing, et al. "Inferring regulatory element landscapes and transcription 
#' factor networks from cancer methylomes." Genome biology 16.1 (2015): 1.
#' Method section (Linking enhancer probes with methylation changes to target genes with expression changes).
#' @param data A MultiAssayExperiment with a DNA methylation martrix or a DNA methylation matrix
#' @param K Cut off to consider probes as methylated or unmethylated. Default: 0.3
#' @param percentage The percentage of samples we should have at least considered as methylated and unmethylated
#' @return An object with the same class, but with the probes removed.
#' @importFrom MultiAssayExperiment experiments<-
#' @export
#' @examples
#'  random.probe <- runif(100, 0, 1)
#'  bias_l.probe <- runif(100, 0, 0.3)
#'  bias_g.probe <- runif(100, 0.3, 1)
#'  met <- rbind(random.probe,bias_l.probe,bias_g.probe)
#'  met <- preAssociationProbeFiltering(data = met,  K = 0.3, percentage = 0.05)
#'  met <- rbind(random.probe,random.probe,random.probe)
#'  met <- preAssociationProbeFiltering(met,  K = 0.3, percentage = 0.05)
#'  data <- ELMER:::getdata("elmer.data.example") # Get data from ELMER.data
#'  data <- preAssociationProbeFiltering(data,  K = 0.3, percentage = 0.05)
#'  
#'  cg24741609 <- runif(100, 0, 1)
#'  cg17468663 <- runif(100, 0, 0.3)
#'  cg14036402 <- runif(100, 0.3, 1)
#'  met <- rbind(cg24741609,cg14036402,cg17468663)
#'  colnames(met) <- paste("sample",1:100)
#'  exp <- met
#'  rownames(exp) <- c("ENSG00000141510","ENSG00000171862","ENSG00000171863")
#'  sample.info <- S4Vectors::DataFrame(primary = paste("sample",1:100),
#'                                      sample.type = rep(c("Normal", "Tumor"),50))
#'  rownames(sample.info) <- colnames(exp)
#'  mae <- createMAE(exp = exp, met = met, colData = sample.info, genome = "hg38") 
#'  mae <- preAssociationProbeFiltering(mae,  K = 0.3, percentage = 0.05)
preAssociationProbeFiltering <- function(data, K = 0.3, percentage = 0.05){
  print.header("Filtering probes", type = "section")
  message("For more information see function preAssociationProbeFiltering")
  
  if(is(data,"MultiAssayExperiment")) { 
    met <- assay(getMet(data))
  } else {
    met <- data
  }
  # In percentage how many probes are bigger than K ?
  Meth_B <- rowMeans(met > K, na.rm = TRUE)
  # We should  have at least 5% methylation value < K or at least 5% methylation value > K
  keep <- Meth_B < (1 - percentage) & Meth_B > percentage
  keep[is.na(keep)] <- FALSE
  message("Making sure we have at least ", percentage * 100, "% of beta values lesser than ", K," and ", 
          percentage * 100, "% of beta values greater ",K,".") 
  if(length(keep) - sum(keep) != 0) {
    message("Removing ", length(keep) - sum(keep), " probes out of ", length(keep))
  } else {
    message("There were no probes to be removed")
  }
  if(is(data,"MultiAssayExperiment")) {
    experiments(data)["DNA methylation"][[1]] <- experiments(data)["DNA methylation"][[1]][keep,]
  } else {
    data <- data[keep,,drop = FALSE]
  }
  return(data)
}

#' @importFrom xml2 read_html
#' @importFrom rvest html_table
getHocomocoTable <- function(){
  hocomoco <- tryCatch({
    hocomoco <- "https://hocomoco14.autosome.org/human/mono?full=true" %>% read_html()  %>%  html_table()
    hocomoco <- hocomoco[[1]]
  }, error = function(e) {
    getdata("hocomoco.table")
  })
  
  TF.family <-  createMotifRelevantTfs()
  TF.subfamily <-  createMotifRelevantTfs("subfamily")
  
  x <- do.call(rbind.data.frame,lapply(TF.family, function(x) paste(x,collapse = ";")))
  x$Model <- names(TF.family)
  colnames(x) <- c("TF.family.member","Model")
  hocomoco <- merge(hocomoco,x, by = "Model")
  x <- do.call(rbind.data.frame,lapply(TF.subfamily, function(x) paste(x,collapse = ";")))
  x$Model <- names(TF.subfamily)
  colnames(x) <-  c("TF.subfamily.member","Model")
  hocomoco <- merge(hocomoco,x, by = "Model")
  return(hocomoco)
}

#' @title Get random pairs
#' @description 
#' This function will receive a pair gene probes and will return a 
#' random object with the following pattern, if a probe is linked to R1 and L3 genes
#' the random pairs will be a random probes (a distal probe not in the input pairs) 
#' also linked to its R1 and L3 gene.
#' @param pairs A data frame with probe, gene and side information. See example below.
#' @param met.platform DNA methyaltion platform to retrieve data from: EPIC or 450K (default)
#' @param genome Which genome build will be used: hg38 (default) or hg19.
#' @param cores A interger which defines the number of cores to be used in parallel 
#' process. Default is 1: no parallel process.
#' @return A data frame with the random linkages
#' @export
#' @importFrom dplyr pull filter
#' @importFrom TCGAbiolinks colDataPrepare
#' @examples
#' \dontrun{
#'  data <- ELMER:::getdata("elmer.data.example")
#'  nearGenes <- GetNearGenes(TRange=getMet(data)[c("cg00329272","cg10097755"),],
#'                             geneAnnot=getExp(data))
#'                             
#'  pair <- get.pair(data = data,
#'                   group.col = "definition", 
#'                   group1 = "Primary solid Tumor", 
#'                   group2 = "Solid Tissue Normal",
#'                   mode = "supervised",
#'                   diff.dir = "hypo",
#'                   nearGenes = nearGenes,
#'                   permu.size = 5,
#'                   raw.pvalue =  0.001,
#'                   Pe = 0.2,
#'                   dir.out="./",
#'                   permu.dir = "permu_test",
#'                   label = "hypo")
#' }
#'  pair <- data.frame(Probe = rep("cg00329272",3), 
#'                     GeneID = c("ENSG00000116213","ENSG00000130762","ENSG00000149527"),
#'                     Sides = c("R5","R2","L4"))                    
#'  getRandomPairs(pair)
getRandomPairs <- function(pairs, 
                           genome = "hg38",
                           met.platform = "450K",
                           cores = 1) {
  
  if(missing(pairs)) stop("Please set pairs argument")
  if(is.data.frame(pairs)) pairs <- as.data.frame(pairs)
  
  # Rename column
  if("Side" %in% colnames(pairs)) colnames(pairs)[grep("Side",colnames(pairs))] <- "Sides"
  if(!"Sides" %in% colnames(pairs)) stop("No column Sides in the object")
  
  if("Target" %in% colnames(pairs)) colnames(pairs)[grep("Target",colnames(pairs))] <- "Probe"
  if("ID" %in% colnames(pairs)) colnames(pairs)[colnames(pairs) == "ID"] <- "Probe"
  
  parallel <- FALSE
  if (cores > 1){
    if (cores > detectCores()) cores <- detectCores()
    registerDoParallel(cores)
    parallel = TRUE
  }
  
  # Get Probe information
  met.info <- getInfiniumAnnotation(plat = met.platform, genome = genome) 
  probes.ranges <- as.data.frame(met.info)[,1:5]
  colnames(probes.ranges) <- paste0("probe_",colnames(probes.ranges))
  
  # Get distal probes not in the pairs
  distal.probe <- get.feature.probe(
    genome = genome,
    met.platform = met.platform,
    feature = NULL
  ) # get distal probes
  distal.probe <- distal.probe[!names(distal.probe) %in% pairs$Probe,] # Select probes were not used
  
  nb.pairs <- nrow(pairs)
  nb.probes <- length(unique(pairs$Probe))
  # get gene information
  genes <- getTSS(genome = genome)
  genes$ensembl_transcript_id <- NULL
  genes <- genes[!duplicated(genes$ensembl_gene_id)]
  
  df.random <- NULL
  near.genes.linked <- NULL
  near.genes.df <- NULL
  # We will get the double of random probes, because some will not be used in case it does not matches the position
  # Example: real probe + gene R10 and random probe does not have R10. Discart and get next random
  not.matched <- 1:nb.probes
  numFlankingGenes <- 20
  while(length(not.matched) > 0){
    random.probes <- distal.probe[sample(1:length(distal.probe), length(not.matched)),]
    near.genes <- GetNearGenes(TRange = random.probes, 
                               geneAnnot = genes, 
                               numFlankingGenes = numFlankingGenes)
    
    near.genes.linked <- plyr::alply(1:length(not.matched),
                                     .margins = 1,
                                     .fun = function(x){
                                       side <- pairs %>% 
                                         filter(pairs$Probe == unique(pairs$Probe)[not.matched[x]]) %>% pull('Sides')
                                       ret <- near.genes[near.genes$ID == unique(near.genes$ID)[x] & near.genes$Side %in% side,]
                                       if(!all(side %in% ret$Side)) return(NULL)
                                       ret
                                     },.progress = "time", .parallel = parallel)
    
    not.matched <- not.matched[grep("NULL",near.genes.linked)]
    if (length(not.matched) > 0) {
      aux <- pairs %>% filter(pairs$Probe == unique(pairs$Probe)[not.matched[1]]) %>% pull('Sides') 
      numFlankingGenes <- max(as.numeric(gsub("L|R","",aux))) * 2
    }
    near.genes.df <- rbind(near.genes.df,data.table::rbindlist(near.genes.linked))
    distal.probe <- distal.probe[!names(distal.probe) %in% near.genes.df$Target,] # Remove probes already used
  }
  colnames(near.genes.df)[1] <- "Probe" 
  
  # Add probe metadata to output
  probes.ranges$Probe <- rownames(probes.ranges)
  near.genes.df <- merge(near.genes.df, probes.ranges, by = "Probe")
  return(near.genes.df)
}


# Reading homer output. For each reagion (rows)
# homer will try to find if a given motif was found in it.
# This will read this homer file and create a sparce matrix
# in which 1 means the motif was found and 0 not found.
# this is used to calculate the motif enrichement compared 
# to a background signal.
getMatrix <- function(filename) {
  motifs <- readr::read_tsv(filename)
  # From 1 to 21 we have annotations
  matrix <- Matrix::Matrix(0, nrow = nrow(motifs), ncol = ncol(motifs) - 21 ,sparse = TRUE)
  colnames(matrix) <- gsub(" Distance From Peak\\(sequence,strand,conservation\\)","",colnames(motifs)[-c(1:21)])
  rownames(matrix) <- motifs$PeakID
  matrix[!is.na(motifs[,-c(1:21)])] <- 1
  matrix <- as(matrix, "nsparseMatrix")
  return(matrix)
}

#' @title Calculate motif Erichment 
#' @description Calculates fisher exact test
#' @param foreground A nsparseMatrix object in each 1 means the motif is found in a region, 0 not.
#' @param background A nsparseMatrix object in each 1 means the motif is found in a region, 0 not.
#' @export
#' @examples 
#' foreground <- Matrix::Matrix(sample(0:1,size = 100,replace = TRUE), 
#'                              nrow = 10, ncol = 10,sparse = TRUE)
#' rownames(foreground) <- paste0("region",1:10)
#' colnames(foreground) <- paste0("motif",1:10)
#' background <- Matrix::Matrix(sample(0:1,size = 100,replace = TRUE), 
#'                              nrow = 10, ncol = 10,sparse = TRUE)
#' rownames(background) <- paste0("region",1:10)
#' colnames(background) <- paste0("motif",1:10)
#' calculateEnrichement(foreground,background)
calculateEnrichement <- function(foreground, 
                                 background){
  if(missing(foreground)) stop("foreground argument is missing")
  if(missing(background)) stop("background argument is missing")
  
  # a is the number of probes within the selected probe set that contain one or more motif occurrences; 
  # b is the number of probes within the selected probe set that do not contain a motif occurrence; 
  # c and d are the same counts within the entire enhancer probe set (background)
  # lower boundary of 95% conf idence interval = exp (ln OR - SD)
  a <- Matrix::colSums(foreground)
  b <- nrow(foreground) - Matrix::colSums(foreground)
  c <- Matrix::colSums(background)
  d <- nrow(background) - Matrix::colSums(background)
  fisher <- plyr::adply(seq_len(length(a)),.margins = 1, .fun = function(i)  { 
    x <- fisher.test(matrix(c(a[i],b[i],c[i],d[i]),nrow = 2,ncol = 2))
    ret <- data.frame(x$conf.int[1],x$conf.int[2],x$estimate,x$p.value)
    colnames(ret) <- c("lowerOR","upperOR","OR","p.value")
    ret
  },.id = NULL,.progress = "text")
  rownames(fisher) <- names(a)
  Summary <- data.frame(motif  =  names(a),
                        NumOfRegions = Matrix::colSums(foreground, na.rm=TRUE),
                        fisher,
                        FDR = p.adjust(fisher$p.value,method = "BH"),
                        stringsAsFactors = FALSE)
  Summary <- Summary[order(-Summary$lowerOR),]
  return(Summary)
}


#' @title Use Hocomoco motif and homer to identify motifs in a given region
#' @param regions A GRanges object. Names will be used as the identifier.
#' @param output.filename Final file name
#' @param region.size If NULL the motif will be mapped to the region. If set a window around its center will be considered.
#' For example if region.size is 500, then +-250bp round it will be searched.
#' @param cores A interger which defines the number of cores to be used in parallel 
#' process. Default is 1: no parallel process.
#' @param genome Homer genome (hg38, hg19)
#' @param nstep Number of regions to evaluate in homer, the bigger, more memory it will use at each step.
#' @description 
#' To find for each probe the know motif we will use HOMER software (http://homer.salk.edu/homer/).
#' Homer and genome should be installed before this function is executed
#' Step:
#' 1 - get DNA methylation probes annotation with the regions
#' 2 - Make a bed file from it
#' 3 - Execute section: Finding Instance of Specific Motifs 
#' from http://homer.salk.edu/homer/ngs/peakMotifs.html to the HOCOMOCO TF motifs
#' Also, As HOMER is using more RAM than the available we will split the files in to 100k probes.
#' Obs: for each probe we create a winddow of 500 bp (-size 500) around it. 
#' This might lead to false positives, but will not have false negatives.
#' The false posives will be removed latter with some statistical tests.
#' @importFrom utils write.table
#' @examples 
#' \dontrun{
#'  # use the center of the region and +-250bp around it
#'  gr0 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), 
#'                     c(1, 3, 2, 4)
#'                     ), 
#'                IRanges(1:10, width=10:1)
#'                )
#'  names(gr0) <- paste0("ID",c(1:10))
#'  findMotifRegion(regions = gr0, region.size = 500, genome = "hg38", cores = 1)
#'  
#'  # use the region size itself
#'  gr1 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), 
#'                 IRanges(1:10, width=sample(200:1000,10)))
#'  names(gr1) <- paste0("ID",c(1:10))
#'  findMotifRegion(regions = gr0, genome = "hg38", cores = 1)
#' }
findMotifRegion <- function(regions, 
                            output.filename = "mapped_motifs_regions.txt",
                            region.size = NULL,
                            genome = "hg38",
                            nstep = 10000,
                            cores = 1){
  
  if(!is(regions, "GRanges")) stop("Regions must be a Genomic Ranges object")
  # get all hocomoco 11 motifs
  TFBS.motif <- "http://hocomoco11.autosome.org/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_homer_format_0.0001.motif"
  if(!file.exists(basename(TFBS.motif))) downloader::download(TFBS.motif,basename(TFBS.motif))
  
  
  if(is.null(names(regions))){
    names(regions) <- tibble::as_tibble(regions) %>% tidyr::unite(col = "names","seqnames","start","end") %>% pull(names)
  }
  df <- data.frame(seqnames = seqnames(regions),
                   starts = as.integer(start(regions)),
                   ends = end(regions),
                   names = names(regions),
                   scores = c(rep(".", length(regions))),
                   strands = strand(regions))
  n <- nrow(df)
  step <- ifelse(n > nstep, nstep, n)
  
  if(!file.exists(output.filename)){
    pb <- txtProgressBar(max = floor(n/step), style = 3)
    
    for(j in 0:floor(n/step)){
      setTxtProgressBar(pb, j)
      # STEP 2
      file.aux <- paste0(gsub(".txt","",output.filename),"_",j,".bed")
      if(!file.exists(gsub(".bed",".txt",file.aux))){
        end <- ifelse(((j + 1) * step) > n, n,((j + 1) * step))
        write.table(df[((j * step) + 1):end,], file = file.aux, col.names = F, quote = F,row.names = F,sep = "\t")
        
        # STEP 3 use -mscore to get scores
        # we need to check if annotatePeaks.pl is available!
        cmd <- paste("annotatePeaks.pl" ,file.aux, genome, "-m", basename(TFBS.motif), 
                     ifelse(is.null(region.size),"",paste0("-size ", region.size)),
                     "-cpu", cores,
                     ">", gsub(".bed",".txt",file.aux))
        error <- system(cmd)
        if(error == 127) {
          unlink(file.aux)
          unlink(gsub(".bed",".txt",file.aux))
          stop(paste0("annotatePeaks.pl had an error. Please check homer install.",
                      "\nIf already installed be sure R can see it.",
                      "You can change PATH evn variable with ",
                      "\nSys.setenv(PATH=paste(Sys.getenv('PATH'), '/the/bin/folder/of/bedtools', sep=':'))")
          )
        }
      }
    }
    close(pb)
  }
  # We will merge the results from each file into one
  peaks <- NULL
  pb <- txtProgressBar(max = floor(n/step), style = 3)
  for(j in 0:floor(n/step)){
    setTxtProgressBar(pb, j)
    f <- paste0(gsub(".txt","",output.filename),"_",j,".txt")
    if(file.exists(f)){
      aux <-  readr::read_tsv(f)
      colnames(aux)[1] <- "PeakID"
      if(is.null(peaks)) {
        peaks <- aux
      } else {
        peaks <- rbind(peaks, aux)
      }
    }
    gc()
  }
  close(pb)
  print(paste0("Writing file: ",output.filename))
  if(!is.null(peaks)) {
    readr::write_tsv(peaks, path = output.filename,col_names = TRUE)
  }
  print("DONE!")
}


#' @title Get MR TF binding regions infered by ELMER 
#' @description Saves a bed file with the unmethylated probes (+-250bp) regions that was infered
#' to be bound by a given TF
#' @param tf TF name
#' @param results.dir path to the directory with the results 
#' (i.e. analysis/unsupervised/definition-Primary.solid.Tumor_vs_Solid.Tissue.Normal/hypo/)
#' @param genome Human genome (hg38, hg19)
#' @param met.platform DNA Methylation  Array platform (EPIC, 450K)
#' @importFrom readr read_csv
#' @importFrom dplyr %>% mutate
#' @importFrom tidyr unnest 
#' @importFrom GenomicRanges resize
#' @examples 
#' \dontrun{
#'   getTFBindingSites("HNF1A",
#'                     results.dir = "analysis/unsupervised/group-Tumor_vs_Normal/hypo/")
#' }
#' @export
getTFBindingSites <- function(tf = NULL, 
                              results.dir = NULL,
                              genome = "hg38",
                              met.platform =  "450K"){
  
  if(is.null(tf)) stop("Please set a tf to be searched")
  
  tf.file <- dir(path = results.dir,
                 pattern = "*.significant.TFs.with.motif.summary.csv",
                 recursive = T,
                 full.names = T,
                 all.files = T)
  if(length(tf.file) == 0) stop("No TF results file found")
  tf.tab <- readr::read_csv(tf.file,col_types = readr::cols()) %>% na.omit
  
  pair.file <- dir(path = results.dir,
                   pattern = "*.pairs.significant.withmotif.csv",
                   recursive = T,
                   full.names = T,
                   all.files = T)
  if(length(pair.file) == 0) stop("No pair results file found")
  pair.tab <- readr::read_csv(pair.file,col_types = readr::cols())
  
  
  # for each enriched motif find the one with the TF in the classification (family of subfamily within the 5%)
  for(classification in c("family","subfamily")){
    if(classification == "family"){
      tf.tab <- tf.tab %>% 
        mutate(tf.target = strsplit(as.character(tf.tab$potential.TF.family), ";")) %>% 
        unnest(cols = "tf.target") 
    } else {
      tf.tab <- tf.tab %>% 
        mutate(tf.target = strsplit(as.character(tf.tab$potential.TF.subfamily), ";")) %>% 
        unnest(cols = "tf.target") 
    }
    
    if( !tf %in% tf.tab$tf.target) stop("TF not found")
    motif <- tf.tab %>% filter(tf.tab$tf.target == tf) %>% pull('motif')
    
    # For each enriched motif find in each paired probes it appers
    pair.tab <- pair.tab %>% 
      mutate(motif.target = strsplit(as.character(pair.tab$enriched_motifs), ";")) %>% 
      unnest(cols = "motif.target") 
    probes <- pair.tab %>% filter(pair.tab$motif.target %in% motif) %>% pull('Probe') %>% unique
    
    # for each probe get region and write bed file
    metadata <- getInfiniumAnnotation(plat = met.platform, genome = genome)
    metadata <- metadata[probes,]
    metadata <- resize(metadata,  width = 500, fix = 'center')
    file.out <- file.path(results.dir, paste0(tf, "_",classification,".bed"))
    message("Saving as ", file.out)
    rtracklayer::export.bed(metadata,con = file.out)
  }
  
}
