if(!require("devtools")){
  install.packages("devtools", repos="http://cran.us.r-project.org")
  library("devtools")
}

if(!require("dplyr")){
  install.packages("dplyr", repos="http://cran.us.r-project.org")
  library("dplyr")
}

if(!require("BioThingsClient.R")){
  install_github("biothings/BioThingsClient.R")
  library("BioThingsClient.R")
}

VariantAnnotator = function(variant_file = "FHS_EA_MRS_5e8_snplist.txt",
                            chr = "CHR",
                            pos = "POS",
                            ref = "NEA",
                            alt = "EA",
                            delimiter = " "){
  
  cat("Loading Variant File...\n")
  data = read.delim(variant_file, sep = delimiter)
  
  names(data)[which(names(data)==chr)] = "chr"
  names(data)[which(names(data)==pos)] = "pos"
  names(data)[which(names(data)==ref)] = "ref"
  names(data)[which(names(data)==alt)] = "alt"
  
  data = data %>%
    mutate(chr = ifelse(chr == "X", "23", chr)) %>%
    mutate(chr = as.numeric(chr)) %>%
    mutate(pos = as.numeric(pos))
  
  data = data %>%
    mutate(query = paste0("chr", chr, ":g.", pos, ref, ">", alt))
  
  cat("Generating Variant Annotation...\n")
  annotation = vector()
  for(q in data$query){
    try({
      v = btGet(variant_client, q)[[1]]
      
      gene_symbol = NA
      gene_name = NA
      variant_type = NA
      variant_consequence = NA
      variant_CADD_phred = NA
      variant_phylop_vertebrate = NA
      variant_polyphen = NA
      variant_sift = NA
      variant_rs = NA
      
      gene_symbol = paste(v[["dbsnp"]][["gene"]][["symbol"]], collapse = "; ")
      gene_name = gsub(pattern = ",", ";", v[["dbsnp"]][["gene"]][["name"]])
      variant_type = gsub(",","; ", paste(v[["cadd"]][["consdetail"]],
                                          collapse = "; "))
      variant_consequence = gsub(",","; ", paste(v[["cadd"]][["consequence"]],
                                                 collapse = "; "))
      variant_CADD_phred = v[["cadd"]][["phred"]]
      variant_phylop_vertebrate = v[["cadd"]][["phylop"]][["vertebrate"]]
      variant_polyphen = v[["cadd"]][["polyphen"]][["val"]]
      variant_sift = v[["cadd"]][["sift"]][["val"]]
      variant_rs = v[["dbsnp"]][["rsid"]]
      
      gene_symbol = ifelse(is.null(gene_symbol), NA, gene_symbol)
      gene_name = ifelse(is.null(gene_name), NA, gene_name)
      variant_type = ifelse(is.null(variant_type), NA, variant_type)
      variant_consequence = ifelse(is.null(variant_consequence), NA, variant_consequence)
      
      variant_CADD_phred = ifelse(is.null(variant_CADD_phred), NA, variant_CADD_phred)
      variant_phylop_vertebrate = ifelse(is.null(variant_phylop_vertebrate),
                                         NA, variant_phylop_vertebrate)
      variant_polyphen = ifelse(is.null(variant_polyphen), NA, variant_polyphen)
      variant_sift = ifelse(is.null(variant_sift), NA, variant_sift)
      variant_rs = ifelse(is.null(variant_rs), NA, variant_rs)
      
      anno = tibble(query = q, variant_rs, gene_symbol, gene_name, variant_type, 
                    variant_consequence, variant_CADD_phred, 
                    variant_phylop_vertebrate, variant_polyphen,
                    variant_sift)
      
      annotation = rbind.data.frame(annotation, anno)
      cat(paste0("Annotating ", q," ..\n"))
    }, silent = TRUE)
  }
  
  data_annotated = data %>%
    left_join(annotation, by = "query")
  
  plot_name = gsub(pattern = ".txt", 
                   replacement = "_Annotated.csv",
                   variant_file)
  write.csv(data_annotated, file = plot_name,
            quote = FALSE, row.names = FALSE)
  
}

### NOT RUN
# VariantAnnotator()