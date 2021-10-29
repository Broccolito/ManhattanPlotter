if(!require("dplyr")){
  install.packages("dplyr", repos="http://cran.us.r-project.org")
  library("dplyr")
}

if(!require("ggplot2")){
  install.packages("ggplot2", repos="http://cran.us.r-project.org")
  library("ggplot2")
}

if(!require("ggpubr")){
  install.packages("ggpubr", repos="http://cran.us.r-project.org")
  library("ggpubr")
}

if(!require("latex2exp")){
  install.packages("latex2exp", repos="http://cran.us.r-project.org")
  library("latex2exp")
}

ManhattanPlotter = function(gwas_result_file = "FHS_EA_MRS_merged.txt",
                            chr = "CHR",
                            pos = "POS",
                            pvalue = "p.value",
                            delimiter = " ",
                            color1 = "cornflowerblue",
                            color2 = "blue4",
                            threshold_linecolor = "black"){
  
  cat("Loading GWAS Statistics...\n")
  data = read.delim(gwas_result_file, sep = delimiter)
  cat("GWAS Statistics Loaded...\n")
  
  names(data)[which(names(data)==chr)] = "chr"
  names(data)[which(names(data)==pos)] = "pos"
  names(data)[which(names(data)==pvalue)] = "pvalue"
  
  data = data %>%
    select(chr, pos, pvalue) %>%
    mutate(chr = ifelse(chr == "X", "23", chr)) %>%
    mutate(chr = as.numeric(chr)) %>%
    mutate(pos = as.numeric(pos)) %>%
    mutate(pvalue = as.numeric(pvalue))
  
  data_chromosome_position = data %>%
    group_by(chr) %>%
    summarise(chromosome_length = max(pos)) %>%
    mutate(tot = cumsum(as.numeric(chromosome_length)) - chromosome_length) %>%
    select(-chromosome_length)
  
  data_cleaned = left_join(data, data_chromosome_position, by = "chr") %>%
    arrange(chr, pos) %>%
    mutate(bp_cum = pos + tot) %>%
    mutate(nlogp = -log(pvalue, base = 10)) %>%
    select(chr, pos, pvalue, nlogp, tot, bp_cum)
  
  axisdf = data_cleaned %>% 
    group_by(chr) %>% summarize(center = (max(bp_cum) + min(bp_cum))/2)
  
  plt = ggplot(data_cleaned, aes(x = bp_cum, y = nlogp)) +
    geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c(color1, color2), 22)) +
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
    scale_y_continuous(limits = c(0, max(max(data_cleaned$nlogp)+1, 
                                         -log(5e-8,base = 10)+1)), 
                       expand = c(0, 0)) +
    geom_hline(yintercept = -log(5e-8,base = 10), 
               linetype = "solid", color = threshold_linecolor) +
    geom_hline(yintercept = -log(1e-5,base = 10), 
               linetype = "dashed", color = threshold_linecolor) +
    xlab("Chromosome") + 
    ylab(TeX("-log_{10}(p-value)")) +
    theme_pubr() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20)
    )
  
  cat("Saving the Manhattan Plot as a PNG file...\n")
  plot_name = gsub(pattern = ".txt", replacement = "_ManhattanPlot.png", gwas_result_file)
  ggsave(filename = plot_name, plot = plt, device = "png", dpi = 1200,
         height = 6, width = 15)
  cat(paste0("Plot saved as ", plot_name, "\n"))
  
}

### NOT RUN
# ManhattanPlotter()
