library("GetoptLong")
library("ggplot2")
library("cowplot")
library("dplyr")
library("readr")

GetoptLong(
  "infile=s", "A TSV file of translated variants to plot.",
  "outbase=s", "A basename to use for outputs."
)

#' Base colors for making the single_gene_plots
#' There is one color for each of A, C, T, G, DEL, and INS
#base_colors <- c()
#label_base_colors <- function(x){
#  
#}

#' Amino acid colors for making the single_gene_plots
#' There is one color per amino acid.
#aa_colors <- c()
#label_aa_colors <- function(x){
#  
#}

#' A plot of DNA base / Amino Acid
#' over the length of a gene/genome/feature.
#' 
#' @param x a dataframe containing the variables (pos, base)
#' @param xlimits a vector (length two) of start/end positions to restrict the plot the plotting area to.
#' @param displayAsProportions display a proportional ("filled") bar plot rather than one of counts
#' @return A ggplot2 object
single_gene_plot <- function(x,
                             xlimits=NULL,
                             displayAsProportion=FALSE,
                             displayUniqueBases =FALSE,
                             hideInvariableSites =FALSE){
  require(ggplot2)
  require(cowplot)
  require(dplyr)
  require(scales)
  
  if (displayUniqueBases){
    x <- x %>% distinct(pos, base)
  }
  if (hideInvariableSites){
    keep_sites <- x %>% distinct(pos, base) %>% group_by(pos) %>% mutate(n = n()) %>% filter(n > 1)
    x <- x %>% filter(pos %in% keep_sites$pos)
  }
  geneplot <- ggplot(x)
  if (displayAsProportion){
    geneplot <- geneplot +
      geom_bar(aes(x = pos, fill = base), position = "fill", width = 1)
  }

  else{
    geneplot <- geneplot + geom_bar(aes(x = pos, fill = base), stat="count", width = 1)
  }
    
  geneplot <- geneplot +
    theme_half_open(14) + 
    coord_cartesian(expand=FALSE) +
    xlab("Position") +
    ylab("Count")
  
  
  if (!is.null(xlimits)){
      geneplot <- geneplot +
        scale_x_continuous(breaks=unique(x$pos), limits=c(xlimits[1] -1, xlimits[2]+1))
  }
  else if (hideInvariableSites){
    geneplot <- geneplot + scale_x_continuous(breaks=pretty_breaks())
  }
  else{
    geneplot <- geneplot + scale_x_continuous(breaks=unique(x$pos))
  }
  
  if (is.null(xlimits)){
    geneplot <- geneplot +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  return (geneplot)
} 

#' Plots variants over the positions of every feature,
#' with each unique feature receiving its own facet.
#' 
all_features_plot <- function(x,
                  xlimits=NULL,
                  displayAsProportion=FALSE,
                  displayUniqueBases =FALSE,
                  hideInvariableSites =FALSE){
  geneplot <- single_gene_plot(x,xlimits,displayAsProportion,displayUniqueBases,hideInvariableSites)
  feat_plot <- geneplot + facet_wrap(facets = vars(feature))
  return (feat_plot)
}

#' Plots a simple bar plot of the number of variants within each
#' feature.
per_feature_count_plot <- function(x){
  count_plot <- ggplot(x) +
    geom_bar(aes(x = feature)) +
    theme_half_open(14) +
    guides(fill=F)
  return (count_plot)
}


# test_df <- tibble::tibble(base=c("A", "A", "C", "A", "T","I", "G", "C", "A", "A", "A"), pos=c(1,1,1,1,2,1,3,4,4,4,4))
# single_gene_plot(test_df)
# single_gene_plot(test_df, displayAsProportion = TRUE)
# single_gene_plot(test_df, displayAsProportion = TRUE, hideInvariableSites = TRUE)
# single_gene_plot(test_df, displayUniqueBases = TRUE)
# single_gene_plot(test_df, displayUniqueBases = TRUE, hideInvariableSites = TRUE)
# single_gene_plot(test_df, hideInvariableSites = TRUE)
# single_gene_plot(test_df, xlimits = c(1,3))


fi <- read_tsv(infile)
sgp <- single_gene_plot(fi)
save_plot(paste(outbase, ".genome_variants_plot.pdf", sep=""), sgp, nrow = 1, base_asp = 3)
afp <- all_features_plot(fi)
save_plot(paste(outbase, ".feature_variants_plot.pdf", sep=""), afp, nrow = 1, base_asp = 3)
pfcp <- per_feature_count_plot(fi)
save_plot(paste(outbase, ".per_feature_count_plot.pdf", sep=""), pfcp, nrow = 1, base_asp = 3)





