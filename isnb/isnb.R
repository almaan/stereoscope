#!/usr/bin/Rscript

#' Code used to generate Supplementary Figure 20.
#' Provide as input: (1) a .tsv file with the expression data, rownames
#' should be [x_coordinate]x[y_coordinate]; and (2) a list of genes to
#' be visualized. Output directory can be specified, default is wd.


library(argparse)
library(ggplot2)
library(gridExtra)

TEXT_SIZE = 25

gene.plot <- function(X,
                      col.name,
                      size = 5
                      ) {

  x <- as.numeric(sapply(rownames(X), function(x){strsplit(x,"x")[[1]][1]}))
  y <- as.numeric(sapply(rownames(X), function(x){strsplit(x,"x")[[1]][2]}))
  v <- X[,col.name]
  plotter <- data.frame(x = y, y = -x, v = v)
  g <- ggplot(plotter,aes(x,y)) +
    theme(text = element_text(size = TEXT_SIZE),
          panel.background = element_rect(fill ="white",color ="white"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()

          ) +

    labs(title = col.name) +
    geom_point(aes(fill = v),color = "gray",size = size,shape = 21) +
    scale_fill_gradient(low="white", high="blue",guide = F)+
    coord_fixed()

  return(g)
}


parser <- ArgumentParser()
parser$add_argument("-c",
                    "--count_data",
                    type = "character")
parser$add_argument("-g",
                    "--genes",
                    type = "character")
parser$add_argument("-o",
                    "--out_dir",
                    default = NULL,
                    type = "character")

args <- parser$parse_args()

if is.null(args$out_dir) {
  args$out_dir <- getwd()
}

data <- read.table(args$count_data,
                   sep = '\t',
                   header = T,
                   row.names = 1)

gene.names <- as.character(unlist(read.delim(args$genes,
                                             header = F,
                                             sep = "\n")))

elist <- list()
for (gene in gene.names) {
  print(gene)
  elist[[gene]] <- gene.plot(data,
                             gene,
                             size = 1.8)
}

n.cols <- 3
n.rows <- ceiling(length(genes) / n.cols)


width <- 400*n.cols
height <- 400*n.rows

png(file.path(args$out_dir,
              "expression-visualization.png"),
    width = width,
    height = height,
    units = "px")

E <- grid.arrange(grobs = elist,
                  ncol = n.cols)
print(E)
dev.off()
