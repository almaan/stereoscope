#!/usr/bin/Rscript

#' ----Goodness of Fit comparison across distributions----
#' Code to generate the results presented in Supplementary Figure 21 - 37.
#' Provide (either via the CLI or with hardcoded paths) the count data and genes
#' that are of interest to assess. The code will then :
#' 1. Normalize and Cluster the data using the Seurat suite
#' 2. Fit each of the three distributions : Negative Binomial, Normal and Poisson, to
#'    the expression of each gene within each cluster.
#' 3. Compute BIC for each fit
#' 4. Generate the images displayed in Supplementary Figure 21-37


sh <- suppressPackageStartupMessages
sh(library(Seurat))
sh(library(sctransform))
sh(library(ggplot2))
sh(library(dplyr))
sh(library(gridExtra))
sh(library(grid))
sh(library(argparse))
sh(library(fitdistrplus))
sh(library(wesanderson))
sh(library(ggsci))


TEXT_SIZE <- 20

make.bars <- function(g,crit,n){

  g <- g  + geom_bar(stat = "identity",
                     position ="dodge",
                     color = "black")+

  scale_fill_manual(values = wes_palette(n = n,
                                         name = "Cavalcanti1")) +

  theme(axis.text.x = element_text(color = "black",size =15,angle = 45),
        legend.text = element_text(size = 20),
        axis.title.x = element_text(color = "black",size = 20),
        plot.title = element_text(size = TEXT_SIZE),
        axis.title.y = element_text(color = "black",size = 20))+
  xlab("Gene") +
  ylab(crit)
  return(g)
}


VisualizeGof <- function(gof.stat,cluster){

  g <- list()
  n.dists <- length(unique(gof.stat$Distribution))
  for (metric in unique(gof.stat$Metric)) {

    g.s <- ggplot(data = gof.stat[gof.stat$Metric == metric,],
                  aes(x = Gene,
                      y = Estimate,
                      fill = Distribution))

    g.s <- make.bars(g.s,crit = metric,n = n.dists) +
      labs(title=(paste0("Cluster : ",cluster)))

    g[[metric]] <- g.s
  }
  return(g)
}

FitGene <- function(gene.expr,gene.name){

  fitD.nb <- fitdist(gene.expr,
                     distr ="nbinom")
  fitD.poi <- fitdist(gene.expr,
                      distr ="pois")

  fitD.norm <- fitdist(gene.expr,
                       distr ="norm")

  vals <- seq(0,max(gene.expr + 1))

  ys.nb <- dnbinom(vals,
                   mu = fitD.nb$estimate["mu"],
                   size = fitD.nb$estimate["size"])

  ys.poi <-  dpois(vals,
                   lambda = fitD.poi$estimate["lambda"]
                   )

  ys.norm <- dnorm(vals,
                   mean = fitD.norm$estimate["mean"],
                   sd = fitD.norm$estimate["sd"]
                   )

  fitted <- data.frame(xs = vals,
                       NegativeBinomial = ys.nb,
                       Poisson = ys.poi,
                       Normal = ys.norm,
                       NB.id = rep("NB",length(vals)),
                       Poi.id = rep("Pois",length(vals)),
                       Norm.id = rep("Norm",length(vals))
                       )
  
  comp.scores <- data.frame(Distribution = rep(c("Negative Binomial",
                                                 "Poisson",
                                                 "Normal"),
                                               2),
                            Estimate = c(fitD.nb$bic,
                                         fitD.poi$bic,
                                         fitD.norm$bic,
                                         fitD.nb$aic,
                                         fitD.poi$aic,
                                         fitD.norm$aic
                                         ),
                            Metric = c(rep("BIC",3),rep("AIC",3)),
                            Gene = rep(gene.name,6)
                       )
  return(list(fitted = fitted, comp.scores = comp.scores))

}


VisualizeFitDist <- function(plot.data,fitted,gene.name) {
  g <- ggplot(plot.data,aes_string(x = "vals")) +

    theme_minimal(base_size = 25) +

    geom_density(fill = "gray",alpha =0.7) +

    geom_line(data = fitted,
              inherit.aes = F,
              linetype = "dashed",
              size = 1.5,
              mapping = aes(x=xs,y=NegativeBinomial,color = NB.id)) +

    geom_line(data = fitted,
              inherit.aes = F,
              linetype = "dashed",
              size = 1.5,
              mapping = aes(x=xs,y=Normal,color = Norm.id)) +

    geom_line(data = fitted,
              inherit.aes = F,
              linetype = "dashed",
              size = 1.5,
              mapping = aes(x=xs,y=Poisson,color = Poi.id)) +

  scale_color_manual(
    values = c(NB = "red", Norm = "blue", Pois = "green"),
    labels = c("NegBin.","Norm.","Pois."),
    name = "Dist.",
    ) +

  labs(title = gene.name) +

  xlab("Counts") + 
  ylab("Density") +

  theme(text = element_text(size = TEXT_SIZE),
        axis.ticks = element_blank(),
        ## axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
        )

  return(g)

}


SeuratAnalysis <- function(data) {

  se <- CreateSeuratObject(data,min.cells = 1, min.features = 1 )
  se <- PercentageFeatureSet(se, pattern = "^MT-", col.name = "percent.mt")
  se <- SCTransform(se, vars.to.regress = "percent.mt")

  se <- RunPCA(se, verbose = FALSE)
  se <- RunUMAP(se,
                dims = 1:30,
                verbose = FALSE)

  se <- FindNeighbors(se,
                      dims = 1:30,
                      verbose = FALSE)

  se <- FindClusters(se,
                     verbose = FALSE)

  xcrd <- as.numeric(sapply(colnames(data),function(x){gsub("X"," ",strsplit(x,"x")[[1]][1])}))
  ycrd <- as.numeric(sapply(colnames(data),function(x){gsub("X"," ",strsplit(x,"x")[[1]][2])}))

  se@meta.data$X <- xcrd
  se@meta.data$Y <- ycrd

  res <- data.frame(x = xcrd,
                    y = ycrd ,
                    cluster.id = se@meta.data$SCT_snn_res.0.8,
                    cluster.name = paste("cluster",
                                    se@meta.data$SCT_snn_res.0.8))

  rownames(res) <- colnames(data)

  return(list(cluster.res =  res))
}


VisualizeClusters <- function(plot.data,
                              size = 2,
                              fontsize = 10) {

  clrs <- c(pal_jco()(10),pal_aaas()(10))
  clrs <- clrs[1:length(unique(plot.data$cluster.name))]
  print(clrs)

  g <- ggplot(plot.data, aes(x=y,y=-x)) +
    geom_point(aes(color = cluster.name),
               pch = 19,
               size = size) +
    coord_fixed() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          legend.text = element_text(size = fontsize)) +
  scale_color_manual(values = clrs)

  return(g)
}


VisualizeExpression <- function(X,
                                row.name,
                                size = 5
                                ) {

  x <- as.numeric(sapply(colnames(X),
                         function(x){strsplit(x,"x")[[1]][1]}))
  y <- as.numeric(sapply(colnames(X),
                         function(x){strsplit(x,"x")[[1]][2]}))

  v <- as.numeric(X[row.name,])

  plotter <- data.frame(x = y, y = -x, v = v)
  print(head(plotter))

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

    labs(title = row.name) +
    geom_point(aes(fill = v),color = "gray",size = size,shape = 21) +
  scale_fill_gradient(low="white",
                      high="blue",
                      guide = F)+
    coord_fixed()

  return(g)
}

main <- function(data.pth,
                 gene.list.pth,
                 odir
                 ){

  # General vars
  n.cols <- 3
  w.size <- 400
  h.size <- 400

  # Read Expression data
  data <-  read.table(data.pth,
                      header = 1,
                      row.names = 1)

  data <- as.data.frame(t(data))

  ## Read gene names to plot
  gene.names <- as.character(unlist(read.delim(gene.list.pth,
                                               header = F,
                                               sep = "\n")))

  print("Will work with:")
  print(gene.names)

  n.genes <- length(gene.names)
  n.rows <- ceiling(n.genes / n.cols)

  ## Plot genes

  grob.expr <- list()
  for (gn in gene.names){
    print(paste0("Plot : ", gn))
    g.expr <- VisualizeExpression(data,
                                  row.name = gn,
                                  size = 1.8
                                  )
    grob.expr[[gn]] <- g.expr
  }

  png(file.path(odir,
                "genes-spatial.png"),
      width = w.size * n.cols,
      heigh = h.size * n.rows,
      units = "px")

  grid.arrange(grobs = grob.expr,ncol = n.cols)
  dev.off()

  ## Cluster data using Seurat
  set.seed(1337)
  seu.res <- SeuratAnalysis(data)
  seu.clusters <- seu.res[["cluster.res"]]

  cluster.labels <- unique(seu.clusters$cluster.id)
  n.clusters <- length(cluster.labels)

  ## Plot cluster
  png(file.path(odir,
                "cluster-spatial.png"
              ),
      width = 1000,
      height =1000,
      units = "px"
      )

  g.clu <- VisualizeClusters(plot.data = seu.clusters,
                            size = 6,
                            fontsize =20
                            )
  print(g.clu)
  dev.off()

  ## Visualize distribution for each genes
  ## within every cluster

  grob.bic <- list()
  grob.aic <- list()
  for (lab in cluster.labels){

    spots.idx <- colnames(data)[seu.clusters$cluster.id == lab]

    grobs.fit <- list()
    gof.stat <- data.frame()

    for(gene.name in gene.names){
      print(sprintf("Cluster : %s |Gene : %s",lab,gene.name))

      gene.expr <- as.numeric(data[gene.name,spots.idx])

      fit.res <- FitGene(gene.expr = gene.expr,gene.name)
      fitted <- fit.res$fitted
      gof.stat <- rbind(gof.stat,
                        fit.res$comp.scores)

      plot.data <- data.frame(vals = gene.expr)

      g.fit <- VisualizeFitDist(plot.data,fitted,gene.name)

      grobs.fit[[gene.name]] <- g.fit
    }

    # Visualize count distribution with each cluster
    g.gof <- VisualizeGof(gof.stat,cluster = lab)

    n.dists <- length(colnames(fitted)) - 1
    n.genes <- length(gene.names)

    png(file.path(odir,
                  paste("cluster-",
                        lab,
                        "-fitdist.png",
                        sep = "")),
        width = 400*n.cols,
        height =400*n.rows,
        units = "px")

    Fd <- grid.arrange(grobs = grobs.fit,
                       ncol = n.cols,
                       top = textGrob(paste0("Cluster : ",lab), gp=gpar(fontsize=2*TEXT_SIZE)
                                      )
                       )
    print(Fd)
    dev.off()

    grob.bic[[lab]] <- g.gof$BIC
    grob.aic[[lab]] <- g.gof$AIC
  }

  grob.bic <- grob.bic[order(as.numeric(names(grob.bic)))]
  grob.aic <- grob.bic[order(as.numeric(names(grob.aic)))]

  n.cols.score <- 2
  n.rows.score <- ceiling(n.clusters/n.cols.score)

  png(file.path(odir,
                paste("score-",
                      "BIC.png",
                      sep = "")
                ),
      width = 20*n.genes * n.dists * n.cols.score,
      height = 200 * n.rows.score,
      units = "px")

  grid.arrange(grobs = grob.bic,
               ncol = n.cols.score
               )
  dev.off()

  png(file.path(odir,
                paste("score-",
                      "AIC.png"
                     ,sep = "")),
      width = 20*n.genes * n.dists * n.cols.score,
      height = 200 * n.rows.score,
      units = "px")

  grid.arrange(grobs = grob.aic,
               ncol = n.cols.score
               )

  dev.off()



}

if (!(interactive())) {
  parser <- ArgumentParser()
  parser$add_argument("-c","--count_data",type = "character")
  parser$add_argument("-g","--genes",type = "character")
  parser$add_argument("-o","--out_dir",type = "character")

  args <- parser$parse_args()

  main(data.pth = args$count_data,
       gene.list.pth = args$genes,
       odir = args$out_dir
       )

} else {
  data.pth <- "PATH_TO_DATA"
  gene.list.pth <- "GENE_LIST_PATH"
  odir <- "OUTPUT_DIRECTORY"

  main(data.pth = data.pth,
       gene.list.pth = gene.list.pth,
       odir = odir
       )
}
