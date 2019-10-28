#!/usr/bin/Rscript

library(deconvSeq)
library(SingleCellExperiment)
library(argparse)

parser <- ArgumentParser()

parser$add_argument('-sc','--sc_data',
                    type = 'character',
                    help = ''
                    )

parser$add_argument('-st','--st_data',
                    type = 'character',
                    help = ''
                    )

parser$add_argument('-mt','--meta_data',
                    type = 'character',
                    help = ''
                    )

parser$add_argument('-o','--out_dir',
                    type = 'character',
                    default = NULL,
                    help =''
                    )

args <- parser$parse_args()


sc_cnt_pth <- args$sc_data
sc_mta_pth <- args$meta_data
st_cnt_pth <- args$st_data 

dir.create(args$out_dir,
           showWarnings = F)

lbl_vec <- read.table(sc_mta_pth,
                      row.names = 1,
                      sep = '\t',
                      header = T,
                      stringsAsFactors = F)

lbl_vec <- lbl_vec[,"bio_celltype"]

print("read sc count data")

cnts.sc <- read.table(sc_cnt_pth,
                      row.names = 1,
                      sep = '\t',
                      header = T,
                      stringsAsFactors = F)

cnts.sc <- t(cnts.sc)

print("read st count data")

cnts.st <- read.table(st_cnt_pth,
                      row.names = 1,
                      header = T,
                      sep = '\t',
                      stringsAsFactors = F)

cnts.st <- t(cnts.st)

keep.cells <- colSums(cnts.sc) > 300
keep.genes <- rowMeans(cnts.sc) > 0.05

cnts.sc <- cnts.sc[keep.genes,keep.cells]
lbl_vec <- lbl_vec[keep.cells]

ori_lnames <- levels(as.factor(lbl_vec))
lbl_vec <- gsub(",","",lbl_vec)
lbl_vec <- as.factor(lbl_vec)
design.sc = model.matrix(~-1+lbl_vec)

colnames(design.sc) <- levels(lbl_vec)
rownames(design.sc) <- colnames(cnts.sc)

inter <- intersect(rownames(cnts.st),rownames(cnts.sc))
cnts.st <- as.matrix(cnts.st[inter,])
cnts.sc <- as.matrix(cnts.sc[inter,])

print(dim(cnts.st))
print(dim(cnts.sc))


print("Prepping sc data")
#cnts.sc <- prep_scrnaseq(cnts.sc,
#                         genenametype = "hgnc_symbol",
#                         cellcycle =NULL,
#                         count.threshold = 0.05)
#cnts.sc <-getcellcycle(cnts.sc,"G1")
#
#
#cnts.st <- prep_scrnaseq(cnts.st,
#                         genenametype = "hgnc_symbol",
#                         cellcycle =NULL,
#                         count.threshold = 0.05)
#
#cnts.st <-getcellcycle(cnts.st,"G1")


cnts.sc <- SingleCellExperiment(assays = list(counts=cnts.sc))  #counts must be matrix
cnts.sc <- assay(cnts.sc, i ="counts")
#
cnts.st <- SingleCellExperiment(assays = list(counts=cnts.st))  #counts must be matrix
cnts.st <- assay(cnts.st, i ="counts")

print("hej")

dge.sc = getdge(cnts.sc,
                design.sc,
                ncpm.min=1,
                nsamp.min=4,
                method="bin.loess")

print("fit b0")

b0.sc = getb0.rnaseq(dge.sc,
                     design.sc,
                     ncpm.min=1,
                     nsamp.min=4)

print("fir st")
dge.st = getdge(cnts.st,
                NULL,
                ncpm.min=1,
                nsamp.min=4,
                method="bin.loess")

res = getx1.rnaseq(NB0=200,
                   b0.sc,
                   dge.st)

print(colnames(res$x1))
print(ori_lnames)
colnames(res$x1) <- ori_lnames
save(res,file = file.path(args$out_dir,"results.R"))

write.table(res$x1,
            file = file.path(args$out_dir,
                             "deconvSeq-proportions.tsv"),
            sep = '\t',
            quote = F,
            col.names = T,
            row.names = T
            )

