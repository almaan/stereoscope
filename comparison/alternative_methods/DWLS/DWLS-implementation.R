#!/usr/bin/Rscript

library(zeallot)
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

parser$add_argument('-wd','--workdir',
                    type = 'character',
                    default = NULL,
                    help =''
                    )

args <- parser$parse_args()

workdir <-  ifelse(is.null(args$workdir),
                 'project',
                 args$workdir)

sc_cnt_pth <- args$sc_data
sc_mta_pth <- args$meta_data
st_cnt_pth <- args$st_data 




dir.create(file.path(workdir,"results"), showWarnings = F)

source("Modded_Deconvolution_functions.R")
print("Loading SC count data")
dataSC <-  read.table(sc_cnt_pth,
                      sep = '\t',
                      header = T,
                      row.names = 1,
                      stringsAsFactors = F
                      )

setwd(workdir)

print(dataSC[1:10,1:10])
print("Loading ST count data")

dataST <- read.table(st_cnt_pth,
                     sep = '\t',
                     header = T,
                     row.names = 1,
                     stringsAsFactors = F
                     )
print(dataST[1:10,1:10])
print("Loading SC labels")

labels <- read.table(sc_mta_pth,
                     sep = '\t',
                     header = T,
                     row.names = 1,
                     stringsAsFactors = F
                     )
labels <- labels['bio_celltype']
print(labels[1:10,])
print("Prepare Data for analysis")
intermeta <- intersect(rownames(dataSC),rownames(labels))
dataSC <- dataSC[intermeta,]
labels <- labels[intermeta,]
old_labels <- labels
labels <- gsub(',| |\\.|-','_',labels, perl = T)
types <- unique(labels)
print(types)
n_types <- length(types)

interst <- intersect(colnames(dataSC),colnames(dataST))
dataSC <- dataSC[,interst]
dataST <- dataST[,interst]

dataSC <- t(dataSC) # transposition converts to matrix
labels <- as.vector(unlist(labels))

Signatures <- buildSignatureMatrixMAST(scdata=dataSC,
                                       id=labels,
                                       path="results",
                                       fitDEA = T,
                                       diff.cutoff=0.5,
                                       pval.cutoff=0.01
                                       ) 
c(n_spots,n_genes) %<-% dim(dataST)

prop_mat <- as.data.frame(matrix(0,
                                 nrow = n_spots,
                                 ncol = n_types
                                )
                         )
rownames(prop_mat) <- rownames(dataST)
colnames(prop_mat) <- unique(old_labels)

print("Estimate propotions in each spot")
dataST <- as.matrix(dataST)
for (s in 1:n_spots) {

    print(sprintf("Estimating proportion for spot : %d / %d",
                  s,n_spots)
          )

    spot <- dataST[s,]
    tr <- trimData(Signatures,spot)

    tr$sig <- tr$sig[,colSums(tr$sig) > 0]
    is_pd <- eigen(t(tr$sig)%*%tr$sig)$values
    print(is_pd)
    is_pd <- all(is_pd > 10e-8)

    if (!(is_pd)) { 
        next
    }

    solDWLS <- solveDampenedWLS(tr$sig,tr$bulk)
    print("Proportions >> ")
    print(solDWLS)
    prop_mat[s,names(solDWLS)] <- solDWLS
}

write.table(prop_mat,
            file = 'results/DWLS-proportions.tsv',
            sep = '\t',
            quote = F,
            col.names = T,
            row.names = T
            )

