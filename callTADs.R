library(tidyverse)
library(ggplot2)
if (!require("rGMAP")) {
  library(devtools)
  install_github("ningbioinfostruggling/rGMAP")
}
library(rGMAP)

options(warn=-1)

args = commandArgs(trailingOnly=TRUE)

## arguments
## 1: matrix file
## 2: bin file
## 3: resolution
## 4: chr
## 5: output bed name
## 6: output plot name

if(length(args) < 6){
  stop("6 arguments are expected, including 1.matrix file 2.bin file 3.resolution 4.chr 5.output bed name 6.output plot name.")
}

mfile <- args[1]
bfile <- args[2]
res <- as.integer(args[3])
chr <- args[4]

if (args[4] != "all") {
  idx <- read_delim(bfile, delim = '\t', col_names = F) %>%
    filter(X1 == chr)

  startbin <- idx[1,4] %>% as.integer()

  endbin <- idx[nrow(idx),4] %>% as.integer()

  hic <- read_delim(mfile, delim = "\t", col_names = F) %>%
    filter(X1 >= startbin) %>%
    filter(X1 <= endbin)
} else if (args[4] == "all"){
  hic <- read_delim(mfile, delim = "\t", col_names = F)
}

hic_tad <- rGMAP(hic, index_file = bfile, resl = res)

hic_tad$tads %>%
  mutate(start = as.character(start), end = as.character(end)) %>%
  write_delim(args[5], delim = '\t', col_names = F)

pp <- plotdom(hic, hic_tad$hierTads, start_bin = 1, end_bin = 500, cthr = 20, resl = res)

png(filename = args[6])

pp$p2 + ggtitle(paste("TADs and sub-TADs with", res, "resolution")) +
  theme_bw() +
  xlab("Coordinate (Mb)") +ylab("Normalised interaction count")

dev.off()
