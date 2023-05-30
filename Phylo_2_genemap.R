
### Check package availability ###

list.of.packages <- c("ape", "ggtree","tidyverse","gggenes","tidytree")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if("ggtree" %in% new.packages){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ggtree")
}

### Load packages ###

suppressMessages(suppressWarnings(library("ape")))
suppressMessages(suppressWarnings(library("ggtree")))
suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(suppressWarnings(library("gggenes")))
suppressMessages(suppressWarnings(library("tidytree")))

### Parse Arguments ###

args <- commandArgs(trailingOnly = TRUE)

#args <- c("--tree","/Users/muratbuyukyorukmsu/Desktop/ANK/final_tree.nexus","--genemap","/Users/muratbuyukyorukmsu/Desktop/ANK/Pro_ANK_containing_hits_genemap_inuse.txt","--anchor","ANK_1", "--clade", "/Users/muratbuyukyorukmsu/Desktop/ANK/test_for_R.txt","--output","outfile.txt")

filename <- args[which(args=="--tree")+1]
genemap_in <-  args[which(args=="--genemap")+1]
anc_point <-  args[which(args=="--anchor")+1]
clade_in <-  args[which(args=="--clade")+1]
file.out <- args[which(args=="--output")+1]

### Import Tree ###

treein <- read.tree(filename)

x <- treein
p1 <- ggtree(x)

### Color branches if Clade info provided ###

if(length(clade_in)!=0){
  quants <- list()
  
  clade_df <- read.delim(clade_in, header = TRUE)
  clade_append <- NULL
  
  for(i in 1:ncol(clade_df)){
    clade_info <- paste(as.list(clade_df[,i][clade_df[,i]!=""]),collapse = "|")
    assign(paste("clade",i,sep = "_"),grep(clade_info,treein$tip.label))
    if(i==1){
      clade_append <- clade_info
    }else{
      clade_append <- paste(clade_append,clade_info,sep = "|")
    }
    if(i==ncol(clade_df)){
      uncladed <- grep(clade_append,treein$tip.label,invert = TRUE)
    }
  }
  for(i in 1:ncol(clade_df)){
    if (length(uncladed)!=0){
      quants[["Rest"]] <- uncladed
    }
    quants[[paste0(colnames(clade_df)[i])]] <- get(paste("clade",i,sep = "_"))
  }
  
  p1 <- groupOTU(p1, quants, 'Clade') + aes(color=Clade)
}

### Define genemap function ###

get_genes <- function(data, genome) {
  filter(data, molecule == genome) %>% pull(gene)
}

### Import genemap dataframe ###

genemap <- read.delim(genemap_in,header = TRUE)

### Generate genemap matrix (takes a long time when the df is long) ###

g <- unique(genemap[,1])
n <- length(g)
d <- matrix(nrow = n, ncol = n)
rownames(d) <- colnames(d) <- g
genes <- lapply(g, get_genes, data = genemap)

for (i in 1:n) {
  for (j in 1:i) {
    jaccard_sim <- length(intersect(genes[[i]], genes[[j]])) / 
      length(union(genes[[i]], genes[[j]]))
    d[j, i] <- d[i, j] <- 1 - jaccard_sim
  }
}

### Map domains to ggtree plot ###

p <- p1 + 
  geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene),
             data = genemap, geom = geom_motif, panel = 'Alignment',
             on = anc_point ,arrowhead_height = unit(0, "mm"), arrowhead_width = unit(0, "mm"), arrow_body_height = unit(0.1,"mm"),size = 0) + 
  scale_x_continuous(expand=c(0,0)) #+ theme_genes()

p <- facet_widths(p, widths=c(1,3))

setEPS()
postscript(paste0(tools::file_path_sans_ext(file.out),".eps"))
Sys.sleep(2)
p
dev.off()

pdf(paste0(tools::file_path_sans_ext(file.out),".pdf"))
Sys.sleep(2)
p
dev.off()

print(paste0("Raw plot is exported as PDF and SVG files and can be found in ",getwd()))


