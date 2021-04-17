library(ape)

args = commandArgs(trailingOnly=TRUE)
tree = read.tree(args[1])

#png(paste0(args[1], ".png"), height=4800, width=600)
png(paste0(args[1], ".png"), height=1200, width=600)
plot(tree,no.margin=TRUE,edge.width=2)
dev.off()
