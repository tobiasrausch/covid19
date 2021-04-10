library(ape)

args = commandArgs(trailingOnly=TRUE)
tree = read.tree(args[1])

png(paste0(args[1], ".png"))
plot(tree,no.margin=TRUE,edge.width=2)
dev.off()
