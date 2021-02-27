library(ggplot2)
library(reshape2)
library(grid)
library(scales)

# Split by plate first
# cat contam.tsv | sed 's/\([A-H]\)\([0-9]*\)\t/\t\1\t\2\t/' > byplate.tsv
# Rscript ../scripts/contam.R byplate.tsv

args = commandArgs(trailingOnly=TRUE)
df = read.table(args[1], header=T)
colnames(df) = c("Plate", "Row", "Column", "Freemix")

for (plate in unique(df$Plate)) {
    sub = df[df$Plate == plate,]
    p = ggplot(data=sub)
    p =  p + geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=1, fill=Freemix))
    p = p + facet_grid(Row~Column)
    p = p + theme(legend.position="bottom", axis.text = element_blank(), axis.ticks = element_blank(), panel.grid  = element_blank())
    p = p + labs(fill="Contamination") + xlab("") + ylab("")
    p = p + geom_text(aes(x=0.5, y=0.5, label=Freemix), vjust=-1)
    p = p + ggtitle(paste0(plate, " contamination estimates"))
    p = p + scale_fill_gradient2(limits=c(0, 1), mid="#f03b20", high="#bd0026", low="#f7f7f7", midpoint=0.5)
    ggsave(p, file=paste0(args[1], ".", plate, ".png"), width=14, height=7)
}
