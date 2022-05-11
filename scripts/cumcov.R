library(ggplot2)
library(reshape2)

txtFontSize=16
axisFontSize=16
axisTtlFontSize=18
lgdTtlFontSize=18
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(angle=45, size=axisFontSize, hjust=1), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize), strip.text.x=element_text(size=axisTtlFontSize))

# column1: chrom, column2: pos, column3...: depth
args = commandArgs(trailingOnly=TRUE)
x=read.table(args[1], header=F)
x[,3:ncol(x)] = x[,3:ncol(x)] / colSums(x[,3:ncol(x)])
x$medcov = apply(x[,3:ncol(x)], 1, median)
png(paste0(args[1], ".png"), width=1200, height=400)
p = ggplot(data=x, aes(x=V2, y=medcov))
p = p + geom_line()
p = p + xlab("SARS-CoV-2 genome") + ylab("Median normalized coverage")
p = p + scienceTheme + theme(legend.position="top")
p = p + geom_vline(xintercept = 22850, linetype="dashed")
p = p + geom_vline(xintercept = 23150, linetype="dashed")
p = p + geom_hline(yintercept = 0, linetype="dashed")
p
dev.off()
