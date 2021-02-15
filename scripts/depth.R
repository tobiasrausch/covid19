library(ggplot2)

args=commandArgs(trailingOnly=TRUE)
x = read.table(args[1], header=F)
x$V2 = x$V2 - 1
print(summary(x))

png(paste0(args[1], ".png"), width=1200, height=600)
p = ggplot(data=x, aes(x=V1, y=V2))
p = p + geom_line() + ylab("#Samples <10x") + xlab("Position")
p = p + scale_x_continuous(breaks=(1:30*1000))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
dev.off()
