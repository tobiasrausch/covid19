# To compare multiple library prep strategies (requires Setup column)

library(ggplot2)


x=read.table("summary.tsv", header=T, sep="\t", comment.char="$")

print(summary(x))

ggplot(data=x, aes(x=Setup, y=RKI)) + geom_jitter(aes(color=Setup), width=0.1, height=0.1)
ggplot(data=x, aes(x=Setup, y=MedianCoverage)) + geom_boxplot(aes(color=Setup), outlier.shape=NA) + geom_jitter(aes(color=Setup), width=0.05, height=0)
ggplot(data=x, aes(x=MedianCoverage, y=X.ConsensusNs)) + geom_point(aes(color=Setup, fill=Setup, shape=RKI))
ggplot(data=x, aes(x=MedianCoverage, y=X.ConsensusAmbiguous)) + geom_point(aes(color=Setup, fill=Setup, shape=RKI))
ggplot(data=x, aes(x=Setup, y=MedianInsertSize)) + geom_boxplot(outlier.shape=NA) + geom_jitter(aes(color=RKI), width=0.05, height=0)
ggplot(data=x, aes(x=MedianCoverage, y=SDCoverage)) + geom_point(aes(color=Setup, fill=Setup, shape=RKI))
