library(seriation)

S <- pc_loadings%*%t(pc_loadings)

max = 2 * sd(S)

S[which(S>=max)] = max

S[which(S<=-max)] = -max


pimage(S)

#S1 <- abs(S[1:618, 1:618])

pdf(file='Covplot.pdf')

summary(c(S))

dev.off()

hist(c(S), breaks = 100)




pimage(S)
