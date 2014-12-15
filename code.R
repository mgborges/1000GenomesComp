setwd("/home/murilo/exomePopulation/dep/")

# files <- list.files("../density/", full.names = T)
# 
# leitor = function (fin)
# {
#   b <- strsplit(basename(fin), split = "[.]")
#   sample_phase <- b[[1]][7]
#   sample_pop <- b[[1]][5]
#   id <- b[[1]][1]
#   
#   a <- read.delim(fin, header = F, sep = " ")
#   a <- cbind (a, sample_phase, sample_pop, id)
#   names(a) <- c("Depth", "Counts", "V3", "P(Depth>d)", "Phase", "Population", "id")
#   return (a[-1,])
# }
# 
# dadosCobertura <- do.call(rbind, lapply (files, leitor))
# 
# ggplot(dadosCobertura, aes(log10(Depth+1), log10(Counts), colour = id)) + 
#   geom_line(alpha = .4) + facet_grid(Phase~Population) 


files <- list.files(path = "cat/", full.names = T, pattern = ".cat")

library(data.table)

Data <- fread("/home/murilo/exomePopulation/dep/ALL.table", sep="\t", stringsAsFactors=FALSE, data.table=TRUE)
setnames(Data, 1:121, c("Location", as.character(files)))
# 
# save(Data, file = "/home/murilo/exomePopulation/dep/ALL.table.Rdata")
# source("/home/murilo/exomePopulation/dep/ALL.table.Rdata")
#Data <- as.data.frame(Data)

idx <- !(substr(Data$Location, 1, 1) %in% c("X", "Y", "MT"))

pop <- sapply(strsplit(colnames(Data), "\\."), "[", 5)[-1]
phase <- sapply(strsplit(colnames(Data), "\\."), "[", 7)[-1]


res=vector("list", length(unique(pop))*length(unique(phase)))
group <- vector("character", length(unique(pop))*length(unique(phase)))

#Data <- as.data.table(Data)
i=1
for (p in unique(pop))
{
  for (ph in unique(phase))
  {
    group[i] <- paste(p, ph, sep = "-")
    idpop <- pop==p & phase==ph
    res[[i]] <- rowMeans(log10(Data[idx, which(idpop)+1L, with=FALSE]+1))
    i <- i+1
    print (i)
  }
}

res <- do.call(cbind, res)
colnames(res) <- group

densities <- apply(res, 2, density)
MATX <- do.call(cbind, sapply(densities, "[", "x"))
MATY <- do.call(cbind, sapply(densities, "[", "y"))


library(RColorBrewer)
cols <- brewer.pal(12, "Paired")


#pdf(file="densities.pdf")
matplot(MATX, MATY, type = "l", lty = 1, col = cols, lwd = 2)
legend ("topright", group, col = cols, lwd = 2)

#dev.off()

add=3

#G <- c(1,4,7,10)+add

#pdf(file=paste0(add, "popdensities.pdf"))

G <- (1:3)+3*add

matplot(MATX[,G], MATY[,G], type = "l", lty = 1, col = cols[G], lwd = 2, ylim = c(0,3), xlim = c(0,3))
legend ("topright", group[G], col = cols[G], lwd = 2)

#dev.off()

setwd("/home/murilo/exomePopulation/dep/")

library(Rcpp)

cpp<-sourceCpp("dfdst.cpp")

#Data <-read.delim("ALL.table", header=T, row.names=1)

vector <- dfDst(Data)
dim_data <- dim (Data)

mat <- matrix(0,nc=dim_data[1],nr=dim_data[2])
mat[lower.tri(mat)] = vector

mds<-as.dist(mat)

gc()

fitmds <- cmdscale(mds, eig=TRUE, k=2)

sample_phase <- NULL
sample_pop <- NULL

for (i in seq(1:length(files)))
{
  b <- strsplit(files[i], split = "[.]")
  sample_phase[i] <- b[[1]][7]
  sample_pop[i] <- b[[1]][5]
}

factor <- as.factor(sample_pop)

x <- fitmds$points[,1]
y <- fitmds$points[,2]

pdf(file="MDS.pdf")

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric Multidimensional Scaling for the Average Coverage of the Whole Exome", pch=unclass(factor), col=unclass(as.factor(sample_phase)))

color <- unclass(factor)
legend("topleft", c(unique(as.character(factor)), levels(unclass(as.factor(sample_phase)))), col = c(rep (c(1), times = length(unique(unclass(factor))) ),unclass(as.factor(levels(color)))), pch = c(unique(unclass(factor)),rep(c(15), times = length(levels(color)))), merge = FALSE, border = F, horiz = F, box.lwd = 0, box.col = NULL, bg = NULL)

dev.off()

phaseI<-"20120522"
phaseII<-"20121211"
phaseIII<-"20130415"

fitmds <- cmdscale(mds, eig=TRUE, k=2)

M_pI <- fitmds$points[sample_phase == phaseI,]
M_pII <- fitmds$points[sample_phase == phaseII,]
M_pIII <- fitmds$points[sample_phase == phaseIII,]

library(ICSNP)

pI_versus_pII <- HotellingsT2(M_pI, M_pII)
pI_versus_pIII <- HotellingsT2(M_pI, M_pIII)
pIII_versus_pII <- HotellingsT2(M_pIII, M_pII)

pI_versus_pII
pI_versus_pIII
pIII_versus_pII

##########################

ACB="ACB"
GBR="GBR"
JPT="JPT"
YRI="YRI"

M_ACB <- fitmds$points[sample_pop == ACB,]
M_GBR <- fitmds$points[sample_pop == GBR,]
M_JPT <- fitmds$points[sample_pop == JPT,]
M_YRI <- fitmds$points[sample_pop == YRI,]

GBR_versus_ACB <- HotellingsT2(M_GBR, M_ACB)
GBR_versus_JPT <- HotellingsT2(M_GBR, M_JPT)
GBR_versus_YRI <- HotellingsT2(M_GBR, M_YRI)
ACB_versus_YRI <- HotellingsT2(M_ACB, M_YRI)
ACB_versus_JPT <- HotellingsT2(M_ACB, M_JPT)
YRI_versus_JPT <- HotellingsT2(M_YRI, M_JPT)


GBR_versus_ACB 
GBR_versus_JPT 
GBR_versus_YRI 
ACB_versus_YRI 
ACB_versus_JPT 
YRI_versus_JPT


################## PHASE 

for (i in seq(1:length(files)))
{
  b <- strsplit(files[i], split = "[.]")
  sample_phase[i] <- b[[1]][7]
  sample_pop[i] <- b[[1]][5]
}

phase=phaseI

factor <- as.factor(sample_pop[sample_phase==phase])



x <- fitmds$points[sample_phase == phase,1]
y <- fitmds$points[sample_phase == phase,2]

sample_phase=phase

pdf(file=paste0(phase, "MDS.pdf"))

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric Multidimensional Scaling for the Average Coverage of the Whole Exome", pch=unclass(factor), col=unclass(factor))

color <- unclass(factor)
legend("topleft", c(unique(as.character(factor)), levels(unclass(as.factor(sample_phase)))), col = c(rep (c(1), times = length(unique(unclass(factor))) ),unclass(as.factor(levels(color)))), pch = c(unique(unclass(factor)),rep(c(15), times = length(levels(color)))), merge = FALSE, border = F, horiz = F, box.lwd = 0, box.col = NULL, bg = NULL)

dev.off()

