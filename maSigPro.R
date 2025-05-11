library(maSigPro)

setwd("data/location")
male_sampletable = read.csv("Example_sampletable.csv", header = T, row.names = 1)
male_datatable = read.csv("Example_normalised counts.csv", header = T, row.names = 1)

str(male_datatable)

#Two step polynomial regression model

design <- make.design.matrix(male_sampletable, degree=4)
fit <- p.vector(male_datatable, design, Q = 0.05, MT.adjust = "BH", min.obs = 15)
write.table(fit$SELEC, file="Example_Masigpro_Significant_Genes.csv") #First step of regression model 

tstep <- T.fit(fit, alfa=0.05)
sigs <- get.siggenes(tstep, rsq = 0.6, vars="all") #Second step of regression model 
write.table(sigs$sig.genes$sig.profiles, file="Example_SigProfiles.csv") 
write.table(sigs$sig.genes$sig.pvalues, file="Example_SigPvalues.csv") 

#Venn Diagram of Significant Genes 

sigs <- get.siggenes(tstep, rsq =0.6, vars = "groups")
names(sigs$sig.genes)
write.table(sigs$summary, file="Example_SigGenes_Treatment.csv")
suma2Venn(sigs$summary[,c(1:2)])

#Clustering of Significant Genes 

cluster<-see.genes(sigs$sig.genes$AdlibvsFasting, show.fit=T, dis=design$dis,cluster.method='hclust', cluster.data =1, k=6)

write.table(cluster$cut, file="Example_clusters.csv")
write.table(sigs$sig.genes$AdlibvsFasting$sig.profiles, file="Example_AdlibvsFastingsigprofiles.csv")
write.table(sigs$sig.genes$AdlibvsFasting$sig.pvalues, file="Example_AdlibvsFastingsigpvalues.csv")
