gili_kaks <- read.csv("/Users/jeffcole/Git_repo/gili_dnds.csv",header=T, stringsAsFactors = F)

gili_sperm <-read.csv("/Users/jeffcole/Git_repo/D.plex_sperm_proteome.csv", header= T, stringsAsFactors = F)

gili_chrom <-read.table("/Users/jeffcole/Git_repo/D.plex_chromosomes.txt", header=T)
gili_zlinked1 <-subset(gili_chrom, dpChr==1)
gili_zlinked16 <-subset(gili_chrom, dpChr==16)
gili_zlinked <- rbind(gili_zlinked1, gili_zlinked16)

#this removes NA from Ka row
gili_no_NA <- subset(gili_kaks, !(is.na(gili_kaks$Ka)))
#there are still rows from Ks that have NA, this resolves that
gili_no_NA1<- subset(no_NA, !(is.na(gili_no_NA$Ks)))
#restore NA to kaks object if needed
gili_kaks <-gili_no_NA1

#####this is an ok spot to skip the nodup step and proceed to splitting into type subsets, merging, and run analysis to see
#if it made a difference filtering it to have on hit per gene instead of one hit per contig

spermy <- gili_kaks[gili_kaks$gene %in% gili_sperm$OGS2 ,]
autosomy <- gili_kaks[!(gili_kaks$gene %in% gili_zlinked$gene.dp) ,]
autosomy$KaKs <- as.numeric(autosomy$KaKs)
spermy$KaKs <- as.numeric(spermy$KaKs)
wilcox.test(spermy$KaKs,autosomy$KaKs)

######kwallis for non duplicated##############################

gili_kaksy <- gili_kaks[,c("gene","Ka","Ks", "KaKs")] 
gili_kaksy$Ka <- as.numeric(gili_kaksy$Ka)
gili_kaksy$Ks <- as.numeric(gili_kaksy$Ks)
gili_kaksy$KaKs <- as.numeric(gili_kaksy$KaKs)
gili_kaksy[gili_kaksy>2]=NA
gili_kaksy$gene <- gili_kaks$gene



gili_apyreney <- gili_kaksy[gili_kaksy$gene %in% gili_apyrene$OGS2 ,]
gili_eupyreney <- gili_kaksy[gili_kaksy$gene %in% gili_eupyrene$OGS2 ,]
gili_sharedy <- gili_kaksy[gili_kaksy$gene %in% gili_apy.eupy$OGS2 ,]
gili_dnds_apyreney <- data.frame(x=1:length(gili_apyreney$Ka), gili_apyrene_Ka=gili_apyreney$Ka,gili_apyrene_Ks=gili_apyreney$Ks, gili_apyrene_KaKs=gili_apyreney$KaKs)
gili_dnds_eupyreney <- data.frame(x=1:length(gili_eupyreney$Ka), gili_eupyrene_Ka=gili_eupyreney$Ka,gili_eupyrene_Ks=gili_eupyreney$Ks, gili_eupyrene_KaKs=gili_eupyreney$KaKs)
gili_dnds_sharedy <- data.frame(x=1:length(gili_sharedy$Ka), gili_apy.eupy_Ka=gili_sharedy$Ka,gili_apy.eupy_Ks=gili_sharedy$Ks, gili_apy.eupy_KaKs=gili_sharedy$KaKs)
gili_dnds_listy <- list(gili_dnds_apyreney, gili_dnds_eupyreney, gili_dnds_sharedy)
gili_dndsy <- Reduce(function(...) merge(..., all=T), gili_dnds_listy)

kay<- data.frame(Ka=c(gili_dndsy$gili_apyrene_Ka,gili_dndsy$gili_eupyrene_Ka,gili_dndsy$gili_apy.eupy_Ka),type= c(rep("apyrene" , length(gili_dndsy$gili_apyrene_Ka)),rep("eupyrene",length(gili_dndsy$gili_eupyrene_Ka)),rep("apy.eupy",length(gili_dndsy$gili_apyrene_Ka))))
kruskal.test(kay$Ka~kay$type)
ksy<- data.frame(Ks=c(gili_dndsy$gili_apyrene_Ks,gili_dndsy$gili_eupyrene_Ks,gili_dndsy$gili_apy.eupy_Ks),type= c(rep("apyrene" , length(gili_dndsy$gili_apyrene_Ks)),rep("eupyrene",length(gili_dndsy$gili_eupyrene_Ks)),rep("apy.eupy",length(gili_dndsy$gili_apyrene_Ks))))
kruskal.test(ksy$Ks~ksy$type)
kaksy<- data.frame(KaKs=c(gili_dndsy$gili_apyrene_KaKs,gili_dndsy$gili_eupyrene_KaKs,gili_dndsy$gili_apy.eupy_KaKs),type= c(rep("apyrene" , length(gili_dndsy$gili_apyrene_KaKs)),rep("eupyrene",length(gili_dndsy$gili_eupyrene_KaKs)),rep("apy.eupy",length(gili_dndsy$gili_apyrene_KaKs))))
kruskal.test(kaksy$KaKs~kaksy$type)
##############################################################

# this sorts gene ID in descending order, for duplicate gene IDs its length will also be listed in descending order 
gili_kaks_order <- gili_kaks[order(gili_kaks$gene, gili_kaks$length),]
#this chops all duplicated gene ID leaving highest value length
gili_kaks_nodup <- gili_kaks_order[!duplicated(gili_kaks_order$gene),]
gili_kaks_nodup <- gili_kaks_nodup[-1,]
gili_kaks_nodup$Ka <- as.numeric(gili_kaks_nodup$Ka)
gili_kaks_nodup$Ks <- as.numeric(gili_kaks_nodup$Ks)
gili_kaks_nodup$KaKs <- as.numeric(gili_kaks_nodup$KaKs)
#from here we need to know what percentage of the genome we sampled
#percent genome sampled = length(genes_with_dnds)/15,130
gili_percent_genome_sampled<-length(unique(gili_kaks_nodup$gene))/15130

#a small example to demonstrate duplicate lengths for duplicate gene IDs is not a problem
#d <- data.frame(id=c(4,1,2,3,3,3,4),values=c(1,3,5,8,8,5,9))
#d_order <- d[order(d$id, -d$values),]
#d_nodup <- d_order[!duplicated(d_trimmed$id),]


gili_autosome <- gili_kaks_nodup[!(gili_kaks_nodup$gene %in% gili_zlinked$gene.dp) ,]
#we now have all of the autosome genes that have dnds ratios calculated for them
#we need to know percentage of autosome sampled
#total autosome = 15,130-zlinked
#percent autosome sampled = gili_autosome/total autosome
gili_percent_autosome_sampled<-length(unique(gili_autosome$gene))/(15130-length(unique(gili_zlinked$gene.dp)))


#uncomment if necessary to remove sperm proteins from autosomal data for analysis purposes
#autosome_kaks <- autosome[!(autosome$gene %in% sperm$OGS2) ,]

gili_sperm_kaks <- gili_kaks_nodup[gili_kaks_nodup$gene %in% gili_sperm$OGS2 ,]
#from here, find out what percent of sperm proteome we sampled
#total sperm = length(unique(gili_sperm$OGS2))
#percent sperm sampled =  length(gili_sperm_kaks$gene)/total sperm
gili_percent_sperm_sampled<-length(unique(gili_sperm_kaks$gene))/length(unique(gili_sperm$OGS2))


gili_zlinked_kaks <- gili_kaks_nodup[gili_kaks_nodup$gene %in% gili_zlinked$gene.dp ,]
#from here, find out what percent of zlinked we sampled
#total zilnked = length(zlinked$dpChr)
#percent zlinked sampled = length(gili_zlinked_kaks$gene)/total zlinked
gili_percent_zlinked_sampled<-length(unique(gili_zlinked_kaks$gene))/length(unique(gili_zlinked$gene.dp))
#let's find out how many sperm genes are zlinked
sperm_in_z<-(gili_zlinked_kaks[gili_zlinked_kaks$gene %in% sperm$OGS2,])
length(sperm_in_z$gene)

gili_apyrene <- subset(gili_sperm, Apyrene==TRUE)
gili_apyrene_kaks <- gili_kaks_nodup[gili_kaks_nodup$gene %in% gili_apyrene$OGS2 ,]
#from here, find out what percentage of apyrene was sampled
#total apyrene = length(unique(gili_apyrene$genome))
#percent apyrene sampled = gili_apyrene_kaks$gene/total apyrene
gili_percent_apyrene_sampled<-length(unique(gili_apyrene_kaks$gene))/length(unique(gili_apyrene$OGS2))

gili_eupyrene <- subset(gili_sperm, Eupyrene==TRUE)
gili_eupyrene_kaks <- gili_kaks_nodup[gili_kaks_nodup$gene %in% gili_eupyrene$OGS2 ,]
#from here, find out what percentage of eupyrene was sampled
#total eupyrene = length(unique(gili_eupyrene$genome))
#percent eupyrene sampled = gili_eupyrene_kaks$gene/total eupyrene
gili_percent_eupyrene_sampled<-length(unique(gili_eupyrene_kaks$gene))/length(unique(gili_eupyrene$OGS2))

gili_apy.eupy <- subset(gili_sperm, AE==TRUE)
gili_apy.eupy_kaks <- gili_kaks_nodup[gili_kaks_nodup$gene %in% gili_apy.eupy$OGS2 ,]
#from here, find out what percentage of shared was sampled
#total shared = length(unique(gili_apy.eupy$genome))
#percent shared sampled = gili_apy.eupy$gene/total shared
gili_percent_shared_sampled<-length(unique(gili_apy.eupy_kaks$gene))/length(unique(gili_apy.eupy$OGS2))

#put percent sampled into one dataframe
gili_percent_sampled<-data.frame(gili_genome=gili_percent_genome_sampled,gili_autosome=gili_percent_autosome_sampled,gili_zlinked=gili_percent_zlinked_sampled,gili_sperm=gili_percent_sperm_sampled,gili_apyrene=gili_percent_apyrene_sampled,gili_eupyrene=gili_percent_eupyrene_sampled,gili_shared=gili_percent_shared_sampled)


#from here plots and analysis can begin
#dN_boxplot <- boxplot(kaks_nodup$Ka, autosome$Ka,zlinked_kaks$Ka, sperm_kaks$Ka, apyrene_kaks$Ka, eupyrene_kaks$Ka, apy.eupy_kaks$Ka, main="dN", names=c("genome","autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
#dS_boxplot <- boxplot(kaks_nodup$Ks, autosome$Ks,zlinked_kaks$Ks, sperm_kaks$Ks, apyrene_kaks$Ks, eupyrene_kaks$Ks, apy.eupy_kaks$Ks, main="dS", names=c("genome","autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
#dNdS_boxplot <- boxplot(kaks_nodup$KaKs, autosome$KaKs,zlinked_kaks$KaKs, sperm_kaks$KaKs, apyrene_kaks$KaKs, eupyrene_kaks$KaKs, apy.eupy_kaks$KaKs, main="dN/dS", names=c("genome","autosome","Z-linked","sperm","apyrene","eupyrene","shared"))


#put into one dataframe
gili_dnds_autosome <- data.frame(x=1:length(gili_autosome$Ka), gili_autosome_Ka=gili_autosome$Ka,gili_autosome_Ks=gili_autosome$Ks, gili_autosome_KaKs=gili_autosome$KaKs)
gili_dnds_zlinked <- data.frame(x=1:length(gili_zlinked_kaks$Ka), gili_zlinked_Ka=gili_zlinked_kaks$Ka,gili_zlinked_Ks=gili_zlinked_kaks$Ks, gili_zlinked_KaKs=gili_zlinked_kaks$KaKs)
gili_dnds_sperm <- data.frame(x=1:length(gili_sperm_kaks$Ka), gili_sperm_Ka=gili_sperm_kaks$Ka,gili_sperm_Ks=gili_sperm_kaks$Ks, gili_sperm_KaKs=gili_sperm_kaks$KaKs)
gili_dnds_apyrene <- data.frame(x=1:length(gili_apyrene_kaks$Ka), gili_apyrene_Ka=gili_apyrene_kaks$Ka,gili_apyrene_Ks=gili_apyrene_kaks$Ks, gili_apyrene_KaKs=gili_apyrene_kaks$KaKs)
gili_dnds_eupyrene <- data.frame(x=1:length(gili_eupyrene_kaks$Ka), gili_eupyrene_Ka=gili_eupyrene_kaks$Ka,gili_eupyrene_Ks=gili_eupyrene_kaks$Ks, gili_eupyrene_KaKs=gili_eupyrene_kaks$KaKs)
gili_dnds_apy.eupy <- data.frame(x=1:length(gili_apy.eupy_kaks$Ka), gili_apy.eupy_Ka=gili_apy.eupy_kaks$Ka,gili_apy.eupy_Ks=gili_apy.eupy_kaks$Ks, gili_apy.eupy_KaKs=gili_apy.eupy_kaks$KaKs)
gili_dnds_list <- list(gili_dnds_autosome, gili_dnds_zlinked, gili_dnds_sperm, gili_dnds_apyrene, gili_dnds_eupyrene, gili_dnds_apy.eupy)
gili_dnds <- Reduce(function(...) merge(..., all=T), gili_dnds_list)
#see if sperm are diff before removing outliers
wilcox.test(gili_dnds$gili_sperm_KaKs, gili_dnds$gili_autosome_KaKs)
#filter values>2
gili_dnds[gili_dnds>2]=NA




##############THINGS I NEED TO DO#################
######1. using final dnds dataframe#######
#create boxplot for Ka, Ks, and KaKs

boxplot(gili_dnds$gili_autosome_Ka, gili_dnds$gili_zlinked_Ka,gili_dnds$gili_sperm_Ka, gili_dnds$gili_apyrene_Ka, gili_dnds$gili_eupyrene_Ka, gili_dnds$gili_apy.eupy_Ka, outline=F, main="dN", names=c("autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
boxplot(gili_dnds$gili_autosome_Ks, gili_dnds$gili_zlinked_Ks,gili_dnds$gili_sperm_Ks, gili_dnds$gili_apyrene_Ks, gili_dnds$gili_eupyrene_Ks, gili_dnds$gili_apy.eupy_Ks, outline=F, main="dS", names=c("autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
boxplot(gili_dnds$gili_autosome_KaKs, gili_dnds$gili_zlinked_Ka,gili_dnds$gili_sperm_KaKs, gili_dnds$gili_apyrene_KaKs, gili_dnds$gili_eupyrene_KaKs, gili_dnds$gili_apy.eupy_KaKs, outline=F, main="dNdS", names=c("autosome","Z-linked","sperm","apyrene","eupyrene","shared"))

#create histograms for all Ka, Ks and KaKs
hist(gili_dnds$gili_autosome_Ka)
hist(gili_dnds$gili_autosome_Ks)
hist(gili_dnds$gili_autosome_KaKs)
hist(gili_dnds$gili_zlinked_Ka)
hist(gili_dnds$gili_zlinked_Ks)
hist(gili_dnds$gili_zlinked_KaKs)
hist(gili_dnds$gili_sperm_Ka)
hist(gili_dnds$gili_sperm_Ks)
hist(gili_dnds$gili_sperm_KaKs)
hist(gili_dnds$gili_apyrene_Ka)
hist(gili_dnds$gili_apyrene_Ks)
hist(gili_dnds$gili_apyrene_KaKs)
hist(gili_dnds$gili_eupyrene_Ka)
hist(gili_dnds$gili_eupyrene_Ks)
hist(gili_dnds$gili_eupyrene_KaKs)
hist(gili_dnds$gili_apy.eupy_Ka)
hist(gili_dnds$gili_apy.eupy_Ks)
hist(gili_dnds$gili_apy.eupy_KaKs)



#run wilcox test for sperm Ka, Ks, and KaKs vs autosome, 
#also findout if sperm is all sperm or just leftovers, may need new column
wilcox.test(gili_dnds$gili_autosome_Ka, gili_dnds$gili_sperm_Ka)
wilcox.test(gili_dnds$gili_autosome_Ks, gili_dnds$gili_sperm_Ks)
wilcox.test(gili_dnds$gili_autosome_KaKs, gili_dnds$gili_sperm_KaKs)


#run kWallis for apyrene, eupyrene, shared with Ka, Ks, and KaKs (maybe extend to sperm general and zlinked)
#try this, if not it may need a separate column denoting categorical factor
#kruskal.test(dnds$apyrene_Ka, dnds$eupyrene_Ka, dnds$apy.eupy_Ka)
#kruskal.test(dnds$apyrene_Ks, dnds$eupyrene_Ks, dnds$apy.eupy_Ks)
#kruskal.test(dnds$apyrene_KaKs, dnds$eupyrene_KaKs, dnds$apy.eupy_KaKs)

#dataframe for kwallis


#Ka_apy <- data.frame(ka=dnds$apyrene_Ka, type=rep("apyrene", length(dnds$apyrene_Ka))
#Ka_eupy <- data.frame(ka=dnds$eupyrene_Ka, type=rep("eupyrene", length(dnds$apyrene_Ka))

gili_Ka_2col<- data.frame(Ka=c(gili_dnds$gili_apyrene_Ka,gili_dnds$gili_eupyrene_Ka,gili_dnds$gili_apy.eupy_Ka),type= c(rep("apyrene" , length(gili_dnds$gili_apyrene_Ka)),rep("eupyrene",length(gili_dnds$gili_eupyrene_Ka)),rep("apy.eupy",length(gili_dnds$gili_apyrene_Ka))))
kruskal.test(gili_Ka_2col$Ka~gili_Ka_2col$type)

gili_Ks_2col<- data.frame(Ks=c(gili_dnds$gili_apyrene_Ks,gili_dnds$gili_eupyrene_Ks,gili_dnds$gili_apy.eupy_Ks),type= c(rep("apyrene" , length(gili_dnds$gili_apyrene_Ks)),rep("eupyrene",length(gili_dnds$gili_eupyrene_Ks)),rep("apy.eupy",length(gili_dnds$gili_apyrene_Ks))))
kruskal.test(gili_Ks_2col$Ks~gili_Ks_2col$type)

gili_KaKs_2col<- data.frame(KaKs=c(gili_dnds$gili_apyrene_KaKs,gili_dnds$gili_eupyrene_KaKs,gili_dnds$gili_apy.eupy_KaKs),type= c(rep("apyrene" , length(gili_dnds$gili_apyrene_KaKs)),rep("eupyrene",length(gili_dnds$gili_eupyrene_KaKs)),rep("apy.eupy",length(gili_dnds$gili_apyrene_KaKs))))
kruskal.test(gili_KaKs_2col$KaKs~gili_KaKs_2col$type)


##############2. compare outlier removed from outlier remained#################
#determine how many outliers were removed from each column, if > 100, no bueno
#also maybe investigate whether or not trimming the duplicates affected the significance by performing tests on untrimmed data

##################3. determine percentage of genes sampled######################
#there are ~15,000 genes in D.plex, how many of those genes did we get ratios for
#there are ~800 genes in the sperm proteome, how many of those did we get ratios for
#try repeating this for all subsets, find proportion of genes with measurements to possible genes
#I did this piece by piece as I was constructing each individual element



##################4. ANCOVA dN~dS,sperm,zlinked#####################
#read ch 10 in Jamie's book
plot(gili_kaks_nodup$Ks,gili_kaks_nodup$Ka)
#super messy plot, maybe there are predictors? There may not be
anco_dnds <- kaks[,c("gene","Ka","Ks","KaKs")]
anco_dnds$sperm <- anco_dnds$gene %in% sperm$OGS2
anco_dnds<-subset(anco_dnds, Ks<2)
anco_dnds<-subset(anco_dnds, !(is.na(anco_dnds$Ka)))
mod_dnds <- aov(Ka~Ks*sperm,data=anco_dnds)
summary(mod_dnds)
anco_dnds$zlinked <- anco_dnds$gene %in% gili_zlinked$gene.dp
mod2_dnds<-aov(Ka~Ks*sperm*zlinked, data=anco_dnds) #dN as a function of dS given sperm or zlinked

anco_sperm <- subset(anco_dnds, sperm==TRUE)
anco_sperm$apyrene <- anco_sperm$gene %in% gili_apyrene$OGS2
anco_sperm$eupyrene <- anco_sperm$gene %in% gili_eupyrene$OGS2
anco_sperm$shared <- anco_sperm$gene %in% gili_apy.eupy$OGS2
mod3_dnds<-aov(Ka~Ks*apyrene*eupyrene*shared,data=anco_sperm)

yes_sperm <- subset(anco_dnds, sperm==TRUE)
no_sperm <- subset(anco_dnds, sperm==FALSE)

yes_zlinked <- subset(anco_dnds, zlinked==TRUE)
no_zlinked <- subset(anco_dnds, zlinked==FALSE)

plot(Ka~Ks, data=anco_dnds)
points(yes_sperm$Ka, yes_sperm$Ks,pch=20)
points(yes_zlinked$Ka, yes_zlinked$Ks,pch=1)

#plot(sperm, Ka,data=anco_dnds)


################5. BOOTSTRAP, insert here because many of the subsets and files loaded make sense for it to be here
divs <- kaks[,c("gene","s.sites", "n.sites", "s.substitions", "n.substitutions" )] # extract relevant columns, just to make life easier

names(divs) <- c("gene","s.sites", "n.sites", "s.subs", "n.subs") # rename columns, just to make life easier.  also note typo: s.substitions
apply(divs, 2, function(x){sum(is.na(x))})  # count how many NAs are in each column, just to know what the data look like

# convert NA to 0 for substitution counts.
divs$s.subs[is.na(divs$s.subs)] <- 0
divs$n.subs[is.na(divs$n.subs)] <- 0
#force site counts to numeric.  I don't know why this is reading in as character, but we need to force it to numeric.
divs[[2]] <- as.numeric(divs[[2]])
divs[[3]] <- as.numeric(divs[[3]])
divs <- head(divs,-3) #takes care of random 3 lines of NA
divs$s.sites[is.na(divs$s.sites)] <- 0
divs$n.sites[is.na(divs$n.sites)] <- 0

bulk.kaks <- function(divframe) {  # function to get ka, ks, and ka/ks (ie w, or "omega") from the "divs" dataframe, or some subset of it. Function returns a vector with names.
	ka <- sum(divframe$n.subs) / sum(divframe$n.sites)
	ks <- sum(divframe$s.subs) / sum(divframe$s.sites)
	w <- ka/ks
	return(c("ka"=ka, "ks"=ks, "w"= w))
} 


bulk.kaks(divs[1:500,])  # test function on a smaller subset of data, as an example
#use this function on subsets that we care about, zlinked, sperm, etc
#sperm_div<-divs[divs$gene %in% sperm$OGS2,]
library(boot)  # load the boot package


# create the function required as input for bootstrapping the statistics we want
boot.kaks <- function(divdata, indices) {
	divs2 <- divdata[indices, ]
	kaks.out <- bulk.kaks(divs2)
}


# generate a bootstrap object for various subsets of the data (could be Z vs Autosome, or sperm vs genomic background...)
div.boot.A <- boot(data=divs[1:500,], statistic=boot.kaks, R=1000)
div.boot.B <- boot(data=divs[501:1000,], statistic=boot.kaks, R=1000)

div.boot.genome <-boot(data=divs,statistic=boot.kaks, R=1000)
div.boot.autosome <-boot(data=divs[!(divs$gene %in% gili_zlinked$gene.dp),],statistic=boot.kaks, R=1000)
div.boot.zlinked <-boot(data=divs[divs$gene %in% gili_zlinked$gene.dp,],statistic=boot.kaks, R=1000)
div.boot.sperm <-boot(data=divs[divs$gene %in% sperm$OGS2,],statistic=boot.kaks, R=1000)
div.boot.apyrene <-boot(data=divs[divs$gene %in% gili_apyrene$OGS2,],statistic=boot.kaks, R=1000)
div.boot.eupyrene <-boot(data=divs[divs$gene %in% gili_eupyrene$OGS2,],statistic=boot.kaks, R=1000)
div.boot.apy.eupy <-boot(data=divs[divs$gene %in% gili_apy.eupy$OGS2,],statistic=boot.kaks, R=1000)



#
#
# get confidence intervals for each point estimate.  call help on boot.ci to understand what gets returned...
ka.ci.A <- boot.ci(div.boot.A, index = 1, type = "basic")
kw.ci.A <- boot.ci(div.boot.A, index = 2, type = "basic")
w.ci.A  <- boot.ci(div.boot.A, index = 3, type = "basic")

ka.ci.B <- boot.ci(div.boot.B, index = 1, type = "basic")
kw.ci.B <- boot.ci(div.boot.B, index = 2, type = "basic")
w.ci.B  <- boot.ci(div.boot.B, index = 3, type = "basic")

ka.ci.genome <- boot.ci(div.boot.genome, index = 1, type = "basic")
kw.ci.genome <- boot.ci(div.boot.genome, index = 2, type = "basic")
w.ci.genome <- boot.ci(div.boot.genome, index = 3, type = "basic")

ka.ci.autosome <- boot.ci(div.boot.autosome, index = 1, type = "basic")
kw.ci.autosome <- boot.ci(div.boot.autosome, index = 2, type = "basic")
w.ci.autosome <- boot.ci(div.boot.autosome, index = 3, type = "basic")

ka.ci.zlinked <- boot.ci(div.boot.zlinked, index = 1, type = "basic")
kw.ci.zlinked <- boot.ci(div.boot.zlinked, index = 2, type = "basic")
w.ci.zlinked <- boot.ci(div.boot.zlinked, index = 3, type = "basic")

ka.ci.sperm <- boot.ci(div.boot.sperm, index = 1, type = "basic")
kw.ci.sperm <- boot.ci(div.boot.sperm, index = 2, type = "basic")
w.ci.sperm <- boot.ci(div.boot.sperm, index = 3, type = "basic")

ka.ci.apyrene <- boot.ci(div.boot.apyrene, index = 1, type = "basic")
kw.ci.apyrene <- boot.ci(div.boot.apyrene, index = 2, type = "basic")
w.ci.apyrene <- boot.ci(div.boot.apyrene, index = 3, type = "basic")

ka.ci.eupyrene <- boot.ci(div.boot.eupyrene, index = 1, type = "basic")
kw.ci.eupyrene <- boot.ci(div.boot.eupyrene, index = 2, type = "basic")
w.ci.eupyrene <- boot.ci(div.boot.eupyrene, index = 3, type = "basic")

ka.ci.apy.eupy <- boot.ci(div.boot.apy.eupy, index = 1, type = "basic")
kw.ci.apy.eupy <- boot.ci(div.boot.apy.eupy, index = 2, type = "basic")
w.ci.apy.eupy <- boot.ci(div.boot.apy.eupy, index = 3, type = "basic")



boot.ci(div.boot.A, type = "basic")

boot.ci(div.boot.genome, type = "basic")
boot.ci(div.boot.autosome, type = "basic")
boot.ci(div.boot.zlinked, type = "basic")
boot.ci(div.boot.sperm, type = "basic")
boot.ci(div.boot.apyrene, type = "basic")
boot.ci(div.boot.eupyrene, type = "basic")
boot.ci(div.boot.apy.eupy, type = "basic")


#Plot version 1: simple lines for error bars
plot(y=c(w.ci.A$t0, w.ci.B$t0), x = 1:2, xlim=c(0.5, 2.5), ylim = c(0.05,.15), pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka/Ks")
axis(side = 1, labels = c("groupA", "groupB"), at = 1:2)
segments(x0 = 1, y0=w.ci.A$basic[4], y1=w.ci.A$basic[5], col = "blue" )  # simple plotting of errorbars, without end ticks
segments(x0 = 2, y0=w.ci.B$basic[4], y1=w.ci.B$basic[5], col = "blue" )

plot(y=c(w.ci.autosome$t0, w.ci.sperm$t0), x = 1:2, xlim=c(0.5, 2.5), ylim = c(0.05,.15), pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka/Ks")
axis(side = 1, labels = c("autosome", "sperm"), at = 1:2)
segments(x0 = 1, y0=w.ci.autosome$basic[4], y1=w.ci.autosome$basic[5], col = "blue" ) 
segments(x0 = 2, y0=w.ci.sperm$basic[4], y1=w.ci.sperm$basic[5], col = "blue" )


# Plot version 2: using errorbars with end ticks

#first make a function for generating errorbars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) {
		stop("vectors must be same length")
		}
	arrows(x0 = x, y0=lower, y1=upper, angle=90, code=3, length=length, ...)
} 

plot(y=c(w.ci.A$t0, w.ci.B$t0), x = 1:2, xlim=c(0.5, 2.5), ylim = c(0.05,.15), pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka/Ks")
axis(side = 1, labels = c("groupA", "groupB"), at = 1:2)
error.bar(x = 1:2, y =c(w.ci.A$t0, w.ci.B$t0), lower = c(w.ci.A$basic[4], w.ci.B$basic[4]), upper =c(w.ci.A$basic[5], w.ci.B$basic[5]), length = 0.05, col="blue"  )


plot(y=c(w.ci.genome$t0, w.ci.autosome$t0,w.ci.zlinked$t0,w.ci.sperm$t0,w.ci.apyrene$t0,w.ci.eupyrene$t0,w.ci.apy.eupy$t0), x = 1:7,  ylim = c(0.05,.12), pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka/Ks")
axis(side = 1, labels = c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared"), at = 1:7)
error.bar(x = 1:7, y =c(w.ci.genome$t0, w.ci.autosome$t0,w.ci.zlinked$t0,w.ci.sperm$t0,w.ci.apyrene$t0,w.ci.eupyrene$t0,w.ci.apy.eupy$t0), lower = c(w.ci.genome$basic[4], w.ci.autosome$basic[4],w.ci.zlinked$basic[4],w.ci.sperm$basic[4],w.ci.apyrene$basic[4],w.ci.eupyrene$basic[4],w.ci.apy.eupy$basic[4]), upper =c(w.ci.genome$basic[5], w.ci.autosome$basic[5],w.ci.zlinked$basic[5],w.ci.sperm$basic[5],w.ci.apyrene$basic[5],w.ci.eupyrene$basic[5],w.ci.apy.eupy$basic[5]), length = 0.05, col="blue"  )

plot(y=c(ka.ci.genome$t0, ka.ci.autosome$t0,ka.ci.zlinked$t0,ka.ci.sperm$t0,ka.ci.apyrene$t0,ka.ci.eupyrene$t0,ka.ci.apy.eupy$t0), ylim=c(0.01,0.02),x = 1:7, pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka")
axis(side = 1, labels = c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared"), at = 1:7)
error.bar(x = 1:7, y =c(ka.ci.genome$t0, ka.ci.autosome$t0,ka.ci.zlinked$t0,ka.ci.sperm$t0,ka.ci.apyrene$t0,ka.ci.eupyrene$t0,ka.ci.apy.eupy$t0), lower = c(ka.ci.genome$basic[4], ka.ci.autosome$basic[4],ka.ci.zlinked$basic[4],ka.ci.sperm$basic[4],ka.ci.apyrene$basic[4],ka.ci.eupyrene$basic[4],ka.ci.apy.eupy$basic[4]), upper =c(ka.ci.genome$basic[5], ka.ci.autosome$basic[5],ka.ci.zlinked$basic[5],ka.ci.sperm$basic[5],ka.ci.apyrene$basic[5],ka.ci.eupyrene$basic[5],ka.ci.apy.eupy$basic[5]), length = 0.05, col="blue"  )

plot(y=c(kw.ci.genome$t0, kw.ci.autosome$t0,kw.ci.zlinked$t0,kw.ci.sperm$t0,kw.ci.apyrene$t0,kw.ci.eupyrene$t0,kw.ci.apy.eupy$t0), ylim=c(0.15,0.25),x = 1:7, pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ks")
axis(side = 1, labels = c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared"), at = 1:7)
error.bar(x = 1:7, y =c(kw.ci.genome$t0, kw.ci.autosome$t0,kw.ci.zlinked$t0,kw.ci.sperm$t0,kw.ci.apyrene$t0,kw.ci.eupyrene$t0,kw.ci.apy.eupy$t0), lower = c(kw.ci.genome$basic[4], kw.ci.autosome$basic[4],kw.ci.zlinked$basic[4],kw.ci.sperm$basic[4],kw.ci.apyrene$basic[4],kw.ci.eupyrene$basic[4],kw.ci.apy.eupy$basic[4]), upper =c(kw.ci.genome$basic[5], kw.ci.autosome$basic[5],kw.ci.zlinked$basic[5],kw.ci.sperm$basic[5],kw.ci.apyrene$basic[5],kw.ci.eupyrene$basic[5],kw.ci.apy.eupy$basic[5]), length = 0.05, col="blue"  )





#
#
#
#
###################for all subsets, <thing>_TF needs to be copied into kaks_nodup, possibly reneamed to whatever is as <stuff> =TRUE
#sperm_TF <- kaks_nodup$gene %in% sperm$OSG2
#sperm_kaks <- kaks_nodup[kaks_nodup$gene %in% kaks_nodup$sperm,]
#apyrene_TF <- kaks_nodup$gene %in% sperm$Apyrene 
#apyrene_kaks <- subset(kaks_nodup, apyrene=TRUE)
#eupyrene_TF <- kaks$gene %in% sperm$Eupyrene
#eupyrene_kaks <-subset(kaks_nodup, eupyrene=TRUE)
#apy.eupy_TF <- kaks$gene %in% sperm$AE
#apy.eupy_kaks <- subset(kaks_nodup, apy.eupy=TRUE) #
#zlinked_TF <- kaks$gene %in% zlinked$OSG2
#zlinked <- subset(

#converts all columns to numeric
#con[cols.num] <-sapply(con[cols.num],as.numeric)
#boxplot of dN, dS, and dN/dS
#boxplot(con_num$Ka.genome,con_num$Ka.zlinked,con_num$Ka.sperm,con_num$Ka.apyrene,con_num$Ka.eupyrene,con_num$Ka.apy.eupy,main="dN",outline=F, names=c("genome","Z-linked","sperm","apyrene","eupyrene","shared"))
#boxplot(con_num$Ks.genome,con_num$Ks.zlinked,con_num$Ks.sperm,con_num$Ks.apyrene,con_num$Ks.eupyrene,con_num$Ks.apy.eupy,main="dS",outline=F, names=c("genome","Z-linked","sperm","apyrene","eupyrene","shared"))
#boxplot(con_num$KaKs.genome,con_num$KaKs.zlinked,con_num$KaKs.sperm,con_num$KaKs.apyrene,con_num$KaKs.eupyrene,con_num$KaKs.apy.eupy,main="dN/dS",outline=F, names=c("genome","Z-linked","sperm","apyrene","eupyrene","shared"))
