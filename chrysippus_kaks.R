chrys_kaks <- read.csv("/Users/jeffcole/Git_repo/chrys_dnds.csv",header=T, stringsAsFactors = F)

chrys_sperm <-read.csv("/Users/jeffcole/Git_repo/D.plex_sperm_proteome.csv", header= T, stringsAsFactors = F)

chrys_chrom <-read.table("/Users/jeffcole/Git_repo/D.plex_chromosomes.txt", header=T)
chrys_zlinked1 <-subset(chrys_chrom, dpChr==1)
chrys_zlinked16 <-subset(chrys_chrom, dpChr==16)
chrys_zlinked <- rbind(chrys_zlinked1, chrys_zlinked16)

#this removes NA from Ka row
chrys_no_NA <- subset(chrys_kaks, !(is.na(chrys_kaks$Ka)))
#there are still rows from Ks that have NA, this resolves that
chrys_no_NA1<- subset(no_NA, !(is.na(chrys_no_NA$Ks)))
#restore NA to kaks object if needed
chrys_kaks <-chrys_no_NA1

# this sorts gene ID in descending order, for duplicate gene IDs its length will also be listed in descending order 
chrys_kaks_order <- chrys_kaks[order(chrys_kaks$gene, chrys_kaks$length),]
#this chops all duplicated gene ID leaving highest value length
chrys_kaks_nodup <- chrys_kaks_order[!duplicated(chrys_kaks_order$gene),]
chrys_kaks_nodup <- chrys_kaks_nodup[-1,]
chrys_kaks_nodup$Ka <- as.numeric(chrys_kaks_nodup$Ka)
chrys_kaks_nodup$Ks <- as.numeric(chrys_kaks_nodup$Ks)
chrys_kaks_nodup$KaKs <- as.numeric(chrys_kaks_nodup$KaKs)


#a small example to demonstrate duplicate lengths for duplicate gene IDs is not a problem
#d <- data.frame(id=c(4,1,2,3,3,3,4),values=c(1,3,5,8,8,5,9))
#d_order <- d[order(d$id, -d$values),]
#d_nodup <- d_order[!duplicated(d_trimmed$id),]


chrys_autosome <- chrys_kaks_nodup[!(chrys_kaks_nodup$gene %in% chrys_zlinked$gene.dp) ,]
#uncomment if necessary to remove sperm proteins from autosomal data for analysis purposes
#autosome_kaks <- autosome[!(autosome$gene %in% sperm$OGS2) ,]

chrys_sperm_kaks <- chrys_kaks_nodup[chrys_kaks_nodup$gene %in% chrys_sperm$OGS2 ,]
chrys_zlinked_kaks <- chrys_kaks_nodup[chrys_kaks_nodup$gene %in% chrys_zlinked$gene.dp ,]

chrys_apyrene <- subset(chrys_sperm, Apyrene==TRUE)
chrys_apyrene_kaks <- chrys_kaks_nodup[chrys_kaks_nodup$gene %in% chrys_apyrene$OGS2 ,]

chrys_eupyrene <- subset(chrys_sperm, Eupyrene==TRUE)
chrys_eupyrene_kaks <- chrys_kaks_nodup[chrys_kaks_nodup$gene %in% chrys_eupyrene$OGS2 ,]
 
chrys_apy.eupy <- subset(chrys_sperm, AE==TRUE)
chrys_apy.eupy_kaks <- chrys_kaks_nodup[chrys_kaks_nodup$gene %in% chrys_apy.eupy$OGS2 ,]
#sometimes at this point the Ka Ks and KaKs are stored as chr, it may be necessary to convert them to num
#these dataframes may also need to be merged for conveniencem thought that is not necessary

################create dataframe to determine percent sampled
chrys_percent_genome_sampled<-length(unique(chrys_kaks_nodup$gene))/15130
chrys_percent_autosome_sampled<-length(unique(chrys_autosome$gene))/(15130-length(unique(chrys_zlinked$gene.dp)))
chrys_percent_zlinked_sampled<-length(unique(chrys_zlinked_kaks$gene))/length(unique(chrys_zlinked$gene.dp))
chrys_percent_sperm_sampled<-length(unique(chrys_sperm_kaks$gene))/length(unique(chrys_sperm$OGS2))
chrys_percent_apyrene_sampled<-length(unique(chrys_apyrene_kaks$gene))/length(unique(chrys_apyrene$OGS2))
chrys_percent_eupyrene_sampled<-length(unique(chrys_eupyrene_kaks$gene))/length(unique(chrys_eupyrene$OGS2))
chrys_percent_shared_sampled<-length(unique(chrys_apy.eupy_kaks$gene))/length(unique(chrys_apy.eupy$OGS2))

chrys_percent_sampled<-data.frame(chrys_genome=chrys_percent_genome_sampled,chrys_autosome=chrys_percent_autosome_sampled,chrys_zlinked=chrys_percent_zlinked_sampled,chrys_sperm=chrys_percent_sperm_sampled,chrys_apyrene=chrys_percent_apyrene_sampled,chrys_eupyrene=chrys_percent_eupyrene_sampled,chrys_shared=chrys_percent_shared_sampled)


#from here plots and analysis can begin
#dN_boxplot <- boxplot(kaks_nodup$Ka, autosome$Ka,zlinked_kaks$Ka, sperm_kaks$Ka, apyrene_kaks$Ka, eupyrene_kaks$Ka, apy.eupy_kaks$Ka, main="dN", names=c("genome","autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
#dS_boxplot <- boxplot(kaks_nodup$Ks, autosome$Ks,zlinked_kaks$Ks, sperm_kaks$Ks, apyrene_kaks$Ks, eupyrene_kaks$Ks, apy.eupy_kaks$Ks, main="dS", names=c("genome","autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
#dNdS_boxplot <- boxplot(kaks_nodup$KaKs, autosome$KaKs,zlinked_kaks$KaKs, sperm_kaks$KaKs, apyrene_kaks$KaKs, eupyrene_kaks$KaKs, apy.eupy_kaks$KaKs, main="dN/dS", names=c("genome","autosome","Z-linked","sperm","apyrene","eupyrene","shared"))


#put into one dataframe
chrys_dnds_autosome <- data.frame(x=1:length(chrys_autosome$Ka), chrys_autosome_Ka=chrys_autosome$Ka,chrys_autosome_Ks=chrys_autosome$Ks, chrys_autosome_KaKs=chrys_autosome$KaKs)
chrys_dnds_zlinked <- data.frame(x=1:length(chrys_zlinked_kaks$Ka), chrys_zlinked_Ka=chrys_zlinked_kaks$Ka,chrys_zlinked_Ks=chrys_zlinked_kaks$Ks, chrys_zlinked_KaKs=chrys_zlinked_kaks$KaKs)
chrys_dnds_sperm <- data.frame(x=1:length(chrys_sperm_kaks$Ka), chrys_sperm_Ka=chrys_sperm_kaks$Ka,chrys_sperm_Ks=chrys_sperm_kaks$Ks, chrys_sperm_KaKs=chrys_sperm_kaks$KaKs)
chrys_dnds_apyrene <- data.frame(x=1:length(chrys_apyrene_kaks$Ka), chrys_apyrene_Ka=chrys_apyrene_kaks$Ka,chrys_apyrene_Ks=chrys_apyrene_kaks$Ks, chrys_apyrene_KaKs=chrys_apyrene_kaks$KaKs)
chrys_dnds_eupyrene <- data.frame(x=1:length(chrys_eupyrene_kaks$Ka), chrys_eupyrene_Ka=chrys_eupyrene_kaks$Ka,chrys_eupyrene_Ks=chrys_eupyrene_kaks$Ks, chrys_eupyrene_KaKs=chrys_eupyrene_kaks$KaKs)
chrys_dnds_apy.eupy <- data.frame(x=1:length(chrys_apy.eupy_kaks$Ka), chrys_apy.eupy_Ka=chrys_apy.eupy_kaks$Ka,chrys_apy.eupy_Ks=chrys_apy.eupy_kaks$Ks, chrys_apy.eupy_KaKs=chrys_apy.eupy_kaks$KaKs)
chrys_dnds_list <- list(chrys_dnds_autosome, chrys_dnds_zlinked, chrys_dnds_sperm, chrys_dnds_apyrene, chrys_dnds_eupyrene, chrys_dnds_apy.eupy)
chrys_dnds <- Reduce(function(...) merge(..., all=T), chrys_dnds_list)

#filter values>2
chrys_dnds[chrys_dnds>2]=NA




##############THINGS I NEED TO DO#################
######1. using final dnds dataframe#######
#create boxplot for Ka, Ks, and KaKs

boxplot(chrys_dnds$chrys_autosome_Ka, chrys_dnds$chrys_zlinked_Ka,chrys_dnds$chrys_sperm_Ka, chrys_dnds$chrys_apyrene_Ka, chrys_dnds$chrys_eupyrene_Ka, chrys_dnds$chrys_apy.eupy_Ka, outline=F, main="dN", names=c("autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
boxplot(chrys_dnds$chrys_autosome_Ks, chrys_dnds$chrys_zlinked_Ks,chrys_dnds$chrys_sperm_Ks, chrys_dnds$chrys_apyrene_Ks, chrys_dnds$chrys_eupyrene_Ks, chrys_dnds$chrys_apy.eupy_Ks, outline=F, main="dS", names=c("autosome","Z-linked","sperm","apyrene","eupyrene","shared"))
boxplot(chrys_dnds$chrys_autosome_KaKs, chrys_dnds$chrys_zlinked_Ka,chrys_dnds$chrys_sperm_KaKs, chrys_dnds$chrys_apyrene_KaKs, chrys_dnds$chrys_eupyrene_KaKs, chrys_dnds$chrys_apy.eupy_KaKs, outline=F, main="dNdS", names=c("autosome","Z-linked","sperm","apyrene","eupyrene","shared"))

#create histograms for all Ka, Ks and KaKs
hist(chrys_dnds$chrys_autosome_Ka)
hist(chrys_dnds$chrys_autosome_Ks)
hist(chrys_dnds$chrys_autosome_KaKs)
hist(chrys_dnds$chrys_zlinked_Ka)
hist(chrys_dnds$chrys_zlinked_Ks)
hist(chrys_dnds$chrys_zlinked_KaKs)
hist(chrys_dnds$chrys_sperm_Ka)
hist(chrys_dnds$chrys_sperm_Ks)
hist(chrys_dnds$chrys_sperm_KaKs)
hist(chrys_dnds$chrys_apyrene_Ka)
hist(chrys_dnds$chrys_apyrene_Ks)
hist(chrys_dnds$chrys_apyrene_KaKs)
hist(chrys_dnds$chrys_eupyrene_Ka)
hist(chrys_dnds$chrys_eupyrene_Ks)
hist(chrys_dnds$chrys_eupyrene_KaKs)
hist(chrys_dnds$chrys_apy.eupy_Ka)
hist(chrys_dnds$chrys_apy.eupy_Ks)
hist(chrys_dnds$chrys_apy.eupy_KaKs)



#run wilcox test for sperm Ka, Ks, and KaKs vs autosome, 
#also findout if sperm is all sperm or just leftovers, may need new column
wilcox.test(chrys_dnds$chrys_autosome_Ka, chrys_dnds$chrys_sperm_Ka)
wilcox.test(chrys_dnds$chrys_autosome_Ks, chrys_dnds$chrys_sperm_Ks)
wilcox.test(chrys_dnds$chrys_autosome_KaKs, chrys_dnds$chrys_sperm_KaKs)


#run kWallis for apyrene, eupyrene, shared with Ka, Ks, and KaKs (maybe extend to sperm general and zlinked)
#try this, if not it may need a separate column denoting categorical factor
#kruskal.test(dnds$apyrene_Ka, dnds$eupyrene_Ka, dnds$apy.eupy_Ka)
#kruskal.test(dnds$apyrene_Ks, dnds$eupyrene_Ks, dnds$apy.eupy_Ks)
#kruskal.test(dnds$apyrene_KaKs, dnds$eupyrene_KaKs, dnds$apy.eupy_KaKs)

#dataframe for kwallis


#Ka_apy <- data.frame(ka=dnds$apyrene_Ka, type=rep("apyrene", length(dnds$apyrene_Ka))
#Ka_eupy <- data.frame(ka=dnds$eupyrene_Ka, type=rep("eupyrene", length(dnds$apyrene_Ka))

chrys_Ka_2col<- data.frame(Ka=c(chrys_dnds$chrys_apyrene_Ka,chrys_dnds$chrys_eupyrene_Ka,chrys_dnds$chrys_apy.eupy_Ka),type= c(rep("apyrene" , length(chrys_dnds$chrys_apyrene_Ka)),rep("eupyrene",length(chrys_dnds$chrys_eupyrene_Ka)),rep("apy.eupy",length(chrys_dnds$chrys_apyrene_Ka))))
kruskal.test(chrys_Ka_2col$Ka~chrys_Ka_2col$type)

chrys_Ks_2col<- data.frame(Ks=c(chrys_dnds$chrys_apyrene_Ks,chrys_dnds$chrys_eupyrene_Ks,chrys_dnds$chrys_apy.eupy_Ks),type= c(rep("apyrene" , length(chrys_dnds$chrys_apyrene_Ks)),rep("eupyrene",length(chrys_dnds$chrys_eupyrene_Ks)),rep("apy.eupy",length(chrys_dnds$chrys_apyrene_Ks))))
kruskal.test(chrys_Ks_2col$Ks~chrys_Ks_2col$type)

chrys_KaKs_2col<- data.frame(KaKs=c(chrys_dnds$chrys_apyrene_KaKs,chrys_dnds$chrys_eupyrene_KaKs,chrys_dnds$chrys_apy.eupy_KaKs),type= c(rep("apyrene" , length(chrys_dnds$chrys_apyrene_KaKs)),rep("eupyrene",length(chrys_dnds$chrys_eupyrene_KaKs)),rep("apy.eupy",length(chrys_dnds$chrys_apyrene_KaKs))))
kruskal.test(chrys_KaKs_2col$KaKs~chrys_KaKs_2col$type)



################5. BOOTSTRAP, insert here because many of the subsets and files loaded make sense for it to be here
chrys_divs <- chrys_kaks[,c("gene","s.sites", "n.sites", "s.substitions", "n.substitutions" )] # extract relevant columns, just to make life easier

names(chrys_divs) <- c("gene","s.sites", "n.sites", "s.subs", "n.subs") # rename columns, just to make life easier.  also note typo: s.substitions
apply(chrys_divs, 2, function(x){sum(is.na(x))})  # count how many NAs are in each column, just to know what the data look like

# convert NA to 0 for substitution counts.
chrys_divs$s.subs[is.na(chrys_divs$s.subs)] <- 0
chrys_divs$n.subs[is.na(chrys_divs$n.subs)] <- 0
#force site counts to numeric.  I don't know why this is reading in as character, but we need to force it to numeric.
chrys_divs[[2]] <- as.numeric(chrys_divs[[2]])
chrys_divs[[3]] <- as.numeric(chrys_divs[[3]])
chrys_divs <- head(chrys_divs,-3) #takes care of random 3 lines of NA
chrys_divs$s.sites[is.na(chrys_divs$s.sites)] <- 0
chrys_divs$n.sites[is.na(chrys_divs$n.sites)] <- 0

chrys_bulk.kaks <- function(divframe) {  # function to get ka, ks, and ka/ks (ie w, or "omega") from the "divs" dataframe, or some subset of it. Function returns a vector with names.
	chrys_ka <- sum(divframe$n.subs) / sum(divframe$n.sites)
	chrys_ks <- sum(divframe$s.subs) / sum(divframe$s.sites)
	chrys_w <- chrys_ka/chrys_ks
	return(c("ka"=chrys_ka, "ks"=chrys_ks, "w"= chrys_w))
} 


chrys_bulk.kaks(divs[1:500,])  # test function on a smaller subset of data, as an example
#use this function on subsets that we care about, zlinked, sperm, etc
#sperm_div<-divs[divs$gene %in% sperm$OGS2,]
library(boot)  # load the boot package


# create the function required as input for bootstrapping the statistics we want
chrys_boot.kaks <- function(divdata, indices) {
	chrys_divs2 <- divdata[indices, ]
	chrys_kaks.out <- chrys_bulk.kaks(chrys_divs2)
}


# generate a bootstrap object for various subsets of the data (could be Z vs Autosome, or sperm vs genomic background...)
div.boot.A <- boot(data=divs[1:500,], statistic=boot.kaks, R=1000)
div.boot.B <- boot(data=divs[501:1000,], statistic=boot.kaks, R=1000)

chrys_div.boot.genome <-boot(data=chrys_divs,statistic=chrys_boot.kaks, R=1000)
chrys_div.boot.autosome <-boot(data=chrys_divs[!(chrys_divs$gene %in% chrys_zlinked$gene.dp),],statistic=chrys_boot.kaks, R=1000)
chrys_div.boot.zlinked <-boot(data=chrys_divs[chrys_divs$gene %in% chrys_zlinked$gene.dp,],statistic=chrys_boot.kaks, R=1000)
chrys_div.boot.sperm <-boot(data=chrys_divs[chrys_divs$gene %in% chrys_sperm$OGS2,],statistic=chrys_boot.kaks, R=1000)
chrys_div.boot.apyrene <-boot(data=chrys_divs[chrys_divs$gene %in% chrys_apyrene$OGS2,],statistic=chrys_boot.kaks, R=1000)
chrys_div.boot.eupyrene <-boot(data=chrys_divs[chrys_divs$gene %in% chrys_eupyrene$OGS2,],statistic=chrys_boot.kaks, R=1000)
chrys_div.boot.apy.eupy <-boot(data=chrys_divs[chrys_divs$gene %in% chrys_apy.eupy$OGS2,],statistic=chrys_boot.kaks, R=1000)



#
#
# get confidence intervals for each point estimate.  call help on boot.ci to understand what gets returned...
ka.ci.A <- boot.ci(div.boot.A, index = 1, type = "basic")
kw.ci.A <- boot.ci(div.boot.A, index = 2, type = "basic")
w.ci.A  <- boot.ci(div.boot.A, index = 3, type = "basic")

ka.ci.B <- boot.ci(div.boot.B, index = 1, type = "basic")
kw.ci.B <- boot.ci(div.boot.B, index = 2, type = "basic")
w.ci.B  <- boot.ci(div.boot.B, index = 3, type = "basic")

chrys_ka.ci.genome <- boot.ci(chrys_div.boot.genome, index = 1, type = "basic")
chrys_kw.ci.genome <- boot.ci(chrys_div.boot.genome, index = 2, type = "basic")
chrys_w.ci.genome <- boot.ci(chrys_div.boot.genome, index = 3, type = "basic")

chrys_ka.ci.autosome <- boot.ci(chrys_div.boot.autosome, index = 1, type = "basic")
chrys_kw.ci.autosome <- boot.ci(chrys_div.boot.autosome, index = 2, type = "basic")
chrys_w.ci.autosome <- boot.ci(chrys_div.boot.autosome, index = 3, type = "basic")

chrys_ka.ci.zlinked <- boot.ci(chrys_div.boot.zlinked, index = 1, type = "basic")
chrys_kw.ci.zlinked <- boot.ci(chrys_div.boot.zlinked, index = 2, type = "basic")
chrys_w.ci.zlinked <- boot.ci(chrys_div.boot.zlinked, index = 3, type = "basic")

chrys_ka.ci.sperm <- boot.ci(chrys_div.boot.sperm, index = 1, type = "basic")
chrys_kw.ci.sperm <- boot.ci(chrys_div.boot.sperm, index = 2, type = "basic")
chrys_w.ci.sperm <- boot.ci(chrys_div.boot.sperm, index = 3, type = "basic")

chrys_ka.ci.apyrene <- boot.ci(chrys_div.boot.apyrene, index = 1, type = "basic")
chrys_kw.ci.apyrene <- boot.ci(chrys_div.boot.apyrene, index = 2, type = "basic")
chrys_w.ci.apyrene <- boot.ci(chrys_div.boot.apyrene, index = 3, type = "basic")

chrys_ka.ci.eupyrene <- boot.ci(chrys_div.boot.eupyrene, index = 1, type = "basic")
chrys_kw.ci.eupyrene <- boot.ci(chrys_div.boot.eupyrene, index = 2, type = "basic")
chrys_w.ci.eupyrene <- boot.ci(chrys_div.boot.eupyrene, index = 3, type = "basic")

chrys_ka.ci.apy.eupy <- boot.ci(chrys_div.boot.apy.eupy, index = 1, type = "basic")
chrys_kw.ci.apy.eupy <- boot.ci(chrys_div.boot.apy.eupy, index = 2, type = "basic")
chrys_w.ci.apy.eupy <- boot.ci(chrys_div.boot.apy.eupy, index = 3, type = "basic")



boot.ci(div.boot.A, type = "basic")

boot.ci(chrys_div.boot.genome, type = "basic")
boot.ci(chrys_div.boot.autosome, type = "basic")
boot.ci(chrys_div.boot.zlinked, type = "basic")
boot.ci(chrys_div.boot.sperm, type = "basic")
boot.ci(chrys_div.boot.apyrene, type = "basic")
boot.ci(chrys_div.boot.eupyrene, type = "basic")
boot.ci(chrys_div.boot.apy.eupy, type = "basic")


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
#remove ylim
plot(y=c(w.ci.A$t0, w.ci.B$t0), x = 1:2, xlim=c(0.5, 2.5), ylim = c(0.05,.15), pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka/Ks")
axis(side = 1, labels = c("groupA", "groupB"), at = 1:2)
error.bar(x = 1:2, y =c(w.ci.A$t0, w.ci.B$t0), lower = c(w.ci.A$basic[4], w.ci.B$basic[4]), upper =c(w.ci.A$basic[5], w.ci.B$basic[5]), length = 0.05, col="blue"  )


#plot boot bit by bit 
plot(y=(chrys_w.ci.genome$t0) #should be 0.09094763




labels <- c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared")
par(las=2)
plot(y=c(chrys_w.ci.genome$t0, chrys_w.ci.autosome$t0,chrys_w.ci.zlinked$t0,chrys_w.ci.sperm$t0,chrys_w.ci.apyrene$t0,chrys_w.ci.eupyrene$t0,chrys_w.ci.apy.eupy$t0), x = 1:7,  pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "",ylim = c(0.06,.11), ylab = "Ka/Ks")
axis(side = 1, labels = c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared"), at = 1:7)
error.bar(x = 1:7, y =c(chrys_w.ci.genome$t0, chrys_w.ci.autosome$t0,chrys_w.ci.zlinked$t0,chrys_w.ci.sperm$t0,chrys_w.ci.apyrene$t0,chrys_w.ci.eupyrene$t0,chrys_w.ci.apy.eupy$t0), lower = c(chrys_w.ci.genome$basic[4], chrys_w.ci.autosome$basic[4],chrys_w.ci.zlinked$basic[4],chrys_w.ci.sperm$basic[4],chrys_w.ci.apyrene$basic[4],chrys_w.ci.eupyrene$basic[4],chrys_w.ci.apy.eupy$basic[4]), upper =c(chrys_w.ci.genome$basic[5], chrys_w.ci.autosome$basic[5],chrys_w.ci.zlinked$basic[5],chrys_w.ci.sperm$basic[5],chrys_w.ci.apyrene$basic[5],chrys_w.ci.eupyrene$basic[5],chrys_w.ci.apy.eupy$basic[5]), length = 0.05, col="blue"  )

plot(y=c(chrys_ka.ci.genome$t0, chrys_ka.ci.autosome$t0,chrys_ka.ci.zlinked$t0,chrys_ka.ci.sperm$t0,chrys_ka.ci.apyrene$t0,chrys_ka.ci.eupyrene$t0,chrys_ka.ci.apy.eupy$t0),ylim=c(0.0135,0.02), x = 1:7, pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ka")
axis(side = 1, labels = c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared"), at = 1:7)
error.bar(x = 1:7, y =c(chrys_ka.ci.genome$t0, chrys_ka.ci.autosome$t0,chrys_ka.ci.zlinked$t0,chrys_ka.ci.sperm$t0,chrys_ka.ci.apyrene$t0,chrys_ka.ci.eupyrene$t0,chrys_ka.ci.apy.eupy$t0), lower = c(chrys_ka.ci.genome$basic[4], chrys_ka.ci.autosome$basic[4],chrys_ka.ci.zlinked$basic[4],chrys_ka.ci.sperm$basic[4],chrys_ka.ci.apyrene$basic[4],chrys_ka.ci.eupyrene$basic[4],chrys_ka.ci.apy.eupy$basic[4]), upper =c(chrys_ka.ci.genome$basic[5], chrys_ka.ci.autosome$basic[5],chrys_ka.ci.zlinked$basic[5],chrys_ka.ci.sperm$basic[5],chrys_ka.ci.apyrene$basic[5],chrys_ka.ci.eupyrene$basic[5],chrys_ka.ci.apy.eupy$basic[5]), length = 0.05, col="blue"  )

plot(y=c(chrys_kw.ci.genome$t0, chrys_kw.ci.autosome$t0,chrys_kw.ci.zlinked$t0,chrys_kw.ci.sperm$t0,chrys_kw.ci.apyrene$t0,chrys_kw.ci.eupyrene$t0,chrys_kw.ci.apy.eupy$t0),ylim=c(0.165,0.215), x = 1:7, pch = 19, col = ("blue"), cex = 1.5, xaxt="n", xlab = "", ylab = "Ks")
axis(side = 1, labels = c("genome", "autosome","zlinked","sperm","apyrene","eupyrene","shared"), at = 1:7)
error.bar(x = 1:7, y =c(chrys_kw.ci.genome$t0, chrys_kw.ci.autosome$t0,chrys_kw.ci.zlinked$t0,chrys_kw.ci.sperm$t0,chrys_kw.ci.apyrene$t0,chrys_kw.ci.eupyrene$t0,chrys_kw.ci.apy.eupy$t0), lower = c(chrys_kw.ci.genome$basic[4], chrys_kw.ci.autosome$basic[4],chrys_kw.ci.zlinked$basic[4],chrys_kw.ci.sperm$basic[4],chrys_kw.ci.apyrene$basic[4],chrys_kw.ci.eupyrene$basic[4],chrys_kw.ci.apy.eupy$basic[4]), upper =c(chrys_kw.ci.genome$basic[5], chrys_kw.ci.autosome$basic[5],chrys_kw.ci.zlinked$basic[5],chrys_kw.ci.sperm$basic[5],chrys_kw.ci.apyrene$basic[5],chrys_kw.ci.eupyrene$basic[5],chrys_kw.ci.apy.eupy$basic[5]), length = 0.05, col="blue"  )









########################################################





######2. compare outlier removed from outlier remained
#determine how many outliers were removed from each column, if > 100, no bueno
#also maybe investigate whether or not trimming the duplicates affected the significance by performing tests on untrimmed data

######3. determine percentage of genes sampled########
#there are ~15,000 genes in D.plex, how many of those genes did we get ratios for
#there are ~800 genes in the sperm proteome, how many of those did we get ratios for
#try repeating this for all subsets, find proportion of genes with measurements to possible genes







#outlierReplace <- function(dataframe, cols,rows,newValue = NA){
#	if(any(rows)){
#		set(dataframe,rows,cols,newValue)
#		}
#	}
#outlierReplace(dnds,"dnds$autosome_Ka", which(dnds$autosome>2),NA)





#for all subsets, <thing>_TF needs to be copied into kaks_nodup, possibly reneamed to whatever is as <stuff> =TRUE
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
