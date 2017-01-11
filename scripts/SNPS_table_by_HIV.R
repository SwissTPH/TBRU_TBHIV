####R version: 3.2.2 (2015-08-14) -- "Fire Safety"####
####Copyright (C) 2015 The R Foundation for Statistical Computing###
####Platform: x86_64-pc-linux-gnu (64-bit)###
####Daniela Brites###
####Date:1.12.2016####
####content:This scripts extracts SFS of Mtb/HIV+ and Mtb/HIV- mutations and performs chi-square tests per Tuberculist categories.

#read snps table, this table contains all position after filtering IS elements, transposases,phages and PPE/PE-PGRS genes
table_snps <- read.table("~/data/table.HIV_UGII.genotypes",na.strings="",stringsAsFactors=F,sep="\t",header=T)

dim(table_snps)

#read annotation
annotation <- read.table("~/data/HIV_UgandaII.annovar_in.finalout_updated",na.strings="",stringsAsFactors=F,sep="\t",header=F,fill=T)
dim(annotation)

colnames(annotation) <- c("position","ref","mut","NS/S/I","Rv","in_vitro","macrophage","tuberculist_cat","size_aa","DR")
table_snps_annotation <-cbind(table_snps,annotation)

##remove drug resistance positions from the table
DR_mutations <- table_snps_annotation$DR=="DR_position"
to.select_DR <- which(DR_mutations)
table_snps_annotation_filtered<- table_snps_annotation[-c(to.select_DR),]

dim(table_snps_annotation_filtered)

#remove ESX family recently duplicated genes (J (Rv1038c), W (Rv3620c), K (Rv1197), P (Rv2347c), M (Rv1792), I (Rv1037c), V (Rv3619c), N (Rv1793), L (Rv1198), O (Rv2346c), G (Rv0287), S (Rv3020c), H (Rv0288), R (Rv3019c), Q (Rv3017c) )

esx_duplicated <- c("Rv1038c", "Rv3620c", "Rv1197", "Rv2347c","Rv1792","Rv1037c", "Rv3619c", "Rv1793", "Rv1198", "Rv2346c", "Rv0287", "Rv3020c", "Rv0288","Rv3019c", "Rv3017c")
no_esx_duplicated <- apply(table_snps_annotation_filtered, 1, function(x) !any(x %in% esx_duplicated ))
table_snps_annotation_filtered <- table_snps_annotation_filtered[c(no_esx_duplicated),]
dim(table_snps_annotation_filtered)

##the position of rrs have been removed alredy 
DR_genes <- c("Rv2428", "Rv1483", "Rv1484", "Rv1908c","Rv0667","Rv3795","Rv2043c", "Rv1630", "Rv0006", "Rv0682","Rv3919c","Rv1694","Rv2416c")
no_DR_genes <- apply(table_snps_annotation_filtered, 1, function(x) !any(x %in% DR_genes))
table_snps_annotation_filtered <- table_snps_annotation_filtered[c(no_DR_genes),]
dim(table_snps_annotation_filtered)

genes_repetions <- c("Rv1572c","Rv1574","Rv1575","Rv0336","Rv0515","IG1195_Rv1174c-Rv1175c","IG127_Rv0126-Rv0127","IG1711_Rv1682-Rv1683","IG18_Rv0018c-Rv0019c","IG3012_Rv2965c-Rv2966c","IG3013_Rv2966c-Rv2967c","IG533_Rv0525-Rv0526","IG559_Rv0551c-Rv0552","IG622_Rv0612-Rv0613c","IG71_Rv0071-Rv0072","IG784_Rv0769-Rv0770","IG877_Rv0861c-Rv0862c","Rv0031","Rv0094c","Rv0095c","Rv0096","Rv0109","Rv0124","Rv0151c","Rv0152c","Rv0159c","Rv0160c","Rv0256c","Rv0257","Rv0277c","Rv0278c","Rv0279c","Rv0280","Rv0285","Rv0286","Rv0297","Rv0304c","Rv0305c","Rv0335c","Rv0353","Rv0354c","Rv0355c","Rv0387c","Rv0388c","Rv0393","Rv0397","Rv0442c","Rv0453","Rv0487","Rv0490","Rv0532","Rv0538","Rv0578c","Rv0605","Rv0605","Rv0740","Rv0741","Rv0742","Rv0746","Rv0747","Rv0750","Rv0754","Rv0755A","Rv0755c","Rv0795","Rv0796","Rv0797","Rv0814c","Rv0823c","Rv0829","Rv0832","Rv0833","Rv0834c","Rv0850","Rv0867c","Rv0872c","Rv0878c","Rv0915c","Rv0916c","Rv0920c","Rv0921","Rv0922","Rv0977","Rv0978c","Rv0980c","Rv1034c","Rv1035c","Rv1036c","Rv1037c","Rv1038c","Rv1039c","Rv1040c","Rv1041c","Rv1042c","Rv1047","Rv1067c","Rv1068c","Rv1087","Rv1088","Rv1089","Rv1091","Rv1128c","Rv1135c","Rv1148c","Rv1149","Rv1150","Rv1168c","Rv1169c","Rv1172c","Rv1173","Rv1195","Rv1196","Rv1197","Rv1198","Rv1199c","Rv1214c","Rv1243c","Rv1288","Rv1295","Rv1313c","Rv1318c","Rv1319c","Rv1325c","Rv1361c","Rv1369c","Rv1370c","Rv1386","Rv1387","Rv1396c","Rv1430","Rv1441c","Rv1450c","Rv1452c","Rv1458c","Rv1468c","Rv1489A","Rv1493","Rv1548c","Rv1557","Rv1558","Rv1573","Rv1574","Rv1575","Rv1576c","Rv1577c","Rv1578c","Rv1579c","Rv1580c","Rv1581c","Rv1582c","Rv1583c","Rv1584c","Rv1585c","Rv1586c","Rv1587c","Rv1588c","Rv1646","Rv1651c","Rv1702c","Rv1705c","Rv1706c","Rv1753c","Rv1756c","Rv1757","Rv1758","Rv1759c","Rv1763","Rv1764","Rv1765A","Rv1765c","Rv1768","Rv1787","Rv1788","Rv1789","Rv1790","Rv1791","Rv1793","Rv1800","Rv1801","Rv1802","Rv1803c","Rv1806","Rv1807","Rv1808","Rv1809","Rv1818c","Rv1829","Rv1840c","Rv1910c","Rv1911c","Rv1917c","Rv1918c","Rv1945","Rv1983","Rv2013","Rv2014","Rv2015c","Rv2048c","Rv2082","Rv2085","Rv2090","Rv2105","Rv2106","Rv2107","Rv2108","Rv2112c","Rv2123","Rv2126c","Rv2162c","Rv2167c","Rv2168c","Rv2177c","Rv2196","Rv2258c","Rv2277c","Rv2278","Rv2279","Rv2328","Rv2340c","Rv2346c","Rv2347c","Rv2352c","Rv2353c","Rv2354","Rv2355","Rv2356c","Rv2371","Rv2396","Rv2408","Rv2424c","Rv2430c","Rv2431c","Rv2460c","Rv2461c","Rv2479c","Rv2480c","Rv2487c","Rv2489","Rv2490c","Rv2512c","Rv2519","Rv2543","Rv2544","Rv2591","Rv2608","Rv2615c","Rv2634c","Rv2648","Rv2649","Rv2650c","Rv2651c","Rv2652c","Rv2653c","Rv2654c","Rv2655c","Rv2656c","Rv2657c","Rv2659c","Rv2665","Rv2666","Rv2673","Rv2680","Rv2689c","Rv2690c","Rv2741","Rv2768c","Rv2769c","Rv2770c","Rv2774c","Rv2791c","Rv2792c","Rv2805","Rv2807","Rv2810c","Rv2812","Rv2814c","Rv2815c","Rv2825c","Rv2828c","Rv2853","Rv2859c","Rv2882c","Rv2885c","Rv2886c","Rv2892c","Rv2931","Rv2932","Rv2943","Rv2943A","Rv2944","Rv2961","Rv2977c","Rv2978c","Rv2979c","Rv2980","Rv3018A","Rv3018c","Rv3021c","Rv3022A","Rv3022c","Rv3023c","Rv3097c","Rv3115","Rv3125c","Rv3135","Rv3136","Rv3144c","Rv3159c","Rv3184","Rv3185","Rv3186","Rv3187","Rv3191c","Rv3281","Rv3325","Rv3326","Rv3327","Rv3343c","Rv3344c","Rv3345c","Rv3346c","Rv3347c","Rv3348","Rv3349c","Rv3350c","Rv3355c","Rv3367","Rv3380c","Rv3381c","Rv3386","Rv3387","Rv3388","Rv3424c","Rv3425","Rv3426","Rv3427c","Rv3428c","Rv3429","Rv3430c","Rv3431c","Rv3466","Rv3467","Rv3474","Rv3475","Rv3477","Rv3478","Rv3507","Rv3508","Rv3511","Rv3512","Rv3513c","Rv3514","Rv3515c","Rv3532","Rv3533c","Rv3539","Rv3558","Rv3590c","Rv3595c","Rv3611","Rv3619c","Rv3620c","Rv3621c","Rv3622c","Rv3636","Rv3637","Rv3638","Rv3639c","Rv3640c","Rv3650","Rv3680","Rv3710","Rv3738c","Rv3739c","Rv3746c","Rv3798","Rv3812","Rv3826","Rv3827c","Rv3828c","Rv3844","Rv3873","Rv3876","Rv3892c","Rv3893c")
no_genes_repetions <- apply(table_snps_annotation_filtered, 1, function(x) !any(x %in% genes_repetions))
table_snps_annotation_filtered <- table_snps_annotation_filtered[ c(no_genes_repetions),]
dim(table_snps_annotation_filtered)

hiv_status <- read.table("~/data/list_name_status",na.strings="",stringsAsFactors=F,sep="\t",header=F,fill=T)

###calculations of SFS

#d_hiv contains all strains (rows) and positons (5312) plus the two last collums are hiv status and year

dim(table_snps_annotation_filtered)

d <- table_snps_annotation_filtered[,1:180]
dim(d)

d_hiv <- cbind(as.data.frame(t(d), stringsAsFactors = F),hiv_status[,2:3])
dim(d_hiv)
colnames(d_hiv)[5313] <- "hiv"
colnames(d_hiv)[5314] <- "year"

annotation_filtered <-  as.data.frame(t(table_snps_annotation_filtered[,181:190]))
dim(annotation_filtered)

#define mutation (derived allele) 
mut <- annotation_filtered[3,]
hiv <- d_hiv$hiv
d <- d_hiv[,-5313:-5314]

##define function to detect invariable sites. In this way invariable sites will be those where the mutation is fixed in all strains, i.e if there is an intermination or gap, that site will not be consired as fixed anymore.

is.not.variable <- function(x) length(table(x[x %in% c("A","G","C","T","N")]))==1

not.variable <- apply(d,2,is.not.variable)
table(not.variable)

to.remove0 <- which(not.variable)
d <- d[, -c(to.remove0)]
dim(d)
mut <-  mut[,-c(to.remove0)]

annotation_filtered <- annotation_filtered[, -c(to.remove0)]

dim(annotation_filtered)
table(unlist(d))

d[d=="K"] <- NA
d[d=="M"] <- NA
d[d=="R"] <- NA
d[d=="S"] <- NA
d[d=="W"] <- NA
d[d=="Y"] <- NA

############## dataset prepared ###########

#create a dataset the same size as d but with the information of mut.
dmut <- mut[rep(1,dim(d)[1]),]
dim(dmut)
ismut <- (d == dmut)
mutg <- colSums(ismut, na.rm=T)
table(mutg, useNA="always")

table(mutg==1)

#define factor singleton and non-singleton 
glob.singl <- factor(mutg==1, labels = c("non-singleton","singleton"))
table(glob.singl, useNA="always")

##counts for HIV+ and HIV-
ismuthiv <- (d[hiv==1,] == dmut[hiv==1,])
mut1 <- colSums(ismuthiv, na.rm=T)
table(mut1, useNA="always")
table(mut1==1)

hiv.factor <- cut(mut1, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
table(hiv.factor, useNA="ifany")

table(hiv.factor[glob.singl=="singleton"], useNA="always")
table(hiv.factor[glob.singl=="non-singleton"], useNA="always")

ismutneg <- (d[hiv==0,] == dmut[hiv==0,])
mut0 <- colSums(ismutneg, na.rm=T)

neg.factor <- cut(mut0, breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
table(neg.factor, useNA="always")
table(neg.factor[glob.singl=="singleton"], useNA="always")

#mutation which are present in MTb/HIV+ but not in Mtb/HIV-
d_hiv_specific <- d[neg.factor=="no-mut"]
dmut_hiv_specific <-dmut[neg.factor=="no-mut"]
ismutHIV_specific <-(d_hiv_specific[hiv==1,]==dmut_hiv_specific[hiv==1,])
mutHIV_specific <-colSums(ismutHIV_specific, na.rm=T)
annotation_HIV <- annotation_filtered[neg.factor=="no-mut"]

##split Mtb/HIV+ specific mutations in singletons and nonsingletons
hiv.specific.factor<-cut(mutHIV_specific,breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
table(hiv.specific.factor,useNA = "always")
table(hiv.specific.factor, t(annotation_HIV[4,]))
table(hiv.specific.factor,t(annotation_HIV[8,]))


#mutation which are present in MTb/HIV- but not in Mtb/HIV+
d_hivneg_specific <- d[hiv.factor=="no-mut"]
dmut_hivneg_specific <-dmut[hiv.factor=="no-mut"]
ismutHIVneg_specific <-(d_hivneg_specific[hiv==0,]==dmut_hivneg_specific[hiv==0,])
mutHIVneg_specific <-colSums(ismutHIVneg_specific, na.rm=T)
annotation_HIVneg <-  annotation_filtered[hiv.factor=="no-mut"]

##split Mtb/HIV+ specific mutations in singletons and nonsingletons
hivneg.specific.factor<-cut(mutHIVneg_specific,breaks=c(0,1,2,10000), labels=c("no-mut","singel","non-sing"), right = F)
table(hivneg.specific.factor,t(annotation_HIVneg[4,]))
table(hivneg.specific.factor,t(annotation_HIVneg[8,]))

##gene enrichment for the specific mutations
##cell wall
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,1],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,1]))

##conserved hypotheticals
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,2],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,2]))

#information pathways
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,3],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,3]))

#intermediary metabolism and respiration
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,5],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,5]))

#lipid metabolism
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,6],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,6]))

#regulatory proteins
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,7],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,7]))

#virulence, detoxification, adaptation
chisq.test(cbind(table(hivneg.specific.factor,t(annotation_HIVneg[8,]))[2:3,9],table(hiv.specific.factor,t(annotation_HIV[8,]))[2:3,9]))

counts_HIV <-table(mutHIV_specific)
counts_HIVneg <-table(mutHIVneg_specific)
counts <-rbind(counts_HIV,counts_HIVneg)

barplot(counts,beside=T,col=c("red","darkblue"),xlab="# derived mutations",ylab="N°sites",main="SFS of specific Mtb/HIV+ and Mtb/HIV- mutations",space=c(0.1,0,0.1,0,0.1,0,0.1,0))

legend(6,1200,legend=c( "Mtb/HIV-","Mtb/HIV+"),col=c("darkblue", "red"), lty=1,cex=1.5,bty='n')

