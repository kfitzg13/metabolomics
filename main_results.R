library(WGCNA)
library(gee)
# set directory to location of github_metabolite_info and metab_data_template_final.csv files

mets <- read.csv(file="github_metabolite_info.csv", header=T, sep=",")
all_metab <- read.csv(file="metab_data_template_final.csv")

pheno3<- all_metab[, colnames(all_metab) %in% c("record_id", "status", "age", "gender", "race", "batch", "subtype", "dmt_clean", "edss")]
combat_edata3 <- all_metab[, !(colnames(all_metab) %in% c("record_id", "status", "age", "gender", "race", "batch", "subtype", "dmt_clean", "edss"))]

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(as.matrix(combat_edata3), powerVector = powers, verbose = 5)
sft

## select 4

net = blockwiseModules(as.matrix(combat_edata3), power=4, TOMType = "signed", minModuleSize = 5,  
                       reassignThreshold = 0, mergeCutHeight = 0.4, numericLabels = TRUE, pamRespectsDendro = FALSE, 
                       saveTOMs = F, verbose = 3)
table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
modulesA1=net$colors
moduleColorsToExtend = labels2colors(net$colors)

# calculate module scores from modules derived in the controls and project on the cases. 
MEso=moduleEigengenes(as.matrix(combat_edata3), colors= moduleColors)
MEso_results= MEso$eigengenes

status=pheno3$status
age= pheno3$age
gender= pheno3$gender
race= factor(pheno3 $race)
brid=factor(pheno3$record_id, labels=1:length(unique(pheno3$record_id)))
race<- relevel(race, ref="White")

wgcna_results <- NULL
for ( i in 1:ncol(MEso_results)){
  temp= MEso_results[,i]
  td=data.frame(status, age, gender, temp, race, brid=brid)
  td=td[order(td$brid),]
  m1=gee(temp~age+gender+status+race, data=td, id=brid, corstr="exchangeable")
  beta<-coef(summary(m1))[4,1]
  se<-coef(summary(m1))[4,4]
  pval<-2*pnorm(abs(coef(summary(m1))[4,5]), lower.tail=FALSE)
  mod<-c(beta, se, pval)
  wgcna_results<-rbind(wgcna_results, mod)
}

wgcna_results<-data.frame(wgcna_results)
colnames(wgcna_results)<-c("beta", "se", "pvalue")
wgcna_results$module<- colnames(MEso_results)
wgcna_results$fdrpv<- p.adjust(wgcna_results$pval)
wgcna_results$beta<- wgcna_results$beta*100 
wgcna_results$se <- wgcna_results$se* 100 
wgcna_results$lower <- wgcna_results$beta-1.96*wgcna_results$se
wgcna_results$upper <- wgcna_results$beta+1.96*wgcna_results$se
t3<- paste0(sprintf('%.2f', round(wgcna_results$beta,2)), " (", sprintf('%.2f', round(wgcna_results$lower,2)),", ",  sprintf('%.2f', round(wgcna_results$upper,2)), ")")

# main wgcna findings
wgcna_results$t3 <- t3

individual_results <- NULL
mean_diff=NULL
for ( i in 1:dim(combat_edata3)[2]){
  temp= combat_edata3[,i]/sd(combat_edata3[,i])
  temp2= combat_edata3[,i]
  td=data.frame(status, age, gender, temp, race, brid=brid, temp2)
  td=td[order(td$brid),]
  m1=gee(temp~age+gender+status+race, data=td, id=brid, corstr="exchangeable")
  beta<-coef(summary(m1))[4,1]
  se<-coef(summary(m1))[4,4]
  pval<-2*pnorm(abs(coef(summary(m1))[4,5]), lower.tail=FALSE)
  mod<-c(beta, se, pval)
  individual_results <- rbind(individual_results, mod)
  temp_mean_diff=mean(td$temp[which(td$status=="MS")])-mean(td$temp[which(td$status=="HC")])
  mean_diff=c(mean_diff, temp_mean_diff)
}
individual_results <- data.frame(individual_results)
individual_results$mean_diff= mean_diff
colnames(individual_results)<-c("beta", "se", "pvalue", "mean_diff")
individual_results $metabolite=colnames(combat_edata3)


individual_results $fdrpv<- p.adjust(individual_results $pval)
individual_results$n=1:nrow(individual_results)
individual_results $beta<- individual_results $beta* 10 
individual_results $se <- individual_results $se* 10
individual_results $lower <- individual_results $beta-1.96* individual_results $se
individual_results $upper <- individual_results $beta+1.96* individual_results $se
individual_results $ t3<- paste0(sprintf('%.2f', round(individual_results $beta,2)), " (", sprintf('%.2f', round(individual_results $lower,2)),", ",  sprintf('%.2f', round(individual_results $upper,2)), ")")


individual_results =merge(individual_results, mets, by.y="kfcomp", by.x="metabolite")
individual_results <-individual_results[order(individual_results$pvalue),]

individual_results_with_module <- merge(individual_results, module_results, by.x="metabolite",  by.y="kfcomp" )
individual_results_with_module <- individual_results_with_module[order(individual_results_with_module$pvalue),]
individual_results_with_module_significant <- individual_results_with_module[which(individual_results_with_module $fdrpv<(0.05)),]


edss=pheno3$edss
###############
wgcna_results_edss <- NULL
for ( i in 1:dim(MEso_results)[2]){
  temp= MEso_results[,i]
  td=data.frame(status, age, gender, temp, race, brid=brid, edss=edss, dmtclass=pheno3$dmtclass)
  td=td[order(td$brid),]
  m1=gee(temp~edss+age+gender+race, data=td, id=brid, corstr="exchangeable", na.action="na.omit", maxiter=1000)
  beta<-coef(summary(m1))[2,1]
  se<-coef(summary(m1))[2,4]
  pval<-2*pnorm(abs(coef(summary(m1))[2,5]), lower.tail=FALSE)
  mod<-c(beta, se, pval)
  wgcna_results_edss <-rbind(wgcna_results_edss, mod)
}

wgcna_results_edss <-data.frame(wgcna_results_edss)
colnames(wgcna_results_edss)<-c("beta", "se", "pvalue")
wgcna_results_edss $module<- colnames(MEso_results)
wgcna_results_edss $fdrpv<- p.adjust(wgcna_results_edss $pval)
wgcna_results_edss[which(wgcna_results_edss $fdrpv<0.05),]
wgcna_results_edss[which(wgcna_results_edss $pval<0.00313),]

individual_results_edss <- NULL
for ( i in 1:dim(combat_edata3)[2]){
  temp= combat_edata3[,i]/sd(combat_edata3[,i])
  td=data.frame(status, age, gender, temp, race, brid=brid, edss=edss, dmtclass=pheno3$dmtclass)
  td=td[order(td$brid),]
  m1=gee(temp~edss+age+gender+race, data=td, id=brid, corstr="exchangeable", na.action="na.omit", maxiter=1000)
  beta<-coef(summary(m1))[2,1]
  se<-coef(summary(m1))[2,4]
  pval<-2*pnorm(abs(coef(summary(m1))[2,5]), lower.tail=FALSE)
  mod<-c(beta, se, pval)
  individual_results_edss <- rbind(individual_results_edss, mod)
}
individual_results_edss <- data.frame(individual_results_edss)
colnames(individual_results_edss)<-c("beta", "se", "pvalue")
individual_results_edss $metabolite=colnames(combat_edata3)
individual_results_edss $fdrpv<- p.adjust(individual_results_edss $pval)
individual_results_edss =merge(individual_results_edss, mets, by.y="kfcomp", by.x="metabolite")
individual_results_edss <-individual_results_edss[order(individual_results_edss $pvalue),]
individual_results_edss[which(individual_results_edss $pvalue<(0.001)),]
individual_results_edss $lower <- individual_results_edss $beta-1.96* individual_results_edss $se
individual_results_edss $upper <- individual_results_edss $beta+1.96* individual_results_edss $se

individual_results_edss $ t3<- paste0(sprintf('%.2f', round(individual_results_edss $beta,2)), " (", 
                                      sprintf('%.2f', round(individual_results_edss $lower,2)),", ",  
                                      sprintf('%.2f', round(individual_results_edss $upper,2)), ")")

library(ggrepel)
  
individual_results $Significant <- ifelse(individual_results $fdrpv < 0.05, "FDR < 0.05", "Not Sig")
individual_results$logpv<-log10(individual_results$pvalue)
ggplot(individual_results, aes(x = mean_diff, y =-logpv)) +
    geom_point(aes(color = Significant), size=2) + xlim(c(-0.7, 0.4))+
    scale_color_manual(values = c("red", "grey")) +
    theme_bw(base_size = 12) + xlab("Mean difference: MS vs. HC")+ ylab("-log(p-value)") +theme(legend.position="none", axis.title.y=element_text(size=14, face="bold"), axis.title.x=element_text(size=14, face="bold"), panel.background = element_rect(fill = "white", colour = "black", size =1, linetype = "solid"), panel.grid.major = element_line(size = 0.2, 	linetype = 'solid', colour = "gray90"),panel.grid.minor = element_line(size = 0.1, linetype = 'solid',colour = "gray90"))+ geom_text_repel(
      data = subset(individual_results, fdrpv < 0.05),
      aes(label = BIOCHEMICAL_cr),
      size=2.2,
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.2, "lines")
    )


disease_colors_ms=c("blue", "darkorange")
combat_edata_resid <- matrix(NA, nrow=nrow(combat_edata3), ncol=ncol(combat_edata3))

for ( i in 1:dim(combat_edata3)[2]){
temp= combat_edata3[,i]/sd(combat_edata3[,i])
td=data.frame(status=pheno3$status, age=pheno3$age, gender=pheno3$gender,  temp, race=pheno3$race)
m2<- lm(temp~age+gender, data=td)
combat_edata_resid [,i]<- resid(m2)
}
#mahalanobis
library(vegan)
D <- as.matrix(vegdist(combat_edata_resid, method="mahalanobis"))
#mahalanobis
ref_set<-which(pheno3$status=="HC")
activity_index<- NULL
for (i in 1:nrow(D)){
	refset_no_i <- D[i, ref_set]
	activity_index <- c(activity_index , median(refset_no_i))
	}
t.test(activity_index[which(pheno3$status=="HC")], activity_index[which(pheno3$status=="MS")] )	
summary(glm(activity_index~pheno3$subtype, family="gaussian" ))

m1=lm(activity_index~pheno3$subtype)
m2=lm(activity_index~1)
anova(m2, m1)
metabolomics_diff <- data.frame(activity_index= activity_index, subtype=pheno3$subtype, status=pheno3$status)

def_plot <- ggplot(data= metabolomics_diff, aes(x=activity_index, color=status, fill=status)) +
    geom_density(alpha=0.85, size=0.6) +
    geom_vline(xintercept=disease_activity_threshold, size=.8, linetype="dashed") +
    scale_color_manual(values= disease_colors_ms, name="Diagnosis") +
    scale_fill_manual(values= disease_colors_ms, name="Diagnosis") + 
    xlab("Metabolic Dysfunction Score") + ylab("Density")+
    annotate("rect", xmin= disease_activity_threshold , xmax=Inf, ymin=0, ymax=Inf, alpha=0.15)+theme(legend.position="bottom",legend.title=element_text(size=12, face="bold"), axis.title.y=element_text(size=13, face="bold"), axis.title.x=element_text(size=13, face="bold"), panel.background = element_rect(fill = "white", colour = "black", size =1, linetype = "solid"), panel.grid.major = element_line(size = 0.2, 	linetype = 'solid', colour = "gray90"),panel.grid.minor = element_line(size = 0.1, linetype = 'solid',colour = "gray90"))+xlim(18,30)

plot(def_plot)

