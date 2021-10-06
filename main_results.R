library(WGCNA)
library(gee)
# set directory to location of github_metabolite_info and metab_data_template_final.csv files

mets <- read.csv(file="github_metabolite_info.csv", header=T, sep=",")
all_metab <- read.csv(file="metab_data_template_final.csv")

pheno3<- all_metab[, colnames(all_metab) %in% c("record_id", "status", "age", "gender", "race", "batch")]
combat_edata3 <- all_metab[, !(colnames(all_metab) %in% c("record_id", "status", "age", "gender", "race", "batch"))]


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
