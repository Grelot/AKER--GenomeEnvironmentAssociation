
############################
#fenetre 1kB
test = read.csv("tables_w1Kb/train/all_stats_propre.csv",header=F)
#fenetre 100kB avec id
test = read.csv("tables_w100Kb/avec_id/all_stats_propre.csv")
##fenetre 20kB avec id
test = read.csv("tables_w20Kb/avec_id/all_stats_propre.csv")
test = read.csv("tables_w20Kb/avec_id/all_stats.csv")

colnames(test)= c("id","chrom","position", "is_outlier", "nt_diversity","nt_diversity_per_site","D_tajima","Fu_li_D")

#snp density
snpdensity = test$nt_diversity/test$nt_diversity_per_site
snpdensity=replace(snpdensity,is.na(snpdensity),0)
test = cbind(test,snpdensity)


outlier = test[which(test$is_outlier==1),]
inlier = test[which(test$is_outlier==0),]

zin = inlier[which(inlier$snpdensity > 2),]


#wilcoxon tester liaison outlier/inlier (quali) versus nt_divers (quanti)
wilcox.test(outlier$nt_diversity,inlier$nt_diversity)
wilcox.test(outlier$nt_diversity_per_site,inlier$nt_diversity_per_site)
wilcox.test(outlier$D_tajima,inlier$D_tajima)



res = t.test(outlier$nt_diversity,inlier$nt_diversity)










