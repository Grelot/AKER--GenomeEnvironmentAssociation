library(PopGenome)
library(ggplot2)
library(grid)



# Plot nt diversity and tajima D among sliding windows from a VCF
#
# author = Pierre-Edouard Guerin
#
# using ggplot2, grid and Popgenome R libraries, plot nt diversity and tajima D
# calculated on sliding windows on whole chromosome
# outliers SNP can be added as a list
# outputs calculation into tablesheets
#
# "monchrom" argument is a genome class object from package Popgenome
# "taille_fenetre" is size of sliding windows
# "saut_fenetre" is space between each "beginning point" of windows
# "min_pos" are first SNP position in the chromosome VCF
# "max_pos" are last SNP position in the chromosome VCF
# "monchrom_outliers" is a list of outliers SNP
# outliers SNP are display on nt diversity and tajima calculated on windows
#
calcul_chr = function(monchrom, taille_fenetre,saut_fenetre,min_pos,max_pos,monchrom_outliers,tid) {

##outliers SNP position
pos_out_d = data.frame(poutx=monchrom_outliers$V1)


##fenetre de taille_fenetre avec un saut de saut_fenetre
monchrom.slide10 = sliding.window.transform(monchrom,jump=saut_fenetre,width=taille_fenetre,type=2)
##calcule diversité nucleotidique
monchrom.slide.div = diversity.stats(monchrom.slide10)
##nombre de SNP par fenetre
n_w_snp=(lapply(monchrom.slide.div@region.data@biallelic.sites,length))
n_w_snp1=n_w_snp
n_w_snp1[n_w_snp1 ==0] =1
#calcul le tajima
monchrom.slide.neut = neutrality.stats(monchrom.slide10)


#min_pos = 156825
#max_pos = 42152870

#diversite nuc totale par fenetre
nuc_y = monchrom.slide.div@nuc.diversity.within
#position des fenetres sur lesquelles est calculée diversité nuc
decale_fenetre = taille_fenetre/2
nuc_x=seq(min_pos,max_pos,length.out=length(nuc_y))+decale_fenetre
nuc_d = data.frame(nuc_x=nuc_x,nuc_y=nuc_y)
nuc_d = replace(nuc_d, is.na(nuc_d), 0)


#diversite nuc par site par fenetre
nucs_y=nuc_y/as.numeric(n_w_snp1)
nucs_d=data.frame(nucs_x=nuc_x,nucs_y=nucs_y)


#tajima D par fenetre
D_y=monchrom.slide.neut@Tajima.D
D_y= replace(D_y, is.na(D_y), 0)
D_d=data.frame(D_x=nuc_x,D_y=D_y)

#Fu and Li D par fenetre
FL_y=monchrom.slide.neut@Fu.Li.D
FL_y= replace(FL_y, is.na(FL_y), 0)
FL_d=data.frame(FL_x=nuc_x,FL_y=FL_y)


#position des SNP sur le chromosome
pos_snp_x = monchrom@region.data@biallelic.sites[[1]]
pos_snp_d = data.frame(psx=pos_snp_x)


#SNP density
snpd_y=as.numeric(n_w_snp)
snpd_d=data.frame(snpd_x=nuc_x,snpd_y=snpd_y)


#######TABLE


pos_w_start=seq(min_pos,max_pos,length.out=length(nuc_y))
monchrom_outliers_id_w = findInterval(monchrom_outliers$V1, pos_w_start)
monchrom_outliers_is = rep(1,times=length(monchrom_outliers$V1))
total_snp = monchrom@region.data@biallelic.sites[[1]]
monchrom_inliers =setdiff(total_snp,monchrom_outliers)
monchrom_inliers_id_w = findInterval(monchrom_inliers, pos_w_start)
monchrom_inliers_is = rep(0,times=length(monchrom_inliers))
monchrom_liers= monchrom_outliers$V1
monchrom_liers = c(monchrom_liers,monchrom_inliers)
monchrom_liers_id_w = monchrom_outliers_id_w
monchrom_liers_id_w= c(monchrom_liers_id_w,monchrom_inliers_id_w)
monchrom_liers_is = monchrom_outliers_is
monchrom_liers_is = c(monchrom_liers_is,monchrom_inliers_is)
chr_is = rep(tid, times=length(monchrom_liers))
table_total = cbind(chr_is,monchrom_liers,monchrom_liers_is,nuc_y[monchrom_liers_id_w],nucs_y[monchrom_liers_id_w],D_y[monchrom_liers_id_w],FL_y[monchrom_liers_id_w])
colnames(table_total) = c("chrom","position", "is_outlier", "nt_diversity","nt_diversity_per_site","D_tajima","D_fu_li") 




#######PLOT

##plot SNP positions
#p_snp = ggplot()+geom_vline(size=0.1,data=pos_snp_d,mapping=aes(xintercept=psx))+
#theme_minimal()+
#theme(aspect.ratio=1/12,panel.grid.major = element_blank(),axis.text.x=element_blank())+
#xlim(1,max_pos+saut_fenetre)+
#xlab("")+
#ylab("SNP positions")

##plot SNP density
p_snp = ggplot()+
theme_minimal()+
theme(aspect.ratio=1/12,axis.text.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor= element_blank())+
geom_line(data=snpd_d,aes(x=snpd_x,y=snpd_y), color="violetred2",size=0.3)+
xlim(1,max_pos+saut_fenetre)+
xlab("")+
ylab("SNP density")+
geom_vline(size=0.2,data=pos_out_d,mapping=aes(xintercept=poutx))


##plot nucleotide diversity
p_nd = ggplot()+
theme_minimal()+
theme(aspect.ratio=1/12,axis.text.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor= element_blank())+
geom_line(data=nuc_d,aes(x=nuc_x,y=nuc_y), color="tomato1",size=0.3)+
xlim(1,max_pos+saut_fenetre)+
xlab("")+
ylab(expression(pi))+
geom_vline(size=0.2,data=pos_out_d,mapping=aes(xintercept=poutx))
#ggsave("monchrom_nt_div.pdf",dpi=300,device="pdf")

##plot nucleotide diversity PER SITE
p_ndps = ggplot()+
theme_minimal()+
theme(aspect.ratio=1/12,axis.text.x=element_blank(),panel.grid.major = element_blank(),panel.grid.minor= element_blank())+
geom_line(data=nucs_d,aes(x=nucs_x,y=nucs_y), color="dodgerblue1",size=0.3)+
xlim(1,max_pos+saut_fenetre)+
xlab("")+
ylab( expression(pi~per~site))+
geom_vline(size=0.2,data=pos_out_d,mapping=aes(xintercept=poutx))+
geom_hline(yintercept=mean(nucs_y[which(nucs_y != 0)]), color="red")
#ggsave("monchrom_nt_div_PER_SITE.pdf",dpi=300,device="pdf")


##plot Tajima D
p_tD = ggplot()+
theme_minimal()+
theme(aspect.ratio=1/12,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank())+
geom_line(data=D_d,aes(x=D_x,y=D_y), color="orange2",size=0.3)+
xlim(1,max_pos+saut_fenetre)+
ylim(min(D_y)-6,max(D_y)+1)+
annotate("segment", x = 250000, xend =10250000, y = -5.5, yend = -5.5, arrow=arrow(ends="both", angle=90, length=unit(.2,"cm")))+
annotate("text",x=5000000,y=-4.5,label="10Mb",size=4)+
geom_hline(yintercept=0, color="red")+
xlab("")+
ylab("D")+
geom_vline(size=0.2,data=pos_out_d,mapping=aes(xintercept=poutx))
#ggsave("monchrom_tajima_d.pdf",dpi=300,device="pdf")


##plot Fu&Li D
p_FL = ggplot()+
theme_minimal()+
theme(aspect.ratio=1/12,panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_blank())+
geom_line(data=FL_d,aes(x=FL_x,y=FL_y), color="green2",size=0.3)+
xlim(1,max_pos+saut_fenetre)+
geom_hline(yintercept=0, color="red")+
xlab("")+
ylab("D*")+
geom_vline(size=0.2,data=pos_out_d,mapping=aes(xintercept=poutx))




grid.newpage()
grid.draw(rbind(ggplotGrob(p_snp),ggplotGrob(p_nd), ggplotGrob(p_ndps),ggplotGrob(p_FL), ggplotGrob(p_tD), size = "last"))


return(table_total)
}



metacalcul_chr = function(pref, chr_min_pos,chr_max_pos,montid,taille_sw,jump_sw) {

vcfpath= paste("mes_vcf/",tolower(pref),sep="")
vcfpath = paste(vcfpath,".vcf.sorted.gz",sep="")
current_chrom = readVCF(vcfpath,tid=montid, frompos=chr_min_pos,topos=chr_max_pos,numcols=100000)
outlierspath = paste("mes_outliers/",tolower(pref),sep="")
outlierspath = paste(outlierspath, "_outliers.list",sep="")
chrom_outliers = read.table(outlierspath)
figurepath = paste("figures/",tolower(pref),sep="")
figurepath = paste(figurepath,"_stats.pdf",sep="")
tiff(figurepath, units="in", width=20, height=20, res=300)
chrom_table = calcul_chr(current_chrom, taille_sw,jump_sw,chr_min_pos,chr_max_pos,chrom_outliers,as.integer(montid))
dev.off()
#dev.print(pdf, figurepath)
tablepath = paste("tables/",tolower(pref),sep="")
tablepath = paste(tablepath,"_stats.csv",sep="")
write.csv(chrom_table,file=tablepath)

}


###################

taille_sw =20000
jump_sw=5000

###CHR1
chr1_min_pos=156825
chr1_max_pos=42152870
metacalcul_chr("CHR1", chr1_min_pos,chr1_max_pos,"1",taille_sw,jump_sw)

###CHR2
chr2_min_pos=154312
chr2_max_pos=45789363
metacalcul_chr("CHR2", chr2_min_pos,chr2_max_pos,"2",taille_sw,jump_sw)

###CHR3
chr3_min_pos=153905
chr3_max_pos=32334931
metacalcul_chr("CHR3", chr3_min_pos,chr3_max_pos,"3",taille_sw,jump_sw)

###CHR4
chr4_min_pos=151482
chr4_max_pos=47312317
metacalcul_chr("CHR4", chr4_min_pos,chr4_max_pos,"4",taille_sw,jump_sw)

###CHR5
chr5_min_pos=164491
chr5_max_pos=59641741
metacalcul_chr("CHR5", chr5_min_pos,chr5_max_pos,"5",taille_sw,jump_sw)

###CHR6
chr6_min_pos=156535
chr6_max_pos=63538038
metacalcul_chr("CHR6", chr6_min_pos,chr6_max_pos,"6",taille_sw,jump_sw)

###CHR7
chr7_min_pos=150155
chr7_max_pos=51878146
metacalcul_chr("CHR7", chr7_min_pos,chr7_max_pos,"7",taille_sw,jump_sw)
###CHR8
chr8_min_pos=150642
chr8_max_pos=41504428
metacalcul_chr("CHR8", chr8_min_pos,chr8_max_pos,"8",taille_sw,jump_sw)
###CHR9
chr9_min_pos=161334
chr9_max_pos=47785164
metacalcul_chr("CHR9", chr9_min_pos,chr9_max_pos,"9",taille_sw,jump_sw)


################



