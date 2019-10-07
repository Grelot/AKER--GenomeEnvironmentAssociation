a = read.csv("data/Noms.marqueurs.LFMM.gINLAnd.csv",h=T,sep=";")

lffm = a[which(a$LFMM==1),]
gin = a[which(a$gINLAnd==1),]
sub84=merge(lffm, gin, by.x="probeset_id",by.y="probeset_id")

merge(lffm, gin, by.x="probeset_id",by.y="probeset_id")
write.table(file="nom_84_outliers.txt",sub84[,1],row.names=F,col.names=F,quote=F)


