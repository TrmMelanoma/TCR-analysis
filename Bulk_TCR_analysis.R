##########Venn diagram number of clonal matches using DNA sequences################

bulk_TCR_601<-read.csv("../bulkTCR/match_nt/601_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_601
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_601_all_three_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.triple.venn(area1=X1, area2=X2, area3=X3,  n12=X12, n13=X13, 
                 n23=X23,  n123=X123, category = c("skin","blood","lung_met"),
                 lwd = rep(2, 3), lty = rep("solid", 3), col =
                   c("#54278f","#fc4e2a","#fff7bc"), fill = c("#54278f","#fc4e2a","#fff7bc"), alpha = rep(0.5, 3))

dev.off()
colnames(bulk_TCR_601)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","Metastasis","Skin")
xx<-bulk_TCR_601[bulk_TCR_601$Present.In>1,]
xx[,"PT"]<-rep("PT601",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==3
tag2<-xx$PBMC>0 & xx$Metastasis>0 & xx$Skin==0
tag3<-xx$PBMC>0 & xx$Metastasis==0 & xx$Skin>0
tag4<-xx$PBMC==0 & xx$Metastasis>0 & xx$Skin>0
xx[tag1,"match_level"]<-"Skin Metastasis Blood Match"
xx[tag2,"match_level"]<-"Metastasis Blood Match Only"
xx[tag3,"match_level"]<-"Skin Blood Match Only"
xx[tag4,"match_level"]<-"Skin Metastasis Match Only"
myoutf <- "../bulkTCR/match_nt/PT601_matches_summary.xls"
write.table(xx,myoutf,sep='\t')





bulk_TCR_604<-read.csv("../bulkTCR/match_nt/604_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_604
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X5","X4","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
X4<-sum(venn_clonotype$X4>0)
X5<-sum(venn_clonotype$X5>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X14<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0)
X15<-sum(venn_clonotype$X1>0&venn_clonotype$X5>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X24<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0)
X25<-sum(venn_clonotype$X2>0&venn_clonotype$X5>0)
X34<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0)
X35<-sum(venn_clonotype$X3>0&venn_clonotype$X5>0)
X45<-sum(venn_clonotype$X4>0&venn_clonotype$X5>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
X124<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0)
X125<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X5>0)
X134<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X135<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X5>0)
X145<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
X234<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X235<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X5>0)
X245<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
X345<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
X1234<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X1235<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X5>0)
X1245<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
X1345<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
X2345<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
X12345<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0&venn_clonotype$X5>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_604_all_five_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.quintuple.venn(area1=X1, area2=X2, area3=X3, area4=X4, area5=X5, n12=X12, n13=X13, n14=X14, n15=X15,
                    n23=X23, n24=X24, n25=X25, n34=X34, n35=X35, n45=X45, n123=X123, n124=X124, n125=X125, n134=X134,
                    n135=X135, n145=X145, n234=X234, n235=X235, n245=X245, n345=X345, n1234=X1234, n1235=X1235,
                    n1245=X1245, n1345=X1345, n2345=X2345, n12345=X12345, category = c("skin","blood","LN_met","primary","bladder_met"),
                    lwd = rep(2, 5), lty = rep("solid", 5), col =
                      c("#54278f","#fc4e2a","#fff7bc","#f768a1","#78c679"), fill = c("#54278f","#fc4e2a","#fff7bc","#f768a1","#78c679"), alpha = rep(0.5, 5))

dev.off()
colnames(bulk_TCR_604)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","Metastasis","Primary","LN_Met","Skin")
xx<-bulk_TCR_604[bulk_TCR_604$Present.In>1,]
xx[,"PT"]<-rep("PT604",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==5
xx[tag1,"match_level"]<-"Skin Metastasis Primary LN_Met Blood Match"
tar<-colnames(xx)[4:8]
cob2<-combn(tar,2)
cob3<-combn(tar,3)
cob4<-combn(tar,4)
for(i in 1:ncol(cob2))
{
  target<-cob2[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0&xx[,non_target[3]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ","Match Only")
}

for(i in 1:ncol(cob3))
{
  target<-cob3[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ","Match")
}
for(i in 1:ncol(cob4))
{
  target<-cob4[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,target[4]]>0&xx[,non_target[1]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ",target[4]," ","Match")
}
myoutf <- "../bulkTCR/match_nt/PT604_matches_summary.xls"
write.table(xx,myoutf,sep='\t')



bulk_TCR_612<-read.csv("../bulkTCR/match_nt/612_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_612
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X4","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
X4<-sum(venn_clonotype$X4>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X14<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X24<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0)
X34<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
X124<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0)
X134<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X234<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X1234<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_612_all_four_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.quad.venn(area1=X1, area2=X2, area3=X3, area4=X4, n12=X12, n13=X13, n14=X14,
               n23=X23, n24=X24, n34=X34, n123=X123, n124=X124,  n134=X134,
               n234=X234, n1234=X1234, category = c("skin","blood","Primary","LN_met"),
               lwd = rep(2, 4), lty = rep("solid", 4), col =
                 c("#54278f","#fc4e2a","#fff7bc","#f768a1"), fill = c("#54278f","#fc4e2a","#fff7bc","#f768a1"), alpha = rep(0.5, 4))

dev.off()
colnames(bulk_TCR_612)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","Metastasis","Primary","Skin")
xx<-bulk_TCR_612[bulk_TCR_612$Present.In>1,]
xx[,"PT"]<-rep("PT612",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==4
xx[tag1,"match_level"]<-"Skin Metastasis Primary Blood Match"
tar<-colnames(xx)[4:7]
cob2<-combn(tar,2)
cob3<-combn(tar,3)
for(i in 1:ncol(cob2))
{
  target<-cob2[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ","Match Only")
}
for(i in 1:ncol(cob3))
{
  target<-cob3[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,non_target[1]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ","Match")
}
myoutf <- "../bulkTCR/match_nt/PT612_matches_summary.xls"
write.table(xx,myoutf,sep='\t')







bulk_TCR_625<-read.csv("../bulkTCR/match_nt/625_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_625
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X4","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
X4<-sum(venn_clonotype$X4>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X14<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X24<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0)
X34<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
X124<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0)
X134<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X234<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X1234<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_625_all_four_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.quad.venn(area1=X1, area2=X2, area3=X3, area4=X4, n12=X12, n13=X13, n14=X14,
               n23=X23, n24=X24, n34=X34, n123=X123, n124=X124,  n134=X134,
               n234=X234, n1234=X1234, category = c("skin","blood","Metastasis","Primary"),
               lwd = rep(2, 4), lty = rep("solid", 4), col =
                 c("#54278f","#fc4e2a","#fff7bc","#f768a1"), fill = c("#54278f","#fc4e2a","#fff7bc","#f768a1"), alpha = rep(0.5, 4))

dev.off()
colnames(bulk_TCR_625)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","Primary","Metastasis","Skin")
xx<-bulk_TCR_625[bulk_TCR_625$Present.In>1,]
xx[,"PT"]<-rep("PT625",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==4
xx[tag1,"match_level"]<-"Skin Metastasis Primary Blood Match"
tar<-colnames(xx)[4:7]
cob2<-combn(tar,2)
cob3<-combn(tar,3)
for(i in 1:ncol(cob2))
{
  target<-cob2[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ","Match Only")
}
for(i in 1:ncol(cob3))
{
  target<-cob3[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,non_target[1]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ","Match")
}
myoutf <- "../bulkTCR/match_nt/PT625_matches_summary.xls"
write.table(xx,myoutf,sep='\t')








bulk_TCR_626<-read.csv("../bulkTCR/match_nt/626_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_626
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X4","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
X4<-sum(venn_clonotype$X4>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X14<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X24<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0)
X34<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
X124<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0)
X134<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X234<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X1234<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_626_all_four_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.quad.venn(area1=X1, area2=X2, area3=X3, area4=X4, n12=X12, n13=X13, n14=X14,
               n23=X23, n24=X24, n34=X34, n123=X123, n124=X124,  n134=X134,
               n234=X234, n1234=X1234, category = c("skin","blood","Metastasis","LN_Met"),
               lwd = rep(2, 4), lty = rep("solid", 4), col =
                 c("#54278f","#fc4e2a","#fff7bc","#f768a1"), fill = c("#54278f","#fc4e2a","#fff7bc","#f768a1"), alpha = rep(0.5, 4))

dev.off()





colnames(bulk_TCR_626)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","LN_Met","Metastasis","Skin")
xx<-bulk_TCR_626[bulk_TCR_626$Present.In>1,]
xx[,"PT"]<-rep("PT626",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==4
xx[tag1,"match_level"]<-"Skin LN_Met Metastasis Blood Match"
tar<-colnames(xx)[4:7]
cob2<-combn(tar,2)
cob3<-combn(tar,3)
for(i in 1:ncol(cob2))
{
  target<-cob2[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ","Match Only")
}
for(i in 1:ncol(cob3))
{
  target<-cob3[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,non_target[1]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ","Match")
}
myoutf <- "../bulkTCR/match_nt/PT626_matches_summary.xls"
write.table(xx,myoutf,sep='\t')









bulk_TCR_628<-read.csv("../bulkTCR/match_nt/628_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_628
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X4","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
X4<-sum(venn_clonotype$X4>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X14<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X24<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0)
X34<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
X124<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0)
X134<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X234<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X1234<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_628_all_four_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.quad.venn(area1=X1, area2=X2, area3=X3, area4=X4, n12=X12, n13=X13, n14=X14,
               n23=X23, n24=X24, n34=X34, n123=X123, n124=X124,  n134=X134,
               n234=X234, n1234=X1234, category = c("skin","blood","Metastasis1","Metastasis2"),
               lwd = rep(2, 4), lty = rep("solid", 4), col =
                 c("#54278f","#fc4e2a","#fff7bc","#f768a1"), fill = c("#54278f","#fc4e2a","#fff7bc","#f768a1"), alpha = rep(0.5, 4))

dev.off()







colnames(bulk_TCR_628)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","Metastasis1","Metastasis2","Skin")
xx<-bulk_TCR_628[bulk_TCR_628$Present.In>1,]
xx[,"PT"]<-rep("PT628",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==4
xx[tag1,"match_level"]<-"Skin Metastasis1 Metastasis2 Blood Match"
tar<-colnames(xx)[4:7]
cob2<-combn(tar,2)
cob3<-combn(tar,3)
for(i in 1:ncol(cob2))
{
  target<-cob2[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ","Match Only")
}
for(i in 1:ncol(cob3))
{
  target<-cob3[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,non_target[1]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ","Match")
}
myoutf <- "../bulkTCR/match_nt/PT628_matches_summary.xls"
write.table(xx,myoutf,sep='\t')








bulk_TCR_635<-read.csv("../bulkTCR/match_nt/635_matches_nt.tsv",sep='\t')
venn_clonotype<-bulk_TCR_635
colnames(venn_clonotype)<- c("Amino.Acid","Sum..Templates.","Present.In","X2","X4","X3","X1")
X1<-sum(venn_clonotype$X1>0)
X2<-sum(venn_clonotype$X2>0)
X3<-sum(venn_clonotype$X3>0)
X4<-sum(venn_clonotype$X4>0)
venn_clonotype<-venn_clonotype[order(venn_clonotype$Present.In),]
venn_clonotype<-venn_clonotype[venn_clonotype$Present.In>1,]
X12<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0)
X13<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0)
X14<-sum(venn_clonotype$X1>0&venn_clonotype$X4>0)
X23<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0)
X24<-sum(venn_clonotype$X2>0&venn_clonotype$X4>0)
X34<-sum(venn_clonotype$X3>0&venn_clonotype$X4>0)
X123<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0)
X124<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X4>0)
X134<-sum(venn_clonotype$X1>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X234<-sum(venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
X1234<-sum(venn_clonotype$X1>0&venn_clonotype$X2>0&venn_clonotype$X3>0&venn_clonotype$X4>0)
myoutf <- paste0("../bulkTCR/match_nt/Venn_635_all_four_tissue.pdf")
pdf(myoutf,width=10,height=10)
draw.quad.venn(area1=X1, area2=X2, area3=X3, area4=X4, n12=X12, n13=X13, n14=X14,
               n23=X23, n24=X24, n34=X34, n123=X123, n124=X124,  n134=X134,
               n234=X234, n1234=X1234, category = c("skin","blood","Metastasis1","LN_Met"),
               lwd = rep(2, 4), lty = rep("solid", 4), col =
                 c("#54278f","#fc4e2a","#fff7bc","#f768a1"), fill = c("#54278f","#fc4e2a","#fff7bc","#f768a1"), alpha = rep(0.5, 4))

dev.off()
colnames(bulk_TCR_635)<-c("Amino.Acid","Sum..Templates.","Present.In","PBMC","LN_Met","Metastasis","Skin")
xx<-bulk_TCR_635[bulk_TCR_635$Present.In>1,]
xx[,"PT"]<-rep("PT635",nrow(xx))
xx[,"match_level"]<-rep("Not matched",nrow(xx))
tag1<-xx$Present.In==4
xx[tag1,"match_level"]<-"Skin LN_Met Metastasis Blood Match"
tar<-colnames(xx)[4:7]
cob2<-combn(tar,2)
cob3<-combn(tar,3)
for(i in 1:ncol(cob2))
{
  target<-cob2[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,non_target[1]]==0&xx[,non_target[2]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ","Match Only")
}
for(i in 1:ncol(cob3))
{
  target<-cob3[,i]
  non_target<-tar[!tar%in%target]
  tag<-xx[,target[1]]>0&xx[,target[2]]>0&xx[,target[3]]>0&xx[,non_target[1]]==0
  xx[tag,"match_level"]<-paste0(target[1]," ",target[2]," ",target[3]," ","Match")
}
myoutf <- "../bulkTCR/match_nt/PT635_matches_summary.xls"
write.table(xx,myoutf,sep='\t')




############Calculate the Match_Index##########

rm(list=ls())


myinf<-"../bulkTCR/all_arrangements/"
files<-list.files(myinf)
files<-files[grep(".tsv",files)]
tar<-gsub(".tsv","",files)


library(gtools)

for (i in 1:length(files))
{
  cat("\r",i)
  data<-read.csv(paste0(myinf,files[i]),sep='\t')
  matches<-matrix(0,length(unique(data$sample_name)),length(unique(data$sample_name)))
  row.names(matches)<-unique(data$sample_name)
  colnames(matches)<-unique(data$sample_name)
  
  
  all_number<-matrix(0,length(unique(data$sample_name)),length(unique(data$sample_name)))
  row.names( all_number)<-unique(data$sample_name)
  colnames( all_number)<-unique(data$sample_name)
  
  all_comb<-permutations(length(unique(data$sample_name)), 2, v=as.character(unique(data$sample_name)),
                         set=F, repeats.allowed=T)
  
  for (j in 1:nrow(all_comb))
  {
    
    cat("\r",j)
    pair<-all_comb[j,]
    tag1<-pair[1]
    tag2<-pair[2]
    xx1<-data[data$sample_name==tag1,]
    xx2<-data[data$sample_name==tag2,]
    
    match_num<-length(unique((merge(xx1,xx2,by="rearrangement")$rearrangement)))
    matches[tag1,tag2]<-match_num
    matches[tag2,tag1]<-match_num
    
    timed_number<-nrow(xx1)*nrow(xx2)
    all_number[tag1,tag2]<-timed_number
    all_number[tag2,tag1]<-timed_number
  }
  myoutf<-paste0("../bulkTCR/match/different_tissue_each_PT_match_number/PT_",tar[i],"_matched_clonotypes_number.xls")
  write.table(matches,myoutf,quote=F,sep='\t')
  
  match_level_index<-matches/all_number*100000
  
  myoutf<-paste0("../bulkTCR/match/different_tissue_each_PT_match_index/PT_",tar[i],"_matched_clonotypes_level_index.xls")
  write.table(match_level_index,myoutf,quote=F,sep='\t')
  
  
}




###########Plot the Match_Index##########



library(ggcorrplot)
library(reshape2)
library(ggpubr)

rm(list=ls())

myinf<-"../bulkTCR/match/different_tissue_each_PT_match_index/"
files<-list.files(myinf)
files<-files[grep(".xls",files)]
tar<-strsplit(files,"_")
tar<-paste0(lapply(tar,`[[`,1),lapply(tar,`[[`,2))

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}


reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


myinf2<-"../bulkTCR/match/different_tissue_each_PT_match_number/"
files1<-list.files(myinf2)
files1<-files1[grep(".xls",files)]



breaklist<-seq(0,3,0.1)

for (i in 1:length(files))
{
  cat("\r",i)
  data<-read.csv(paste0(myinf,files[i]),sep='\t')
  diag(data)<-0
  tag1<-grep("skin",row.names(data))
  tag2<-grep("PBMC",row.names(data))
  tag3<-grep("melanoma",row.names(data))
  
  data<-data[c(tag3,tag2,tag1),c(tag3,tag2,tag1)]
  
  data_number<-read.csv(paste0(myinf2,files1[i]),sep='\t')
  diag(data_number)<-0
  data_number<-data_number[c(tag3,tag2,tag1),c(tag3,tag2,tag1)]
  
  
  data<-get_lower_tri(data)
  data_number<-get_lower_tri(data_number)
  
  melted_data <- melt(as.matrix(data), na.rm = TRUE)
  melted_data_number <- melt(as.matrix(data_number), na.rm = TRUE)
  
  
  myoutf<-paste0("../bulkTCR/match/different_tissue_each_PT_match_index/plots/restructured/lower_tri_PT_",tar[i],"_matched_level_index_corrplot.pdf")
  pdf(myoutf,width=6,height=5)
  print(ggplot(melted_data, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
      geom_text(aes(Var2, Var1, label = melted_data_number$value), color = "black", size = 6)+
    scale_fill_gradient2(low = "grey90", high = "red1", mid = "tomato", 
                         midpoint = round(max(melted_data$value+0.5),0)/2, limit = c(0,round(max(melted_data$value+0.5),0)), space = "Lab", 
                         name="Match level index") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 6, hjust = 1,color="black"),
          axis.text.y = element_text(angle = 0, vjust = 0, 
                                     size = 6, hjust = 1,color="black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank())+
    coord_fixed())
  
  dev.off()
  
}











###########scatter plot show the frequency of skin tumor matches#############

rm(list=ls())

library(DescTools)
library(scales) 
myinf<-"../bulkTCR/all_arrangements/"
files<-list.files(myinf)
files<-files[grep(".tsv",files)]
tar<-gsub(".tsv","",files)




library(gtools)
ST_all<-data.frame()
STB_all<-data.frame()


ST_freq_all<-data.frame()
STB_freq_all<-data.frame()


for (i in 1:length(files))
{
  cat("\r",i)
  data<-read.csv(paste0(myinf,files[i]),sep='\t')
  
  all_tumor<-as.character(unique(data$sample_name[grep("melanoma",data$sample_name)]))
  all<-matrix(0,length(unique(data$rearrangement)),length(all_tumor)+2)
  row.names(all)<-unique(data$rearrangement)
  colnames(all)<-combine("skin","blood",all_tumor)
  
  
  all_freq<-matrix(0,length(unique(data$rearrangement)),length(all_tumor)+2)
  row.names(all_freq)<-unique(data$rearrangement)
  colnames(all_freq)<-combine("skin","blood",all_tumor)
  
  
  #all<-as.data.frame(all)
  
  tag<-grep("PBMC",data$sample_name)
  data_blood<-data[tag,]
  
  tag<-grep("skin",data$sample_name)
  data_skin<-data[tag,]
  
  
  all[as.character(data_skin$rearrangement),"skin"]<-data_skin$templates
  all[as.character(data_blood$rearrangement),"blood"]<-data_blood$templates
  
  
  all_freq[as.character(data_skin$rearrangement),"skin"]<-data_skin$productive_frequency
  all_freq[as.character(data_blood$rearrangement),"blood"]<-data_blood$productive_frequency
  
  
  
  for(k in 1:length(all_tumor))
  {
    cat("\r",k)
    tag<-grep(all_tumor[k],data$sample_name)
    data_tumor<-data[tag,]
    all[as.character(data_tumor$rearrangement),all_tumor[k]]<-data_tumor$templates
  }
  
  
  for(k in 1:length(all_tumor))
  {
    cat("\r",k)
    tag<-grep(all_tumor[k],data$sample_name)
    data_tumor<-data[tag,]
    all_freq[as.character(data_tumor$rearrangement),all_tumor[k]]<-data_tumor$productive_frequency
  }
  
  
  
  
  myoutf<-paste0("../bulkTCR/scatter_ST_frequency/",tar[i],"/all_clonotypes_",tar[i],".xls")
  write.table(all,myoutf,sep='\t')
  
  
  
  myoutf<-paste0("../bulkTCR/scatter_ST_frequency/",tar[i],"/all_clonotypes_inFreq_",tar[i],".xls")
  write.table(all_freq,myoutf,sep='\t')
  
  
  
  for(m in 1:length(all_tumor))
  {
    ST<-as.data.frame(all[,c("skin",all_tumor[m])])
    ST<-ST[!(ST[,1]==0&ST[,2]==0),]
    ST[,"group"]<-rep("both",nrow(ST))
    ST$group[ST[,1]==0]<-"Only tumor"
    ST$group[ST[,2]==0]<-"Only skin"
    colnames(ST)<-c("skin","tumor","type")
    ST[,"ID"]<-rep(all_tumor[m],nrow(ST))
    ST_all<-rbind(ST_all,ST)
    
    
    
    
    
    ST_freq<-as.data.frame(all_freq[,c("skin",all_tumor[m])])
    ST_freq<-ST_freq[!(ST_freq[,1]==0&ST_freq[,2]==0),]
    ST_freq[,"group"]<-rep("both",nrow(ST_freq))
    ST_freq$group[ST_freq[,1]==0]<-"Only tumor"
    ST_freq$group[ST_freq[,2]==0]<-"Only skin"
    colnames(ST_freq)<-c("skin","tumor","type")
    ST_freq[,"ID"]<-rep(all_tumor[m],nrow(ST_freq))
    ST_freq_all<-rbind(ST_freq_all,ST_freq)
    
    
    
    
    
    myoutf<-paste0("../bulkTCR/scatter_ST_frequency/",tar[i],"/scatterplot_skin",all_tumor[m],".pdf")
    pdf(myoutf,width=6,height=5)
    print(ggplot(ST, aes(x=skin, y=tumor, color=type, shape=type,size=type)) +
            geom_jitter()+
            scale_shape_manual(values=c(17, 16, 16))+ 
            scale_size_manual(values=c(1.5, 1, 1))+
            scale_x_continuous(breaks=seq(0,50,5))+
            scale_y_continuous(breaks=seq(0,50,5))+
            coord_cartesian(xlim = c(0,50), ylim = c(0,50), expand = TRUE,
                            default = FALSE, clip = "on")+
            scale_color_manual(values=c('#de2d26','#d9d9d9', '#d9d9d9'))+
            geom_hline(yintercept=3, linetype="dashed", color = "blue")+
            geom_vline(xintercept=3, linetype="dashed", color = "blue")+
            theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
    dev.off()
    
    
    
    
    
    
    
    
    myoutf<-paste0("../bulkTCR/scatter_ST_frequency/",tar[i],"/Freq_scatterplot_skin",all_tumor[m],".pdf")
    pdf(myoutf,width=6,height=5)
    print(ggplot(ST_freq, aes(x=skin, y=tumor, color=type, shape=type,size=type)) +
            geom_jitter()+
            scale_shape_manual(values=c(17, 16, 16))+ 
            scale_x_log10(breaks=c(.001,.01,.1))+
            scale_y_log10(breaks=c(.001,.01,.1))+
            
            scale_size_manual(values=c(1.5, 1, 1))+
            scale_color_manual(values=c('#de2d26','#d9d9d9', '#d9d9d9'))+
            theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
    dev.off()
    
    
    
    
    
    
    
    
    
    
    
    
    ST<-as.data.frame(all[,c("skin",all_tumor[m])])
    ST<-ST[(ST[,2]>=10),]
    STB<-cbind(ST,all[row.names(ST),"blood"])
    colnames(STB)[3]<-"blood"
    colnames(STB)[2]<-"tumor"
    
    STB[,"type"]<-rep("Also in Blood",nrow(STB))
    STB$type[STB[,3]==0]<-"Only in Skin and Tumor"
    STB[,"ID"]<-rep(all_tumor[m],nrow(STB))
    STB_all<-rbind(STB_all,STB)
    
    
    
    
    
    
    
    ST_freq<-as.data.frame(all_freq[,c("skin",all_tumor[m])])
    #ST_freq<-ST_freq[(ST_freq[,2]>0.002),]
    STB_freq<-cbind(ST_freq,all_freq[row.names(ST_freq),"blood"])
    colnames(STB_freq)[3]<-"blood"
    colnames(STB_freq)[2]<-"tumor"
    
    STB_freq[,"type"]<-rep("Also in Blood",nrow(STB_freq))
    STB_freq$type[STB_freq[,3]==0]<-"Only in Skin and Tumor"
    STB_freq[,"ID"]<-rep(all_tumor[m],nrow(STB_freq))
    STB_freq_all<-rbind(STB_freq_all,STB_freq)
    
    
  }
  
}


tar<-c("601_melanoma_lung_met_FFPE_TCRB" ,"604_melanoma_primary_FFPE_TCRB", "604_melanoma_bladder_met_FFPE_TCRB", "604_melanoma_sentinel_LN_FFPE_TCRB",
       "612_melanoma_primary_FFPE_TCRB","625_melanoma_primary_FFPE_TCRB","635_melanoma_inguinal_LN_FFPE_TCRB")


STB_freq_all_3yreas<-STB_freq_all[STB_freq_all$ID%in%tar,]
tag1<-STB_freq_all_3yreas$skin>0&STB_freq_all_3yreas$tumor>0&STB_freq_all_3yreas$blood==0
tag2<-STB_freq_all_3yreas$skin==0&STB_freq_all_3yreas$tumor>0&STB_freq_all_3yreas$blood>0
tag3<-STB_freq_all_3yreas$skin>0&STB_freq_all_3yreas$tumor>0&STB_freq_all_3yreas$blood>0
STB_freq_all_3yreas$type<-rep("Other",nrow(STB_freq_all_3yreas))
STB_freq_all_3yreas$type[tag1]<-"ST"
STB_freq_all_3yreas$type[tag2]<-"BT"
STB_freq_all_3yreas$type[tag3]<-"STB"

tag<-STB_freq_all_3yreas$tumor>0.002
STB_freq_all_3yreas_0.002<-STB_freq_all_3yreas[tag,]
colnames(STB_freq_all_3yreas_0.002)[4:5]<-c("type","ID")

ST_freq_all_3yreas<-ST_freq_all[ST_freq_all$ID%in%tar,]

colnames(STB_freq_all_3yreas)[4:5]<-c("type","ID")
colnames(ST_freq_all_3yreas)[3:4]<-c("type","ID")

color = brewer.pal(n = 7, name ="Set1")
#color=c(color,"#F0F0F0")

options(scipen=999)
myoutf<-paste0("../bulkTCR/scatter_ST_frequency/All_three_years_tumor/In_frequency_tumor_expanded_0.002_ST_frequency.pdf")
pdf(myoutf,width=10,height=7)
print(ggplot(STB_freq_all_3yreas_0.002, aes(x=skin, y=tumor, color=ID)) +
        geom_point()+
        scale_size_manual(values=c(1.2))+
        scale_x_continuous(trans = pseudo_log_trans(0.0001, 10),
                           breaks=c(0.0001, 0.001, 0.01,0.05))+
        scale_y_continuous(trans = pseudo_log_trans(0.01, 10),
                           breaks=c(0.001, 0.002, 0.01,0.05))+
        expand_limits(x = 0, y = 0)+
        scale_color_manual(values=color)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
dev.off()




tag1<-STB_freq_all_3yreas_0.002$tumor>0&STB_freq_all_3yreas_0.002$skin>0
tag2<-STB_freq_all_3yreas_0.002$tumor>0&STB_freq_all_3yreas_0.002$skin==0

tt1<-STB_freq_all_3yreas_0.002[tag1,]
tt2<-STB_freq_all_3yreas_0.002[tag2,]



myoutf<-paste0("../bulkTCR/scatter_ST_frequency/All_three_years_tumor/Freq_tumor_expanded_0.002_ST_bloodfrequency.pdf")
pdf(myoutf,width=12,height=7)
print(ggplot(tt1, aes(x=blood, y=tumor, color=ID)) +
        geom_point()+
        scale_size_manual(values=c(1.2))+
        scale_x_continuous(trans = pseudo_log_trans(0.0001, 10),
                           breaks=c(0.0001, 0.001, 0.01,0.05))+
        scale_y_continuous(trans = pseudo_log_trans(0.01, 10),
                           breaks=c(0.001, 0.002, 0.01,0.05))+
        expand_limits(x = 0, y = 0)+
        scale_color_manual(values=color)+
        geom_vline(xintercept=0.00001, linetype="dashed", color = "blue",size=0.6)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
dev.off()




myoutf<-paste0("../bulkTCR/scatter_ST_frequency/All_three_years_tumor/Freq_tumor_expanded_0.002_tumor_only_bloodfrequency.pdf")
pdf(myoutf,width=12,height=7)
print(ggplot(tt2, aes(x=blood, y=tumor, color=ID)) +
        geom_point()+
        #geom_jitter(group = "ID",width = 0.05, height = 0.05)+
        scale_size_manual(values=c(1.2))+
        scale_x_continuous(trans = pseudo_log_trans(0.0001, 10),
                           breaks=c(0.0001, 0.001, 0.01,0.05))+
        scale_y_continuous(trans = pseudo_log_trans(0.01, 10),
                           breaks=c(0.001, 0.002, 0.01,0.05))+
        expand_limits(x = 0, y = 0)+
        #coord_cartesian(xlim = c(0,50), expand = TRUE,
        #                default = FALSE, clip = "on")+
        scale_color_manual(values=color)+
        geom_vline(xintercept=0.00001, linetype="dashed", color = "blue",size=0.6)+
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
dev.off()















###########Match bulkTCR data with single cell TCR data##########


##########Process single cell TCR data

tar<-c(1,3,5,6)
barcode_info_all<-data.frame()

for(i in 1:length(tar))
{
  myinf1<-paste0("../scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_filtered_contig_annotations.csv")
  myinf2<-paste0("../scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_filtered_contig_annotations.csv")
  myinf3<-paste0("../scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_filtered_contig_annotations.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$raw_clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$raw_clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$raw_clonotype_id,"_",data_tumor$label,tar[i])
  
  barcode_info<-rbind(data_skin,data_blood,data_tumor)
  
  
  tag<-barcode_info$is_cell=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$high_confidence=="True"
  barcode_info<-barcode_info[tag,]
  tag<-barcode_info$productive=="True"
  barcode_info<-barcode_info[tag,]
  
  barcode_info$barcode<-gsub("-1","",barcode_info$barcode)
  barcode_info$barcode<-paste0(barcode_info$label,"_",barcode_info$barcode)
  barcode_info$barcode<-gsub("_",paste0("_",tar[i],"_"),barcode_info$barcode)
  barcode_info_all<-rbind(barcode_info,barcode_info_all)
}





tar<-c(1,3,5,6)
clonotype_STB_all_skin<-character()
clonotype_STB_all_tumor<-character()
clonotype_STB_all_blood<-character()

clonotype_ST_all_tumor<-character()
clonotype_ST_all_skin<-character()

clonotype_BT_all_tumor<-character()
clonotype_BT_all_blood<-character()



all_info_STB<-data.frame()

all_info_ST<-data.frame()

all_info_BT<-data.frame()

all_info_blood_only<-data.frame()


clonotype_label_STB_all<-character()
clonotype_label_ST_all<-character()


for(i in 1:length(tar))
{
  cat('\r',i)
  myinf1<-paste0("../scTCRseq/Run_",tar[i],"/TCR/data/skin/skin_clonotypes.csv")
  myinf2<-paste0("../scTCRseq/Run_",tar[i],"/TCR/data/blood/blood_clonotypes.csv")
  myinf3<-paste0("../scTCRseq/Run_",tar[i],"/TCR/data/tumor/tumor_clonotypes.csv")
  data_skin<-read.csv(myinf1)
  data_blood<-read.csv(myinf2)
  data_tumor<-read.csv(myinf3)
  data_skin[,"label"]<-rep("skin",nrow(data_skin))
  data_blood[,"label"]<-rep("blood",nrow(data_blood))
  data_tumor[,"label"]<-rep("tumor",nrow(data_tumor))
  
  data_skin[,"labeled_clonotype_id"]<-paste0(data_skin$clonotype_id,"_",data_skin$label,tar[i])
  data_blood[,"labeled_clonotype_id"]<-paste0(data_blood$clonotype_id,"_",data_blood$label,tar[i])
  data_tumor[,"labeled_clonotype_id"]<-paste0(data_tumor$clonotype_id,"_",data_tumor$label,tar[i])
  
  
  tag<-data_blood$cdr3s_nt%in%data_skin$cdr3s_nt|data_blood$cdr3s_nt%in%data_tumor$cdr3s_nt
  data_blood_1<-data_blood[!tag,]
  all_info_blood_only<-rbind(all_info_blood_only,data_blood_1)
  
  
  ST_match<-merge(data_skin,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
  STB_match<-merge( ST_match,data_blood,by="cdr3s_nt")
  clonotype_STB_all_skin<-combine(clonotype_STB_all_skin,STB_match$labeled_clonotype_id.x)
  clonotype_STB_all_tumor<-combine(clonotype_STB_all_tumor,STB_match$labeled_clonotype_id.y)
  clonotype_STB_all_blood<-combine(clonotype_STB_all_blood,STB_match$labeled_clonotype_id)
  all_info_STB<-rbind(all_info_STB,STB_match)
  
  
  tag<-ST_match$cdr3s_nt%in%data_blood$cdr3s_nt
  ST_spec_match<-ST_match[!tag,]
  clonotype_ST_all_tumor<-combine(clonotype_ST_all_tumor,ST_spec_match$labeled_clonotype_id.y)
  clonotype_ST_all_skin<-combine(clonotype_ST_all_skin,ST_spec_match$labeled_clonotype_id.x)
  all_info_ST<-rbind( all_info_ST,ST_spec_match)
  
  
  BT_match<-merge(data_blood,data_tumor,by.x="cdr3s_nt",by.y="cdr3s_nt")
  tag<-BT_match$cdr3s_nt%in%STB_match$cdr3s_nt
  BT_spec_match<-BT_match[!tag,]
  clonotype_BT_all_tumor<-combine(clonotype_BT_all_tumor,BT_spec_match$labeled_clonotype_id.y)
  clonotype_BT_all_blood<-combine(clonotype_BT_all_blood,BT_spec_match$labeled_clonotype_id.x)
  all_info_BT<-rbind(all_info_BT, BT_spec_match)
  
  
}

all_info_STB<-data.frame(lapply(all_info_STB, as.character), stringsAsFactors=FALSE)







##########Process bulk TCR data and match with single cell TCR data


myinf<-"../bulkTCR/all_arrangements/"
files<-list.files(myinf)
files<-files[grep(".tsv",files)]

tar<-gsub(".tsv","",files)

data_all<-data.frame()

for (i in 1:length(files))
{
  xx<-read.csv(paste0(myinf,files[i]),sep='\t')
  data_all<-rbind(data_all,xx)
}

zzz<-data.frame()
tag<-unique(names(unlist(sapply(yy$cdr3_nt, grep, row.names(STB_freq_all_3yreas_0.002)))))
for(i in 1:length(tag))
{
  xxx<-(grep(tag[i],row.names(tt1)))
  yyy<-tt1[xxx,]
  zzz<-rbind(zzz,yyy)
  
}

data_635<-data_all[grep("635",data_all$sample_name),]

long_term<-rbind(data_628[!grepl("scRNA",data_628$sample_name),],data_635[!grepl("scRNA",data_635$sample_name),])
tar<-c("635")

for(i in 1:length(tar))
{
  data<-long_term[grep(tar[i],long_term$sample_name),]
  data_skin<-data[grep("skin",data$sample_name),]
  data_tumor<-data[grep("melanoma",data$sample_name),]
  data_blood<-data[grep("PBMC",data$sample_name),]
  data_skin<-data_skin[data_skin$templates>0,]
  data_tumor<-data_tumor[data_tumor$templates>10,]
  data_blood<-data_blood[data_blood$templates>0,]
  
  
  
  ST<-merge(data_skin,data_tumor,by="rearrangement")
  BT<-merge(data_blood,data_tumor,by="rearrangement")
  STB<-merge(ST,data_blood,by="rearrangement")
  ST_spec<-ST[!ST$rearrangement%in%STB$rearrangement,]
  BT_spec<-BT[!BT$rearrangement%in%STB$rearrangement,]
  
  
  STB_freq_all_3yreas_0.002_635<-STB_freq_all_3yreas_0.002[grep("635",STB_freq_all_3yreas_0.002$ID),]
  
  
  STB_long_term_tumor_expanded_635<-STB_freq_all_3yreas_0.002_635[STB_freq_all_3yreas_0.002_635$type=='STB',]
  ST_long_term_tumor_expanded_635<-STB_freq_all_3yreas_0.002_635[STB_freq_all_3yreas_0.002_635$type=='ST',]
  BT_long_term_tumor_expanded_635<-STB_freq_all_3yreas_0.002_635[STB_freq_all_3yreas_0.002_635$type=='BT',]
  
  
  
  
  
  barcode_info_all<- data.frame(lapply(barcode_info_all, as.character), stringsAsFactors=FALSE)
  STB<- data.frame(lapply(STB, as.character), stringsAsFactors=FALSE)
  ST_spec<-data.frame(lapply(ST_spec, as.character), stringsAsFactors=FALSE)
  BT_spec<-data.frame(lapply(BT_spec, as.character), stringsAsFactors=FALSE)
  
  tag1<-unique(names(unlist( sapply(barcode_info_all$cdr3_nt, grep, row.names(STB_long_term_tumor_expanded_635) ) )))
  STB_sc_info<-barcode_info_all[barcode_info_all$cdr3_nt%in%tag1,]
  
  tag2<-unique(names(unlist( sapply(barcode_info_all$cdr3_nt, grep, row.names(ST_long_term_tumor_expanded_635) ) )))
  ST_spec_sc_info<-barcode_info_all[barcode_info_all$cdr3_nt%in%tag2,]
  
  tag3<-unique(names(unlist( sapply(barcode_info_all$cdr3_nt, grep, row.names(BT_long_term_tumor_expanded_635) ) )))
  BT_spec_sc_info<-barcode_info_all[barcode_info_all$cdr3_nt%in%tag3,]
  
  STB_sc_info_brief<-unique(STB_sc_info[,c("barcode","cdr3_nt","label","labeled_clonotype_id")])
  ST_spec_sc_info<-unique(ST_spec_sc_info[,c("barcode","cdr3_nt","label","labeled_clonotype_id")])
  BT_spec_sc_info<-unique(BT_spec_sc_info[,c("barcode","cdr3_nt","label","labeled_clonotype_id")])
  
  myinf<-"../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
  CD8_STB1356<-readRDS(myinf)
  
  pair_ST_clonotype<-rbind(cbind(all_info_STB$labeled_clonotype_id.x,all_info_STB$labeled_clonotype_id.y),cbind(all_info_ST$labeled_clonotype_id.x,all_info_ST$labeled_clonotype_id.y))
  pair_BT_clonotype<-rbind(cbind(all_info_STB$labeled_clonotype_id,all_info_STB$labeled_clonotype_id.y),cbind(all_info_BT$labeled_clonotype_id.x,all_info_BT$labeled_clonotype_id.y))
  col<-c('#2ca25f',"#de2d26","#9ecae1")
  
  
  
  
  
  
  cluster_number<-matrix(0,4,10)
  row.names(cluster_number)<-paste0("PT_",c(1,3,5,6))
  colnames(cluster_number)<-seq(0,9,1)
  for(w in 1:ncol(cluster_number))
  {
    tag<-CD8_STB1356@meta.data$res.0.4_combined==colnames(cluster_number)[w]
    xx<-CD8_STB1356@meta.data[tag,]
    cluster_number[1,w]<-sum(xx$PT_number=="PT_1")
    cluster_number[2,w]<-sum(xx$PT_number=="PT_3")
    cluster_number[3,w]<-sum(xx$PT_number=="PT_5")
    cluster_number[4,w]<-sum(xx$PT_number=="PT_6")
  }
  
  
  
  
  
  number_ST_skin_all_clonotype<-matrix(NA,length(unique(ST_spec_sc_info$cdr3_nt)),10)
  row.names(number_ST_skin_all_clonotype)<-paste0("clonotype_",1:length(unique(ST_spec_sc_info$cdr3_nt)))
  colnames(number_ST_skin_all_clonotype)<-c(0:9)
  
  
  number_ST_tumor_all_clonotype<-matrix(NA,length(unique(ST_spec_sc_info$cdr3_nt)),10)
  row.names(number_ST_tumor_all_clonotype)<-paste0("clonotype_",1:length(unique(ST_spec_sc_info$cdr3_nt)))
  colnames(number_ST_tumor_all_clonotype)<-c(0:9)
  
  
  number_ST_blood_all_clonotype<-matrix(NA,length(unique(ST_spec_sc_info$cdr3_nt)),10)
  row.names(number_ST_blood_all_clonotype)<-paste0("clonotype_",1:length(unique(ST_spec_sc_info$cdr3_nt)))
  colnames(number_ST_blood_all_clonotype)<-c(0:9)
  
  for(k in 1:length(unique(ST_spec_sc_info$cdr3_nt)))
  {
    xx<-ST_spec_sc_info[ST_spec_sc_info$cdr3_nt==unique(ST_spec_sc_info$cdr3_nt)[k],]
    
    CD8_STB1356@meta.data[,"bulk_TCR_match"]<-rep("other",nrow(CD8_STB1356@meta.data))
    tag1<-row.names(CD8_STB1356@meta.data)%in%xx$barcode
    CD8_STB1356@meta.data$bulk_TCR_match[tag1]<-paste0("Bulk_STB","_",CD8_STB1356@meta.data$label[tag1])
    if(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))>1){
      
      myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/ST/bulk_TCR_ST_spec_matched_",k,".pdf")
      pdf(myoutf,width=7.5,height=5)
      print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match',shape.by="PT_number",reduction="umap",
                    cols = c(col[1:(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))-1)],"#f0f0f0"),pt.size = 0.3))
      dev.off()
      qq<-table(CD8_STB1356@meta.data$bulk_TCR_match,CD8_STB1356@meta.data$res.0.4_combined)
      if(length(qq[grep("skin",row.names(qq)),])>0){number_ST_skin_all_clonotype[k,]<-qq[grep("skin",row.names(qq)),]}
      if(length(qq[grep("tumor",row.names(qq)),])>0){number_ST_tumor_all_clonotype[k,]<-qq[grep("tumor",row.names(qq)),]}
      if(length(qq[grep("blood",row.names(qq)),])>0){number_ST_blood_all_clonotype[k,]<-qq[grep("blood",row.names(qq)),]}
      
      
      
    }
    
  }
  
  
  
  
  
  
  number_BT_skin_all_clonotype<-matrix(NA,length(unique(BT_spec_sc_info$cdr3_nt)),10)
  row.names(number_BT_skin_all_clonotype)<-paste0("clonotype_",1:length(unique(BT_spec_sc_info$cdr3_nt)))
  colnames(number_BT_skin_all_clonotype)<-c(0:9)
  
  
  number_BT_tumor_all_clonotype<-matrix(NA,length(unique(BT_spec_sc_info$cdr3_nt)),10)
  row.names(number_BT_tumor_all_clonotype)<-paste0("clonotype_",1:length(unique(BT_spec_sc_info$cdr3_nt)))
  colnames(number_BT_tumor_all_clonotype)<-c(0:9)
  
  
  number_BT_blood_all_clonotype<-matrix(NA,length(unique(BT_spec_sc_info$cdr3_nt)),10)
  row.names(number_BT_blood_all_clonotype)<-paste0("clonotype_",1:length(unique(BT_spec_sc_info$cdr3_nt)))
  colnames(number_BT_blood_all_clonotype)<-c(0:9)
  
  for(k in 1:length(unique(BT_spec_sc_info$cdr3_nt)))
  {
    xx<-BT_spec_sc_info[BT_spec_sc_info$cdr3_nt==unique(BT_spec_sc_info$cdr3_nt)[k],]
    
    CD8_STB1356@meta.data[,"bulk_TCR_match"]<-rep("other",nrow(CD8_STB1356@meta.data))
    tag1<-row.names(CD8_STB1356@meta.data)%in%xx$barcode
    CD8_STB1356@meta.data$bulk_TCR_match[tag1]<-paste0("Bulk_STB","_",CD8_STB1356@meta.data$label[tag1])
    if(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))>1){
      myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/BT/bulk_TCR_BT_spec_matched_",k,".pdf")
      pdf(myoutf,width=7.5,height=5)
      print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match',shape.by="PT_number",reduction="umap",
                    cols = c(col[1:(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))-1)],"#f0f0f0"),pt.size = 0.3))
      dev.off()
      
      qq<-table(CD8_STB1356@meta.data$bulk_TCR_match,CD8_STB1356@meta.data$res.0.4_combined)
      if(length(qq[grep("skin",row.names(qq)),])>0){number_BT_skin_all_clonotype[k,]<-qq[grep("skin",row.names(qq)),]}
      if(length(qq[grep("tumor",row.names(qq)),])>0){number_BT_tumor_all_clonotype[k,]<-qq[grep("tumor",row.names(qq)),]}
      if(length(qq[grep("blood",row.names(qq)),])>0){number_BT_blood_all_clonotype[k,]<-qq[grep("blood",row.names(qq)),]}
      
    }
    
    
  }
  
  
  
  
  
  number_STB_skin_all_clonotype<-matrix(NA,length(unique(STB_sc_info_brief$cdr3_nt)),10)
  row.names(number_STB_skin_all_clonotype)<-paste0("clonotype_",1:length(unique(STB_sc_info_brief$cdr3_nt)))
  colnames(number_STB_skin_all_clonotype)<-c(0:9)
  
  
  number_STB_tumor_all_clonotype<-matrix(NA,length(unique(STB_sc_info_brief$cdr3_nt)),10)
  row.names(number_STB_tumor_all_clonotype)<-paste0("clonotype_",1:length(unique(STB_sc_info_brief$cdr3_nt)))
  colnames(number_STB_tumor_all_clonotype)<-c(0:9)
  
  
  number_STB_blood_all_clonotype<-matrix(NA,length(unique(STB_sc_info_brief$cdr3_nt)),10)
  row.names(number_STB_blood_all_clonotype)<-paste0("clonotype_",1:length(unique(STB_sc_info_brief$cdr3_nt)))
  colnames(number_STB_blood_all_clonotype)<-c(0:9)
  
  all_scST_match_bulk_STB<-data.frame()
  all_scSTB_match_bulk_STB<-data.frame()
  
  for(k in 1:length(unique(STB_sc_info_brief$cdr3_nt)))
  {
    cat('\r',k)
    mm<-STB_sc_info_brief$cdr3_nt[k]
    bulk_sequence<-data[grep(mm,data$rearrangement),"rearrangement"]
    xx<-STB_sc_info_brief[STB_sc_info_brief$cdr3_nt==unique(STB_sc_info_brief$cdr3_nt)[k],]
    
    skin_TCRAB<-unique(xx[xx$label=="skin",]$labeled_clonotype_id)
    if( length(skin_TCRAB)>0){
      for(m in 1:length(skin_TCRAB))
      {
        if(sum(grepl(unique(xx[xx$label=="skin",]$labeled_clonotype_id)[m],all_info_STB$ labeled_clonotype_id.x))!=0)
        {tag<-grep(unique(xx[xx$label=="skin",]$labeled_clonotype_id)[m],all_info_STB$ labeled_clonotype_id.x)
        all_TCRAB<-all_info_STB[tag,c("labeled_clonotype_id.x","labeled_clonotype_id.y","labeled_clonotype_id")]
        all_scSTB_match_bulk_STB<-rbind(all_scSTB_match_bulk_STB,all_TCRAB)
        barcode_TCRAB<-subset(xx,labeled_clonotype_id%in%all_TCRAB)$barcode
        
        CD8_STB1356@meta.data[,"bulk_TCR_match"]<-rep("other",nrow(CD8_STB1356@meta.data))
        tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_TCRAB
        CD8_STB1356@meta.data$bulk_TCR_match[tag1]<-paste0("Bulk_STB","_",CD8_STB1356@meta.data$label[tag1])
        if(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))>1){
          
          
          
          qq<-table(CD8_STB1356@meta.data$bulk_TCR_match,CD8_STB1356@meta.data$res.0.4_combined)
          if(length(qq[grep("skin",row.names(qq)),])>0){number_STB_skin_all_clonotype[k,1:10]<-qq[grep("skin",row.names(qq)),]}
          if(length(qq[grep("tumor",row.names(qq)),])>0){number_STB_tumor_all_clonotype[k,1:10]<-qq[grep("tumor",row.names(qq)),]}
          if(length(qq[grep("blood",row.names(qq)),])>0){number_STB_blood_all_clonotype[k,1:10]<-qq[grep("blood",row.names(qq)),]}
          
          myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/STB/TCRABmatch/bulk_TCR_STBmatched_with_singlecellSTB_",k,".",m,".pdf")
          pdf(myoutf,width=7.5,height=5)
          print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match',shape.by="PT_number",reduction="umap",
                        cols = c(col[1:(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))-1)],"#f0f0f0") ,pt.size = 0.3)+
                  labs(title=paste0("sc_",unique(all_info_STB[tag,"cdr3s_nt"]),"bulk_",unique(as.character(bulk_sequence))))
                +theme(plot.title = element_text(size=4)))
          dev.off()
        }
        
        
        
        }
      }}
    
    
    if( length(skin_TCRAB)>0){
      for(m in 1:length(skin_TCRAB))
      {
        if(sum(grepl(unique(xx[xx$label=="skin",]$labeled_clonotype_id)[m],all_info_ST$ labeled_clonotype_id.x))!=0)
        {tag<-grep(unique(xx[xx$label=="skin",]$labeled_clonotype_id)[m],all_info_ST$ labeled_clonotype_id.x)
        all_TCRAB<-all_info_ST[tag,c("labeled_clonotype_id.x","labeled_clonotype_id.y")]
        all_scST_match_bulk_STB<-rbind(all_scST_match_bulk_STB,all_TCRAB)
        barcode_TCRAB<-subset(xx,labeled_clonotype_id%in%all_TCRAB)$barcode
        
        CD8_STB1356@meta.data[,"bulk_TCR_match"]<-rep("other",nrow(CD8_STB1356@meta.data))
        tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_TCRAB
        CD8_STB1356@meta.data$bulk_TCR_match[tag1]<-paste0("Bulk_STB","_",CD8_STB1356@meta.data$label[tag1])
        if(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))>1){
          
          
          
          qq<-table(CD8_STB1356@meta.data$bulk_TCR_match,CD8_STB1356@meta.data$res.0.4_combined)
          if(length(qq[grep("skin",row.names(qq)),])>0){number_STB_skin_all_clonotype[k,1:10]<-qq[grep("skin",row.names(qq)),]}
          if(length(qq[grep("tumor",row.names(qq)),])>0){number_STB_tumor_all_clonotype[k,1:10]<-qq[grep("tumor",row.names(qq)),]}
          if(length(qq[grep("blood",row.names(qq)),])>0){number_STB_blood_all_clonotype[k,1:10]<-qq[grep("blood",row.names(qq)),]}
          
          myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/STB/TCRABmatch/bulk_TCR_STBmatched_with_singlecellST_",k,".",m,".pdf")
          pdf(myoutf,width=7.5,height=5)
          print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match',shape.by="PT_number",reduction="umap",
                        cols = c(col[1:(length(unique(CD8_STB1356@meta.data$bulk_TCR_match))-1)],"#f0f0f0") ,pt.size = 0.3)+
                  labs(title=paste0("sc_",unique(all_info_ST[tag,"cdr3s_nt"]),"bulk_",unique(as.character(bulk_sequence))))
                +theme(plot.title = element_text(size=4)))
          dev.off()
        }
        
        
        
        }
      }
      
      
      
      
      
      
    }}
  
  
}

tar2<-c("PT_6")
fraction_STB_skin_all_clonotypes<-number_STB_skin_all_clonotype<-na.omit(as.data.frame(number_STB_skin_all_clonotype))
fraction_STB_tumor_all_clonotypes<-number_STB_tumor_all_clonotype<-na.omit(as.data.frame(number_STB_tumor_all_clonotype))
fraction_STB_blood_all_clonotypes<-number_STB_blood_all_clonotype<-na.omit(as.data.frame(number_STB_blood_all_clonotype))

for (m in 1:nrow(fraction_STB_skin_all_clonotypes))
{
  c<-fraction_STB_skin_all_clonotypes[m,]
  cluster_number_use<-cluster_number[tar2[i],]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_STB_skin_all_clonotypes[m,1:10]<-fraction
}

for (m in 1:nrow(number_STB_tumor_all_clonotype))
{
  c<-fraction_STB_tumor_all_clonotypes[m,]
  cluster_number_use<-cluster_number[tar2[i],]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_STB_tumor_all_clonotypes[m,1:10]<-fraction
}

for (m in 1:nrow(number_STB_blood_all_clonotype))
{
  c<-fraction_STB_blood_all_clonotypes[m,]
  cluster_number_use<-cluster_number[tar2[i],]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_STB_blood_all_clonotypes[m,1:10]<-fraction
}




fraction_ST_skin_all_clonotypes<-number_ST_skin_all_clonotype<-na.omit(as.data.frame(number_ST_skin_all_clonotype))
fraction_ST_tumor_all_clonotypes<-number_ST_tumor_all_clonotype<-na.omit(as.data.frame(number_ST_tumor_all_clonotype))
fraction_ST_blood_all_clonotypes<-number_ST_blood_all_clonotype<-na.omit(as.data.frame(number_ST_blood_all_clonotype))

for (m in 1:nrow(number_ST_skin_all_clonotype))
{
  c<-fraction_ST_skin_all_clonotypes[m,]
  cluster_number_use<-cluster_number[tar2[i],]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_ST_skin_all_clonotypes[m,1:10]<-fraction
}

for (m in 1:nrow(number_ST_tumor_all_clonotype))
{
  c<-fraction_ST_tumor_all_clonotypes[m,]
  cluster_number_use<-cluster_number[tar2[i],]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_ST_tumor_all_clonotypes[m,1:10]<-fraction
}

if(nrow(fraction_ST_blood_all_clonotypes)>0){
  for (m in 1:nrow(number_ST_blood_all_clonotype))
  {
    c<-fraction_ST_blood_all_clonotypes[m,]
    cluster_number_use<-cluster_number[tar2[i],]
    fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
    fraction_ST_blood_all_clonotypes[m,1:10]<-fraction
  }}



fraction_BT_skin_all_clonotypes<-number_BT_skin_all_clonotype<-na.omit(as.data.frame(number_BT_skin_all_clonotype))
fraction_BT_tumor_all_clonotypes<-number_BT_tumor_all_clonotype<-na.omit(as.data.frame(number_BT_tumor_all_clonotype))
fraction_BT_blood_all_clonotypes<-number_BT_blood_all_clonotype<-na.omit(as.data.frame(number_BT_blood_all_clonotype))

if(nrow(fraction_BT_skin_all_clonotypes)>0){
  for (m in 1:nrow(number_BT_skin_all_clonotype))
  {
    c<-fraction_BT_skin_all_clonotypes[m,]
    cluster_number_use<-cluster_number[tar2[i],]
    fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
    fraction_BT_skin_all_clonotypes[m,1:10]<-fraction
  }}


for (m in 1:nrow(number_BT_tumor_all_clonotype))
{
  c<-fraction_BT_tumor_all_clonotypes[m,]
  cluster_number_use<-cluster_number[tar2[i],]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_BT_tumor_all_clonotypes[m,1:10]<-fraction
}

if(nrow(fraction_BT_blood_all_clonotypes)>0){
  for (m in 1:nrow(number_BT_blood_all_clonotype))
  {
    c<-fraction_BT_blood_all_clonotypes[m,]
    cluster_number_use<-cluster_number[tar2[i],]
    fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
    fraction_BT_blood_all_clonotypes[m,1:10]<-fraction
  }}



myoutf1 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/STB/TCRABmatch/hypergeometric_enrichment_STB_skin_all_clonotypes.xls")
myoutf2 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/STB/TCRABmatch/hypergeometric_enrichment_STB_tumor_all_clonotypes.xls")
myoutf3 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/STB/TCRABmatch/hypergeometric_enrichment_STB_blood_all_clonotypes.xls")


myoutf4 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/ST/hypergeometric_enrichment_ST_skin_all_clonotypes.xls")
myoutf5 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/ST/hypergeometric_enrichment_ST_tumor_all_clonotypes.xls")
myoutf6 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/ST/hypergeometric_enrichment_ST_blood_all_clonotypes.xls")


myoutf7 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/BT/hypergeometric_enrichment_BT_skin_all_clonotypes.xls")
myoutf8 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/BT/hypergeometric_enrichment_BT_tumor_all_clonotypes.xls")
myoutf9 <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/BT/hypergeometric_enrichment_BT_blood_all_clonotypes.xls")


write.table(na.omit (fraction_STB_skin_all_clonotypes),myoutf1,sep='\t')
write.table(na.omit (fraction_STB_tumor_all_clonotypes),myoutf2,sep='\t')
write.table(na.omit (fraction_STB_blood_all_clonotypes),myoutf3,sep='\t')

write.table(na.omit (fraction_ST_skin_all_clonotypes),myoutf4,sep='\t')
write.table(na.omit (fraction_ST_tumor_all_clonotypes),myoutf5,sep='\t')
write.table(na.omit (fraction_ST_blood_all_clonotypes),myoutf6,sep='\t')


write.table(na.omit (fraction_BT_skin_all_clonotypes),myoutf7,sep='\t')
write.table(na.omit (fraction_BT_tumor_all_clonotypes),myoutf8,sep='\t')
write.table(na.omit (fraction_BT_blood_all_clonotypes),myoutf9,sep='\t')





STB_bulk_match_sc_all<-data.frame()    
ST_bulk_match_sc_all<-data.frame()    
BT_bulk_match_sc_all<-data.frame()    


for(i in 1:length(tar))
{
  data<-long_term[grep(tar[i],long_term$sample_name),]
  data_skin<-data[grep("skin",data$sample_name),]
  data_tumor<-data[grep("melanoma",data$sample_name),]
  data_blood<-data[grep("PBMC",data$sample_name),]
  data_skin<-data_skin[data_skin$templates>3,]
  data_tumor<-data_tumor[data_tumor$templates>3,]
  data_blood<-data_blood[data_blood$templates>0,]
  
  
  
  ST<-merge(data_skin,data_tumor,by="rearrangement")
  BT<-merge(data_blood,data_tumor,by="rearrangement")
  STB<-merge(ST,data_blood,by="rearrangement")
  ST_spec<-ST[!ST$rearrangement%in%STB$rearrangement,]
  BT_spec<-BT[!BT$rearrangement%in%STB$rearrangement,]
  
  
  barcode_info_all<- data.frame(lapply(barcode_info_all, as.character), stringsAsFactors=FALSE)
  STB<- data.frame(lapply(STB, as.character), stringsAsFactors=FALSE)
  ST_spec<-data.frame(lapply(ST_spec, as.character), stringsAsFactors=FALSE)
  BT_spec<-data.frame(lapply(BT_spec, as.character), stringsAsFactors=FALSE)
  
  tag1<-unique(names(unlist( sapply(barcode_info_all$cdr3_nt, grep, STB$rearrangement ) )))
  STB_sc_info<-barcode_info_all[barcode_info_all$cdr3_nt%in%tag1,]
  
  tag2<-unique(names(unlist( sapply(barcode_info_all$cdr3_nt, grep, ST_spec$rearrangement ) )))
  ST_spec_sc_info<-barcode_info_all[barcode_info_all$cdr3_nt%in%tag2,]
  
  tag3<-unique(names(unlist( sapply(barcode_info_all$cdr3_nt, grep, BT_spec$rearrangement ) )))
  BT_spec_sc_info<-barcode_info_all[barcode_info_all$cdr3_nt%in%tag3,]
  
  STB_sc_info_brief<-unique(STB_sc_info[,c("barcode","cdr3_nt","label","labeled_clonotype_id")])
  ST_spec_sc_info_brief<-unique(ST_spec_sc_info[,c("barcode","cdr3_nt","label","labeled_clonotype_id")])
  BT_spec_sc_info_brief<-unique(BT_spec_sc_info[,c("barcode","cdr3_nt","label","labeled_clonotype_id")])
  
  
  myinf<-"../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
  CD8_STB1356<-readRDS(myinf)
  
  tag<-row.names(CD8_STB1356@meta.data)%in%STB_sc_info_brief$barcode
  CD8_STB1356@meta.data[,"bulk_TCR_match_all"]<-rep("other",nrow(CD8_STB1356@meta.data))
  CD8_STB1356@meta.data$bulk_TCR_match_all[tag]<-paste0("Bulk_STB","_",CD8_STB1356@meta.data$label[tag])
  myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/all_clonotype_bulk_TCR_STBmatched_.pdf")
  pdf(myoutf,width=7.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match_all',shape.by="PT_number",reduction="umap",
                cols = c('#2ca25f',"#de2d26","#9ecae1","#f0f0f0") ,pt.size = 0.3))
  dev.off()
  
  
  tag<-row.names(CD8_STB1356@meta.data)%in%ST_spec_sc_info_brief$barcode
  CD8_STB1356@meta.data[,"bulk_TCR_match_all"]<-rep("other",nrow(CD8_STB1356@meta.data))
  CD8_STB1356@meta.data$bulk_TCR_match_all[tag]<-paste0("Bulk_ST","_",CD8_STB1356@meta.data$label[tag])
  myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/all_clonotype_bulk_TCR_STmatched_.pdf")
  pdf(myoutf,width=7.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match_all',shape.by="PT_number",reduction="umap",
                cols = c('#2ca25f',"#de2d26","#9ecae1","#f0f0f0") ,pt.size = 0.3))
  dev.off()
  
  tag<-row.names(CD8_STB1356@meta.data)%in%BT_spec_sc_info_brief$barcode
  CD8_STB1356@meta.data[,"bulk_TCR_match_all"]<-rep("other",nrow(CD8_STB1356@meta.data))
  CD8_STB1356@meta.data$bulk_TCR_match_all[tag]<-paste0("Bulk_BT","_",CD8_STB1356@meta.data$label[tag])
  myoutf <- paste0("../scRNAseq/combine/RNA/vitiligo/STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_scTCR_new/",tar[i],"/all_clonotype_bulk_TCR_BTmatched_.pdf")
  pdf(myoutf,width=7.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='bulk_TCR_match_all',shape.by="PT_number",reduction="umap",
                cols = c('#2ca25f',"#de2d26","#9ecae1","#f0f0f0") ,pt.size = 0.3))
  dev.off()
  
}





