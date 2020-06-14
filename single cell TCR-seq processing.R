###############Extract the cell barcodes of single cell TCR sequencing###############
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

barcode_info_all_brief<-unique(barcode_info_all[,c("barcode","label","labeled_clonotype_id")])


#############Clonal match between tissues for all 4 patients and for each patient##############

barcode_info_all_brief_RNAseqed<-subset(barcode_info_all_brief,barcode%in%row.names(CD8_STB1356@meta.data))
barcode_info_all_brief_RNAseqed[,"frequency"]<-rep(0,nrow(barcode_info_all_brief_RNAseqed))
tar<-unique(barcode_info_all_brief_RNAseqed$labeled_clonotype_id)
tar<-tar[!is.na(tar)]

for (i in 1:length(tar))
{
  cat("\r",i)
  tag<-grep(tar[i],barcode_info_all_brief_RNAseqed$labeled_clonotype_id)
  barcode_info_all_brief_RNAseqed[tag,"frequency"]<-as.numeric(length(tag))
}

xx<-subset(barcode_info_all_brief_RNAseqed,label=="skin"|label=="tumor")
xx_expanded<-subset(xx,frequency>2)

barcode_info_all_brief_RNAseqed_expand<-rbind(xx_expanded,subset(barcode_info_all_brief_RNAseqed,label=="blood"))

xx<-strsplit(barcode_info_all_brief_RNAseqed_expand$barcode,"_")
barcode_info_all_brief_RNAseqed_expand[,"PT"]<-paste0("PT_",sapply(xx, "[[", 2))

unique_clonotypes<-unique(barcode_info_all_brief_RNAseqed_expand[,2:5])

table(unique_clonotypes$label,unique_clonotypes$PT)

#      PT_1 PT_3 PT_5 PT_6
#blood  780  636  253  123
#skin     7    3   89   61
#tumor   20   72   31   54


myoutf <- "../scRNAseq/combine/RNA/STB12356/CD8_STB12356/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/figure3A.pdf"
pdf(myoutf,width=10,height=10)
draw.triple.venn(area1 = 242, area2 = 429, area3 = 150, n12 = 33, n23 = 26, n13 = 19, 
                 n123 = 15,category = c("skin","tumor","blood"),
                 lwd = rep(2, 3), lty = rep("solid", 3), col =
                   c("#54278f","#fc4e2a","#fff7bc"), fill = c("#54278f","#fc4e2a","#fff7bc"), alpha = rep(0.5, 3))

dev.off()


myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/supplementary_figure6.pdf"
pdf(myoutf,width=15,height=15)

plot1<-draw.triple.venn(area1 = 780, area2 = 7, area3 = 20, n12 = 2, n23 = 1, n13 = 3, 
                        n123 = 1,category = c("blood","skin","tumor"),
                        lwd = rep(2, 3), lty = rep("solid", 3), col =
                          c("#54278f","#fc4e2a","#fff7bc"), fill = c("#54278f","#fc4e2a","#fff7bc"), alpha = rep(0.5, 3))

plot2<-draw.triple.venn(area1 = 636, area2 = 3, area3 = 72, n12 = 0, n23 = 0, n13 = 2, 
                        n123 = 0,category = c("blood","skin","tumor"),
                        lwd = rep(2, 3), lty = rep("solid", 3), col =
                          c("#54278f","#fc4e2a","#fff7bc"), fill = c("#54278f","#fc4e2a","#fff7bc"), alpha = rep(0.5, 3))

plot3<-draw.triple.venn(area1 = 253, area2 = 89, area3 = 31, n12 = 14, n23 = 13, n13 = 6, 
                        n123 = 6,category = c("blood","skin","tumor"),
                        lwd = rep(2, 3), lty = rep("solid", 3), col =
                          c("#54278f","#fc4e2a","#fff7bc"), fill = c("#54278f","#fc4e2a","#fff7bc"), alpha = rep(0.5, 3))


plot4<-draw.triple.venn(area1 = 123, area2 = 61, area3 = 54, n12 = 9, n23 = 19, n13 = 15, 
                        n123 = 8,category = c("blood","skin","tumor"),
                        lwd = rep(2, 3), lty = rep("solid", 3), col =
                          c("#54278f","#fc4e2a","#fff7bc"), fill = c("#54278f","#fc4e2a","#fff7bc"), alpha = rep(0.5, 3))

CombinePlots(plots = list(plot1, plot2,plot3,plot4))
dev.off()


##############Identify the matched clonotypes between different tissues############
tar<-c(1,3,5,6)
clonotype_STB_all_skin<-character()
clonotype_STB_all_tumor<-character()
clonotype_STB_all_blood<-character()

clonotype_ST_all_tumor<-character()
clonotype_ST_all_skin<-character()

clonotype_BT_all_tumor<-character()
clonotype_BT_all_blood<-character()

clonotype_BS_all_skin<-character()
clonotype_BS_all_blood<-character()

all_info_STB<-data.frame()

all_info_ST<-data.frame()

all_info_BT<-data.frame()

all_info_BS<-data.frame()

all_info_blood_only<-data.frame()

data_skin_all<-data.frame()
data_tumor_all<-data.frame()
data_blood_all<-data.frame()
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
  
  data_skin_all<-rbind(data_skin_all,data_skin)
  data_tumor_all<-rbind(data_tumor_all,data_tumor)
  data_blood_all<-rbind(data_blood_all,data_blood)
  
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
  
  BS_match<-merge(data_blood,data_skin,by.x="cdr3s_nt",by.y="cdr3s_nt")
  tag<-BS_match$cdr3s_nt%in%STB_match$cdr3s_nt
  BS_spec_match<-BS_match[!tag,]
  
  clonotype_BS_all_skin<-combine(clonotype_BS_all_skin,BS_spec_match$labeled_clonotype_id.y)
  clonotype_BS_all_blood<-combine( clonotype_BS_all_blood,BS_spec_match$labeled_clonotype_id.x)
  
  all_info_BS<-rbind(all_info_BS, BS_spec_match)
}

clonotype_data_all<-rbind(data_skin_all,data_blood_all,data_tumor_all)


tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB_all_skin
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB_all_tumor
tag3<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB_all_blood

barcode_STB_skin<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_STB_tumor<-unique(barcode_info_all_brief[tag2,"barcode"])
barcode_STB_blood<-unique(barcode_info_all_brief[tag3,"barcode"])


tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST_all_skin
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST_all_tumor

barcode_ST_skin<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_ST_tumor<-unique(barcode_info_all_brief[tag2,"barcode"])

tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT_all_blood
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT_all_tumor

barcode_BT_blood<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_BT_tumor<-unique(barcode_info_all_brief[tag2,"barcode"])

tag1<-barcode_info_all_brief$labeled_clonotype_id%in%all_info_blood_only$labeled_clonotype_id
barcode_blood_only<-unique(barcode_info_all_brief[tag1,"barcode"])


tag1<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BS_all_skin
tag2<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BS_all_blood

barcode_BS_skin<-unique(barcode_info_all_brief[tag1,"barcode"])
barcode_BS_blood<-unique(barcode_info_all_brief[tag2,"barcode"])



myinf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/Seruat_CD8_STB1356__finished.RDS"
CD8_STB1356<-readRDS(myinf)


tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_STB_skin
tag2<-row.names(CD8_STB1356@meta.data)%in%barcode_STB_tumor
tag3<-row.names(CD8_STB1356@meta.data)%in%barcode_STB_blood

CD8_STB1356@meta.data[,"STB_match"]<-rep("Not Matched",nrow(CD8_STB1356@meta.data))

CD8_STB1356@meta.data$STB_match[tag1]<-"STB_skin"
CD8_STB1356@meta.data$STB_match[tag2]<-"STB_tumor"
CD8_STB1356@meta.data$STB_match[tag3]<-"STB_blood"

table(CD8_STB1356@meta.data$STB_match,CD8_STB1356@meta.data$PT_number)
#              PT_1 PT_3 PT_5 PT_6
#Not Matched 1714 2572 2585 1730
#STB_blood     83   16   57  768
#STB_skin      25    7  185  159
#STB_tumor     41    9  428  279

tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_ST_skin
tag2<-row.names(CD8_STB1356@meta.data)%in%barcode_ST_tumor

CD8_STB1356@meta.data[,"ST_match"]<-rep("Not Matched",nrow(CD8_STB1356@meta.data))

CD8_STB1356@meta.data$ST_match[tag1]<-"ST_skin"
CD8_STB1356@meta.data$ST_match[tag2]<-"ST_tumor"

table(CD8_STB1356@meta.data$ST_match,CD8_STB1356@meta.data$PT_number)
#             PT_1 PT_3 PT_5 PT_6
#Not Matched 1767 2560 2837 2575
#ST_skin       31   18  225  210
#ST_tumor      65   26  193  151



tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_BT_blood
tag2<-row.names(CD8_STB1356@meta.data)%in%barcode_BT_tumor

CD8_STB1356@meta.data[,"BT_match"]<-rep("Not Matched",nrow(CD8_STB1356@meta.data))

CD8_STB1356@meta.data$BT_match[tag1]<-"BT_blood"
CD8_STB1356@meta.data$BT_match[tag2]<-"BT_tumor"

table(CD8_STB1356@meta.data$BT_match,CD8_STB1356@meta.data$PT_number)
#             PT_1 PT_3 PT_5 PT_6
#BT_blood      24   26    6  169
#BT_tumor      23   49    6   66
#Not Matched 1816 2529 3243 2701

tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_BS_blood
tag2<-row.names(CD8_STB1356@meta.data)%in%barcode_BS_skin

CD8_STB1356@meta.data[,"BS_match"]<-rep("Not Matched",nrow(CD8_STB1356@meta.data))

CD8_STB1356@meta.data$BS_match[tag1]<-"BS_blood"
CD8_STB1356@meta.data$BS_match[tag2]<-"BS_skin"

table(CD8_STB1356@meta.data$BS_match,CD8_STB1356@meta.data$PT_number)
#            PT_1 PT_3 PT_5 PT_6
#BS_blood     137  139   62   25
#BS_skin       26   20   96   23
#Not Matched 1700 2445 3097 2888

tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_blood_only

CD8_STB1356@meta.data[,"blood_only_nomatch"]<-rep("Not Matched",nrow(CD8_STB1356@meta.data))

CD8_STB1356@meta.data$blood_only_nomatch[tag1]<-"Blood only"


table(CD8_STB1356@meta.data$blood_only_nomatch,CD8_STB1356@meta.data$PT_number)

PT_1 PT_3 PT_5 PT_6
#Blood only   760  845  216   86
#Not Matched 1103 1759 3039 2850


myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_STB.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = CD8_STB1356,group.by='STB_match',shape.by="PT_number",reduction="umap",
        cols = c('#f0f0f0','#2ca25f',"#de2d26","#9ecae1"),pt.size = 0.3)
dev.off()

myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_ST.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = CD8_STB1356,group.by='ST_match',shape.by="PT_number",reduction="umap",
        cols = c('#f0f0f0',"#de2d26","#9ecae1"),pt.size = 0.3)
dev.off()


myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_BT.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = CD8_STB1356,group.by='BT_match',shape.by="PT_number",reduction="umap",
        cols = c('#2ca25f',"#9ecae1",'#f0f0f0'),pt.size = 0.3)
dev.off()


myoutf <- "../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/seperate_STB_ST_BT/UMAP_clonotypes_Blood_only_no_match.pdf"
pdf(myoutf,width=7.5,height=5)
DimPlot(object = CD8_STB1356,group.by='blood_only_nomatch',shape.by="PT_number",reduction="umap",
        cols = c('#2ca25f','#f0f0f0'),pt.size = 0.3)
dev.off()


tar<-unique(CD8_STB1356@meta.data$PT_number)
for(i in 1:length(tar))
{
  cat("\t",i)
  Idents(CD8_STB1356)<-"PT_number"
  CD8_STB1356_PT<-subset(CD8_STB1356,idents=tar[i])
  myoutf <- paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/seperate_STB_ST_BT/UMAP_matched_clonotypes_STB_",tar[i],"_.pdf")
  pdf(myoutf,width=7.5,height=5)
  
  print(DimPlot(object = CD8_STB1356_PT,group.by='STB_match',shape.by="label",reduction="umap",
                cols = c('#f0f0f0','#2ca25f',"#de2d26","#9ecae1"),pt.size = 0.3))
  
  dev.off()
  
}


tar<-unique(CD8_STB1356@meta.data$PT_number)
for (i in 1:length(tar))
{
  data<-CD8_STB1356@meta.data[CD8_STB1356@meta.data$PT_number==tar[i],]
  blood_barcode<-row.names(data[data$STB_match=="STB_blood",])
  tumor_barcode<-row.names(data[data$STB_match=="STB_tumor",])
  skin_barcode<-row.names(data[data$STB_match=="STB_skin",])
  
  tag1<-barcode_info_all_brief$barcode%in%skin_barcode
  tag2<-barcode_info_all_brief$barcode%in%tumor_barcode
  tag3<-barcode_info_all_brief$barcode%in%blood_barcode
  
  skin_clone<-unique(barcode_info_all_brief[tag1,"labeled_clonotype_id"])
  tumor_clone<-unique(barcode_info_all_brief[tag2,"labeled_clonotype_id"])
  blood_clone<-unique(barcode_info_all_brief[tag3,"labeled_clonotype_id"])
  
  STB_match_number<-sum(clonotype_STB_all_skin%in%skin_clone&clonotype_STB_all_tumor%in%tumor_clone&clonotype_STB_all_skin%in%skin_clone&clonotype_STB_all_blood%in%blood_clone)
  print(paste0(tar[i],' ',"STB_number",' ',STB_match_number))
  
  
  tumor_barcode<-row.names(data[data$ST_match=="ST_tumor",])
  skin_barcode<-row.names(data[data$ST_match=="ST_skin",])
  
  tag1<-barcode_info_all_brief$barcode%in%skin_barcode
  tag2<-barcode_info_all_brief$barcode%in%tumor_barcode
  
  skin_clone<-unique(barcode_info_all_brief[tag1,"labeled_clonotype_id"])
  tumor_clone<-unique(barcode_info_all_brief[tag2,"labeled_clonotype_id"])
  
  ST_match_number<-sum(clonotype_ST_all_skin%in%skin_clone&clonotype_ST_all_tumor%in%tumor_clone)
  print(paste0(tar[i],' ',"ST_number",' ',ST_match_number))
  
  
  
  tumor_barcode<-row.names(data[data$BT_match=="BT_tumor",])
  blood_barcode<-row.names(data[data$BT_match=="BT_blood",])
  
  tag1<-barcode_info_all_brief$barcode%in%blood_barcode
  tag2<-barcode_info_all_brief$barcode%in%tumor_barcode
  
  blood_clone<-unique(barcode_info_all_brief[tag1,"labeled_clonotype_id"])
  tumor_clone<-unique(barcode_info_all_brief[tag2,"labeled_clonotype_id"])
  
  BT_match_number<-sum(clonotype_BT_all_blood%in%blood_clone&clonotype_BT_all_tumor%in%tumor_clone)
  print(paste0(tar[i],' ',"BT_number",' ' ,BT_match_number))
  
}





#############Plot each tumor-associated clonotype##############
CD8_STB1356@meta.data[,"match_TCR_nt"]<-rep("Not_match",nrow(CD8_STB1356@meta.data))
CD8_STB1356@meta.data[,"match_TCR_aa"]<-rep("Not_match",nrow(CD8_STB1356@meta.data))
CD8_STB1356@meta.data[,"TCR_Match_RNAseqed"]<-rep("None",nrow(CD8_STB1356@meta.data))


#Plot single clonotypes that matches between skin and blood and expanded in skin, frequency>2


barcode_BS_blood<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$BS_match=="BS_blood",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_BS_blood
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>0
blood_expanded<-clonotype[tag,]



barcode_BS_skin<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$BS_match=="BS_skin",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_BS_skin
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
skin_expanded<-clonotype[tag,]




tag<-all_info_BS$labeled_clonotype_id.x%in%blood_expanded$labeled_clonotype_id&all_info_BS$labeled_clonotype_id.y%in%skin_expanded$labeled_clonotype_id

all_info_BS_skin_expanded<-all_info_BS[tag,]



col<-c("#de2d26","#9ecae1",'#f0f0f0')

for(i in 1:nrow(all_info_BS_skin_expanded))
{
  cat("\r",i)
  clonotype_SB<-combine(all_info_BS_skin_expanded$labeled_clonotype_id.x[i],
                        all_info_BS_skin_expanded$labeled_clonotype_id.y[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_SB
  barcode_SB<-barcode_info_all_brief$barcode[tag]
  barcode_SB_skin<-barcode_SB[grep("skin",barcode_SB)]
  barcode_SB_blood<-barcode_SB[grep("blood",barcode_SB)]
  
  
  CD8_STB1356@meta.data[,"single_clonotype"]<-rep("Others",nrow( CD8_STB1356@meta.data))
  tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_SB_skin
  tag2<-row.names(CD8_STB1356@meta.data)%in%barcode_SB_blood
  
  CD8_STB1356@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_Skin")
  CD8_STB1356@meta.data$single_clonotype[tag2]<-paste0("Clonotype",i,"_Blood")
  
  CD8_STB1356@meta.data$TCR_Match_RNAseqed[tag1|tag2]<-"Skin Blood Match"
  
  skin_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==all_info_BS_skin_expanded$labeled_clonotype_id.y[i],"cdr3s_nt"])
  skin_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==all_info_BS_skin_expanded$labeled_clonotype_id.y[i],"cdr3s_aa"])
  blood_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==all_info_BS_skin_expanded$labeled_clonotype_id.x[i],"cdr3s_nt"])
  blood_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==all_info_BS_skin_expanded$labeled_clonotype_id.x[i],"cdr3s_aa"])
  
  CD8_STB1356@meta.data$match_TCR_nt[tag1]<-skin_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag1]<-skin_seq_aa
  CD8_STB1356@meta.data$match_TCR_nt[tag2]<-blood_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag2]<-blood_seq_aa
  
  
  
  myoutf <- paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BS_skin_expanded/UMAP_matched_clonotype",i,".pdf")
  pdf(myoutf,width=6.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='single_clonotype',shape.by="PT_number",reduction="umap",
                cols = col,pt.size = 0.2))
  dev.off()
  
  
  
}


#Plot single clonotypes that matches between skin tumor blood and expanded in both skin and tumor, frequency>2

barcode_STB_blood<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$STB_match=="STB_blood",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_blood
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>0
blood_expanded<-clonotype[tag,]



barcode_STB_tumor<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$STB_match=="STB_tumor",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_tumor
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
tumor_expanded<-clonotype[tag,]



barcode_STB_skin<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$STB_match=="STB_skin",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_STB_skin
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))

for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
skin_expanded<-clonotype[tag,]

tag<-all_info_STB$labeled_clonotype_id%in%blood_expanded$labeled_clonotype_id&all_info_STB$labeled_clonotype_id.x%in%skin_expanded$labeled_clonotype_id&all_info_STB$labeled_clonotype_id.y%in%tumor_expanded$labeled_clonotype_id

all_info_STB_ST_expanded<-all_info_STB[tag,]



enrich_skin_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(enrich_skin_all_clonotypes)<-paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(enrich_skin_all_clonotypes)<-seq(0,9,1)

enrich_tumor_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(enrich_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(enrich_tumor_all_clonotypes)<-seq(0,9,1)

enrich_blood_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(enrich_blood_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(enrich_blood_all_clonotypes)<-seq(0,9,1)


number_skin_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(number_skin_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(number_skin_all_clonotypes)<-seq(0,9,1)


number_tumor_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(number_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(number_tumor_all_clonotypes)<-seq(0,9,1)


number_blood_all_clonotypes<-matrix(0,nrow(all_info_STB_ST_expanded),10)
row.names(number_blood_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_STB_ST_expanded),1))
colnames(number_blood_all_clonotypes)<-seq(0,9,1)

PT_all<-character()




col<-c('#F3766E',"#2AB34B","#7094CD",'#C0C0C0')

for(i in 1:nrow(all_info_STB_ST_expanded))
{
  cat("\r",i)
  clonotype_STB<-combine(all_info_STB_ST_expanded$labeled_clonotype_id.x[i],
                         all_info_STB_ST_expanded$labeled_clonotype_id.y[i],
                         all_info_STB_ST_expanded$labeled_clonotype_id[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_STB
  barcode_STB<-barcode_info_all_brief$barcode[tag]
  barcode_STB_skin<-barcode_STB[grep("skin",barcode_STB)]
  barcode_STB_blood<-barcode_STB[grep("blood",barcode_STB)]
  barcode_STB_tumor<-barcode_STB[grep("tumor",barcode_STB)]
  
  
  CD8_STB1356@meta.data[,"single_clonotype"]<-rep("Others",nrow( CD8_STB1356@meta.data))
  tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_STB_skin
  tag2<-row.names(CD8_STB1356@meta.data)%in%barcode_STB_blood
  tag3<-row.names(CD8_STB1356@meta.data)%in%barcode_STB_tumor
  
  CD8_STB1356@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_Skin")
  CD8_STB1356@meta.data$single_clonotype[tag2]<-paste0("Clonotype",i,"_Blood")
  CD8_STB1356@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
  
  xx<-CD8_STB1356@meta.data[tag1,]
  yy<-CD8_STB1356@meta.data[tag3,]
  zz<-CD8_STB1356@meta.data[tag2,]
  
  CD8_STB1356@meta.data$TCR_Match_RNAseqed[tag1|tag2|tag3]<-"Skin Tumor Blood Match"
  
  
  skin_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_STB[1],"cdr3s_nt"])
  skin_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_STB[1],"cdr3s_aa"])
  blood_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_STB[3],"cdr3s_nt"])
  blood_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_STB[3],"cdr3s_aa"])
  tumor_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_STB[2],"cdr3s_nt"])
  tumor_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_STB[2],"cdr3s_aa"])
  
  CD8_STB1356@meta.data$match_TCR_nt[tag1]<-skin_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag1]<-skin_seq_aa
  CD8_STB1356@meta.data$match_TCR_nt[tag2]<-blood_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag2]<-blood_seq_aa
  CD8_STB1356@meta.data$match_TCR_nt[tag3]<-tumor_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag3]<-tumor_seq_aa
  
  
  
  chi_STB_skin<-matrix(0,2,10)
  row.names(chi_STB_skin)=c("match","not match")
  colnames(chi_STB_skin)<-seq(0,9,1)
  
  
  chi_STB_tumor<-matrix(0,2,10)
  row.names(chi_STB_tumor)=c("match","not match")
  colnames(chi_STB_tumor)<-seq(0,9,1)
  
  chi_STB_blood<-matrix(0,2,10)
  row.names(chi_STB_blood)=c("match","not match")
  colnames(chi_STB_blood)<-seq(0,9,1)
  
  
  for (j in 1:length(unique(xx$res.0.4_combined)))
  {
    chi_STB_skin["match",unique(xx$res.0.4_combined)[j]]=sum(xx$res.0.4_combined==unique(xx$res.0.4_combined)[j])
    chi_STB_skin["not match",unique(xx$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="skin",])$res.0.4_combined==unique(xx$res.0.4_combined)[j])-sum(xx$res.0.4_combined==unique(xx$res.0.4_combined)[j])
  }
  
  
  
  for (j in 1:length(unique(yy$res.0.4_combined)))
  {
    chi_STB_tumor["match",unique(yy$res.0.4_combined)[j]]=sum(yy$res.0.4_combined==unique(yy$res.0.4_combined)[j])
    chi_STB_tumor["not match",unique(yy$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="tumor",])$res.0.4_combined==unique(yy$res.0.4_combined)[j])-sum(yy$res.0.4_combined==unique(yy$res.0.4_combined)[j])
  }
  
  
  for (j in 1:length(unique(zz$res.0.4_combined)))
  {
    chi_STB_blood["match",unique(zz$res.0.4_combined)[j]]=sum(zz$res.0.4_combined==unique(zz$res.0.4_combined)[j])
    chi_STB_blood["not match",unique(zz$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="blood",])$res.0.4_combined==unique(zz$res.0.4_combined)[j])-sum(zz$res.0.4_combined==unique(zz$res.0.4_combined)[j])
  }
  
  
  chi_skin<-chisq.test(chi_STB_skin)
  chi_observed<-chi_skin$observed
  chi_expected<-chi_skin$expected
  ROE_skin<-chi_observed/chi_expected
  
  
  
  chi_tumor<-chisq.test(chi_STB_tumor)
  chi_observed<-chi_tumor$observed
  chi_expected<-chi_tumor$expected
  ROE_tumor<-chi_observed/chi_expected
  
  enrich_skin_all_clonotypes[i,1:10]<-ROE_skin[1,]
  enrich_tumor_all_clonotypes[i,1:10]<-ROE_tumor[1,]
  
  number_skin_all_clonotypes[i,1:10]<-chi_STB_skin[1,]
  number_tumor_all_clonotypes[i,1:10]<-chi_STB_tumor[1,]
  number_blood_all_clonotypes[i,1:10]<-chi_STB_blood[1,]
  
  number_skin_all_clonotypes<-as.data.frame(number_skin_all_clonotypes)
  number_tumor_all_clonotypes<-as.data.frame(number_tumor_all_clonotypes)
  number_blood_all_clonotypes<-as.data.frame(number_blood_all_clonotypes)
  
  
  PT_all<-combine(PT_all,unique(xx$PT_number))
  
  myoutf <- paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/UMAP_matched_clonotype",i,".pdf")
  pdf(myoutf,width=6.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='single_clonotype',pt.shape="16",reduction="umap",
                cols = col,pt.size = 0.2))
  dev.off()
  
  if(unique(CD8_STB1356@meta.data$single_clonotype[tag1])!="Others"&CD8_STB1356@meta.data$bulk_STB_match[tag1]!="Not Matched"){
    myoutf <- paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/bulk_TCR_match_with_scTCR/single_clonotype_635_matched_bulkscTCR_UMAP_matched_clonotype",i,".pdf")
    pdf(myoutf,width=6.5,height=5)
    print(DimPlot(object = CD8_STB1356,group.by='single_clonotype',reduction="umap",pt.shape="16",
                  cols = col,pt.size = 0.2))
    dev.off()
  }
  
  
}

number_skin_all_clonotypes<-cbind(number_skin_all_clonotypes,PT_all)
number_tumor_all_clonotypes<-cbind(number_tumor_all_clonotypes,PT_all)
number_blood_all_clonotypes<-cbind(number_blood_all_clonotypes,PT_all)


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/number_of_cells_skin.xls"
write.table(number_skin_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/number_of_cells_tumor.xls"
write.table(number_tumor_all_clonotypes,myoutf,sep='\t')
#Calculate the enrichment of every clonotype in each cluster

cluster_number<-matrix(0,4,10)
row.names(cluster_number)<-paste0("PT_",c(1,3,5,6))
colnames(cluster_number)<-seq(0,9,1)
for(i in 1:ncol(cluster_number))
{
  tag<-CD8_STB1356@meta.data$res.0.4_combined==colnames(cluster_number)[i]
  xx<-CD8_STB1356@meta.data[tag,]
  cluster_number[1,i]<-sum(xx$PT_number=="PT_1")
  cluster_number[2,i]<-sum(xx$PT_number=="PT_3")
  cluster_number[3,i]<-sum(xx$PT_number=="PT_5")
  cluster_number[4,i]<-sum(xx$PT_number=="PT_6")
}


fraction_skin_all_clonotypes<-number_skin_all_clonotypes
fraction_tumor_all_clonotypes<-number_tumor_all_clonotypes
fraction_blood_all_clonotypes<-number_blood_all_clonotypes

for (i in 1:nrow(number_skin_all_clonotypes))
{
  c<-number_skin_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_skin_all_clonotypes[i,1:10]<-fraction
}


for (i in 1:nrow(number_tumor_all_clonotypes))
{
  c<-number_tumor_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_tumor_all_clonotypes[i,1:10]<-fraction
}

for (i in 1:nrow(number_blood_all_clonotypes))
{
  c<-number_blood_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_blood_all_clonotypes[i,1:10]<-fraction
}

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/skin_enrich_hypergeometric_clustersize_eachPT_each_cluster.xls"
write.table(fraction_skin_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/blood_enrich_hypergeometric_clustersize_eachPT_each_cluster.xls"
write.table(fraction_blood_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/tumor_enrich_hypergeometric_clustersize_eachPT_each_cluster.xls"
write.table(fraction_tumor_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/skin_enrich_sig_geometric_clustersize_eachPT_each_cluster.xls"
write.table(fraction_skin_all_clonotypes_sig,myoutf,sep='\t')



ColourCount = 15
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

enrich_skin_plot<-enrich_all_skin[as.character(enrich_all_skin$cluster)%in%c("0","1","3","6","7"),]
enrich_tumor_plot<-enrich_all_tumor[as.character(enrich_all_tumor$cluster)%in%as.character(c(0,1,3,6,7)),]



myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_cluster_in_skin_tumor_skin.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_skin_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_skin_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_cluster_in_skin_tumor_tumor.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_tumor_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_tumor_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()
myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_each_cluster.pdf"
pdf(myoutf,width=15,height=30)
print(ggplot(data=enrich_all, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_each_cluster_skin.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_skin, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


t.test(enrich_all_skin$index[1:15],enrich_all_skin$index[16:30],alternative = "less")


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/STB_tumor_skin_expanded/line_plot_enrichment_each_cluster_tumor.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_tumor, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +geom_hline(yintercept=1, linetype="dashed", 
                  color = "red", size=2)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))
dev.off()


t.test(enrich_all_tumor$index[1:15],enrich_all_tumor$index[16:30],alternative = "less")

#Plot single clonotypes that matches between skin tumor only and expanded in both skin and tumor

barcode_ST_skin<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$ST_match=="ST_skin",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_ST_skin
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
skin_expanded<-clonotype[tag,]


barcode_ST_tumor<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$ST_match=="ST_tumor",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_ST_tumor
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))

for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
tumor_expanded<-clonotype[tag,]

tag<-all_info_ST$labeled_clonotype_id.y%in%tumor_expanded$labeled_clonotype_id & all_info_ST$labeled_clonotype_id.x%in%skin_expanded$labeled_clonotype_id

all_info_ST_ST_expanded<- all_info_ST[tag,]



enrich_skin_all_clonotypes<-matrix(0,nrow(all_info_ST_ST_expanded),10)
row.names(enrich_skin_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1))
colnames(enrich_skin_all_clonotypes)<-seq(0,9,1)

enrich_tumor_all_clonotypes<-matrix(0,nrow(all_info_ST_ST_expanded),10)
row.names(enrich_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1))
colnames(enrich_tumor_all_clonotypes)<-seq(0,9,1)


number_skin_all_clonotypes<-matrix(0,nrow(all_info_ST_ST_expanded),10)
row.names(number_skin_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1))
colnames(number_skin_all_clonotypes)<-seq(0,9,1)


number_tumor_all_clonotypes<-matrix(0,nrow(all_info_ST_ST_expanded),10)
row.names(number_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_ST_ST_expanded),1))
colnames(number_tumor_all_clonotypes)<-seq(0,9,1)




PT_all<-character()


col<-c("#2AB34B","#7094CD",'#C0C0C0')

for(i in 1:nrow(all_info_ST_ST_expanded))
{
  cat("\r",i)
  clonotype_ST<-combine(all_info_ST_ST_expanded$labeled_clonotype_id.x[i],
                        all_info_ST_ST_expanded$labeled_clonotype_id.y[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_ST
  barcode_ST<-barcode_info_all_brief$barcode[tag]
  barcode_ST_skin<-barcode_ST[grep("skin",barcode_ST)]
  
  barcode_ST_tumor<-barcode_ST[grep("tumor",barcode_ST)]
  
  
  CD8_STB1356@meta.data[,"single_clonotype"]<-rep("Others",nrow( CD8_STB1356@meta.data))
  tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_ST_skin
  tag3<-row.names(CD8_STB1356@meta.data)%in%barcode_ST_tumor
  
  CD8_STB1356@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_Skin")
  CD8_STB1356@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
  
  xx<-CD8_STB1356@meta.data[tag1,]
  yy<-CD8_STB1356@meta.data[tag3,]
  
  CD8_STB1356@meta.data$TCR_Match_RNAseqed[tag1|tag3]<-"Skin Tumor Only Match"
  
  
  skin_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_ST[1],"cdr3s_nt"])
  skin_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_ST[1],"cdr3s_aa"])
  
  tumor_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_ST[2],"cdr3s_nt"])
  tumor_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_ST[2],"cdr3s_aa"])
  
  CD8_STB1356@meta.data$match_TCR_nt[tag1]<-skin_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag1]<-skin_seq_aa
  
  CD8_STB1356@meta.data$match_TCR_nt[tag3]<-tumor_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag3]<-tumor_seq_aa

  
  chi_ST_skin<-matrix(0,2,10)
  row.names(chi_ST_skin)=c("match","not match")
  colnames(chi_ST_skin)<-seq(0,9,1)
  
  
  chi_ST_tumor<-matrix(0,2,10)
  row.names(chi_ST_tumor)=c("match","not match")
  colnames(chi_ST_tumor)<-seq(0,9,1)
  
  for (j in 1:length(unique(xx$res.0.4_combined)))
  {
    chi_ST_skin["match",unique(xx$res.0.4_combined)[j]]=sum(xx$res.0.4_combined==unique(xx$res.0.4_combined)[j])
    chi_ST_skin["not match",unique(xx$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="skin",])$res.0.4_combined==unique(xx$res.0.4_combined)[j])-sum(xx$res.0.4_combined==unique(xx$res.0.4_combined)[j])
  }
  
  
  
  for (j in 1:length(unique(yy$res.0.4_combined)))
  {
    chi_ST_tumor["match",unique(yy$res.0.4_combined)[j]]=sum(yy$res.0.4_combined==unique(yy$res.0.4_combined)[j])
    chi_ST_tumor["not match",unique(yy$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="tumor",])$res.0.4_combined==unique(yy$res.0.4_combined)[j])-sum(yy$res.0.4_combined==unique(yy$res.0.4_combined)[j])
  }
  
  
  chi_skin<-chisq.test(chi_ST_skin)
  chi_observed<-chi_skin$observed
  chi_expected<-chi_skin$expected
  ROE_skin<-chi_observed/chi_expected
  
  
  
  chi_tumor<-chisq.test(chi_ST_tumor)
  chi_observed<-chi_tumor$observed
  chi_expected<-chi_tumor$expected
  ROE_tumor<-chi_observed/chi_expected
  
  enrich_skin_all_clonotypes[i,]<-ROE_skin[1,]
  enrich_tumor_all_clonotypes[i,]<-ROE_tumor[1,]
  
  
  number_skin_all_clonotypes[i,1:10]<-chi_ST_skin[1,]
  number_tumor_all_clonotypes[i,1:10]<-chi_ST_tumor[1,]
  
  number_skin_all_clonotypes<-as.data.frame(number_skin_all_clonotypes)
  number_tumor_all_clonotypes<-as.data.frame(number_tumor_all_clonotypes)
  
  PT_all<-combine(PT_all,unique(xx$PT_number))
  
  
  myoutf <- paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/UMAP_matched_clonotype",i,".pdf")
  pdf(myoutf,width=6.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='single_clonotype',reduction="umap",pt.shape="16",
                cols = col,pt.size = 0.2))
  dev.off()
  
}


number_skin_all_clonotypes<-cbind(number_skin_all_clonotypes,PT_all)
number_tumor_all_clonotypes<-cbind(number_tumor_all_clonotypes,PT_all)



myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/number_of_cells_skin.xls"
write.table(number_skin_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/number_of_cells_tumor.xls"
write.table(enrich_tumor_all_clonotypes,myoutf,sep='\t')


#Calculate the enrichment of every clonotype in each cluster

cluster_number<-matrix(0,4,10)
row.names(cluster_number)<-paste0("PT_",c(1,3,5,6))
colnames(cluster_number)<-seq(0,9,1)
for(i in 1:ncol(cluster_number))
{
  tag<-CD8_STB1356@meta.data$res.0.4_combined==colnames(cluster_number)[i]
  xx<-CD8_STB1356@meta.data[tag,]
  cluster_number[1,i]<-sum(xx$PT_number=="PT_1")
  cluster_number[2,i]<-sum(xx$PT_number=="PT_3")
  cluster_number[3,i]<-sum(xx$PT_number=="PT_5")
  cluster_number[4,i]<-sum(xx$PT_number=="PT_6")
}


fraction_skin_all_clonotypes<-number_skin_all_clonotypes
fraction_tumor_all_clonotypes<-number_tumor_all_clonotypes

for (i in 1:nrow(number_skin_all_clonotypes))
{
  c<-number_skin_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_skin_all_clonotypes[i,1:10]<-fraction
}


for (i in 1:nrow(number_tumor_all_clonotypes))
{
  c<-number_tumor_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_tumor_all_clonotypes[i,1:10]<-fraction
}




myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/skin_enrich_hypergeometric_clustersize_eachPT_each_cluster.txt"
write.table(fraction_skin_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/tumor_enrich_hypergeometric_clustersize_eachPT_each_cluster.txt"
write.table(fraction_tumor_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/skin_enrich_sig_geometric_clustersize_eachPT_each_cluster.xls"
write.table(fraction_skin_all_clonotypes_sig,myoutf,sep='\t')




colourCount = 18
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_each_cluster.pdf"
pdf(myoutf,width=15,height=30)
print(ggplot(data=enrich_all, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_each_cluster_skin.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_skin, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


t.test(enrich_all_skin$index[1:15],enrich_all_skin$index[16:30],alternative = "less")


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_each_cluster_tumor.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_tumor, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +geom_hline(yintercept=1, linetype="dashed", 
                  color = "red", size=2)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))
dev.off()


t.test(enrich_all_tumor$index[1:15],enrich_all_tumor$index[16:30],alternative = "less")


enrich_skin_plot<-enrich_all_skin[as.character(enrich_all_skin$cluster)%in%c("0","1","3","6","7"),]
enrich_tumor_plot<-enrich_all_tumor[as.character(enrich_all_tumor$cluster)%in%as.character(c(0,1,3,6,7)),]



myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_cluster_in_skin_tumor_skin.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_skin_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_skin_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/ST_skin_tumor_expanded/line_plot_enrichment_cluster_in_skin_tumor_tumor.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_tumor_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_tumor_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


#Plot single clonotypes that matches between blood and tumor only and expanded in both blood and tumor##

barcode_BT_blood<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$BT_match=="BT_blood",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_BT_blood
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))
for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>0
blood_expanded<-clonotype[tag,]


barcode_BT_tumor<-row.names(CD8_STB1356@meta.data[CD8_STB1356@meta.data$BT_match=="BT_tumor",])
tag<-unique(barcode_info_all_brief)$barcode%in%barcode_BT_tumor
clonotype<-barcode_info_all_brief[tag,]
clonotype[,"frequency"]<-rep(0,nrow(clonotype))

for(i in 1:length(unique(clonotype$labeled_clonotype_id)))
{
  tag<-clonotype$labeled_clonotype_id==unique(clonotype$labeled_clonotype_id)[i]
  xx<-clonotype[tag,]
  clonotype$frequency[tag]=nrow(xx)
  
}

tag<-clonotype$frequency>2
tumor_expanded<-clonotype[tag,]

tag<-all_info_BT$labeled_clonotype_id.y%in%tumor_expanded$labeled_clonotype_id & all_info_BT$labeled_clonotype_id.x%in%blood_expanded$labeled_clonotype_id

all_info_BT_BT_expanded<- all_info_BT[tag,]



enrich_blood_all_clonotypes<-matrix(0,nrow(all_info_BT_BT_expanded),10)
row.names(enrich_blood_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1))
colnames(enrich_blood_all_clonotypes)<-seq(0,9,1)

enrich_tumor_all_clonotypes<-matrix(0,nrow(all_info_BT_BT_expanded),10)
row.names(enrich_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1))
colnames(enrich_tumor_all_clonotypes)<-seq(0,9,1)





number_blood_all_clonotypes<-matrix(0,nrow(all_info_BT_BT_expanded),10)
row.names(number_blood_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1))
colnames(number_blood_all_clonotypes)<-seq(0,9,1)


number_tumor_all_clonotypes<-matrix(0,nrow(all_info_BT_BT_expanded),10)
row.names(number_tumor_all_clonotypes)=paste0("clonotype_",seq(1,nrow(all_info_BT_BT_expanded),1))
colnames(number_tumor_all_clonotypes)<-seq(0,9,1)

PT_all<-character()


col<-c('#F3766E',"#7094CD",'#C0C0C0')

for(i in 1:nrow(all_info_BT_BT_expanded))
{
  cat("\r",i)
  clonotype_BT<-combine(all_info_BT_BT_expanded$labeled_clonotype_id.x[i],
                        all_info_BT_BT_expanded$labeled_clonotype_id.y[i])
  
  
  tag<-barcode_info_all_brief$labeled_clonotype_id%in%clonotype_BT
  barcode_BT<-barcode_info_all_brief$barcode[tag]
  barcode_BT_blood<-barcode_BT[grep("blood",barcode_BT)]
  
  barcode_BT_tumor<-barcode_BT[grep("tumor",barcode_BT)]
  
  
  CD8_STB1356@meta.data[,"single_clonotype"]<-rep("Others",nrow( CD8_STB1356@meta.data))
  tag1<-row.names(CD8_STB1356@meta.data)%in%barcode_BT_blood
  tag3<-row.names(CD8_STB1356@meta.data)%in%barcode_BT_tumor
  
  CD8_STB1356@meta.data$single_clonotype[tag1]<-paste0("Clonotype",i,"_blood")
  CD8_STB1356@meta.data$single_clonotype[tag3]<-paste0("Clonotype",i,"_Tumor")
  
  xx<-CD8_STB1356@meta.data[tag1,]
  yy<-CD8_STB1356@meta.data[tag3,]
  
  CD8_STB1356@meta.data$TCR_Match_RNAseqed[tag1|tag3]<-"Blood Tumor Only Match"
  
  blood_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_BT[1],"cdr3s_nt"])
  blood_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_BT[1],"cdr3s_aa"])
  
  tumor_seq_nt<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_ST[2],"cdr3s_nt"])
  tumor_seq_aa<-as.character(clonotype_data_all[clonotype_data_all$labeled_clonotype_id==clonotype_ST[2],"cdr3s_aa"])
  
  CD8_STB1356@meta.data$match_TCR_nt[tag1]<-blood_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag1]<-blood_seq_aa
  
  CD8_STB1356@meta.data$match_TCR_nt[tag3]<-tumor_seq_nt
  CD8_STB1356@meta.data$match_TCR_aa[tag3]<-tumor_seq_aa
  
  chi_BT_blood<-matrix(0,2,10)
  row.names(chi_BT_blood)=c("match","not match")
  colnames(chi_BT_blood)<-seq(0,9,1)
  
  
  chi_BT_tumor<-matrix(0,2,10)
  row.names(chi_BT_tumor)=c("match","not match")
  colnames(chi_BT_tumor)<-seq(0,9,1)
  
  for (j in 1:length(unique(xx$res.0.4_combined)))
  {
    chi_BT_blood["match",unique(xx$res.0.4_combined)[j]]=sum(xx$res.0.4_combined==unique(xx$res.0.4_combined)[j])
    chi_BT_blood["not match",unique(xx$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="blood",])$res.0.4_combined==unique(xx$res.0.4_combined)[j])-sum(xx$res.0.4_combined==unique(xx$res.0.4_combined)[j])
  }
  
  
  
  for (j in 1:length(unique(yy$res.0.4_combined)))
  {
    chi_BT_tumor["match",unique(yy$res.0.4_combined)[j]]=sum(yy$res.0.4_combined==unique(yy$res.0.4_combined)[j])
    chi_BT_tumor["not match",unique(yy$res.0.4_combined)[j]]=sum((CD8_STB1356@meta.data[CD8_STB1356@meta.data$label=="tumor",])$res.0.4_combined==unique(yy$res.0.4_combined)[j])-sum(yy$res.0.4_combined==unique(yy$res.0.4_combined)[j])
  }
  
  
  chi_blood<-chisq.test(chi_BT_blood)
  chi_observed<-chi_blood$observed
  chi_expected<-chi_blood$expected
  ROE_blood<-chi_observed/chi_expected
  
  
  
  chi_tumor<-chisq.test(chi_BT_tumor)
  chi_observed<-chi_tumor$observed
  chi_expected<-chi_tumor$expected
  ROE_tumor<-chi_observed/chi_expected
  
  enrich_blood_all_clonotypes[i,]<-ROE_blood[1,]
  enrich_tumor_all_clonotypes[i,]<-ROE_tumor[1,]
  
  
  number_blood_all_clonotypes[i,1:10]<-chi_BT_blood[1,]
  number_tumor_all_clonotypes[i,1:10]<-chi_BT_tumor[1,]
  
  number_blood_all_clonotypes<-as.data.frame(number_blood_all_clonotypes)
  number_tumor_all_clonotypes<-as.data.frame(number_tumor_all_clonotypes)
  
  PT_all<-combine(PT_all,unique(xx$PT_number))
  
  
  
  
  myoutf <- paste0("../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/UMAP_matched_clonotype",i,".pdf")
  pdf(myoutf,width=6.5,height=5)
  print(DimPlot(object = CD8_STB1356,group.by='single_clonotype',reduction="umap",pt.shape="16",
                cols = col,pt.size = 0.2))
  dev.off()
  
}

number_blood_all_clonotypes<-cbind(number_blood_all_clonotypes,PT_all)
number_tumor_all_clonotypes<-cbind(number_tumor_all_clonotypes,PT_all)



myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/blood_enrichment_score_each_cluster.xls"
write.table( enrich_blood_all_clonotypes,myoutf,sep='\t')

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/tumor_enrichment_score_each_cluster.xls"
write.table(enrich_tumor_all_clonotypes,myoutf,sep='\t')




#Calculate the enrichment of every clonotype in each cluster

cluster_number<-matrix(0,4,10)
row.names(cluster_number)<-paste0("PT_",c(1,3,5,6))
colnames(cluster_number)<-seq(0,9,1)
for(i in 1:ncol(cluster_number))
{
  tag<-CD8_STB1356@meta.data$res.0.4_combined==colnames(cluster_number)[i]
  xx<-CD8_STB1356@meta.data[tag,]
  cluster_number[1,i]<-sum(xx$PT_number=="PT_1")
  cluster_number[2,i]<-sum(xx$PT_number=="PT_3")
  cluster_number[3,i]<-sum(xx$PT_number=="PT_5")
  cluster_number[4,i]<-sum(xx$PT_number=="PT_6")
}


fraction_blood_all_clonotypes<-number_blood_all_clonotypes
fraction_tumor_all_clonotypes<-number_tumor_all_clonotypes

for (i in 1:nrow(number_blood_all_clonotypes))
{
  c<-number_blood_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_blood_all_clonotypes[i,1:10]<-fraction
}


for (i in 1:nrow(number_tumor_all_clonotypes))
{
  c<-number_tumor_all_clonotypes[i,]
  cluster_number_use<-cluster_number[c$PT_all,]
  fraction<-c[1:10]*sum(cluster_number_use)/sum(c[1:10])/cluster_number_use
  fraction_tumor_all_clonotypes[i,1:10]<-fraction
}





myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/blood_enrich_hypergeometric_clustersize_eachPT_each_cluster.xlsx"
write.xlsx(fraction_blood_all_clonotypes,myoutf)

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/tumor_enrich_hypergeometric_clustersize_eachPT_each_cluster.xlsx"
write.xlsx(fraction_tumor_all_clonotypes,myoutf)

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/blood_enrich_sig_geometric_clustersize_eachPT_each_cluster.xls"
write.table(fraction_blood_all_clonotypes_sig,myoutf,sep='\t')


colourCount = 11
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_each_cluster.pdf"
pdf(myoutf,width=15,height=30)
print(ggplot(data=enrich_all, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_each_cluster_blood.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_blood, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()


t.test(enrich_all_blood$index[1:10],enrich_all_blood$index[11:20],alternative = "less")


myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_each_cluster_tumor.pdf"
pdf(myoutf,width=15,height=20)
print(ggplot(data=enrich_all_tumor, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_all$cluster))+
        geom_line(aes(color=clonotype),size=1.5)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +geom_hline(yintercept=1, linetype="dashed", 
                  color = "red", size=2)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))
dev.off()


t.test(enrich_all_tumor$index[1:15],enrich_all_tumor$index[16:30],alternative = "less")


enrich_blood_plot<-enrich_all_blood[as.character(enrich_all_blood$cluster)%in%c("0","1","3","6","7"),]
enrich_tumor_plot<-enrich_all_tumor[as.character(enrich_all_tumor$cluster)%in%as.character(c(0,1,3,6,7)),]



myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_cluster_in_blood_tumor_blood.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_blood_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_blood_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()

myoutf<-"../STB1356/CD8_STB1356/TCRseqed/clean_CD8_exp_fromraw/sctransform/PCx/res.0.4/matched_TCR/single_clonotype_match/BT_blood_tumor_expanded/line_plot_enrichment_cluster_in_blood_tumor_tumor.pdf"
pdf(myoutf,width=30,height=20)
print(ggplot(data=enrich_tumor_plot, aes(x=cluster, y=index, group=clonotype)) + 
        scale_x_discrete(limits=unique(enrich_tumor_plot$cluster))+
        scale_y_continuous(expand = c(0,0.025))+
        geom_line(aes(color=clonotype),size=1.5)+
        geom_hline(yintercept=1, linetype="dashed", 
                   color = "red", size=2)+
        scale_color_manual(values=getPalette(colourCount))+
        geom_point(aes(color=clonotype,shape=tissue),size=8)
      +theme(axis.text.x = element_text(angle = 60, hjust = 1,size=6),axis.title.x = element_blank(),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black")),
      legend.text=element_text(size=3),
      legend.key.size = unit(0.1, "cm"),
      legend.key.width = unit(0.1,"cm"))

dev.off()

