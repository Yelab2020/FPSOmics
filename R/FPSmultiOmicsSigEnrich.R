
#' Function enrichment
#'
#' Study the differences of different ferroptosis status in cancer hallmarks, KEGG pathway by Gene Set Enrichment Analysis (GSEA) and GO, KEGG enrichment.
#'
#' @param input_PSM data.frame of the PSM result from FPSmultiOmicsSig function
#' @param input_gs_kegg the genelist of two column for GSEA enrichemnt, the first column is term and the second column is gene
#' @param temp_folder the temporal folder for output result
#' @return table and figures of GO, KEGG and GSEA enrichemnt result
#' @export
#'
#' @examples
#' data(m4_input_gs_kegg,package='FPSOmics')
#' data(m4_input_PSM,package='FPSOmics')
#' enrichment_result=FPSmultiOmicsSigEnrich(m4_input_PSM,m4_input_gs_kegg,'temp_folder')
#'
FPSmultiOmicsSigEnrich <- function(input_PSM,input_gs_kegg,temp_folder) {
#mRNA <- read.delim("/work/dy/Collaboration/he/Ferroptosis_website_code/result/mRNAseqlog2.genes.SKCM.txt")
  library(magrittr)
  library(dplyr)
  library(ggplot2)
  if (!file.exists(temp_folder)) { dir.create(temp_folder) }

mRNA=input_PSM
mRNA_sig <- mRNA[mRNA$fdr.sig<0.05&abs(mRNA$coef.sig)>1,]
mRNA_sig$Class <- ifelse(mRNA_sig$coef.sig>0,"Up","Down")
SigGene <-  mRNA_sig
formula_res <- clusterProfiler::compareCluster(
  keyType="SYMBOL",
  feature.sig ~ Class,
  data=SigGene,
  fun="enrichGO",
  OrgDb="org.Hs.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
Simple_ego <- clusterProfiler::simplify(
  formula_res,
  cutoff=0.7,
  by="p.adjust",
  select_fun=min
)

Simple_ego_d <- Simple_ego %>% data.frame()

#输出 GO 功能富集 表 1
write.table(Simple_ego_d,file=paste(temp_folder,"/GO_Up_down.tab",sep=''),sep="\t",row.names = F,quote = F)

#输出GO功能富集 图 1   top 10
pdf(paste(temp_folder,"/GO_Up_down.pdf",sep=''),width=8,height = 4)
GO_Up_down=clusterProfiler::dotplot(Simple_ego,showCategory = 10)+
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(11,'Spectral')[2:11])+
  labs(title = "high_score_vs_low_score_GOplot")
print(GO_Up_down)
dev.off()


genelist <- SigGene$feature.sig
genelist <- clusterProfiler::bitr(genelist,fromType = "SYMBOL",
                                  toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)


## KEGG ###############################################################################################################################################################################
SigSub <- merge( SigGene,genelist,by.x="feature.sig",by.y="SYMBOL")
KEGG_res <- clusterProfiler::compareCluster(
  #keyType="SYMBOL",
  ENTREZID~Class,
  data=SigSub,
  organism="hsa",
  fun="enrichKEGG",
  # OrgDb="org.Mm.eg.db",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

KEGG_Result <- KEGG_res %>% data.frame()
#输出 KEGG 功能富集 表 2
write.table(KEGG_Result,file=paste(temp_folder,"/KEGG_Up_down.tab",sep=''),sep="\t",row.names = F,quote = F)

#输出 KEGG 功能富集 图 2
pdf(paste(temp_folder,"/KEGG_Up_down.pdf",sep=''),width=8,height = 4)
KEGG_Up_down=clusterProfiler::dotplot( KEGG_res,showCategory = 10)+
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(11,'Spectral')[2:11])+
  labs(title ="high_score_vs_low_score_KEGGplot")
print(KEGG_Up_down)
dev.off()

##GSEA 富集分析，GO，KEGG
options(stringsAsFactors = F)
library(org.Hs.eg.db)
library(KEGG.db)
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(magrittr)

res1_order <- mRNA[rev(order(mRNA$coef.sig)),]
subRankList <- res1_order$coef.sig
names(subRankList) <- res1_order$feature.sig
#subRankList  <- rev(sort(subRankList))
fgseaRes_sub <- GSEA(geneList = subRankList,
                     nPerm=1000,TERM2GENE=input_gs_kegg ,
                     minGSSize=5,
                     maxGSSize=1000,
                     pvalueCutoff = 0.8)
fgseaRes <- fgseaRes_sub %>% data.frame()

#rio::export(fgseaRes, file = "/work/dy/Collaboration/he/Ferroptosis_website_code/result/GSEA_mRNASig.xlsx")
write.table(fgseaRes,file=paste(temp_folder,"/SEA_Up_down.tab",sep=''),sep="\t",row.names = F,quote = F)

#
fgseaRes_top_down <- fgseaRes[order(fgseaRes$NES),][c(1:10),]
fgseaRes_top_up <- fgseaRes[rev(order(fgseaRes$NES)),][c(1:10),]
fgseaRes_top=rbind(fgseaRes_top_up,fgseaRes_top_down)

#输出 GSEA 功能富集 图 3
pdf(paste(temp_folder,"/GSEA_Up_down_diff_mRNA.pdf",sep=''),width=6,height = 4)
GSEA_Up_down_diff_mRNA=ggplot(data = fgseaRes_top, aes(x = NES, y = Description, fill = -log10(pvalue))) +
  geom_bar(stat="identity",position = "stack",width = 0.7) +
  #  scale_y_continuous(expand=c(0,0)) +
  scale_y_discrete(limit=c(fgseaRes_top_up$Description,rev(fgseaRes_top_down$Description)),label=gsub("HALLMARK_","",c(fgseaRes_top_up$Description,rev(fgseaRes_top_down$Description))))+
  theme(axis.text.x=element_text(color="black",hjust=1,vjust=0.5,angle=90),axis.text.y = element_text(color="black"),
        axis.title.y = element_blank(),
        panel.background = element_rect(color="black",fill=NA),panel.grid=element_blank(),
        legend.title=element_text(size=8), axis.ticks.x = element_blank(),
        legend.text = element_text(size=10)) + #legend.position = "horizontal"
  #scale_color_manual(limits=c("Sig","Nonsig"),values=c("black",NA),guide=F)+
  scale_fill_gradientn(limit=c(-log10(0.05),3),colors= colorRampPalette(RColorBrewer::brewer.pal(5,"Reds"),space="rgb")(20),
                       breaks = c(-log10(0.05),2,3),labels = c(0.05,0.01,0.001),name="Pvalue")+
  labs(title = "high_score vs low_score mRNA GSEA",x="NES")
print(GSEA_Up_down_diff_mRNA)
dev.off()
}

