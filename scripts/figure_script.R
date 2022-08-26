setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))
 
### load R packages####

library(tidyverse)

### ggplot theme ####

theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=10,  color = "black"),
        legend.text = ggplot2::element_text(size=10,  color = "black"),
        legend.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.text =  ggplot2::element_text(size=10,  color = "black"),
        strip.text = ggplot2::element_text(size=10, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))

ancestry.colours <- c("gold2","plum4", "darkorange1",   "lightskyblue2", 
                      "springgreen4", "lightpink2",  "deepskyblue4", 
                      "yellow3",  "yellow4",  
                      'black','red2', 'cornflowerblue', 'magenta', 'darkolivegreen4', 
                      'indianred1', 'tan4', 'darkblue', 'yellowgreen', "tan1",
                      'darkgray', 'wheat4', '#DDAD4B', 'chartreuse','seagreen1',
                      'moccasin', 'mediumvioletred', 'cadetblue1',"darkolivegreen1" ,"#7CE3D8",
                      "gainsboro","#E69F00","#009E73", "#F0E442", "sienna4", "#0072B2", 
                      "mediumpurple4","#D55E00", "burlywood3","gray51","#CC79A7","gray19", "firebrick") 
 
period_size_color <- c("1"="gold2","2"="plum4","3"="darkorange1",
                       "4"= "lightskyblue2","5" = "springgreen4", "6" = "lightpink2" )


#### func ####


fisher_enrich_func <- function(pt_deg){
  sis <- pt_deg$sig_n  
  silp <- pt_deg$n - sis  
  fis <- pt_deg$sig_total - sis   
  filp <- pt_deg$total - pt_deg$sig_total - silp  
  ftp <- fisher.test(matrix(c(sis,silp,fis,filp), 2, 2), alternative='greater')
  return(ftp$p.value)
}

 

############# Figure  1    ###############
#          pSTR distribution             #
##########################################

table_s1 <- data.table::fread("../processed_data/table_s1_refSTRs_pSTRs.txt")  

table_s1_poly <- table_s1 %>% 
  dplyr::filter(pSTRs=="Yes")
 
###### fig_1a ######
str_dist_polym <- table_s1_poly %>% 
  dplyr::filter(!Chr=="MtDNA")  

domain_count_polym <- str_dist_polym  %>% 
  dplyr::group_by(Chr,domain,domain_start ,domain_end ) %>% 
  dplyr::count() %>%
  dplyr::mutate(npm=n*1e6/(domain_end+1-domain_start),
                Pos=(domain_end-domain_start)/2+domain_start) 

fig_1a <-  ggplot() + 
  geom_histogram(data=str_dist_polym, aes(x=start/1e6), bins = 50,fill="gray69") +
  geom_line(data=domain_count_polym, aes(x=Pos/1e6,y=npm/2),color="red",size=0.5 ) +
  geom_point(data=domain_count_polym, aes(x=Pos/1e6,y=npm/2),fill="red",shape=25,size=2 ) +
  scale_y_continuous( name = "Number of\npSTRs", sec.axis = sec_axis( trans=~.*2, name="Number of\npSTRs / Mb") )+
  facet_grid(.~Chr,scales = "free", space="free") +
  theme_cust +
  xlab("Genomic Position (Mb)") 

###### fig_1b ######

polystr_motif_count  <- table_s1_poly %>% 
  dplyr::mutate( motif_geno=motif_geno_fwd ) %>% 
  dplyr::group_by(motif_geno,motif_length ) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::top_n(10, n)  %>% 
  dplyr::arrange(  n )

polystr_motif_count$motif <- factor(polystr_motif_count$motif_geno, levels = polystr_motif_count$motif_geno )


fig_1b <- ggplot(polystr_motif_count,aes(x=motif,y=n ,fill=factor(motif_length))) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=period_size_color) +
  theme_cust +
  coord_flip() +
  theme(legend.position = "none") +
  labs(y="Number of sites",
       x="STR motif") +
  scale_y_continuous(breaks=c(0, 1000,1800),limits = c(0,1900) )

###### fig_1c ######

PolySTR_region <- table_s1_poly %>% 
  dplyr::group_by(gfeature,motif_length) %>% 
  dplyr::add_count(name = "STRbyPS_REGION")%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::add_count(name = "STRby_REGION") %>% 
  dplyr::mutate(ps2region=100*STRbyPS_REGION/STRby_REGION) %>% 
  dplyr::distinct(gfeature,motif_length,STRbyPS_REGION,STRby_REGION,ps2region) %>% 
  dplyr::mutate(STRby_REGIONs=ifelse(motif_length==6,STRby_REGION,NA))


PolySTR_region$gfeatures<- factor(PolySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

#plot

fig_1c <- ggplot(PolySTR_region,aes(x=gfeatures,y=ps2region,fill=factor(motif_length))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_cust +
  xlab("Genomic features")+
  ylab("Percent of pSTRs (%)")  +
  labs(fill="Motif\nlength")+
  scale_fill_manual(values=period_size_color) +
  geom_text(aes(label=STRby_REGIONs),y=115,size = 10*5/14)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_1d ######

enrich_PolySTR_region_stats <- PolySTR_region %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(sig_n=STRbyPS_REGION,
                sig_total=STRby_REGION) %>% 
  dplyr::group_by(motif_length) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>% 
  dplyr::group_by(motif_length,gfeature)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.1D")%>% 
  dplyr::mutate(group_factor="gfeature, motif_length",
                group_factor_catogory=paste0(gfeature,  ", ", motif_length),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrich_PolySTR_region <- enrich_PolySTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05) %>% 
  dplyr::mutate(logp=-log10(fisherp_adj))  

enrich_PolySTR_region$gfeatures<- factor(enrich_PolySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))


fig_1d <-ggplot(enrich_PolySTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F)+
  theme_cust+
  theme( legend.position = "none")+
  scale_color_manual(values=period_size_color) +
  scale_x_continuous(breaks = c(0,150,300 ),labels = c("0","150", "300" )  )   +
  ylab("Genomic\nfeatures")+
  xlab(expression(-log[10](italic(p)))) 

###### fig_1e ######
 
enrichMotif_polySTR_region_stats <-  table_s1_poly %>% 
  dplyr::mutate( motif_geno=motif_geno_fwd ) %>% 
  dplyr::group_by(gfeature, motif_geno ) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno )  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.1E")%>% 
  dplyr::mutate(group_factor="gfeature, motif_geno",
                group_factor_catogory=paste0(gfeature,  ", ", motif_geno),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrichMotif_polySTR_region <- enrichMotif_polySTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05)   %>%  
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::top_n(3, logp) 


enrichMotif_polySTR_region$gfeatures<- factor(enrichMotif_polySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))



fig_1e <- ggplot(enrichMotif_polySTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 3,segment.linetype=6,
                           max.overlaps=Inf,size=10*5/14,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic\nfeatures")+
  xlab(expression(-log[10](italic(p)))) 




###### fig 1 #####

fig_1bc <-  cowplot::plot_grid(fig_1b, fig_1c,
                              labels = c('', 'C'), 
                              rel_widths =  c(1.3,2),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "t",
                              nrow = 1)

fig_1de <-  cowplot::plot_grid(fig_1d,fig_1e,
                               labels = c('', 'E'), 
                               rel_widths =  c(1.2 ,2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)


fig_1 <-  cowplot::plot_grid(fig_1a, fig_1bc,fig_1de,  
                            labels = c('A', 'B',"D" ), 
                            rel_heights =  c(1,1.5,1 ),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            nrow = 3)

ggsave(fig_1, filename = paste( "../figures/Fig_1.png",sep = ""), units = "mm",height = 160, width = 170)

 
############# Figure  2    ###############
#         expansion contraciton          #
##########################################

###### fig_2a ######  
bp_diff_BPD <- table_s1_poly %>% 
  dplyr::select(ref_STR,Chr,start,BPDIFFS) %>% 
  dplyr::arrange(Chr,start) %>% 
  splitstackshape::cSplit("BPDIFFS",",", direction = "long",sep = ",")  

fig_2a <- ggplot(bp_diff_BPD,aes(BPDIFFS))   + 
  geom_histogram(color="black", fill="white",bins =80,size=0.2) +
  theme_cust +
  ylab("Number of\nalleles") +
  xlab("Base-pair difference")

###### fig_2b ###### 


expansion_contractionS <- table_s1_poly %>% 
   dplyr::select(ref_STR ,expansion_score,contraction_score,motif_length) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length) %>% 
  dplyr::filter(!score==0) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(diff,motif_length) %>% 
  dplyr::mutate(lg_sc=log10(abs(score))) %>% 
  dplyr::mutate(mean_lg_sc=mean(lg_sc),
                mean_sc=mean(score),
                median_sc=median(score)) %>% 
  dplyr::mutate(id=row_number(),
                mean_lg_sc=ifelse(id==1,round(mean_lg_sc,digits = 2),NA))


fig_2b <- ggplot(expansion_contractionS ,aes(score,fill=diff)) +
  geom_histogram(color="black", 
                 bins =100,size=0.2) +
  theme_cust +
  xlab("Contraction / Expansion score")+
  ylab("Number of\nSTRs")   +
  scale_fill_manual(values=c("#E7B800", "#00AFBB")) +
  theme( legend.position = "none") 

###### fig_2c ###### 

 
fig_2c <- ggplot(expansion_contractionS,aes(x=diff,y=lg_sc ,color=diff))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=mean_lg_sc),size=0.5,color="red") +
  theme_cust +
  ylab("Contraction /\nExpansion\ntransformed scores") +
  ggpubr::stat_compare_means(  label = "p.signif", method = "wilcox.test", 
    label.x=1.35,
    symnum.args = list(cutpoints = c(0, 0.00001, 0.0001, 0.001,  1), 
                       symbols = c("****","***","**",  "ns")) ,
    label.y = 0.2) +
  facet_grid(.~motif_length,scales="free") +
  scale_color_manual(values=c("#E7B800", "#00AFBB")) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        legend.position = "none")  + 
   scale_y_continuous(expand = c(0, 0), limits = c(-2, 0.5))
 
###### fig_2d ###### 

expand_frac <- table_s1_poly %>% 
  dplyr::select(ref_STR ,Contraction_frac,Expansion_frac,motif_length) %>% 
  tidyr::gather(diff,frac,-ref_STR,-motif_length) %>% 
  dplyr::mutate(frac=frac/100) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(diff,motif_length)%>% 
  dplyr::mutate(mean_frac=mean(frac),
                median_frac=median(frac)) %>% 
  dplyr::mutate(id=row_number(),
                mean_frac=ifelse(id==1,round(mean_frac,digits = 2),NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate( diff=sub("(.*)(_frac)","\\1",diff))



fig_2d <- ggplot(expand_frac,aes(x=diff,y=frac ,color=diff))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=abs(mean_frac)),size=0.5,color="red") +
  theme_cust +
  ylab("Allele frequency") +
  xlab("Variation to median allele")+
  ggpubr::stat_compare_means(label = "p.signif",
                              label.y = 0.53 ,  label.x=1.35,
                              symnum.args = list(cutpoints = c(0,  0.001,  1), 
                                                 symbols = c("****",   "ns")),
                              method = "wilcox.test" ) +
  facet_grid(.~motif_length,scales="free") +
  scale_color_manual(values=c("#E7B800", "#00AFBB")) +
  theme(axis.text.x = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6)) 



 


###### fig_2e ###### 

expansion_contractionS_di <- table_s1_poly %>% 
  dplyr::filter(!grepl("/",motif_geno)) %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>%
  dplyr::select(ref_STR ,expansion_score,contraction_score,motif_length,motif_geno) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length,-motif_geno) %>% 
  dplyr::filter(motif_length==2) %>% 
  dplyr::mutate(diff=ifelse(score==0,"substitution",diff)) %>% 
  dplyr::group_by(motif_geno,diff,motif_length) %>% 
  dplyr::count()  %>% 
  dplyr::group_by(motif_geno) %>% 
  dplyr::mutate(count_motif=sum(n)) %>% 
  dplyr::mutate(frac=100*n/count_motif) %>% 
  dplyr::mutate(count_motif2=ifelse(diff=="substitution",count_motif,NA)) %>% 
  dplyr::mutate(diff=sub("(.*)(_score)","\\1",diff))


fig_2e <- ggplot(expansion_contractionS_di,aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=105,size = 10*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,110) ) 



###### fig 2 #####

fig2ab <-  cowplot::plot_grid(fig_2a, fig_2b,
                              labels = c('', 'B'), 
                              # rel_widths =  c(1,2),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "tb",
                              # align = "h",
                              nrow = 1)

fig2cd <-  cowplot::plot_grid(fig_2c,fig_2d,
                              labels = c('', 'D'), 
                              rel_heights =  c(1.1,1 ),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              align = "v",
                              nrow = 2)

fig2 <-  cowplot::plot_grid(fig2ab, fig2cd,fig_2e,
                            labels = c('A', 'C' , "E"), 
                            rel_heights =  c(1,2,1.2),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            #   align = "v",
                            nrow = 3)

ggsave(fig2, filename = paste( "../figures/Fig_2.png",sep = ""), units = "mm",height = 180, width = 170)


 
############# Figure  3    ###############
#             popgenet                    #
##########################################

###### fig_3a ######

num_allele <- bp_diff_BPD %>% # fig_2a
  dplyr::group_by(ref_STR) %>% 
  dplyr::count() %>% 
  dplyr::mutate(n=n+1)

fig_3a <- ggplot(num_allele,aes(n))   + 
  geom_histogram(color="black", fill="white",bins =21,size=0.2) +
  theme_cust +
  xlab("Number of alleles\nper pSTR") +
  ylab("Number of pSTRs") + 
  scale_y_continuous(breaks=c(0, 3000, 1500), limits=c(0, 3000))+ 
  scale_x_continuous(breaks=c(2,10,21) )  +
  theme( axis.title.y = ggplot2::element_text(size=10,  color = "black", hjust = 0.2 )) 
 
###### fig_3b ######

majorAF_ExpectedHe <- data.table::fread("../processed_data/Expected_Heterozygosity_swept_div_all.csv")

majorAF_data <- majorAF_ExpectedHe %>% 
  dplyr::group_by(ref_STR,n_st) %>% 
  dplyr::mutate( major_af=max(af)) %>% 
  dplyr::filter(af==major_af)  


fig_3b <- ggplot()   + 
  ggplot2::stat_density(data=majorAF_data,aes(x=major_af,color=n_st ), geom="line",position = "identity",
                         size=0.8 ) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +
  theme_cust +
  theme(legend.position = c(.4,.6) ,
        legend.spacing  = unit(0.1, 'cm'),
        legend.margin=margin(1,1,1,1)) + 
  ylab("Density") +
  xlab("Major\nallele frequency") +
  scale_x_continuous(breaks=c(0, 0.5, 1),limits = c(0,1) )+
  labs(color="Strains")

###### fig_3c ######

strain_ALT_frac <- data.table::fread("../processed_data/strain_ALT_frac.tsv")

fig_3c <- ggplot()+
  geom_point(data=strain_ALT_frac,aes(x=alt_frac,y=hets_frac,color=sweep),size=0.5,alpha=0.8)+
  theme_cust+
  labs(x="Homozygous ALT\npSTRs in each strain (%)",
       y= "Heterozygous pSTRs\nin each strain (%)") +
  scale_color_manual(values = c("#E7B800", "#FC4E07")) +
  theme( legend.background = element_rect(colour = 'gray79', fill = 'white', linetype='solid'),
         legend.spacing  = unit(0.1, 'cm'),
         legend.margin=margin(1,1,1,1),
         plot.margin = unit(c(0, 2, 0, 5), "mm"),
         legend.position="none")

###### fig_3d ######

pca_pSTR_SNV <- data.table::fread("../processed_data/pca_pSTR_SNV.tsv")

pca_pSTR <- pca_pSTR_SNV %>% 
  dplyr::filter(data=="pSTRs")

fig_3d <- ggplot(pca_pSTR, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=cluster ),size=1,alpha=0.8) +
  theme_cust + 
  theme(legend.position="none")+
  labs(x=paste0("PC1: ",unique(pca_pSTR$PC1_var_exp)[1],"%"),
       y=paste0("PC2: ",unique(pca_pSTR$PC2_var_exp)[1],"%")) +
  ggtitle("9,691 pSTRs") +
  scale_color_manual(values=ancestry.colours)




###### fig_3e ######

pca_SNV <- pca_pSTR_SNV %>% 
  dplyr::filter(data=="SNVs")


fig_3e <- ggplot(pca_SNV, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=cluster ),size=1,alpha=0.8) +
  theme_cust + 
  labs(x=paste0("PC1: ",unique(pca_SNV$PC1_var_exp)[1],"%"),
       y=paste0("PC2: ",unique(pca_SNV$PC2_var_exp)[1],"%"),
       color = "Sample locations") +
  ggtitle("13,580 SNVs") +
  scale_color_manual(values=ancestry.colours)



###### fig_3f ######

data_fig_3f <- majorAF_ExpectedHe %>% 
  dplyr::select(Chr,start,ref_STR,hets,n_st) %>% 
  dplyr::distinct()  %>% 
  dplyr::filter(Chr != "MtDNA")


fig_3f <- ggpubr::ggscatter(data_fig_3f, 
                                      x = "start", y = "hets",  point=FALSE, palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                                      add = "loess", conf.int = FALSE,color="n_st" ) + 
  labs(x= "Genomic position (Mb)",
       y= expression(italic(H)[E])) +
  facet_grid(.~Chr,scales="free")+
  theme_cust +
  theme( axis.text.x = element_blank(), 
         legend.position = "none",
         panel.spacing = unit(0.2,"line"))


###### fig_3 ######

fig3abc <-  cowplot::plot_grid(fig_3a, fig_3b,fig_3c,
                               labels = c('', 'B','C'), 
                               rel_widths =  c(0.9,0.9 , 1.2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "tb",
                               align = "h",
                               nrow = 1)

fig3de <- cowplot::plot_grid(fig_3d,  fig_3e,   
                             labels = c('', 'E'  ), 
                             rel_widths = c( 1,1.75  ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             nrow = 1)

fig3 <-  cowplot::plot_grid(fig3abc, fig3de, fig_3f,  
                            labels = c('A' ,'D','F' ), 
                            rel_heights =  c(1.2,1.5,1.1 ),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            align = "v",
                            nrow = 3)

ggsave(fig3, filename = paste( "../figures/Fig_3.png",sep = ""), units = "mm",height = 180, width = 170)

 


############# Figure  4    ###############
#          MA lines                      #
##########################################


data_fig_4 <- data.table::fread("../processed_data/MA174_pSTRs_mutationRate.tsv")

###### fig_4a ######

data_fig_4a <- data_fig_4 %>% 
  dplyr::filter(motif_length=="1-6" &
                  motif_geno=="All" &
                  gfeature=="All" &
                  OMA=="O1MA") %>% 
  dplyr::group_by(strain,Line) %>% 
  dplyr::mutate(mutation_rate=sum(mutation_rate)) %>% 
  dplyr::select(-mutation,-N_mutation ) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(mutation="All mutations")

 

fig_4a <- ggplot(data_fig_4a,aes(x=strain,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8 ,aes(color=strain)) +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=strain)) +
  theme_cust +
  labs(x="Strain",y="ANC-O1MA\nMutation rate" ) +
  theme( legend.position = "none"  ) +
  scale_color_manual(values = c("orange","#007e2f","#ffcd12","#721b3e") ) +
  ggpubr::stat_compare_means(label.y = 0.000135,
    label.x = 1.5,
    size = 4,
    label = "p.signif", 
    method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1.5e-04), breaks = c(0, 7.5e-05, 1.5e-04))   



###### fig_4b ######


data_fig_4b <- data_fig_4 %>% 
  dplyr::filter(motif_length=="1-6" &
                  motif_geno=="All" &
                  gfeature=="All" &
                  OMA=="O1MA") 
 

fig_4b <- ggpubr::ggboxplot(data_fig_4b, x="mutation",y="mutation_rate",outlier.shape = NA,
                            color="strain"  ) +
  geom_point( position = position_jitterdodge(jitter.width = 0.2) ,aes(color=strain), size=0.5, alpha=0.8)+
  theme_cust +
  theme(legend.position = "none")+
  scale_color_manual(values = c("orange","#007e2f","#ffcd12","#721b3e") ) +
  labs(x="Mutations",y="ANC-O1MA\nMutation rate" ) + 
  ggpubr::stat_compare_means( 
    aes(group = strain),
    label = "p.signif",   label.y = 4.8e-05,
    size = 4,
    method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0, 5.2e-05), breaks = c(0, 2.5e-05, 5e-05))   


 
###### fig_4c ######

data_fig_4c <- data_fig_4 %>% 
  dplyr::filter(motif_length=="1-6" &
                  motif_geno=="All" &
                  gfeature!="All" &
                  OMA=="O1MA") %>% 
  dplyr::group_by(strain,Line,gfeature) %>% 
  dplyr::mutate(mutation_rate=sum(mutation_rate)) %>% 
  dplyr::select(-mutation,-N_mutation ) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(mutation="All mutations") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(gfeature=ifelse(gfeature=="pseudogene","pseudogene   ",gfeature))


data_fig_4c$gfeatures<- factor(data_fig_4c$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene   ","RNAs & TEs","intergenic"))
 
fig_4c <-ggpubr::ggboxplot(data_fig_4c, x="gfeatures",y="mutation_rate",outlier.shape = NA  ) +
  geom_point( position = position_jitterdodge(jitter.width = 0.2) ,
              aes(color=strain), size=0.5, alpha=0.8)+
  theme_cust +
  theme(legend.position = "none",
        axis.text.x =  ggplot2::element_text(size=8, 
                                             color = "black") )+
  scale_color_manual(values = c("orange","#007e2f","#ffcd12","#721b3e") ) +
  labs(x="Genomic features",y="ANC-O1MA\nMutation rate" ) + 
  ggpubr::stat_compare_means( 
    label = "p.signif", 
    label.y = 2.3e-04,
    ref.group = "CDS",
    size = 4,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001,  1), 
                       symbols = c("****","***",  "ns")),
    method = "wilcox.test")+ 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0, 2.5e-04), 
                     breaks = c(0, 1e-04,  2e-04 ))   


###### fig_4d ######

data_fig_4d <- data_fig_4 %>% 
  dplyr::filter(motif_length=="1-6" &
                  motif_geno=="All" &
                  gfeature=="All" &
                  OMA=="O2MA",
                comparison=="O1MA-O2MA") %>% 
  dplyr::group_by(strain,Line) %>% 
  dplyr::mutate(mutation_rate=sum(mutation_rate)) %>% 
  dplyr::select(-mutation,-N_mutation ) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(mutation="All mutations")%>% 
  dplyr::ungroup()

 


fig_4d <- ggplot(data_fig_4d,aes(x=Fitness,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8 ,aes(color=Fitness)) +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=Fitness)) +
  theme_cust +
  labs(x="Fitness of O1MA progenitors",y="O1MA-O2MA\nMutation rate" ) +
  theme( legend.position = "none" ,
         axis.title.x =  ggplot2::element_text(size=10,  color = "black",hjust =  0.9) ) +
  scale_color_manual(values = c( "black", 'darkgray') ) +
  ggpubr::stat_compare_means( label.y = 8.8e-05,
    label.x = 1.5,
    size = 4,
    label = "p.signif", 
    method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1e-04), breaks = c(0, 5e-05, 1e-04))   


###### fig_4e ######


data_fig_4e <- data_fig_4 %>% 
  dplyr::filter(motif_length=="1-6" &
                  motif_geno=="All" &
                  gfeature=="All" &
                  OMA=="O2MA",
                comparison=="O1MA-O2MA") 

  
fig_4e <- ggpubr::ggboxplot(data_fig_4e, x="mutation",y="mutation_rate",outlier.shape = NA,
                            color="Fitness"  ) +
  geom_point( position = position_jitterdodge(jitter.width = 0.2) ,aes(color=Fitness), size=0.5, alpha=0.8)+
  theme_cust +
  theme(legend.position = "none")+
  scale_color_manual(values = c( "black", 'darkgray') ) +
  labs(x="Mutations",y="O1MA-O2MA\nMutation rate" ) + 
  ggpubr::stat_compare_means( 
    aes(group = Fitness),
    label = "p.signif",   label.y = 7e-05,
    size = 4,
    method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,8e-05), breaks = c(0,   4e-05,   8e-05))   


###### fig_4f ######

data_fig_4f <- data_fig_4 %>% 
  dplyr::filter(motif_length!="1-6" &
                  motif_geno=="All" &
                  gfeature=="All" &
                  OMA=="O2MA",
                comparison=="O1MA-O2MA")  
 

fig_4f <- ggpubr::ggboxplot(subset(data_fig_4f, motif_length %in% c(1,2)), x="mutation",y="mutation_rate",outlier.shape = NA,
                            color="Fitness"  ) +
  geom_point( position = position_jitterdodge(jitter.width = 0.2) ,aes(color=Fitness), size=0.5, alpha=0.8)+
  theme_cust +
  theme(legend.position = "none")+
  scale_color_manual(values = c( "black", 'darkgray') ) +
  labs(x="Mutations",y="O1MA-O2MA\nMutation rate" ) + 
  facet_grid(.~motif_length,scales="free") +
  ggpubr::stat_compare_means( 
    aes(group = Fitness),
    label = "p.signif",   label.y = 1e-04,
    symnum.args = list(cutpoints = c(0, 0.001, 0.003,  1), 
                       symbols = c("***","*",  "ns")),
    size = 4,
    method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1.2e-04), breaks = c(0,   5e-05,   1e-04))   



### cow fig 4 ####

fig_4ab <-  cowplot::plot_grid(  fig_4a, fig_4b,   
                                 labels = c('', 'B' ), 
                                 rel_widths =  c( 2 ,3 ),
                                 label_size = 12, 
                                 label_fontfamily="Helvetica",
                                 align = "h",
                                 axis = "tb",
                                 nrow =1)




fig_4de <-  cowplot::plot_grid(  fig_4d, fig_4e,   
                                 labels = c('', 'E' ), 
                                 rel_widths =  c( 2 ,3 ),
                                 label_size = 12, 
                                 label_fontfamily="Helvetica",
                                 align = "h",
                                 axis = "tb",
                                 nrow =1)



fig_4 <- cowplot::plot_grid(  fig_4ab, fig_4c,  fig_4de,fig_4f,
                              labels = c('A', 'C' ,'D' ,'F' ), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow =4)


ggsave(fig_4, filename = paste( "../figures/Fig_4.png",sep = ""), units = "mm",height =  180, width = 170)


############# Figure  5    ###############
#          STR and phenotypes           #
##########################################

data_fig_5 <- data.table::fread("../processed_data/table_s4_pSTRonPhenotypes.txt")


traits_name_df <- data.frame(organismal_trait=c("broods","abamectin_norm.n","arsenic_pc1","Albendazole_q90.TOF",
                                                "bleo_medianEXT","dauer_frac","amsacrine_f.L1","etoposide_median.TOF",
                                                "propionicL1survival","telomere_resids","zinc_norm_EXT"),
                             trait_name=c("Lifetime fecundity","Response to abamectin","Response to arsenic","Response to albendazole",
                                          "Response to bleomycin","Dauer formation","Response to amsacrine","Response to etoposide",
                                          "Response to propionate","Telomere length","Response to zinc") )




data_fig_5_pos <- table_s1_poly %>% dplyr::select(pSTR=ref_STR,Chr,start) %>% 
  dplyr::left_join(data_fig_5) %>% 
  na.omit() %>% 
  dplyr::left_join(traits_name_df)


STRmanha_plt_list=list()

for( tr in unique(data_fig_5_pos$organismal_trait)) {
  
  data_fig_5_tr <- data_fig_5_pos %>% 
    dplyr::filter(organismal_trait==tr) %>% 
    dplyr::filter(!Chr=="MtDNA")
  
  STRmanha_plt <- ggplot(data_fig_5_tr,aes(x=start/1e6,y=-log10(lrt_pvalue),color=BF_significant))+
    geom_point(size=0.3)+
    theme_cust+
    facet_grid(.~Chr,scales = "free")+
    scale_color_manual(values = c("gray69","red"))+
    ylab(expression(-log[10](italic(p)))) +
    xlab("Genomic position (Mb)")+
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1,"line"))+
    ggtitle(unique(data_fig_5_tr$trait_name))
  
  STRmanha_plt_list[[tr]] <- STRmanha_plt
}

fig_5 <-  cowplot::plot_grid(STRmanha_plt_list$broods+ theme(axis.title.x = element_blank()),  
                             STRmanha_plt_list$dauer_frac+ theme(axis.title.x = element_blank()), 
                             STRmanha_plt_list$telomere_resids+ theme(strip.text = element_blank(), axis.title.x = element_blank()),
                             STRmanha_plt_list$abamectin_norm.n+ theme(strip.text = element_blank(), axis.title.x = element_blank()), 
                             STRmanha_plt_list$Albendazole_q90.TOF+ theme(strip.text = element_blank(), axis.title.x = element_blank()),
                             STRmanha_plt_list$arsenic_pc1+ theme(strip.text = element_blank(), axis.title.x = element_blank()), 
                             STRmanha_plt_list$bleo_medianEXT+ theme(strip.text = element_blank(), axis.title.x = element_blank()), 
                             STRmanha_plt_list$amsacrine_f.L1+ theme(strip.text = element_blank(), axis.title.x = element_blank()), 
                             STRmanha_plt_list$propionicL1survival+ theme(strip.text = element_blank(), axis.title.x = element_blank()), 
                             STRmanha_plt_list$zinc_norm_EXT+ theme(strip.text = element_blank(), axis.title.x = element_blank()), 
                             STRmanha_plt_list$etoposide_median.TOF + theme(strip.text = element_blank()),
                             labels = c('A', 'B', 'C', 'D', 'E', 'F',
                                        'G', 'H', 'I', 'J', 'K'), 
                             rel_heights =  c(1.2,1,1,1,1,1.2 ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             axis = "lr",
                             align = "v",
                             ncol = 2,
                             nrow = 6)

ggsave(fig_5, filename = paste( "../figures/Fig_5.png",sep = ""), units = "mm",height = 200, width = 170)


##############################################


############# Figure  S1    ###############
#          refSTR distribution             #
##########################################

 str_dist <- table_s1 %>% 
  dplyr::filter(!Chr=="MtDNA")  

domain_count  <- str_dist  %>% 
  dplyr::group_by(Chr,domain,domain_start ,domain_end ) %>% 
  dplyr::count() %>%
  dplyr::mutate(npm=n*1e6/(domain_end+1-domain_start),
                Pos=(domain_end-domain_start)/2+domain_start) 

fig_S1a <-  ggplot() + 
  geom_histogram(data=str_dist, aes(x=start/1e6), bins = 50,fill="gray69") +
  geom_line(data=domain_count, aes(x=Pos/1e6,y=npm/2),color="blue",size=0.5 ) +
  geom_point(data=domain_count, aes(x=Pos/1e6,y=npm/2),fill="blue",shape=25,size=2 ) +
  scale_y_continuous( name = "Number of\nSTRs", sec.axis = sec_axis( trans=~.*2, name="Number of\nSTRs / Mb") )+
  facet_grid(.~Chr,scales = "free", space="free") +
  theme_cust +
  xlab("Genomic Position (Mb)") 


fig_S1b <- ggplot() + 
  geom_histogram(data=str_dist,aes(x=start/1e6,fill=factor(motif_length) ),bins = 50 )  + 
  facet_grid(motif_length~Chr,scales = "free")+
  theme_cust +
  xlab("Genomic Position (Mb)")+
  ylab("Number of STRs")  +
  scale_fill_manual(values=period_size_color) +
  theme(legend.position = "none" ,
        plot.title =   ggplot2::element_text(size=10,  color = "black"))


fig_S1 <-  cowplot::plot_grid(fig_S1a, fig_S1b,  
                              labels = c('A', 'B'  ), 
                              rel_heights =  c(1.3,4  ),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow = 2)

ggsave(fig_S1, filename = paste( "../figures/Supp_fig1_ref_dist.png",sep = ""), units = "mm",height = 180, width = 170)




############# Figure  S2    ###############
#          composition                  #
##########################################
 

perfect <- table_s1 %>% 
  dplyr::filter(!is.na(perfect_str)) %>% 
  dplyr::select(ref_STR,Chr,start,end,polymorphic,perfect_str,pSTRs) %>% 
  dplyr::group_by(perfect_str) %>% 
  dplyr::add_count(name = "All_n") %>% 
  dplyr::group_by(perfect_str,pSTRs) %>% 
  dplyr::add_count(name = "pSTRs_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(perfect_str,All_n,pSTRs_n,pSTRs) %>% 
  tidyr::gather(category,Counts,-perfect_str,-pSTRs) %>% 
  dplyr::filter(!(pSTRs=="No" & category=="pSTRs_n")) %>% 
  dplyr::select(-pSTRs) %>% 
  dplyr::distinct() 



perfect_plt_data <- perfect %>% 
  dplyr::mutate(perfect_str=gsub('_', '-', perfect_str),
                category=ifelse(category=="pSTRs_n", "9,691 pSTRs","27,667 STRs with polymorphisms"))

perfect_plt_data$perfect_str2<- factor(perfect_plt_data$perfect_str,levels = c("compound-interrupted",
                                                                               "compound-center-perfect",
                                                                               "compound-perfect",
                                                                               "simple-interrupted",
                                                                               "simple-center-perfect",
                                                                               "simple-perfect"))

fig_S2 <-ggplot(perfect_plt_data,aes(x=perfect_str2,y=Counts))+
  geom_bar(stat='identity',fill=NA,color="black") +
  theme_cust+
  coord_flip() +
  facet_grid(.~category,scales = "free",space = "free") +
  theme(axis.title.y=element_blank(),
        plot.margin = unit(c(1, 4, 1, 1), "mm"))+
  scale_y_continuous(breaks=seq(0, 20000,5000)  )


ggsave(fig_S2, filename = paste( "../figures/Supp_fig2_perfectSTR.png",sep = ""), units = "mm",height = 100, width = 170)


############# Figure  S3    ###############
#          pSTR distribution             #
##########################################

fig_S3 <- ggplot() + 
  geom_histogram(data=str_dist_polym,aes(x=start/1e6,fill=factor(motif_length) ),bins = 50 )  + 
  facet_grid(motif_length~Chr,scales = "free")+
  theme_cust +
  xlab("Genomic Position (Mb)")+
  ylab("Number of polymorphic STRs")  +
  scale_fill_manual(values=period_size_color) +
  theme(legend.position = "none" ,
        plot.title =   ggplot2::element_text(size=10,  color = "black"))

ggsave(fig_S3, filename = paste( "../figures/Supp_fig3_polystr_dist_size.png",sep = ""), units = "mm",height = 160, width = 170)




############# Figure  S4    ###############
#          refSTR motif enrichment       #
##########################################
 
###### fig_S4a  ######

refstr_motif_count  <- table_s1 %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>%
  dplyr::group_by(motif_geno,motif_length) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::top_n(10, n)  %>% 
  dplyr::arrange(  n )

refstr_motif_count$motif <- factor(refstr_motif_count$motif_geno, levels = refstr_motif_count$motif_geno )


fig_S4a <- ggplot(refstr_motif_count,aes(x=motif,y=n ,fill=factor(motif_length))) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=period_size_color) +
  theme_cust +
  coord_flip() +
  theme(legend.position = "none") +
  labs(y="Number of sites",
       x="STR motif") +
  scale_y_continuous(breaks=c(0, 4000,8000), limits = c(0,9000))

###### fig_S4b ######

refSTR_region <- table_s1 %>% 
  dplyr::group_by(gfeature,motif_length) %>% 
  dplyr::add_count(name = "STRbyPS_REGION")%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::add_count(name = "STRby_REGION") %>% 
  dplyr::mutate(ps2region=100*STRbyPS_REGION/STRby_REGION) %>% 
  dplyr::distinct(gfeature,motif_length,STRbyPS_REGION,STRby_REGION,ps2region) %>% 
  dplyr::mutate(STRby_REGIONs=ifelse(motif_length==6,STRby_REGION,NA))


refSTR_region$gfeatures<- factor(refSTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

#plot

fig_S4b <- ggplot(refSTR_region,aes(x=gfeatures,y=ps2region,fill=factor(motif_length))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_cust +
  xlab("Genomic features")+
  ylab("Percent of pSTRs (%)")  +
  labs(fill="Motif\nlength")+
  scale_fill_manual(values=period_size_color) +
  geom_text(aes(label=STRby_REGIONs),y=115,size = 10*5/14)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_S4c ######

enrich_refSTR_region_stats <- refSTR_region %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(sig_n=STRbyPS_REGION,
                sig_total=STRby_REGION) %>% 
  dplyr::group_by(motif_length) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>% 
  dplyr::group_by(motif_length,gfeature)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.S3C")%>% 
  dplyr::mutate(group_factor="gfeature, motif_length",
                group_factor_catogory=paste0(gfeature,  ", ", motif_length),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrich_refSTR_region <- enrich_refSTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05) %>% 
  dplyr::mutate(logp=-log10(fisherp_adj))  %>% 
  dplyr::mutate(logp=ifelse(logp=="Inf",400,logp))  

enrich_refSTR_region$gfeatures<- factor(enrich_refSTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))


fig_S4c <- ggplot(enrich_refSTR_region,
                aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F)+
  theme_cust+
  theme( legend.position = "none")+
  scale_color_manual(values=period_size_color) +
  scale_x_continuous(breaks = c(0,100,200,300,400),labels = c("0","100","200","300","Inf")  )   +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p)))) 

###### fig_S4d ######


enrichMotif_refSTR_region_stats  <-  table_s1 %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>% 
  dplyr::group_by(gfeature, motif_geno ) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno )  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.S3D") %>% 
  dplyr::mutate(group_factor="gfeature, motif_geno",
                group_factor_catogory=paste0(gfeature,  ", ", motif_geno),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrichMotif_refSTR_region <- enrichMotif_refSTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05)   %>%  
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::top_n(3, logp) 


enrichMotif_refSTR_region$gfeatures<- factor(enrichMotif_refSTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))



fig_S4d <- ggplot(enrichMotif_refSTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 3,segment.linetype=6,size=10*5/14,
                           max.overlaps=Inf,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p)))) 




###### fig S3 #####

fig_S4ab <-  cowplot::plot_grid(fig_S4a, fig_S4b,
                               labels = c('', 'B'), 
                               rel_widths =  c(1.3,2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)

fig_S4cd <-  cowplot::plot_grid(fig_S4c,fig_S4d,
                               labels = c('', 'D'), 
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)


fig_S4 <-  cowplot::plot_grid(fig_S4ab, fig_S4cd,  
                             labels = c('A', 'C'  ), 
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             axis = "lr",
                             nrow = 2)

ggsave(fig_S4, filename = paste( "../figures/Supp_fig4_refmotif.png",sep = ""), units = "mm",height = 160, width = 170)






############# Figure  S5    ###############
#           motif 134                     #
##########################################

expansion_contractionS_134 <- table_s1_poly %>% 
  dplyr::filter(!grepl("/",motif_geno)) %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>%
  dplyr::select(ref_STR ,expansion_score,contraction_score,motif_length,motif_geno) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length,-motif_geno) %>% 
  dplyr::filter(motif_length %in% c(1,3,4)) %>% 
  dplyr::mutate(diff=ifelse(score==0,"substitution",diff)) %>% 
  dplyr::group_by(motif_geno,diff,motif_length) %>% 
  dplyr::count()  %>% 
  dplyr::group_by(motif_geno) %>% 
  dplyr::mutate(count_motif=sum(n)) %>% 
  dplyr::mutate(frac=100*n/count_motif) %>% 
  dplyr::mutate(diff=sub("(.*)(_score)","\\1",diff)) %>% 
  dplyr::group_by(motif_geno) %>% 
  dplyr::mutate(id=row_number(),
                count_motif2=ifelse(id==1,count_motif,NA))


fig_S5a <- ggplot(subset(expansion_contractionS_134,motif_length==1),aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  facet_wrap(.~motif_length,scales = "free") +
  coord_flip() +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=109,size = 10*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,113) ) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) 


fig_S5b <- ggplot(subset(expansion_contractionS_134,motif_length==3),aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  facet_wrap(.~motif_length,scales = "free") +
  coord_flip() +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=109,size = 10*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,113) ) +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 3)) 


fig_S5c <- ggplot(subset(expansion_contractionS_134,motif_length==4),aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  facet_wrap(.~motif_length,scales = "free") +
  coord_flip() +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=109,size = 10*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,113) ) +
  theme(legend.position = "none" ,
        axis.text.y = ggplot2::element_text(size=10,  color = "black"),) 


fig_S5ab <-  cowplot::plot_grid(fig_S5a, fig_S5b,
                                labels = c('', 'B'), 
                                rel_heights  =  c(1,4),
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "lr",
                                align = "v",
                                nrow = 2)



fig_S5 <-  cowplot::plot_grid(fig_S5ab, fig_S5c,
                              labels = c('A', 'C'  ), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "tb",
                              nrow = 1)

ggsave(fig_S5, filename = paste( "../figures/Supp_fig5_frac_expansion.png",sep = ""), units = "mm",height = 200, width = 170)




############# Figure  S6    ###############
#          constrained CDS                #
##########################################

###### fig_S6a  ######

pstr_fea <- table_s1_poly %>% dplyr::select(ref_STR,gfeature)

majorAF_ExpectedHe_fea <- majorAF_ExpectedHe %>% 
  dplyr::left_join(pstr_fea)%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(id=row_number(),
                mean_het=round(mean(hets),digits = 3),
                mean_het2=ifelse(id==1,round(mean(hets),digits = 2),NA)) %>% 
  dplyr::ungroup()

majorAF_ExpectedHe_fea$gfeatures<- factor(majorAF_ExpectedHe_fea$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

fig_S6a <- ggplot(majorAF_ExpectedHe_fea,aes(x=gfeatures,y=hets))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2)+ 
  geom_point(aes(x=gfeatures,y=mean_het2),size=0.5,color="red")+
  theme_cust +
  xlab("Genomic features")+
  ylab(expression(italic(H)[E]))+
  ylim(0,1) +
  theme(axis.text.x =  element_blank(),
        axis.title.x= element_blank() ) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.9,  p.adjust.method = "bonferroni", 
                             ref.group = "CDS")
 
 
###### fig_S6b  ######

repeat_var <- table_s1_poly %>% 
  dplyr::select(ref_STR,Chr,start,BPDIFFS,motif_length,gfeature) %>% 
  dplyr::arrange(Chr,start) %>% 
  splitstackshape::cSplit("BPDIFFS",",", direction = "long",sep = ",") %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate( repeatN_var=abs(BPDIFFS)/motif_length)  %>% 
  dplyr::group_by(gfeature,ref_STR) %>% 
  dplyr::mutate(mean_repeatN_var=mean(repeatN_var)) %>% 
  dplyr::distinct(gfeature,ref_STR,mean_repeatN_var ) %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(id=row_number(),
                gmean_repeatN_var=ifelse(id==1,round(mean(mean_repeatN_var),digits = 2),NA))%>% 
  dplyr::ungroup()



repeat_var$gfeatures<- factor(repeat_var$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

#plot

fig_S6b <- ggplot(repeat_var,aes(x=gfeatures,y=mean_repeatN_var))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=gfeatures,y=gmean_repeatN_var),size=0.5,color="red") +
  theme_cust +
  xlab("Genomic features")+
  ylab("Mean repeat\nnumber variance") +
  theme(axis.text.x =  element_blank(),
        axis.title.x= element_blank() ) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 14,  
                             p.adjust.method = "bonferroni", 
                             symnum.args = list(cutpoints = c(0, 0.0001,0.05,  1), 
                                                symbols = c("****","*",  "ns")),
                             ref.group = "CDS")+
  ylim(0,15)

 


###### fig_S6c  ######


poly_gc <- table_s1_poly  %>% 
  dplyr::mutate(countA=stringr::str_count(motif_geno, "A"),
                countG=stringr::str_count(motif_geno, "G"),
                countC=stringr::str_count(motif_geno, "C"),
                countT=stringr::str_count(motif_geno, "T"),
                contentGC=100*(countG+countC)/(countG+countC+countA+countT)) %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(id=row_number(),
                gmean_gc=ifelse(id==1,round(mean(contentGC),digits = 2),NA)) %>% 
  dplyr::ungroup() 

poly_gc$gfeatures<- factor(poly_gc$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

#plot

fig_S6c <- ggplot(poly_gc,aes(x=gfeatures,y=contentGC))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=gfeatures,y=gmean_gc),size=0.5,color="red") +
  theme_cust +
  xlab("Genomic features")+
  ylab("GC content (%)")+
  theme(axis.text.x =  element_text(size=10,   color = "black", angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 101,
                             p.adjust.method = "bonferroni", 
                             symnum.args = list(cutpoints = c(0, 0.00001,0.0001,0.05,  1), 
                                                symbols = c("****","***","**",  "ns")),
                             ref.group = "CDS")+
  ylim(0,105)

 
fig_S6 <-  cowplot::plot_grid(fig_S6a, fig_S6b,fig_S6c,  
                             labels = c('A' ,'B' ,'C'), 
                             rel_heights =  c(1,1,1.5 ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             align = "v",
                             nrow = 3)

ggsave(fig_S6, filename = paste( "../figures/Supp_fig6_constrainedCDS.png",sep = ""), units = "mm",height = 170, width = 170)



############# Figure  S7    ###############
#          MA_pSTR motif enrichment       #
##########################################

 
data_fig_S7 <- data.table::fread("../processed_data/table_s3_MA_pSTRs.txt")   

###### fig_S7a  ######


MAstr_motif_count  <- data_fig_S7 %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>%
  dplyr::group_by(motif_geno,motif_length) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::top_n(10, n)  %>% 
  dplyr::arrange(  n )

MAstr_motif_count$motif <- factor(MAstr_motif_count$motif_geno, levels = MAstr_motif_count$motif_geno )


fig_S7a <- ggplot(MAstr_motif_count,aes(x=motif,y=n ,fill=factor(motif_length))) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=period_size_color) +
  theme_cust +
  coord_flip() +
  theme(legend.position = "none") +
  labs(y="Number of sites",
       x="STR motif") +
  scale_y_continuous(breaks=c(0, 300,600) )

###### fig_S7b ######

MASTR_region <- data_fig_S7  %>% 
  dplyr::group_by(gfeature,motif_length) %>% 
  dplyr::add_count(name = "STRbyPS_REGION")%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::add_count(name = "STRby_REGION") %>% 
  dplyr::mutate(ps2region=100*STRbyPS_REGION/STRby_REGION) %>% 
  dplyr::distinct(gfeature,motif_length,STRbyPS_REGION,STRby_REGION,ps2region) %>% 
  dplyr::mutate(STRby_REGIONs=ifelse(motif_length==1,STRby_REGION,NA))


MASTR_region$gfeatures<- factor(MASTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

#plot

fig_S7b <- ggplot(MASTR_region,aes(x=gfeatures,y=ps2region,fill=factor(motif_length))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_cust +
  xlab("Genomic features")+
  ylab("Percent of pSTRs (%)")  +
  labs(fill="Motif\nlength")+
  scale_fill_manual(values=period_size_color) +
  geom_text(aes(label=STRby_REGIONs),y=115,size = 10*5/14)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_S7c ######

enrich_MASTR_region_stats <- MASTR_region %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(sig_n=STRbyPS_REGION,
                sig_total=STRby_REGION) %>% 
  dplyr::group_by(motif_length) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>% 
  dplyr::group_by(motif_length,gfeature)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.S6C") %>% 
  dplyr::mutate(group_factor="gfeature, motif_length",
                group_factor_catogory=paste0(gfeature,  ", ", motif_length),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrich_MASTR_region<-enrich_MASTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05) %>% 
  dplyr::mutate(logp=-log10(fisherp_adj)) 

enrich_MASTR_region$gfeatures<- factor(enrich_MASTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))


fig_S7c <- ggplot(enrich_MASTR_region,
                  aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F)+
  theme_cust+
  theme( legend.position = "none")+
  scale_color_manual(values=period_size_color) +
  scale_x_continuous(limits = c(0,100) )   +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p)))) 

###### fig_S7d ######


enrichMotif_MASTR_region_stats <-  data_fig_S7 %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>% 
  dplyr::group_by(gfeature, motif_geno ) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno )  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni"))%>% 
  dplyr::mutate(figure="FIG.S6C") %>% 
  dplyr::mutate(group_factor="gfeature, motif_geno",
                group_factor_catogory=paste0(gfeature,  ", ", motif_geno),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  


enrichMotif_MASTR_region <- enrichMotif_MASTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05)   %>% 
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::top_n(3, logp) 


enrichMotif_MASTR_region$gfeatures<- factor(enrichMotif_MASTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))



fig_S7d <- ggplot(enrichMotif_MASTR_region,
                  aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 2,segment.linetype=6, size=10*5/14,
                           max.overlaps=Inf,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p))))  




###### fig S7 #####

fig_S7ab <-  cowplot::plot_grid(fig_S7a, fig_S7b,
                                labels = c('', 'B'), 
                                rel_widths =  c(1.3,2),
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "t",
                                nrow = 1)

fig_S7cd <-  cowplot::plot_grid(fig_S7c,fig_S7d,
                                labels = c('', 'D'), 
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "t",
                                nrow = 1)


fig_S7 <-  cowplot::plot_grid(fig_S7ab, fig_S7cd,  
                              labels = c('A', 'C'  ), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow = 2)

ggsave(fig_S7, filename = paste( "../figures/Supp_fig7_MAmotif.png",sep = ""), units = "mm",height = 160, width = 170)





 
############# Figure  S8    ###############
#          MA ps u                       #
##########################################

###### fig_S8  ######
 
data_fig_S8 <- data_fig_4 %>% 
  dplyr::filter(motif_length!="1-6" &
                  motif_geno=="All" &
                  gfeature=="All" &
                  OMA=="O1MA") 


data_fig_S8_stats <-ggpubr::compare_means( mutation_rate ~ strain, 
                       data= data_fig_S8 ,
                       group.by = c( "mutation" ,"motif_length"),  
                       p.adjust.method = "bonferroni", 
                       label = "p.signif", 
                       method = "wilcox.test" ) %>% 
  dplyr::mutate(figure="FIG.S7") %>% 
  dplyr::select(-p.format,-p.signif) %>% 
  dplyr::mutate(group_factor="mutation, motif_length",
                group_factor_catogory=paste0(mutation,  ", ", motif_length),
                method="two-sided Wilcoxon test",
                padjustment="BF") %>% 
  dplyr::select(figure,method,group_factor,group_factor_catogory,'.y.', group1,group2,p,padjustment,p.adj)





fig_S8 <- ggpubr::ggboxplot(data_fig_S8, x="mutation",y="mutation_rate",outlier.shape = NA,
                            color="strain"  ) +
  geom_point( position = position_jitterdodge(jitter.width = 0.2) ,aes(color=strain), size=0.5, alpha=0.8)+
  theme_cust +
  theme(legend.position = "bottom")+
  facet_grid( motif_length~.,scales = "free")+
  scale_color_manual(values = c("orange","#007e2f","#ffcd12","#721b3e") ) +
  labs(x="Mutations",y="ANC-O1MA mutation rate",color="Strain" ) + 
  ggpubr::stat_compare_means( 
    aes(group = strain),
    label = "p.signif",  
    size = 3,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 1), 
                       symbols = c("****","**","*",  "ns")),
    method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),expand = c(0.1, 0 ) )   

ggsave(fig_S8, filename = paste( "../figures/Supp_fig8_MA_ps_u.png",sep = ""), units = "mm",height = 200, width = 170)


 
 
 ############# Figure  S9    ###############
 #          pSTR pxg                      #
 ##########################################
 
 data_fig_S9 <- data.table::fread("../processed_data/Organismal_Lrt_STRs_pxg.tsv")
 
 
 
 traits_name_df <- data.frame(organismal_trait=c("broods","abamectin_norm.n","arsenic_pc1","Albendazole_q90.TOF",
                                                 "bleo_medianEXT","dauer_frac","amsacrine_f.L1","etoposide_median.TOF",
                                                 "propionicL1survival","telomere_resids","zinc_norm_EXT"),
                              trait_name=c("Lifetime fecundity","Response to abamectin","Response to arsenic","Response to albendazole",
                                           "Response to bleomycin","Dauer formation","Response to amsacrine","Response to etoposide",
                                           "Response to propionate","Telomere length","Response to zinc") )
 
 
 
 data_fig_S9_N2 <- data_fig_S9 %>% 
   dplyr::left_join(traits_name_df) %>% 
   dplyr::group_by(organismal_trait) %>% 
   dplyr::filter(adjusted_p==min(adjusted_p)) %>% 
   dplyr::ungroup() %>% 
   dplyr::group_by(strain) %>% 
   dplyr::add_count() %>% 
   dplyr::filter( strain=="N2") %>% 
   dplyr::select(organismal_trait, pSTR,ref_strain=strain,ref_str_length=STR_length)
 
 data_fig_S9_top <- data_fig_S9 %>% 
   dplyr::left_join(traits_name_df) %>% 
   dplyr::group_by(organismal_trait) %>% 
   dplyr::filter(adjusted_p==min(adjusted_p)) %>% 
   dplyr::ungroup() %>% 
   dplyr::left_join(data_fig_S9_N2) %>% 
   dplyr::mutate(allele=ifelse(is.na(ref_strain) & STR_length==18, "REF",
                               ifelse(is.na(ref_strain) & STR_length==15, "ALT", 
                                      ifelse(ref_str_length==STR_length,"REF","ALT"))) ,
                 trait_name=paste0(trait_name,"\n",pSTR))
 
 
 
 data_fig_S9_top$trait_name2<- factor(data_fig_S9_top$trait_name,levels = c("Lifetime fecundity\nSTR_25031","Dauer formation\nSTR_30772",
                                                                            "Telomere length\nSTR_10678","Response to abamectin\nSTR_25415",
                                                                            "Response to albendazole\nSTR_6418","Response to arsenic\nSTR_8043",
                                                                            "Response to bleomycin\nSTR_31790","Response to amsacrine\nSTR_1102",
                                                                            "Response to propionate\nSTR_29860"))
 
 
 fig_S9 <- ggplot(data_fig_S9_top,aes(x=as.factor(STR_length),y=trait_value ))+
   geom_jitter(shape=21,position=position_jitter(0.2),aes(fill=as.factor(allele)) ) +
   geom_boxplot( outlier.shape = NA,alpha=0.5) +
   theme_cust+
   ggplot2::scale_fill_manual(values = c("REF"="orange","ALT"="blue") )+
   xlab("Length of STR") + ylab("Trait values") +  
   facet_wrap(.~trait_name2 ,scales = "free",nrow=3) +
   theme( legend.position = "none"  )
 
 
 ggsave(fig_S9, filename = paste( "../figures/Supp_fig9_pxg.png",sep = ""), units = "mm",height =  180, width = 170)
  
 
#######