setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/"))
#setwd("../")

### load R packages####

library(tidyverse)

### ggplot theme ####

theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=12,  color = "black"),
        legend.text = ggplot2::element_text(size=12,  color = "black"),
        legend.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.title =  ggplot2::element_text(size=12,  color = "black"),
        axis.text =  ggplot2::element_text(size=12,  color = "black"),
        strip.text = ggplot2::element_text(size=12, vjust = 1,  color = "black"),
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
  sis <- pt_deg$sig_n # success-in-sample  ## number of positions with variants in arm
  silp <- pt_deg$n - sis # success-in-left-part ## number of positions with variants in other parts of chromosomes
  fis <- pt_deg$sig_total - sis   # failure-in-sample  ## number of positions with no variants in arm
  filp <- pt_deg$total - pt_deg$sig_total - silp # failure-in-left-part ## number of positions with no variants in other parts of chromosomes
  ftp <- fisher.test(matrix(c(sis,silp,fis,filp), 2, 2), alternative='greater')
  return(ftp$p.value)
}

############# Figure  1    ###############
#          pSTR distribution             #
##########################################

table_s1 <- data.table::fread("../processed_data/table_s1_refSTRs_pSTRs.tsv")

table_s1_poly <- table_s1 %>% 
  dplyr::filter(polymorphic=="Yes")

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
  dplyr::group_by(motif_geno,motif_length) %>% 
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
  scale_y_continuous(breaks=c(0, 750,1500),limits = c(0,1600) )

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
  geom_text(aes(label=STRby_REGIONs),y=115,size = 18/5)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_1d ######

enrich_PolySTR_region <- PolySTR_region %>% 
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


enrichMotif_polySTR_region <-  table_s1_poly %>% 
  dplyr::group_by(gfeature, motif_geno) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::filter(fisherp_adj<0.05)   %>% 
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::top_n(10, logp) 


enrichMotif_polySTR_region$gfeatures<- factor(enrichMotif_polySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))



fig_1e <- ggplot(enrichMotif_polySTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 3,segment.linetype=6,
                           max.overlaps=Inf,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic\nfeatures")+
  xlab(expression(-log[10](italic(p)))) 




###### fig 1 #####

fig_1bc <-  cowplot::plot_grid(fig_1b, fig_1c,
                              labels = c('', 'C'), 
                              rel_widths =  c(1.2,2),
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
 # dplyr::left_join(hipstr_reference)  %>% 
 # dplyr::left_join(gfeat) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(diff,motif_length) %>% 
  dplyr::mutate(mean_sc=mean(score)) %>% 
  dplyr::mutate(id=row_number(),
                mean_sc=ifelse(id==1,round(mean(mean_sc),digits = 2),NA))


fig_2b <- ggplot(expansion_contractionS ,aes(score,fill=diff)) +
  geom_histogram(color="black",# fill="lightskyblue2",
                 bins =100,size=0.2) +
  theme_cust +
  xlab("Contraction / Expansion score")+
  ylab("Number of\nSTRs")   +
  scale_fill_manual(values=c("#E7B800", "#00AFBB")) +
  theme( legend.position = "none") 

###### fig_2c ###### 

expansion_contractionS_abs <- expansion_contractionS %>% dplyr::mutate(score=abs(score))


fig_2c <- ggplot(expansion_contractionS_abs,aes(x=diff,y=abs(score),color=diff))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=abs(mean_sc)),size=0.5,color="red") +
  theme_cust +
  ylab("Contraction /\nExpansion\nabsolute scores") +
  ggpubr::stat_compare_means( comparisons = list(c("contraction_score","expansion_score")), 
                              label = "p.signif", method = "wilcox.test", 
                              p.adjust.method = "bonferroni",
                              symnum.args = list(cutpoints = c(0, 0.00001, 0.0001, 0.001,  1), 
                                                 symbols = c("****","***","**",  "ns")),
                              # method.args = list(alternative = "less" ),
                              label.y = 2.3) +
  facet_grid(.~motif_length,scales="free") +
  scale_color_manual(values=c("#E7B800", "#00AFBB")) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        #    strip.text = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.7))

ggpubr::compare_means( score ~ diff, 
                       data= expansion_contractionS_abs ,
                       group.by = c("motif_length"),  
                       # comparisons = list(c("Contraction","Expansion")), 
                       # method.args = list(alternative = "greater" ),
                       p.adjust.method = "bonferroni", 
                       label = "p.signif", 
                       method = "wilcox.test" )

###### fig_2d ###### 

expand_frac <- table_s1_poly %>% 
  dplyr::select(ref_STR ,Contraction_frac,Expansion_frac,motif_length) %>% 
  tidyr::gather(diff,frac,-ref_STR,-motif_length) %>% 
  dplyr::mutate(frac=frac/100) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(diff,motif_length)%>% 
  dplyr::mutate(mean_frac=mean(frac)) %>% 
  dplyr::mutate(id=row_number(),
                mean_frac=ifelse(id==1,round(mean(mean_frac),digits = 2),NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate( diff=sub("(.*)(_frac)","\\1",diff))
# dplyr::filter(!frac==0)



fig_2d <- ggplot(expand_frac,aes(x=diff,y=frac ,color=diff))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=abs(mean_frac)),size=0.5,color="red") +
  theme_cust +
  ylab("Allele frequency") +
  xlab("Variation to median allele")+
  ggpubr::stat_compare_means( comparisons = list(c("Contraction","Expansion")), 
                              label = "p.signif",
                              label.y = 0.53 , 
                              #   hide.ns = TRUE,
                              # method.args = list(alternative = "greater" ),
                              p.adjust.method = "bonferroni",
                              symnum.args = list(cutpoints = c(0,  0.001,  1), 
                                                 symbols = c("****",   "ns")),
                              method = "wilcox.test" 
  ) +
  facet_grid(.~motif_length,scales="free") +
  scale_color_manual(values=c("#E7B800", "#00AFBB")) +
  theme(axis.text.x = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.6)) 




ggpubr::compare_means(frac ~ diff,
                      comparisons = list(c("Contraction","Expansion")),
                      expand_frac,
                      method.args = list(alternative = "greater"), 
                      method = "wilcox.test", 
                      p.adjust.method = "bonferroni", 
                      group.by = c("motif_length"))


###### fig_2e ###### 

expansion_contractionS_di <- table_s1_poly %>% 
  dplyr::select(ref_STR ,expansion_score,contraction_score,motif_length,motif_geno) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length,-motif_geno) %>% 
  dplyr::filter(motif_length==2) %>% 
  dplyr::mutate(diff=ifelse(score==0,"substitution",diff)) %>% 
  dplyr::group_by(motif_geno,diff,motif_length) %>% 
  dplyr::count()  %>% 
  dplyr::filter(!grepl("/",motif_geno)) %>% 
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
#          MA lines                      #
##########################################

data_fig_3 <- data.table::fread("../processed_data/MA174_pSTRs_mutationRate.tsv") %>% 
  dplyr::group_by(comparison,strain) %>% 
  dplyr::mutate(mean_gen=round(mean(N_generation),digits = 0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(motif_length=ifelse(motif_length=="Motif sizes 1-6", "1-6",motif_length)) %>% 
  dplyr::mutate(comparison=ifelse(comparison=="Ancestor-O1MA","ANC-O1MA",
                                  ifelse(comparison=="Ancestor-O2MA","ANC-O2MA",comparison))) %>% 
  dplyr::mutate(strain_com=paste(strain,comparison,mean_gen,sep="\n"),
                comparison2=paste( comparison,mean_gen,sep="\n"))

 

# wilcox test

data_fig_3_sig <- ggpubr::compare_means( mutation_rate ~ mutation, 
                                         data= subset(data_fig_3, mutation %in% c("insertions","deletions")) ,
                                         group.by = c("motif_length","strain","comparison"),  
                                         p.adjust.method = "bonferroni", 
                                         label = "p.signif", 
                                         method = "wilcox.test" ) %>% 
  dplyr::filter(p.adj<0.05)


 
#### fig3a ####

O1MA_ps_all <- subset(data_fig_3, motif_length=="1-6"  ) %>% 
  dplyr::group_by(strain,comparison)  


fig_3a <- ggplot(O1MA_ps_all,aes(x=mutation,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8,color="gray69") +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=mutation)) +
  facet_grid(motif_length~strain_com ,scales = "free") +
  theme_cust+
  labs(x="Mutation",y="Mutation rate", color="Mutation")+
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1, 0, 1, 3), "mm"),
        axis.title.y = element_blank() ,
        axis.title.x = element_blank()  ) +
  scale_color_manual(values = c("#007e2f","#ffcd12","#721b3e") ) +
  ggpubr::stat_compare_means( comparisons = list(c("insertions","deletions")),
                              p.adjust.method = "bonferroni",
                              label.y = 8e-05,
                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001,  1), 
                                                  symbols = c("****","**",  "ns")),
                              label = "p.signif", 
                              method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1e-04), breaks = c(0, 5e-05, 1e-04))   


  


#### fig3bcde ####


O1MA_ps  <- subset(data_fig_3, motif_length!="1-6"  ) %>% 
  dplyr::group_by(strain,comparison)  

O1MA_ps$strain_com2<- factor(O1MA_ps$strain_com,levels = c("N2\nANC-O1MA\n234", "N2\nANC-O2MA\n373","N2\nO1MA-O2MA\n143","PB306\nANC-O1MA\n227" ) )


fig_3_ps1 <- ggplot(subset(O1MA_ps,    motif_length=="1" ),aes(x=mutation,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8,color="gray69") +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=mutation)) +
  facet_grid(motif_length~strain_com2 ,scales = "free") +
  theme_cust+
  labs(x="Mutation",y="Mutation rate", color="Mutation") +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        #  plot.margin = unit(c(1, 0, 1, 1), "mm"),
        axis.title.y = element_blank() ,
        axis.title.x = element_blank() ,
        strip.text.x = element_blank() ) +
  scale_color_manual(values = c("#007e2f","#ffcd12","#721b3e") ) +
  ggpubr::stat_compare_means( comparisons = list(c("insertions","deletions")),
                              p.adjust.method = "bonferroni",
                              symnum.args = list(cutpoints = c(0, 0.000001, 0.0001, 1), 
                                                 symbols = c("****","***",  "ns")),
                              label.y = 8e-05,
                              label = "p.signif", 
                              method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1e-04), breaks = c(0,5e-05, 1e-04))   



fig_3_ps2 <- ggplot(subset(O1MA_ps,    motif_length=="2" ),aes(x=mutation,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8,color="gray69") +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=mutation)) +
  facet_grid(motif_length~strain_com2 ,scales = "free") +
  theme_cust+
  labs(x="Mutation",y="Mutation rate", color="Mutation") +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        #  plot.margin = unit(c(1, 0, 1, 1), "mm"),
        # axis.title.y = element_blank() ,
        axis.title.x = element_blank() ,
        strip.text.x = element_blank() ) +
  scale_color_manual(values = c("#007e2f","#ffcd12","#721b3e") ) +
  ggpubr::stat_compare_means( comparisons = list(c("insertions","deletions")),
                              p.adjust.method = "bonferroni",
                              label.y = 8e-06,
                              symnum.args = list(cutpoints = c(0, 0.00001, 0.01, 1), 
                                                 symbols = c("****","*",  "ns")),
                              label = "p.signif", 
                              method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1e-05), breaks = c(0,5e-06, 1e-05))   



fig_3_ps3 <- ggplot(subset(O1MA_ps,    motif_length=="3" ),aes(x=mutation,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8,color="gray69") +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=mutation)) +
  facet_grid(motif_length~strain_com2 ,scales = "free") +
  theme_cust+
  labs(x="Mutation",y="Mutation rate", color="Mutation") +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        #  plot.margin = unit(c(1, 0, 1, 1), "mm"),
        axis.title.y = element_blank() ,
        axis.title.x = element_blank() ,
        strip.text.x = element_blank() ) +
  scale_color_manual(values = c("#007e2f","#ffcd12","#721b3e") ) +
  ggpubr::stat_compare_means( comparisons = list(c("insertions","deletions")),
                              p.adjust.method = "bonferroni",
                              label.y = 2e-06,
                              symnum.args = list(cutpoints = c(0,   0.0001,0.001, 1), 
                                                 symbols = c("**","*",  "ns")),
                              label = "p.signif", 
                              method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,2.5e-06), breaks = c(0,1e-06,  2e-06))   




fig_3_ps4 <- ggplot(subset(O1MA_ps,    motif_length=="4" ),aes(x=mutation,y=mutation_rate))+
  geom_jitter( shape=20,position=position_jitter(0.4), size=1, alpha=0.8,color="gray69") +
  geom_boxplot(outlier.shape = NA,fill=NA ,aes(color=mutation)) +
  facet_grid(motif_length~strain_com2 ,scales = "free") +
  theme_cust+
  labs(x="Mutation",y="Mutation rate", color="Mutation") +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom"  ,
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0),
        #  plot.margin = unit(c(1, 0, 1, 1), "mm"),
        axis.title.y = element_blank() ,
        axis.title.x = element_blank() ,
        strip.text.x = element_blank() ) +
  scale_color_manual(values = c("#007e2f","#ffcd12","#721b3e") ) +
  ggpubr::stat_compare_means( comparisons = list(c("insertions","deletions")),
                              p.adjust.method = "bonferroni",
                              label.y = 0.8e-06,
                              symnum.args = list(cutpoints = c(0,   0.001, 1), 
                                                 symbols = c("*",  "ns")),
                              label = "p.signif", 
                              method = "wilcox.test") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),limits = c(0,1e-06), breaks = c(0, 5e-07,1e-06 ))   





fig3 <-  cowplot::plot_grid(  fig_3a, fig_3_ps1,  fig_3_ps2,fig_3_ps3,fig_3_ps4 ,
                                     # labels = c('A', 'B','C' ), 
                                     rel_heights =  c( 1.5,1 ,1,1,1.2 ),
                                     label_size = 12, 
                                     label_fontfamily="Helvetica",
                                     align = "v",
                                     axis = "lr",
                                     nrow =5)
#fig3_ps1234

######

 

ggsave(fig3, filename = paste( "../figures/Fig_3.png",sep = ""), units = "mm",height = 180, width = 170)








############# Figure  4    ###############
#             popgenet                    #
##########################################

###### fig_4a ######

num_allele <- bp_diff_BPD %>% # fig_2a
  dplyr::group_by(ref_STR) %>% 
  dplyr::count() %>% 
  dplyr::mutate(n=n+1)

fig_4a <- ggplot(num_allele,aes(n))   + 
  geom_histogram(color="black", fill="white",bins =21,size=0.2) +
  theme_cust +
  xlab("Number of alleles\nper pSTR") +
  ylab("Number of pSTRs") + 
  scale_y_continuous(breaks=c(0, 3000, 1500), limits=c(0, 3000))+ 
  scale_x_continuous(breaks=c(2,10,21) )  +
  theme( axis.title.y = ggplot2::element_text(size=12,  color = "black", hjust = 0.2 )) 
 
###### fig_4b ######

majorAF_ExpectedHe <- data.table::fread("../processed_data/majorAF_ExpectedHe.tsv")

majorAF_data <- majorAF_ExpectedHe %>% 
  dplyr::group_by(ref_STR,n_st) %>% 
  dplyr::mutate( major_af=max(af)) %>% 
  dplyr::filter(af==major_af)  


fig_4b <- ggplot()   + 
  ggplot2::stat_density(data=majorAF_data,aes(x=major_af,color=n_st ), geom="line",position = "identity",#linetype="dashed"
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

###### fig_4c ######

strain_ALT_frac <- data.table::fread("../processed_data/strain_ALT_frac.tsv")

fig_4c <- ggplot()+
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

###### fig_4d ######

pca_pSTR_SNV <- data.table::fread("../processed_data/pca_pSTR_SNV.tsv")

pca_pSTR <- pca_pSTR_SNV %>% 
  dplyr::filter(data=="pSTRs")

fig_4d <- ggplot(pca_pSTR, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=cluster ),size=1,alpha=0.8) +
  theme_cust + 
  theme(legend.position="none")+
  labs(x=paste0("PC1: ",unique(pca_pSTR$PC1_var_exp)[1],"%"),
       y=paste0("PC2: ",unique(pca_pSTR$PC2_var_exp)[1],"%")) +
  ggtitle("9,691 pSTRs") +
  scale_color_manual(values=ancestry.colours)




###### fig_4e ######

pca_SNV <- pca_pSTR_SNV %>% 
  dplyr::filter(data=="SNVs")


fig_4e <- ggplot(pca_SNV, aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=cluster ),size=1,alpha=0.8) +
  theme_cust + 
  labs(x=paste0("PC1: ",unique(pca_SNV$PC1_var_exp)[1],"%"),
       y=paste0("PC2: ",unique(pca_SNV$PC2_var_exp)[1],"%"),
       color = "Sample locations") +
  ggtitle("13,650 SNVs") +
  scale_color_manual(values=ancestry.colours)



###### fig_4f ######

fig_4f <- ggpubr::ggscatter(subset(majorAF_ExpectedHe,Chr != "MtDNA"), 
                                      x = "start", y = "hets",  point=FALSE, palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                                      add = "loess", conf.int = FALSE,color="n_st" ) + 
  labs(x= "Genomic position (Mb)",
       y=#"Expected\nHeterozygosity"
         expression(italic(H)[E])) +
  facet_grid(.~Chr,scales="free")+
  theme_cust +
  theme( axis.text.x = element_blank(), 
         legend.position = "none",
         panel.spacing = unit(0.2,"line"))


###### fig_4 ######

fig4abc <-  cowplot::plot_grid(fig_4a, fig_4b,fig_4c,
                               labels = c('', 'B','C'), 
                               rel_widths =  c(0.9,0.9 , 1.2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "tb",
                               align = "h",
                               nrow = 1)

fig4de <- cowplot::plot_grid(fig_4d,  fig_4e,   
                             labels = c('', 'E'  ), 
                             rel_widths = c( 1,1.75  ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             # axis = "tblr",
                             # align = "v",
                             nrow = 1)

fig4 <-  cowplot::plot_grid(fig4abc, fig4de, fig_4f,  
                            labels = c('A' ,'D','F' ), 
                            rel_heights =  c(1.2,1.5,1.1 ),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            #  axis = "lr",
                            align = "v",
                            nrow = 3)

ggsave(fig4, filename = paste( "../figures/Fig_4.png",sep = ""), units = "mm",height = 180, width = 170)

 

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
        plot.title =   ggplot2::element_text(size=12,  color = "black"))


fig_S1 <-  cowplot::plot_grid(fig_S1a, fig_S1b,  
                              labels = c('A', 'B'  ), 
                              rel_heights =  c(1.3,4  ),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow = 2)

ggsave(fig_S1, filename = paste( "../figures/Supp_fig1_ref_dist.png",sep = ""), units = "mm",height = 180, width = 170)




############# Figure  S2    ###############
#          pSTR distribution             #
##########################################

fig_S2 <- ggplot() + 
  geom_histogram(data=str_dist_polym,aes(x=start/1e6,fill=factor(motif_length) ),bins = 50 )  + 
  facet_grid(motif_length~Chr,scales = "free")+
  theme_cust +
  xlab("Genomic Position (Mb)")+
  ylab("Number of polymorphic STRs")  +
  scale_fill_manual(values=period_size_color) +
  theme(legend.position = "none" ,
        plot.title =   ggplot2::element_text(size=12,  color = "black"))

ggsave(fig_S2, filename = paste( "../figures/Supp_fig2_polystr_dist_size.png",sep = ""), units = "mm",height = 160, width = 170)




############# Figure  S3    ###############
#          refSTR motif enrichment       #
##########################################
 
###### fig_S3a  ######

refstr_motif_count  <- table_s1 %>% 
  dplyr::group_by(motif_geno,motif_length) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::top_n(10, n)  %>% 
  dplyr::arrange(  n )

refstr_motif_count$motif <- factor(refstr_motif_count$motif_geno, levels = refstr_motif_count$motif_geno )


fig_S3a <- ggplot(refstr_motif_count,aes(x=motif,y=n ,fill=factor(motif_length))) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=period_size_color) +
  theme_cust +
  coord_flip() +
  theme(legend.position = "none") +
  labs(y="Number of sites",
       x="STR motif") +
  scale_y_continuous(breaks=c(0, 4000,8000), limits = c(0,9000))

###### fig_S3b ######

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

fig_S3b <- ggplot(refSTR_region,aes(x=gfeatures,y=ps2region,fill=factor(motif_length))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_cust +
  xlab("Genomic features")+
  ylab("Percent of pSTRs (%)")  +
  labs(fill="Motif\nlength")+
  scale_fill_manual(values=period_size_color) +
  geom_text(aes(label=STRby_REGIONs),y=115,size = 18/5)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_S3c ######

enrich_refSTR_region <- refSTR_region %>% 
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
  dplyr::filter(fisherp_adj<0.05) %>% 
  dplyr::mutate(logp=-log10(fisherp_adj))  %>% 
  dplyr::mutate(logp=ifelse(logp=="Inf",400,logp))  

enrich_refSTR_region$gfeatures<- factor(enrich_refSTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))


fig_S3c <- ggplot(enrich_refSTR_region,
                aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F)+
  theme_cust+
  theme( legend.position = "none")+
  scale_color_manual(values=period_size_color) +
  scale_x_continuous(breaks = c(0,100,200,300,400),labels = c("0","100","200","300","Inf")  )   +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p)))) 

###### fig_S3d ######


enrichMotif_refSTR_region <-  table_s1 %>% 
  dplyr::group_by(gfeature, motif_geno) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::filter(fisherp_adj<0.05)   %>% 
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::top_n(10, logp) 


enrichMotif_refSTR_region$gfeatures<- factor(enrichMotif_refSTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))



fig_S3d <- ggplot(enrichMotif_refSTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 3,segment.linetype=6,
                           max.overlaps=Inf,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p)))) 




###### fig S3 #####

fig_S3ab <-  cowplot::plot_grid(fig_S3a, fig_S3b,
                               labels = c('', 'B'), 
                               rel_widths =  c(1.3,2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)

fig_S3cd <-  cowplot::plot_grid(fig_S3c,fig_S3d,
                               labels = c('', 'D'), 
                             #  rel_widths =  c(1.2 ,2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)


fig_S3 <-  cowplot::plot_grid(fig_S3ab, fig_S3cd,  
                             labels = c('A', 'C'  ), 
                           #  rel_heights =  c(1,1.5,1 ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             axis = "lr",
                             nrow = 2)

ggsave(fig_S3, filename = paste( "../figures/Supp_fig3_refmotif.png",sep = ""), units = "mm",height = 160, width = 170)






############# Figure  S4    ###############
#           motif 134                     #
##########################################

expansion_contractionS_134 <- table_s1_poly %>% 
  dplyr::select(ref_STR ,expansion_score,contraction_score,motif_length,motif_geno) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length,-motif_geno) %>% 
  dplyr::filter(motif_length %in% c(1,3,4)) %>% 
  dplyr::mutate(diff=ifelse(score==0,"substitution",diff)) %>% 
  dplyr::group_by(motif_geno,diff,motif_length) %>% 
  dplyr::count()  %>% 
  dplyr::filter(!grepl("/",motif_geno)) %>% 
  dplyr::group_by(motif_geno) %>% 
  dplyr::mutate(count_motif=sum(n)) %>% 
  dplyr::mutate(frac=100*n/count_motif) %>% 
  dplyr::mutate(diff=sub("(.*)(_score)","\\1",diff)) %>% 
  dplyr::group_by(motif_geno) %>% 
  dplyr::mutate(id=row_number(),
                count_motif2=ifelse(id==1,count_motif,NA))


fig_S4a <- ggplot(subset(expansion_contractionS_134,motif_length==1),aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  facet_wrap(.~motif_length,scales = "free") +
  coord_flip() +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=109,size = 12*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,113) ) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) 


fig_S4b <- ggplot(subset(expansion_contractionS_134,motif_length==3),aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  facet_wrap(.~motif_length,scales = "free") +
  coord_flip() +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=109,size = 12*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,113) ) +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 3)) 


fig_S4c <- ggplot(subset(expansion_contractionS_134,motif_length==4),aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  facet_wrap(.~motif_length,scales = "free") +
  coord_flip() +
  scale_fill_manual(values=c("#E7B800", "#00AFBB","gray69")) +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=109,size = 12*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,113) ) +
  theme(legend.position = "none" ) 


fig_S4ab <-  cowplot::plot_grid(fig_S4a, fig_S4b,
                                labels = c('', 'B'), 
                                rel_heights  =  c(1,4),
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "lr",
                                align = "v",
                                nrow = 2)



fig_S4 <-  cowplot::plot_grid(fig_S4ab, fig_S4c,
                              labels = c('A', 'C'  ), 
                              # rel_heights =  c(1,1.5,1.2),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "tb",
                              #   align = "v",
                              nrow = 1)

ggsave(fig_S4, filename = paste( "../figures/Supp_fig4_frac_expansion.png",sep = ""), units = "mm",height = 200, width = 170)





############# Figure  S5    ###############
#          MA_pSTR motif enrichment       #
##########################################

data_fig_s5 <- data.table::fread("../processed_data/MA174_pSTRs.tsv")

###### fig_S5a  ######
 

MAstr_motif_count  <- data_fig_s5 %>% 
  dplyr::group_by(motif_geno,motif_length) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::top_n(10, n)  %>% 
  dplyr::arrange(  n )

MAstr_motif_count$motif <- factor(MAstr_motif_count$motif_geno, levels = MAstr_motif_count$motif_geno )


fig_S5a <- ggplot(MAstr_motif_count,aes(x=motif,y=n ,fill=factor(motif_length))) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=period_size_color) +
  theme_cust +
  coord_flip() +
  theme(legend.position = "none") +
  labs(y="Number of sites",
       x="STR motif") +
  scale_y_continuous(breaks=c(0, 300,600) )

###### fig_S5b ######

MASTR_region <- data_fig_s5  %>% 
  dplyr::group_by(gfeature,motif_length) %>% 
  dplyr::add_count(name = "STRbyPS_REGION")%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::add_count(name = "STRby_REGION") %>% 
  dplyr::mutate(ps2region=100*STRbyPS_REGION/STRby_REGION) %>% 
  dplyr::distinct(gfeature,motif_length,STRbyPS_REGION,STRby_REGION,ps2region) %>% 
  dplyr::mutate(STRby_REGIONs=ifelse(motif_length==1,STRby_REGION,NA))


MASTR_region$gfeatures<- factor(MASTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))

#plot

fig_S5b <- ggplot(MASTR_region,aes(x=gfeatures,y=ps2region,fill=factor(motif_length))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_cust +
  xlab("Genomic features")+
  ylab("Percent of pSTRs (%)")  +
  labs(fill="Motif\nlength")+
  scale_fill_manual(values=period_size_color) +
  geom_text(aes(label=STRby_REGIONs),y=115,size = 18/5)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_S5c ######

enrich_MASTR_region <- MASTR_region %>% 
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
  dplyr::filter(fisherp_adj<0.05) %>% 
  dplyr::mutate(logp=-log10(fisherp_adj))  #%>% 
 # dplyr::mutate(logp=ifelse(logp=="Inf",400,logp))  

enrich_MASTR_region$gfeatures<- factor(enrich_MASTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))


fig_S5c <- ggplot(enrich_MASTR_region,
                  aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F)+
  theme_cust+
  theme( legend.position = "none")+
  scale_color_manual(values=period_size_color) +
#  scale_x_continuous(limits = c(0,30) )   +
  scale_x_continuous(limits = c(0,100) )   +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p)))) 

###### fig_S5d ######


enrichMotif_MASTR_region <-  data_fig_s5 %>% 
  dplyr::group_by(gfeature, motif_geno) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::filter(fisherp_adj<0.05)   %>% 
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::top_n(10, logp) 


enrichMotif_MASTR_region$gfeatures<- factor(enrichMotif_MASTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","pseudogene","RNAs & TEs","intergenic"))



fig_S5d <- ggplot(enrichMotif_MASTR_region,
                  aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 2,segment.linetype=6,
                           max.overlaps=Inf,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic features")+
  xlab(expression(-log[10](italic(p))))  




###### fig S5 #####

fig_S5ab <-  cowplot::plot_grid(fig_S5a, fig_S5b,
                                labels = c('', 'B'), 
                                rel_widths =  c(1.3,2),
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "t",
                                nrow = 1)

fig_S5cd <-  cowplot::plot_grid(fig_S5c,fig_S5d,
                                labels = c('', 'D'), 
                                #  rel_widths =  c(1.2 ,2),
                                label_size = 12, 
                                label_fontfamily="Helvetica",
                                axis = "t",
                                nrow = 1)


fig_S5 <-  cowplot::plot_grid(fig_S5ab, fig_S5cd,  
                              labels = c('A', 'C'  ), 
                              #  rel_heights =  c(1,1.5,1 ),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              nrow = 2)

ggsave(fig_S5, filename = paste( "../figures/Supp_fig5_MAmotif.png",sep = ""), units = "mm",height = 160, width = 170)





############# Figure  S6    ###############
#          constrained CDS                #
##########################################

###### fig_S6a  ######

pstr_fea <- table_s1_poly %>% dplyr::select(ref_STR,gfeature)

majorAF_ExpectedHe <- data.table::fread("../processed_data/majorAF_ExpectedHe.tsv")

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

ggpubr::compare_means(hets~gfeatures,
                      majorAF_ExpectedHe_fea, 
                      p.adjust.method = "bonferroni", 
                      method = "wilcox.test", 
                      ref.group = "CDS" )

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


ggpubr::compare_means(mean_repeatN_var~gfeatures,
                      repeat_var, 
                      p.adjust.method = "bonferroni", 
                      method = "wilcox.test", 
                      ref.group = "CDS" )


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
  theme(axis.text.x =  element_text(size=12,   color = "black", angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 101,
                             p.adjust.method = "bonferroni", 
                             symnum.args = list(cutpoints = c(0, 0.00001,0.0001,0.05,  1), 
                                                symbols = c("****","***","**",  "ns")),
                             ref.group = "CDS")+
  ylim(0,105)



ggpubr::compare_means(contentGC~gfeatures,
                      poly_gc,  p.adjust.method = "bonferroni", 
                      method = "wilcox.test", 
                      ref.group = "CDS" )



figS6 <-  cowplot::plot_grid(fig_S6a, fig_S6b,fig_S6c,  
                             labels = c('A' ,'B' ,'C'), 
                             rel_heights =  c(1,1,1.5 ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             #  axis = "lr",
                             align = "v",
                             nrow = 3)

ggsave(figS6, filename = paste( "../figures/Supp_fig6_constrainedCDS.png",sep = ""), units = "mm",height = 170, width = 170)


 

#######