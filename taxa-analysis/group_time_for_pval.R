#####important notes!!!!!!
# cutoff abundance for this analysis session is 0.05
#set working directory
setwd()

#load needed libraries
library(reshape2) #for melt
library(ggalluvial) #for plot genus changes over time
library(tidyverse) #for readr: read_csv
library(stringr) # for function str_sub
library(ggsci) #for palette color pick
library(ggpubr) #for ggarrange common.legends
library(cowplot) #for plot_grid
library(magick) #for clipping plots 
library(gridExtra) #for grid.arrange
library(FSA) #for dunnTest
library(agricolae) # for lsd test



# Import datasets ---------------------------------------------------------
#import metadata
metadata <- read_csv("metadata.csv", col_types = cols(X1 = col_skip(), 
                                                      group = col_factor(levels = c("control", 
                                                                                    "NEC", "LOS"))))
colnames(metadata)[1] <- "sample"

##import taxacounts and get relative abundance
otu <- read.csv("otu.csv", check.names = FALSE)
#get the correspondent relationships between otu and genus
goname <- otu[ , 2:9]
#change the name without eg g__, f__
for (i in 1: (ncol(goname) -1)) {
  goname[ ,i] <- substring(goname[ ,i], 4)
}
rm(i)
#OTU massage to relative abundance
proptaxa <- cbind(goname, 
                  prop.table(data.matrix(otu[ , 10:204]), 2))

##import matrix
matrix <- read.delim("matrix.csv", check.names = FALSE, sep = ',')
matrix <- matrix[ , -1]
#remove numeric part of time colum in matrix
matrix$time <- gsub('[[:digit:]]+', '', matrix$time)
##transform "control1-time" into "con1-time" 便于后文stri都减去5个character 
matrix$time1 <- gsub(pattern = "control", replacement = "con", matrix$time)
matrix$time1 <- str_sub(matrix$time1, 5)
#add group levels to time1
matrix$time1 <- ordered(matrix$time1, 
                        levels = c("early post partum", "early pre-onset", "late pre-onset", 
                                   "early disease", "middle disease", "late disease", "post disease"))



# subset Genus abundance with sample information -----------------------------------------------------
genus <- proptaxa[ ,c(6, 9:203)] #bind genus name colum with other columns
genus <- genus %>% group_by(Genus) %>% summarise_all(funs(sum)) #sum by genus names
#change row names to genra name --> get ready for transpose and innerjoin
rownames(genus) <- genus$Genus
#transpose: rows as sample , colunms as genus
genus <- t(genus) %>% as.data.frame(stringsAsFactors = FALSE)
colnames(genus) <- genus[1, ] #define column names
genus <- genus[-1, ]
genus$sample <- rownames(genus)
genus[ ,1:115] <- sapply(genus[ ,1:115], as.numeric) #column 116 is sample name


# merge data frames and set Relative abundance cutoff = 0.05 ----------------------------------------
##merge data frames
meta_matrix <- merge(select(metadata, c(sample, group, patID)), select(matrix, c(sample, time1)), by = "sample")
meta_matrix_genus <- merge(meta_matrix, genus, by = "sample")
## cutoff = 0.05 and get which colunms contains more than 1 sample that contains >0.05 abundance
meta_matrix_genus0.05 <- as.numeric()
for (i in 1:ncol(meta_matrix_genus)) {
  meta_matrix_genus0.05[i] <- which(meta_matrix_genus[[i]] >= 0.05) %>% length()
}
## cutoff 0.05
meta_matrix_genus <- cbind(meta_matrix_genus[, 1:4], 
                           meta_matrix_genus[ , meta_matrix_genus0.05 != 0])
meta_matrix_genus <- meta_matrix_genus[ ,-c(5,6)]

## remove proc files of this section
rm(meta_matrix) 
rm(meta_matrix_genus0.05)


# NEC alluvial plot -------------------------------------------------------
# subset by group name
genus_nec <- meta_matrix_genus[which(meta_matrix_genus$group == "NEC"), ]
genus_nec_time <- genus_nec %>% group_by(time1) %>% summarise_all(funs(mean))
genus_nec_time <- genus_nec_time[,-c(2:4)]
# melt data for plot
mgenus_nec_time <- melt(genus_nec_time)
colnames(mgenus_nec_time) <- c("group", "Genera", "value")
# ggalluvial nec_time
pgenus_nec_time <- ggplot(data = mgenus_nec_time, 
                          aes(x = group, y = value, alluvium = Genera)) + 
  geom_alluvium(aes(fill = Genera, color = Genera), 
                alpha = 0.9, decreasing = NA) +
  theme_minimal() + 
  scale_fill_manual(values = pal_ucscgb("default")(21)) +
  scale_color_manual(values = rep("#F2F3F4", 21)) +
  labs(x = "", y = "Trend of Genera% (NEC)") +
  scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
  theme(legend.position = "none")
# check ggalluvial genus_nec_time
pgenus_nec_time

# LOS alluvial plot -------------------------------------------------------
# subset by group name
genus_los <- meta_matrix_genus[which(meta_matrix_genus$group == "LOS"), ]
genus_los_time <- genus_los %>% group_by(time1) %>% summarise_all(funs(mean))
genus_los_time <- genus_los_time[,-c(2:4)]
# melt data for plot
mgenus_los_time <- melt(genus_los_time)
colnames(mgenus_los_time) <- c("group", "Genera", "value")
# ggalluvial los_time
pgenus_los_time <- ggplot(data = mgenus_los_time, 
                          aes(x = group, y = value, alluvium = Genera)) + 
  geom_alluvium(aes(fill = Genera, color = Genera), 
                alpha = 0.9, decreasing = FALSE) +
  theme_minimal() + 
  scale_fill_manual(values = pal_ucscgb("default")(21)) +
  scale_color_manual(values = rep("#F2F3F4", 21)) +
  labs(x = "", y = "Trend of Genera% (LOS)") +
  scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
  theme(legend.position = "none")
# check ggalluvial genus_los_time
pgenus_los_time

# control alluvial plot -------------------------------------------------------
# subset by group name
genus_control <- meta_matrix_genus[which(meta_matrix_genus$group == "control"), ]
genus_control_time <- genus_control %>% group_by(time1) %>% summarise_all(funs(mean))
genus_control_time <- genus_control_time[,-c(2:4)]
# melt data for plot
mgenus_control_time <- melt(genus_control_time)
colnames(mgenus_control_time) <- c("group", "Genera", "value")
# ggalluvial control_time
pgenus_control_time <- ggplot(data = mgenus_control_time, 
                              aes(x = group, y = value, alluvium = Genera)) + 
  geom_alluvium(aes(fill = Genera, color = Genera), 
                alpha = 0.9, decreasing = NA) +
  theme_minimal() + 
  scale_fill_manual(values = pal_ucscgb("default")(21)) +
  scale_color_manual(values = rep("#F2F3F4", 21)) +
  labs(x = "", y = "Trend of Genera% (control)") +
  scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
  theme(legend.position = "none")
# check ggalluvial genus_control_time
pgenus_control_time

# concatenate three plots -------------------------------------------------
pgroup_time <- ggarrange(pgenus_nec_time, 
                         pgenus_los_time, 
                         pgenus_control_time, 
                         common.legend = TRUE, 
                         nrow = 4,
                         ncol = 1, 
                         legend = "bottom", 
                         labels = c('a', 'b', 'c'))
  # check the concatenate plot
  pgroup_time
 
  # save separately for nec and los 4*10, for control 4*6.4
  
# statistical analysis of genus among group----------------------------------------------------
# Utilize "goname" to get nomenclature of genus based on OTU 
  #matrix 1 OTU significant FROM ZIBR 736, 730, 51, 734, 728, 401, 743, 89, 774 
  #matrix 2 OTU significant FROM ZIBR 1243, 234, 37
  goname[which(goname$OTU == "OTU37"), "Genus"]
# duplicate a data frame for escherichia name
  meta_matrix_genus_dupe <- meta_matrix_genus # duplicate 
  colnames(meta_matrix_genus_dupe)[11] <- "Escherichia"
#define compare strategy
#  compare <- list(c("NEC", "LOS"), c("NEC", "control"), c("LOS", "control"))
# compare_means(comparisons = compare, data = epo)
# early post partum 
  epp <- meta_matrix_genus_dupe[which(meta_matrix_genus$time1 == "early post partum"), ]
  compare_means(Escherichia ~ group,  data = epp, method = "anova")
  epp_genus_stat <- list()
  for (i in 5:25) {
    epp_genus_stat[[i]] <- compare_means((paste(names(epp)[i], " ~ group")) %>% as.formula, 
                                         data = epp, 
                                         method = "anova", 
                                         p.adjust.method = "fdr")  
  }
  epp_genus_stat

  #ttest among arbituray two groups
  t_epp <- setNames(as.data.frame(matrix(ncol = 4, nrow = 25)), c("genus", "NECvsLOS", "NECvscontrol", "LOSvscontrol"))
  for (i in 5:25) {
    t_epp$genus[i] <- colnames(epp)[i]
    t_epp$NECvsLOS[i] <- t.test((paste(names(epp)[i], " ~ group")) %>% as.formula, 
                                data = epp[which(epp$group != "control"), ])$p.value
    t_epp$NECvscontrol[i] <- t.test((paste(names(epp)[i], " ~ group")) %>% as.formula, 
                                    data = epp[which(epp$group != "LOS"), ])$p.value
    t_epp$LOSvscontrol[i] <- t.test((paste(names(epp)[i], " ~ group")) %>% as.formula, 
                                    data = epp[which(epp$group != "NEC"), ])$p.value
  }
  

# early pre-onset
  epo <- meta_matrix_genus_dupe[which(meta_matrix_genus_dupe$time1 == "early pre-onset"), ]
  epo_genus_stat <- list() #wilcoxon
  for (i in 5:25) {
    epo_genus_stat[[i]] <- compare_means((paste(names(epo)[i], " ~ group")) %>% as.formula, 
                                         data = epo, 
                                         method = "anova", 
                                         p.adjust.method = "fdr")  
  }
  epo_genus_stat
  
  #ttest among arbituray two groups
  t_epo <- setNames(as.data.frame(matrix(ncol = 4, nrow = 25)), c("genus", "NECvsLOS", "NECvscontrol", "LOSvscontrol"))
  for (i in 5:25) {
    t_epo$genus[i] <- colnames(epo)[i]
    t_epo$NECvsLOS[i] <- t.test((paste(names(epo)[i], " ~ group")) %>% as.formula, 
                                data = epo[which(epo$group != "control"), ])$p.value
    t_epo$NECvscontrol[i] <- t.test((paste(names(epo)[i], " ~ group")) %>% as.formula, 
                                data = epo[which(epo$group != "LOS"), ])$p.value
    t_epo$LOSvscontrol[i] <- t.test((paste(names(epo)[i], " ~ group")) %>% as.formula, 
                                    data = epo[which(epo$group != "NEC"), ])$p.value
  }
# late pre-onset
  lpo <- meta_matrix_genus_dupe[which(meta_matrix_genus_dupe$time1 == "late pre-onset"), ]
  lpo_genus_stat <- list()
  for (i in 5:25) {
    lpo_genus_stat[[i]] <- compare_means((paste(names(lpo)[i], " ~ group")) %>% as.formula, 
                                         data = lpo, 
                                         method = "anova", 
                                         p.adjust.method = "fdr")  
  }
  
  #ttest among arbituray two groups
  t_lpo <- setNames(as.data.frame(matrix(ncol = 4, nrow = 25)), c("genus", "NECvsLOS", "NECvscontrol", "LOSvscontrol"))
  for (i in 5:25) {
    t_lpo$genus[i] <- colnames(lpo)[i]
    t_lpo$NECvsLOS[i] <- t.test((paste(names(lpo)[i], " ~ group")) %>% as.formula, 
                                data = lpo[which(lpo$group != "control"), ])$p.value
    t_lpo$NECvscontrol[i] <- t.test((paste(names(lpo)[i], " ~ group")) %>% as.formula, 
                                    data = lpo[which(lpo$group != "LOS"), ])$p.value
    t_lpo$LOSvscontrol[i] <- t.test((paste(names(lpo)[i], " ~ group")) %>% as.formula, 
                                    data = lpo[which(lpo$group != "NEC"), ])$p.value
  }
  

# early disease
  ed <- meta_matrix_genus_dupe[which(meta_matrix_genus_dupe$time1 == "early disease"), ]
  ed_genus_stat <- list()
  for (i in 5:25) {
    ed_genus_stat[[i]] <- compare_means((paste(names(ed)[i], " ~ group")) %>% as.formula, 
                                         data = ed, 
                                         method = "anova", 
                                         p.adjust.method = "BH")  
  }
  
  #ttest among arbituray two groups
  t_ed <- setNames(as.data.frame(matrix(ncol = 4, nrow = 25)), c("genus", "NECvsLOS", "NECvscontrol", "LOSvscontrol"))
  for (i in 5:25) {
    t_ed$genus[i] <- colnames(ed)[i]
    t_ed$NECvsLOS[i] <- t.test((paste(names(ed)[i], " ~ group")) %>% as.formula, 
                               data = ed[which(ed$group != "control"), ])$p.value
    t_ed$NECvscontrol[i] <- t.test((paste(names(ed)[i], " ~ group")) %>% as.formula, 
                                   data = ed[which(ed$group != "LOS"), ])$p.value
    t_ed$LOSvscontrol[i] <- t.test((paste(names(ed)[i], " ~ group")) %>% as.formula, 
                                   data = ed[which(ed$group != "NEC"), ])$p.value
  }
  
 
# middle disease
  md <- meta_matrix_genus_dupe[which(meta_matrix_genus_dupe$time1 == "middle disease"), ]
  md_genus_stat <- list()
  for (i in 5:25) {
    md_genus_stat[[i]] <- compare_means((paste(names(md)[i], " ~ group")) %>% as.formula, 
                                        data = md, 
                                        method = "anova", 
                                        p.adjust.method = "BH")  
  }
  # t test among two groups
  t_md <- setNames(as.data.frame(matrix(ncol = 2, nrow = 25)), c("genus", "p.val"))
  for (i in 5:25) {
    t_md$genus[i] <- colnames(md)[i]
    t_md$`p.val`[i] <- t.test((paste(names(md)[i], " ~ group")) %>% as.formula, 
                               data = md)$p.value
  }
  

  # late disease
  ld <- meta_matrix_genus_dupe[which(meta_matrix_genus_dupe$time1 == "late disease"), ]
  ld_genus_stat <- list()
  for (i in 5:25) {
    form <- aov((paste(names(ld)[i], " ~ group")) %>% as.formula, ld)
    ld_genus_stat[[i]] <- LSD.test(form, 'group')
    #ld_genus_stat[[i]] <- t.test((paste(names(ld)[i], " ~ group")) %>% as.formula, 
     #                            data = pd)
  }
  names(ld_genus_stat) <- colnames(ld)
  #t.test 
  t_ld <- setNames(as.data.frame(matrix(ncol = 2, nrow = 25)), c("genus", "p.val"))
  for (i in 5:25) {
    t_ld$genus[i] <- colnames(ld)[i]
    t_ld$`p.val`[i] <- t.test((paste(names(ld)[i], " ~ group")) %>% as.formula, 
                              data = ld)$p.value
  }
  
  
  
  # post disease
  pd <- meta_matrix_genus_dupe[which(meta_matrix_genus_dupe$time1 == "post disease"), ]
  pd_genus_stat <- list()
  for (i in 5:25) {
    #pd_genus_stat[[i]] <- t.test(ld[which(pd$group == "NEC"), names(pd)[i]], 
                                 #ld[which(pd$group == "LOS"), names(pd)[i]])
    pd_genus_stat[[i]] <- t.test((paste(names(pd)[i], " ~ group")) %>% as.formula, 
                                 data = pd) 
  }
  names(pd_genus_stat) <- colnames(pd)  
  
