#####important notes!!!!!!
# relative abundance for genus cutoff is 0.1! 
# cutoff abundance for this analysis session is 0.5
#set working directory
setwd()

#load needed libraries
library(reshape2) #for melt
library(ggalluvial) #for plot genus changes over time
library(tidyverse) #for readr: read_csv
library(stringr) # for function str_sub
library(ggsci) #for palette color pick
library(ggpubr) #for ggarrange common.legends



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


# merge data frames and set Relative abundance cutoff = 0.1 ----------------------------------------
##merge data frames
meta_matrix <- merge(select(metadata, c(sample, group, patID)), select(matrix, c(sample, time1)), by = "sample")
meta_matrix_genus <- merge(meta_matrix, genus, by = "sample")
## cutoff = 0.1 and get which colunms contains more than 1 sample that contains >0.05 abundance
meta_matrix_genus0.1 <- as.numeric()
for (i in 1:ncol(meta_matrix_genus)) {
  meta_matrix_genus0.1[i] <- which(meta_matrix_genus[[i]] >= 0.1) %>% length()
}
## cutoff 0.1
meta_matrix_genus <- cbind(meta_matrix_genus[, 1:4], 
                           meta_matrix_genus[ , meta_matrix_genus0.1 != 0])
meta_matrix_genus <- meta_matrix_genus[ ,-c(5,6)]

## remove proc files of this section
rm(meta_matrix) 
rm(meta_matrix_genus0.1)


# postpartum alluvial plot -------------------------------------------------------
# subset by group name
genus_epp <- meta_matrix_genus[which(meta_matrix_genus$time1 == "early post partum"), ]
genus_epp_group <- genus_epp %>% group_by(group) %>% summarise_all(funs(mean))
genus_epp_group <- genus_epp_group[,-c(2:4)]
# melt data for plot
mgenus_epp_group <- melt(genus_epp_group)
colnames(mgenus_epp_group) <- c("group", "Genera", "value")
# ggalluvial epp_groups
pgenus_epp_group <- ggplot (data = mgenus_epp_group, 
                            aes(axis1 = Genera, axis2 = group, 
                                y = value)) + 
  scale_x_discrete(limits = c("Genera", "group"), expand = c(.2, .025)) +
  geom_alluvium(aes(fill = Genera)) +
  scale_fill_manual(values = pal_ucscgb("default")(21)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  labs(x = "", y = "Genera abundance in each group")
# check alluvial epp_groups 
pgenus_epp_group

# stacked barplot for epp_groups
pgenus_epp_group <- ggplot(mgenus_epp_group, aes(x = group, y = value, fill = Genera)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = pal_ucscgb("default")(21))+ 
  theme_minimal()
# check barplot epp_groups 
pgenus_epp_group

#another ggalluvial
ggplot(mgenus_epp_group, 
       aes(x = group, y = value, fill = Genera)) + 
  geom_bar(stat = "identity", width = 0.5, col = 'black') + 
  theme_minimal() + 
  geom_segment(data = )


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
                alpha = 0.9, decreasing = NA) +
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
                         legend = "bottom", 
                         labels = c('a', 'b', 'c'))
  # check the concatenate plot
  pgroup_time

