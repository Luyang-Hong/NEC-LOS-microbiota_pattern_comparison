library(readr)
library(readxl)
library(reshape2)
library(plotly)
library(ggplot2)
library(dplyr)
library(stringr)
library(stringi)
library(ggpubr)

#import alpha diversity
alpha <- read_excel("alpha.xlsx", col_types = c("text", 
                                                "numeric", "numeric", "numeric", "numeric", 
                                                "numeric", "numeric"))
colnames(alpha)[1] <- "sample"
#import beta diversity on otu and genus level
  #genus level beta pc
  genuspcoa <- read_excel("genuspcoa.xlsx")
  colnames(genuspcoa)[1] <- "sample"
  #otu level beta pc
  otupcoa <- read_excel("otupcoa.xlsx")
  colnames(otupcoa)[1] <- "sample"

#impot metadata
metadata <- read_csv("metadata.csv", col_types = cols(X1 = col_skip()))
colnames(metadata)[1] <- "sample"

#import matrixs
  #import matrix: general matrix
  matrix <- read.csv("matrix.csv", check.names = FALSE)
    #remove numeric part of time colum in matrix1
    matrix$time <- gsub('[[:digit:]]+', '', matrix$time)
    ##transform "control1-time" into "con1-time" 便于后文stri都减去5个character 
    matrix$time1 <- gsub(pattern = "control", replacement = "con", matrix$time)
    matrix$time1 <- stri_sub(matrix$time1, 5)
    #add group levels to time1
    matrix$time1 <- ordered(matrix$time1, 
                             levels = c("early post partum", "early pre-onset", "late pre-onset", 
                                        "early disease", "middle disease", "late disease", "post disease"))
  #import matrix1: all groups, from postpartum untill peri-onset
    matrix1 <- read.csv("matrix1.csv", check.names = FALSE)
    #remove numeric part of time colum in matrix1
    matrix1$time <- gsub('[[:digit:]]+', '', matrix1$time)
    ##transform "control1-time" into "con1-time" 便于后文stri都减去5个character 
    matrix1$time1 <- gsub(pattern = "control", replacement = "con", matrix1$time)
    matrix1$time1 <- stri_sub(matrix1$time1, 5)
      #add group levels to time1
      matrix1$time1 <- ordered(matrix1$time1, 
                               levels = c("early post partum", "early pre-onset", "late pre-onset", 
                                          "early disease"))
  #import matrix2: NEC and LOS, from post partum untill post disease
    matrix2 <- read.csv("matrix2.csv", check.names = FALSE)
    #remove numeric part of time colum in matrix1&2
    matrix2$time <- gsub('[[:digit:]]+', '', matrix2$time)
    #同上transform matrix2
    matrix2$time1 <- gsub(pattern = "control", replacement = "con", matrix2$time)
    matrix2$time1 <- stri_sub(matrix2$time1, 5)
    #add group levels to time1
    matrix2$time1 <- ordered(matrix2$time1, 
                             levels = c("early post partum", "early pre-onset", "late pre-onset", 
                                          "early disease", "middle disease", "late disease", "post disease"))
  
##otu based
#merge
meta_otupcoa <- merge(metadata, otupcoa[,1:4])
meta_otupcoa30 <- meta_otupcoa[meta_otupcoa$dol <= 30, ]

#scatter plot
ggplot(meta_otupcoa30,
       aes(x = dol, y = PC1, color = group)) +
  geom_point(shape = 19, size = 2) +
  scale_color_manual(name = NULL,
                     values = c("gray", "dodgerblue", "sienna")) +
  labs(title = "PCoA based on weighted UniFrac distances between bacterial communities plotted over time, OTU level", 
       x = "Age (d)", 
       y = "PCo Axis 1")

#plot with line
ggplot(meta_otupcoa30,
       aes(x = dol, y = PC1, color = PatNo)) +
  geom_point(shape = 19, size = 2) + 
  geom_line() + 
  scale_color_manual(name = NULL, 
                     values = c("gray","gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray",  
                                "dodgerblue", "dodgerblue", "dodgerblue", 
                                "sienna", "sienna", "sienna", "sienna")) +
  labs(title = "PCoA based on weighted UniFrac distances between bacterial communities plotted over time, OTU level, dynamic", 
       x = "Age (d)", 
       y = "PCo Axis 1")
  

#3d scatter using plotly
plot_ly(meta_otupcoa30, x = ~dol, y = ~PC1, z = ~PC2, color = ~group, 
        type = 'scatter3d', mode = 'lines', colors = c("gray", "dodgerblue", "sienna")) %>% 
  add_markers() %>% 
  layout(scene =  list(xaxis = list(title = "Age(d)"), 
                       yaxis = list(title = "PCo Axis 1"), 
                       zaxis = list(title = "PCo Axis 2")))
  
##genus based
#merge
meta_genuspcoa <- merge(metadata, genuspcoa[,1:4])
meta_genuspcoa30 <- meta_otupcoa[meta_otupcoa$dol <= 30, ]

#plot
ggplot(meta_genuspcoa30,
       aes(x = dol, y = PC1, color = group)) +
  geom_point(shape = 19, size = 2) +
  scale_color_manual(name = NULL,
                     values = c("gray", "dodgerblue", "sienna")) +
  labs(title = "PCoA based on weighted UniFrac distances between bacterial communities plotted over time, genus level", 
       x = "Age (d)", 
       y = "PCo Axis 1")

#plot with line
ggplot(meta_genuspcoa,
       aes(x = dol, y = PC1, color = PatNo)) +
  geom_point(shape = 19, size = 2) + 
  geom_line() +
  scale_color_manual(name = NULL, 
                     values = c("gray","gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray",  
                                "dodgerblue", "dodgerblue", "dodgerblue", 
                                "sienna", "sienna", "sienna", "sienna")) +
  labs(title = "PCoA based on weighted UniFrac distances between bacterial communities plotted over time, genus level, dynamic", 
       x = "Age (d)", 
       y = "PCo Axis 1")


#3d scatter using plotly
plot_ly(meta_genuspcoa30, x = ~dol, y = ~PC1, z = ~PC2, color = ~group, 
        type = 'scatter3d', mode = 'lines', colors = c("gray", "dodgerblue", "sienna")) %>% 
  add_markers() %>% 
  layout(scene =  list(xaxis = list(title = "Age(d)"), 
                       yaxis = list(title = "PCo Axis 1"), 
                       zaxis = list(title = "PCo Axis 2")))


####lpha diversity
meta_alpha <- merge(metadata, alpha)

#<= 30d alpha and plot over time
meta_alpha30 <- meta_alpha[meta_alpha$dol <= 30, ]
ggplot(data = meta_alpha, aes(x = dol, y = shannon, color = group)) +
  geom_point(shape = 19, size = 3) +
  geom_line() +
  scale_color_manual(name = 'test',
                     values = c("gray", "dodgerblue", "sienna")) +
  labs(title = "alpha diversity among three groups over time", 
       x= "Age (d)", 
       y = "Shannon index")

#alpha diversity between gender
alphagender <- ggplot(meta_alpha, aes(x = gender, y = shannon, color = gender)) + 
  geom_boxplot() +
  labs(title = "alpha diversity between genders", 
       x = "Gender", 
       y = "Shannon index")
alphagender + geom_jitter(shape = 16, position = position_jitter(0.2))

#alpha between atb usage
alphaatb <- ggplot(meta_alpha, aes(x = atb, y = shannon, color = atb)) +
  geom_boxplot() + 
  labs(x = "Advanced ATB usage")
alphaatb + geom_jitter(shape = 16, position = position_jitter(0.2))

#alpha between group
alphagroup <- ggplot(meta_alpha, aes(x = group, y = shannon, color = group)) +
  geom_boxplot() +
  labs(title = "alpha diversity between groups", 
       x = "Groups", 
       y = "Shannon index")
alphagroup + geom_jitter(shape = 16, position = position_jitter(0.2))

#alpha diversity over time
ggplot(meta_alpha30,
       aes(x = dol, y = shannon, color = PatNo)) +
  geom_point(shape = 19, size = 2) + 
  geom_line() +
  scale_color_manual(name = NULL, 
                     values = c("gray","gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray", "gray", 
                                "gray", "gray",  
                                "dodgerblue", "dodgerblue", "dodgerblue", 
                                "sienna", "sienna", "sienna", "sienna")) +
  labs(title = "alpha diversity between bacterial communities plotted over time", 
       x = "Age (d)", 
       y = "Shannon Index")

##metadat & alphadiversity NEC group
meta_alpha_nec <- meta_alpha[meta_alpha$group == "NEC", ]
ggplot(meta_alpha_nec, 
       aes(x= dol, y = shannon, color = PatNo)) +
  geom_point(shape = 19, size = 2) +
  geom_line() +
  scale_color_manual(name = NULL, 
                     values = c("gray", "dodgerblue","sienna", "red")) +
  labs(title = "alpha diversity", 
       x = "Age", 
       y = "Shannon Index")

##original matrix and alpha diversity
matrix_alpha <- inner_join(select(matrix, c("sample", "time", "time1")), alpha, by = "sample")
meta_matrix_alpha <- inner_join(metadata, matrix_alpha, by = "sample")
  #group by group-time interval and calculate mean shannon
  meta_matrix_alpha %>% group_by(time) %>% 
    summarise_at(vars(-shannon), funs(mean(., na.rm=TRUE)))
  
##matrix1 and alphadiversity
matrix1_alpha <- inner_join(select(matrix1, c("sample", "time", "time1")), alpha, by = "sample")
meta_matrix1_alpha <- inner_join(metadata, matrix1_alpha, by = "sample")
  #specify the comparisons
  shannon_comparisons <- list(c("NEC", "LOS"), c("NEC", "control"), c("control", "LOS"))
  #ggpubr x = time , color = group, y = shannon
  ggboxplot(data = meta_matrix_alpha, x = "time1", y = "shannon", 
            color = "group", palette = c("gray", "dodgerblue","sienna"), 
            add = "jitter",  
            title = "Alpha diversity over post partum time interval",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(aes(group = group), label = "p.format")  #add pairwise comparisons p-value
  
  ggboxplot(data = meta_matrix1_alpha, x = "group", y = "shannon", 
            color = "group", palette = c("gray", "sienna", "dodgerblue"), 
            facet.by = "time1", 
            add = "jitter",  
            title = "Post-partum untill peri-onset comparison of alpha diversity",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(comparisons = shannon_comparisons) + #add pairwise comparisons p-value
    facet_grid(.~time1)


  #specify the comparisons
  shannon_comparisons <- list(c("NEC", "LOS"), c("NEC", "control"), c("control", "LOS"))
  #ggpubr x = time , color = group, y = shannon
  ggboxplot(data = meta_matrix1_alpha, x = "time1", y = "shannon", 
            color = "group", palette = c("gray", "dodgerblue","sienna"), 
            add = "jitter",  
            title = "Post-partum untill peri-onset comparison of alpha diversity",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(aes(group = group), label = "p.format")  #add pairwise comparisons p-value
  
  ggboxplot(data = meta_matrix1_alpha, x = "group", y = "shannon", 
            color = "group", palette = c("gray", "sienna", "dodgerblue"), 
            facet.by = "time1", 
            add = "jitter",  
            title = "Post-partum untill peri-onset comparison of alpha diversity",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(comparisons = shannon_comparisons) + #add pairwise comparisons p-value
    facet_grid(.~time1)

##matrix2 and alphadiversity
matrix2_alpha <- inner_join(select(matrix2, c("sample", "time", "time1")), alpha, by = "sample")
meta_matrix2_alpha <- inner_join(metadata, matrix2_alpha, by = "sample")
  #specify the comparisons
  shannon_comparisons <- list(c("NEC", "LOS"), c("NEC", "control"), c("control", "LOS"))
  #ggpubr x = time , color = group, y = shannon
  ggboxplot(data = meta_matrix2_alpha, x = "time1", y = "shannon", 
            color = "group", palette = c("dodgerblue","sienna"), 
            add = "jitter",  
            title = "Post-partum untill post-disease comparison of alpha diversity",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(aes(group = group), label = "p.format")  #add pairwise comparisons p-value
  
  ggboxplot(data = meta_matrix2_alpha, x = "group", y = "shannon", 
            color = "group", palette = c("sienna", "dodgerblue"), 
            facet.by = "time1", 
            add = "jitter",  
            title = "Post-partum untill peri-onset comparison of alpha diversity",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(aes(group = group), label = "p.format") + #add pairwise comparisons p-value
    facet_grid(.~time1)


