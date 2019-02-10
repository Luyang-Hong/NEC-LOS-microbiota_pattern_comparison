library(readr) #READ CSV
library(readxl)
library(reshape2)
library(plotly)
library(ggplot2)
library(dplyr)
library(stringr)
library(stringi)
library(ggpubr)
library(cowplot)


##################################################################################
######import metadata, alpha diversity and matrix with time information###########
#impot metadata
metadata <- read_csv("metadata.csv", col_types = cols(X1 = col_skip(), 
                                                      group = col_factor(levels = c("control", 
                                                                                    "NEC", "LOS"))))
colnames(metadata)[1] <- "sample"

#import alpha diversity
alpha <- read_csv("alpha.csv", col_types = cols(ace = col_number(), 
                                                chao = col_number(), coverage = col_number(), 
                                                shannon = col_number(), simpson = col_number(), 
                                                sobs = col_number()))
colnames(alpha)[1] <- "sample"
rownames(alpha) <- alpha$sample
  #calculate B0620 B3818 B4025 samples' alphadiversity
  alpha <- rbind(alpha, 
                 #B0620 = B0619 + (B0619-B0622)/3
                 c("B0620", alpha["B0619", 2:7] %>% as.numeric() + (alpha["B0622", 2:7] %>% as.numeric() - alpha["B0619", 2:7] %>% as.numeric())/3), 
                 #B3818 = (B3815 + B3821)/2
                 c("B3818", (alpha["B3815", 2:7] %>% as.numeric() + alpha["B3821", 2:7] %>% as.numeric())/2),
                 #B4025 = B4021 + 4*(B4028-B4021)/7
                 c("B4025", alpha["B4021", 2:7] %>% as.numeric()+ 4*(alpha["B4028", 2:7] %>% as.numeric() - alpha["B4021", 2:7] %>% as.numeric())/7))

  #define rownames of alpha
  rownames(alpha) <- alpha$sample
  #change shannon simpson... as numeric
  alpha[, 2:7] <- sapply(alpha[, 2:7], as.numeric)
  alpha <- as.data.frame(alpha)

#import matrixs
  #import matrix: general matrix
  matrix <- read.delim("matrix.csv", check.names = FALSE, sep = ',')
  matrix <- matrix[ , -1]
    #remove numeric part of time colum in matrix
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


##################################################################################
########################merge matrix, metadata and diversity######################
matrix_alpha <- inner_join(select(matrix, c("sample", "time", "time1")), alpha, by = "sample")
meta_matrix_alpha <- inner_join(metadata, matrix_alpha, by = "sample")

  ##########inter time interval, groups comparisons
  #specify the comparisons
  shannon_comparisons <- list(c("NEC", "LOS"), c("NEC", "control"), c("control", "LOS"))
  #shannon plot: within a time interval, compare shannon diversity of three groups
  shannon_time <- ggboxplot(data = meta_matrix_alpha, x = "group", y = "shannon", 
                          color = "group", 
                          #palette = c("gray", "sienna", "dodgerblue"), 
                          palette = 'lancet', 
                          facet.by = "time1", 
                          add = "jitter",  
                          ylim = c(0,5.5), 
                          xlim = c(1,3),
                          #title = "Alpha diversity over post partum time interval",
                          ylab = "Shannon Diversity Index") +
    stat_compare_means(aes(group =group), label.y = 5.2, label.x = 1.1) +
    stat_compare_means(comparisons = shannon_comparisons, method = "wilcox.test") + #add pairwise comparisons p-value
    facet_wrap(.~time1, ncol = 3) + 
    rremove("xlab")
  shannon_time <- ggpar(shannon_time, legend = "right", legend.title = "groups")
  shannon_time
  #observed plot: within a time interval, compare observed richness of three groups
  sobs_time <- ggboxplot(data = meta_matrix_alpha, x = "group", y = "sobs", 
                         color = "group", 
                         #palette = c("gray", "sienna", "dodgerblue"), 
                         palette = 'lancet', 
                         facet.by = "time1", 
                         add = "jitter",  
                         ylim = c(0,370), 
                         xlim = c(1,3),
                         #title = "Alpha diversity over post partum time interval",
                         ylab = "Sobs") +
    stat_compare_means(aes(group =group), label.y = 350, label.x = 1.1) +
    stat_compare_means(comparisons = shannon_comparisons, method = "wilcox.test") + #add pairwise comparisons p-value
    facet_wrap(.~time1, ncol = 3) + 
    rremove("xlab")
  sobs_time <- ggpar(sobs_time, legend = "right", legend.title = "groups")
  sobs_time
  #ggarrange shannon and sobs index comparison between groups over time
  shannon_sobs_time <- ggarrange(shannon_time, sobs_time, ncol = 2, 
                                 common.legend = TRUE, legend = "bottom", 
                                 labels = c("a", "b"))
  shannon_sobs_time

  ##########inter group, time comparisons
  #specify the comparisons
  shannon_time_comparisons <- list(c("early post partum", "early pre-onset"), 
                                  c("early pre-onset", "late pre-onset"), 
                                  c("late pre-onset", "early disease"), 
                                  c("early disease", "middle disease"),
                                  c("middle disease", "late disease"), 
                                  c("late disease", "post disease"))
  #NEC group shannon diversity changes
  meta_matrix_shannon_nec <- meta_matrix_alpha[meta_matrix_alpha$group == "NEC", ]
  #plot nec group shannon diversity over time
  shannonnec <- ggboxplot(data = meta_matrix_shannon_nec, x = "time1", y = "shannon", 
                        color = "time1", 
                        #palette = c("gray", "sienna", "dodgerblue"), 
                        add = "jitter",  
                        ylim = c(0,4.2), 
                        width = 0.7, 
                        #xlim = c(1,3),
                        #title = "Alpha diversity over post partum time interval",
                        xlab = "Time Interval", ylab = "Shannon Index (NEC)" 
  ) +
    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
    stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
    stat_compare_means(comparisons = shannon_time_comparisons, method = "wilcox.test", 
                       label.y = c(3.6, 3.4, 1.5, 1.2, 1.8, 2.0)) + #adjusting the p values of this comparison at the same )) + #add pairwise comparisons p-value
    rremove("x.ticks") + 
    rremove("legend.title") + rremove("legend")
  shannonnec
  #los group shannon diversity changes
  meta_matrix_shannon_los <- meta_matrix_alpha[meta_matrix_alpha$group == "LOS", ]
  #plot los group shannon diversity over time
  shannonlos <- ggboxplot(data = meta_matrix_shannon_los, x = "time1", y = "shannon", 
                        color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                        add = "jitter",  
                        ylim = c(0,4.2), 
                        palette = "startrek",
                        #width = 0.7, 
                        #xlim = c(1,3),
                        #title = "Post-partum alpha diversity over time",
                        xlab = "Time Interval", ylab = "Shannon Index (LOS)") +
    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
    stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
    stat_compare_means(comparisons = shannon_time_comparisons, method = "wilcox.test", 
                       label.y = c(3.6, 3.4, 1.4, 1.8, 1.56, 2.0)) + #adjusting the p values of this comparison at the same )) + #add pairwise comparisons p-value
    rremove("x.ticks") + 
    rremove("legend.title") + rremove("legend")
  shannonlos
  #control group shannon diversity changes
  meta_matrix_shannon_control <- meta_matrix_alpha[meta_matrix_alpha$group == "control", ]
  #plot control group alpha diversity over time
  shannoncontrol <- ggboxplot(data = meta_matrix_shannon_control, x = "time1", y = "shannon", 
                            color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                            add = "jitter", 
                            palette = "startrek",
                            ylim = c(0,4.2), 
                            width = 0.7,
                            #xlim = c(1,3),
                            #title = "Alpha diversity over post partum time interval",
                            xlab = "Time Interval", ylab = "Shannon Index (control)") +
    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
    stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
    stat_compare_means(comparisons = shannon_time_comparisons[1:3], method = "wilcox.test", 
                       label.y = c(3.6, 3.4, 2.0)) + #adjusting the p values of this comparison at the same )) + #add pairwise comparisons p-value
    rremove("x.ticks") + 
    rremove("legend.title") + rremove("legend")
  shannoncontrol
  #alpha_groups legend for all
  #plot a plot merely for legend
  shannonleg <- ggboxplot(data = meta_matrix_shannon_los, x = "time1", y = "shannon", 
                        color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                        add = "jitter",  
                        palette = "startrek",
                        ylim = c(0,4.2), 
                        width = 0.7, 
                        #xlim = c(1,3),
                        #title = "Post-partum alpha diversity over time",
                        xlab = "Time Interval", ylab = "Shannon Index", 
                        legend.title = "Time Interval") +
    stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
    rremove("x.text") + rremove("x.ticks") 
  shannonleg <- ggpar(shannonleg, legend = "right")
  #get legend of the "merely for legend" plot :)
  shannon_groups_leg <- get_legend(shannonleg) %>% as_ggplot()
  shannon_groups_leg
  #ggarrange alpha diversity change over time plots for all groups
  shannon_groups <- ggdraw() +
    draw_plot(shannonnec, x = 0, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot(shannonlos, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot(shannon_groups_leg, x = 0.32, y = 0, width = 0.3, height = 0.5) +
    draw_plot(shannoncontrol, x = 0, y = 0, width = 0.32, height = 0.5) +
    draw_plot_label(label = c("a", "b", "c"), 
                    size = 15, 
                    x = c(0, 0.5, 0), 
                    y = c(1, 1, 0.5))
  shannon_groups
  
  ##perform tests
    compare_means(data = meta_matrix_shannon_nec, shannon ~ time1, method = "wilcox.test")
  
  ##########inter group, time comparisons for sobs
  
  #NEC group sobs diversity changes
  meta_matrix_sobs_nec <- meta_matrix_alpha[meta_matrix_alpha$group == "NEC", ]
  #plot nec group sobs diversity over time
  sobsnec <- ggboxplot(data = meta_matrix_sobs_nec, x = "time1", y = "sobs", 
                       color = "time1", 
                       palette = "startrek",
                       add = "jitter",  
                       ylim = c(0,238), 
                       width = 0.7, 
                       #xlim = c(1,3),
                       #title = "Alpha diversity over post partum time interval",
                       xlab = "Time Interval", ylab = "sobs (NEC)" 
  ) +
    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
    stat_compare_means(aes(group = time1), label.y = 230, label.x = 1.5) + 
    rremove("x.ticks") + 
    rremove("legend.title") + rremove("legend")
  sobsnec
  #los group sobs diversity changes
  meta_matrix_sobs_los <- meta_matrix_alpha[meta_matrix_alpha$group == "LOS", ]
  #plot los group sobs diversity over time
  sobslos <- ggboxplot(data = meta_matrix_sobs_los, x = "time1", y = "sobs", 
                       color = "time1", 
                       palette = "startrek",
                       add = "jitter",  
                       ylim = c(0, 238), 
                       width = 0.7, 
                       #xlim = c(1,3),
                       #title = "Post-partum alpha diversity over time",
                       xlab = "Time Interval", ylab = "sobs (LOS)") +
    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
    stat_compare_means(aes(group = time1), label.y = 230, label.x = 1.5) + 
    rremove("x.ticks") + 
    rremove("legend.title") + rremove("legend")
  sobslos
  #control group sobs diversity changes
  meta_matrix_sobs_control <- meta_matrix_alpha[meta_matrix_alpha$group == "control", ]
  #plot control group alpha diversity over time
  sobscontrol <- ggboxplot(data = meta_matrix_sobs_control, x = "time1", y = "sobs", 
                      color = "time1", 
                      palette = "startrek", 
                      add = "jitter",  
                      ylim = c(0,238), 
                      width = 0.7,
                      #xlim = c(1,3),
                      #title = "Alpha diversity over post partum time interval",
                      xlab = "Time Interval", ylab = "sobs (control)") +
    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
    stat_compare_means(aes(group = time1), label.y = 230, label.x = 1.5) + 
    rremove("x.ticks") + 
    rremove("legend.title") + rremove("legend")
  sobscontrol
  #alpha_groups legend for all
  #plot a plot merely for legend
  sobsleg <- ggboxplot(data = meta_matrix_sobs_los, x = "time1", y = "sobs", 
                       color = "time1", 
                       palette = "startrek",
                       add = "jitter",  
                       ylim = c(0,4.2), 
                       width = 0.7, 
                       #xlim = c(1,3),
                       #title = "Post-partum alpha diversity over time",
                       xlab = "Time Interval", ylab = "sobs Index", 
                       legend.title = "Time Interval") +
    stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
    rremove("x.text") + rremove("x.ticks") 
  sobsleg <- ggpar(sobsleg, legend = "right")
  #get legend of the "merely for legend" plot :)
  sobs_groups_leg <- get_legend(sobsleg) %>% as_ggplot()
  sobs_groups_leg
  #ggarrange alpha diversity change over time plots for all groups
  sobs_groups <- ggdraw() +
    draw_plot(sobsnec, x = 0, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot(sobslos, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
    draw_plot(sobs_groups_leg, x = 0.32, y = 0, width = 0.3, height = 0.5) +
    draw_plot(sobscontrol, x = 0, y = 0, width = 0.32, height = 0.5) +
    draw_plot_label(label = c("a", "b", "c"), 
                    size = 15, 
                    x = c(0, 0.5, 0), 
                    y = c(1, 1, 0.5))
  sobs_groups
  
  
  


  

  
  

  
  
  
  
####alpha diversity
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
       y = "Shannon Index") + 

 
  #inter time interval, groups comparisons
  #specify the comparisons
  shannon_comparisons <- list(c("NEC", "LOS"), c("NEC", "control"), c("control", "LOS"))
  #ggpubr x = time , color = group, y = shannon
  ggboxplot(data = meta_matrix_alpha, x = "time1", y = "shannon", 
            color = "group", palette = c("gray", "dodgerblue","sienna"), 
            add = "jitter", 
            width = 0.5, 
            title = "Alpha diversity over post partum time interval",
            xlab = "Time Interval", ylab = "Shannon Index") +
    #stat_compare_means(aes(group = group), label = "p.format")
    stat_compare_means(aes(group = group), label = "p.format")  #add pairwise comparisons p-value
  
    alpha_time <- ggboxplot(data = meta_matrix_alpha, x = "group", y = "shannon", 
                            color = "group", palette = c("gray", "sienna", "dodgerblue"), 
                            facet.by = "time1", 
                            add = "jitter",  
                            ylim = c(0,5.5), 
                            xlim = c(1,3),
                            #title = "Alpha diversity over post partum time interval",
                            xlab = "Time Interval", ylab = "Shannon Index") +
                    #stat_compare_means(aes(group = group), label = "p.format")
                    stat_compare_means(aes(group =group), label.y = 5.2, label.x = 1.1) +
                    stat_compare_means(comparisons = c("NEC", "LOS"), method = "wilcox.test") +
                    stat_compare_means(comparisons = shannon_comparisons, method = "wilcox.test") + #add pairwise comparisons p-value
                    facet_wrap(.~time1, ncol = 4)  
    alpha_time <- ggpar(alpha_time, legend = "right", legend.title = "groups")
    alpha_time
    
      #NEC group alpha diversity changes
      meta_matrix_alpha_nec <- meta_matrix_alpha[meta_matrix_alpha$group == "NEC", ]
      #plot nec group alpha diversity over time
      alphanec <- ggboxplot(data = meta_matrix_alpha_nec, x = "time1", y = "shannon", 
                          color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                          add = "jitter",  
                          ylim = c(0,4.2), 
                          width = 0.7, 
                          #xlim = c(1,3),
                          #title = "Alpha diversity over post partum time interval",
                          xlab = "Time Interval", ylab = "Shannon Index (NEC)" 
                          ) +
                    scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
                    stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
                    rremove("x.ticks") + 
                    rremove("legend.title") + rremove("legend")
      alphanec
      #los group alpha diversity changes
      meta_matrix_alpha_los <- meta_matrix_alpha[meta_matrix_alpha$group == "LOS", ]
      #plot los group alpha diversity over time
      alphalos <- ggboxplot(data = meta_matrix_alpha_los, x = "time1", y = "shannon", 
                            color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                            add = "jitter",  
                            ylim = c(0,4.2), 
                            width = 0.7, 
                            #xlim = c(1,3),
                            #title = "Post-partum alpha diversity over time",
                            xlab = "Time Interval", ylab = "Shannon Index (LOS)") +
        scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
        stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
        rremove("x.ticks") + 
        rremove("legend.title") + rremove("legend")
      alphalos
      #control group alpha diversity changes
      meta_matrix_alpha_control <- meta_matrix_alpha[meta_matrix_alpha$group == "control", ]
      #plot control group alpha diversity over time
      alphacontrol <- ggboxplot(data = meta_matrix_alpha_control, x = "time1", y = "shannon", 
                                color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                                add = "jitter",  
                                ylim = c(0,4.2), 
                                width = 0.7,
                                #xlim = c(1,3),
                                #title = "Alpha diversity over post partum time interval",
                                xlab = "Time Interval", ylab = "Shannon Index (control)") +
        scale_x_discrete(labels=c("EPP", "EPO", "LPO", "ED", "MD", "LD", "PD")) +
        stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
        rremove("x.ticks") + 
        rremove("legend.title") + rremove("legend")
      alphacontrol
      #alpha_groups legend for all
        #plot a plot merely for legend
        alphaleg <- ggboxplot(data = meta_matrix_alpha_los, x = "time1", y = "shannon", 
                              color = "time1", #palette = c("gray", "sienna", "dodgerblue"), 
                              add = "jitter",  
                              ylim = c(0,4.2), 
                              width = 0.7, 
                              #xlim = c(1,3),
                              #title = "Post-partum alpha diversity over time",
                              xlab = "Time Interval", ylab = "Shannon Index", 
                              legend.title = "Time Interval") +
          stat_compare_means(aes(group = time1), label.y = 4, label.x = 1.5) + 
          rremove("x.text") + rremove("x.ticks") 
        alphaleg <- ggpar(alphaleg, legend = "right")
          #get legend of the "merely for legend" plot :)
          alpha_groups_leg <- get_legend(alphaleg) %>% as_ggplot()
          alpha_groups_leg
      #ggarrange alpha diversity change over time plots for all groups
      alpha_groups <- ggdraw() +
        draw_plot(alphanec, x = 0, y = 0.5, width = 0.5, height = 0.5) +
        draw_plot(alphalos, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
        draw_plot(alphacontrol, x = 0, y = 0, width = 0.3, height = 0.5) +
        draw_plot(alpha_groups_leg, x = 0.3, y = 0, width = 0.3, height = 0.5) +
        draw_plot_label(label = c("a", "b", "c"), 
                        size = 15, 
                        x = c(0, 0.5, 0), 
                        y = c(1, 1, 0.5))
        
      alpha_groups
      
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
    facet_wrap(.~time1)

