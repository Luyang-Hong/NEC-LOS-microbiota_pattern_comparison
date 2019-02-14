#set working directory
setwd("taxa-analysis")

#load needed libraries
library(reshape2) #for melt
library(ggalluvial) #for plot genus changes over time

##################################################################################
######import metadata, otu with other and matrix with time information###########
  #impot metadata
  metadata <- read_csv("metadata.csv", col_types = cols(X1 = col_skip(), 
                                                        group = col_factor(levels = c("control", 
                                                                                      "NEC", "LOS"))))
  colnames(metadata)[1] <- "sample"
  
  ##import taxacounts
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
    
  #import matrix
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

##################################################################################
######genus analysis###########
    ##massage genus
    genus <- proptaxa[ ,c(6, 9:203)] #bind genus name colum with other columns
    genus <- genus %>% group_by(Genus) %>% summarise_all(funs(sum)) #sum by genus names
    #change row names to genra name --> get ready for transpose
    rownames(genus) <- genus$Genus
    #transpose: rows as sample , colunms as genus
    genus <- t(genus) %>% as.data.frame(stringsAsFactors = FALSE)
    colnames(genus) <- genus[1, ] #define column names
    genus <- genus[-1, ]
    genus$sample <- rownames(genus)
    genus[ ,1:115] <- sapply(genus[ ,1:115], as.numeric)
    
    ##merge data frames
    meta_matrix <- merge(select(metadata, c(sample, group, patID)), select(matrix, c(sample, time1)), by = "sample")
    meta_matrix_genus <- merge(meta_matrix, genus, by = "sample")

      ####intergroup, comparison of time interval
      #### NEC group
      genus_nec <- meta_matrix_genus[which(meta_matrix_genus$group == "NEC"), ]
      # time interval relative abundance
      genus_nec_time <- aggregate(.~time1, genus_nec[,4:119], #this data should include the index"time1"
                                  FUN = mean)
      # select >0.05 genus in time interval 
      genus_nec_time0.05 <- as.numeric()
      for (i in 1:ncol(genus_nec_time)) {
        genus_nec_time0.05[i] <- length(which(genus_nec_time[[i]] >= 0.05))
      }
      # subset genus_nec_time by the 0.05 relative abundance level
      genus_nec_time <- cbind(genus_nec_time$time1, 
                              genus_nec_time[ , which(genus_nec_time0.05 != 0)])
      colnames(genus_nec_time)[1] <- "Time Interval"
      # transpose genus_nec_time for ggalluvial plot
      genus_nec_time <- t(genus_nec_time) %>% as.data.frame(stringsAsFactors = FALSE)
      colnames(genus_nec_time) <- genus_nec_time[1, ]
      genus_nec_time <- genus_nec_time[-1, ]
      genus_nec_time$`Genus` <- rownames(genus_nec_time)
      melt_genus_nec_time <- melt(genus_nec_time, id.vars = "Genus")
      # plot relative abundance during different time intervals(as x axis)
      ggplot()
      