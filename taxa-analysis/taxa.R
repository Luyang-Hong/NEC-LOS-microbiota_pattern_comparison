#set working directory
setwd("taxa-analysis")

#load needed libraries
library()

##################################################################################
######import metadata, otu with other and matrix with time information###########
  #impot metadata
  metadata <- read_csv("metadata.csv", col_types = cols(X1 = col_skip(), 
                                                        group = col_factor(levels = c("control", 
                                                                                      "NEC", "LOS"))))
  colnames(metadata)[1] <- "sample"
  
  #import taxa
  ##otu
  otu <- read.csv("otu.csv", check.names = FALSE)
    #get the correspondent relationships between otu and genus
    goname <- otu[ , c(7,9)]
  
    #OTU massage
    #tranposeOTU and rename colnames by otu name
    genus <- t(otu)
    colnames(genus) <- genus[7, ] #colnames = genus level = 7
    genus <- genus[-c(1:9), ] #remove domain~species taxa information, which occupied first 9 rows of the matrix
    mode(genus) <- "numeric" #change the mode of genus
    genus <- data.frame(genus)
    genus <- genus[-196,]
    genus$sample <- rownames(genus)
    
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
    
    
  

