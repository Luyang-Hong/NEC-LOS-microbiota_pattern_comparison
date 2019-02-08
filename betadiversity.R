library(readr) #READ_CSV
#library(readxl)
#library(reshape2)
#library(plotly)
#library(ggplot2)
#library(dplyr)
#library(stringr)
#library(stringi)
#library(ggpubr)
#library(cowplot)
library(phyloseq)



#####import metadata and beta diversity######
#impot metadata
metadata <- read_csv("metadata.csv", col_types = cols(X1 = col_skip(), 
                                                      group = col_factor(levels = c("control", 
                                                                                    "NEC", "LOS"))))
colnames(metadata)[1] <- "sample"

#import beta diversity
