library(reshape2) #for melt
library(dplyr) #for inner.join
library(stringr) #for finding patterns

sample <- read.csv("sample.csv", check.names = FALSE)
timepoint <- read.csv("paste-time2.csv", check.names = FALSE, comment.char = "#")
otu <- read.csv("otu.csv", check.names = FALSE)

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

#process sample, aim: if the samplename has a real sample
#remove patno of sample
subsample <- sample[,-1]
msubsample <- melt(subsample, id = "patID")
#add sample id
msubsample$sample <- paste(msubsample$patID, msubsample$variable, sep = "")
colnames(msubsample)[3] <- "ifsample"

##process sample, aim: time point和sample id对应
subtimepoint <- timepoint[ , -1]
msubtimepoint <- melt(subtimepoint, id = "patID")
msubtimepoint$sample <- paste(msubtimepoint$patID, msubtimepoint$variable, sep = "")
colnames(msubtimepoint)[3] <- "time"

##merge two matrixs
sampletime <- merge(msubsample, msubtimepoint, by = "sample")
sampletime <- select(sampletime, c("sample", "time", "ifsample"))
sampletime <- sampletime[!is.na(sampletime$ifsample), ]

##otu
#tranposeOTU and rename colnames by otu name
otu <- t(otu)
colnames(otu) <- otu[9, ] #colnames = otunames
otu <- otu[-c(1:9), ] #remove domain~species taxa
mode(otu) <- "numeric"
otu <- data.frame(otu)
otu <- otu[-196,]
otu$sample <- rownames(otu)

##innerjoin out and sampletime matrix
matrix <- inner_join(sampletime, otu, by ="sample")
matrix <- matrix[ , -3]

##delete control-TRUE samples
matrix <- matrix[str_extract(matrix$time, "-TRUE") %>% is.na(), ]


##方案1 舍去middle disease, late disease, post disease
matrix1 <- matrix[str_extract(matrix$time, c("-middle disease")) %>% is.na(), ]
matrix1 <- matrix1[str_extract(matrix1$time, c("-late disease")) %>% is.na(), ]
matrix1 <- matrix1[str_extract(matrix1$time, c("-post disease")) %>% is.na(), ]
  #按time变量分组求均值
  matrix1_average <- aggregate(.~time, matrix1[2:330], mean)
##方案2 比较nec和los
matrix2 <- matrix[str_extract(matrix$time, c("control")) %>% is.na(), ]
matrix2 <- matrix2[str_extract(matrix2$time, c("NEC2-")) %>% is.na(), ]
  #按time变量分组求均值
  matrix2_average <- aggregate(.~time, matrix2[2:330], mean)

#write
write.csv(matrix, "matrix.csv")
write.csv(matrix1, "matrix1.csv")
write.csv(matrix2, "matrix2.csv")
write.csv(matrix1_average, "matrix1_average.csv")
write.csv(matrix2_average, "matrix2_average.csv")


##matrix1——alpha
matrix1_alpha <- inner_join(select(matrix1, c("sample", "time")), alpha, by = "sample")
matrix1_alpha <- matrix1_alpha[,-1]
matrix1_alpha <- aggregate(.~time, matrix1_alpha[2:7], mean)
