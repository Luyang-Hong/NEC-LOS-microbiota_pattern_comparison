library(stringr)

#creating metadata files
metadata <- as.data.frame(sample.names)
metadata$patID <- substr(metadata$sample.names, 1, 3)

#分组
NECID <- c("B04", "B06", "B38", "B40")
LOSID <- c("B14", "B18", "B24")
metadata[metadata$patID == NECID, "group"] <- "NEC"
metadata[metadata$patID == LOSID, "group"] <- "LOS"
metadata[is.na(metadata$group), "group"] <- "control"

#patinfo
metadata[metadata$PatNo == "NEC1", "los"] <- 63
metadata[metadata$PatNo == "NEC2", "los"] <- 22
metadata[metadata$PatNo == "NEC3", "los"] <- 49
metadata[metadata$PatNo == "NEC4", "los"] <- 83
metadata[metadata$PatNo == "LOS1", "los"] <- 65
metadata[metadata$PatNo == "LOS2", "los"] <- 71
metadata[metadata$PatNo == "LOS3", "los"] <- 44
metadata[metadata$PatNo == "control1", "los"] <- 58
metadata[metadata$PatNo == "control2", "los"] <- 56
metadata[metadata$PatNo == "control3", "los"] <- 23
metadata[metadata$PatNo == "control4", "los"] <- 46
metadata[metadata$PatNo == "control5", "los"] <- 22
metadata[metadata$PatNo == "control6", "los"] <- 36
metadata[metadata$PatNo == "control7", "los"] <- 30
metadata[metadata$PatNo == "control8", "los"] <- 53
metadata[metadata$PatNo == "control9", "los"] <- 37
metadata[metadata$PatNo == "control10", "los"] <- 29
metadata[metadata$PatNo == "control11", "los"] <- 29
metadata[metadata$PatNo == "control12", "los"] <- 23
metadata[metadata$PatNo == "control13", "los"] <- 23
metadata[metadata$PatNo == "control14", "los"] <- 19
metadata[metadata$PatNo == "control15", "los"] <- 26
metadata[metadata$PatNo == "control16", "los"] <- 27
metadata[metadata$PatNo == "control17", "los"] <- 22
#收集天数
metadata$dol <- str_sub(metadata$sample.names, -2, -1)
metadata$dol <- as.numeric(metadata$dol)

#source
metadata[, "source"] <- paste(metadata$PatNo, "day", metadata$dol, sep = "-")

write.csv(metadata, file = "metadata.csv")

