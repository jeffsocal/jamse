######################################################################
# AUTH  Jeff Jones -- SoCal Bioinformatics Inc.
# DATE  2018.02.23
# DESC  calculate the peptide descriminate score and estimate the
#       peptide probabilities
#
# install.packages(c("progress","plyr","ggplot2","reshape2","stringr","e1071","rjson"))
######################################################################

rm(list=ls())
library(progress)
library(plyr)
library(ggplot2)
library(reshape2)
library(stringr)
options(warn=-1)

help_text <- "
NAME

	jvlnds.R

REQUIRES

	progress, jsonlite, rjson, plyr, ggplot2, reshape2, stringr, e1071

SYNOPSIS

	jvlnds.R --file=<path to file>


DESCRIPTION

	calculates the discriminant score [psm_dscore] based on PSM metrics.

COMMAND LINE

	--file    : path to and file name

EXAMPLE

	Rscript jvlnds-csv.R --file=/var/jvlnse/dscore/mz5da74e2882b45.csv

"

###############################################################################
# USER INPUT
ui_file                      <- "./"
ui_model_path                <- "/var/jvlnse/bin/dscore_svm_model.rds"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--file", arg) ) ui_file <- arg_value
    if( grepl("--modelpath", arg) ) ui_model_path <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL
if(!file.exists(ui_file)) message <- paste0(message, "  file does not exist\n")
if(!grepl(".csv$", ui_file)) message <- paste0(message, "  csv file (--csv) not a supported format\n")
if(!is.null(message)) stop("ERROR\n", message)

cat("jvln discriminant score\n")

#
# FUNCTION: progress bar
#
progressbar <- function(size, text = "running ..."){
    cat(text, "\n")
    pb <- progress_bar$new(
        format = paste(size, ":percent :elapsed [:bar] :eta"),
        total = size, clear = FALSE, width= 60)
    return(pb)
}

#
# FUNCTION: update the data frame with the modeled discriminant score
#
jvln_dscore <- function(df){
    library(e1071)
    model <- readRDS(ui_model_path)
    fea <- names(model$coefficients)
    if(is.null(fea))
        fea <- colnames(model$SV)
    
    df <- df[!is.na(df$psm_score),]
    df$abs_pre_mma <- abs(df$pre_mma)
    df$pre_mz_log2 <- log2(df$pre_mz)
    df$score <- df$psm_score
    df$psm_dscore <- NA
    df$psm_dscore <- predict(model, newdata = df[,fea])
    df <- ddply(df, 'scan_id', transform,
                hit_rank = rank(-psm_dscore, ties.method = 'first'))
    df <- df[,setdiff(colnames(df), 'score')]
    return(df)
}

#
# read data using jsonlite package
# save data using the rjson package
#
    
df <- c()
d_table <- c()

cat("\n read file                   :", basename(ui_file))
d_table <- read.csv(ui_file)
cat("\n")

cat(" compute discriminant score  :")
d_table <- jvln_dscore(d_table)
cat(" done\n")

cat(" write out file              :")
write.csv(d_table, ui_file, row.names = FALSE)
cat(" done\n\n")
    

