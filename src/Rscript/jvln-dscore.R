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

 jvlnds.R --path=<path>


DESCRIPTION

 calculates the discriminant score [psm_dscore] based on PSM metrics.

COMMAND LINE

 --path    : directory path to json file/files [default: ./]
 --pattern : file pattern to restrict analysis to [default: NULL]

EXAMPLE

 Rscript jvlnds.R --path=./

"

###############################################################################
# USER INPUT
ui_path                      <- "./"
ui_pattern                   <- NULL
ui_model_path                <- "/var/jvlnse/bin/dscore_svm_model.rds"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--path", arg) ) ui_path <- arg_value
    if( grepl("--pattern", arg) ) ui_pattern <- arg_value
    if( grepl("--modelpath", arg) ) ui_model_path <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL
# if(is.null(ui_path)) message <- stop("ERROR\n", "  no mzXML file declared\n")
# if(!grepl("pep.*XML$", ui_path, ignore.case = T)) message <- paste0(message, "  mz file (--xml) not a supported format\n")
# if(is.null(path_csv))
#     path_csv = paste0(ui_path, ".csv")
# if(!grepl(".csv$", path_csv)) message <- paste0(message, "  csv file (--csv) not a supported format\n")

if(!is.null(message)) stop("ERROR\n", message)

cat("jvln discriminant score\n")
cat(" path                        :", ui_path, "\n")

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
# FUNCTION: convert the jvln-json (list) to a tibble
#
jvln_tibble <- function(json, file = NULL){
    d_psms <- do.call(rbind.data.frame, json$jaevalen$psms)
    d_psms$scan_id_sub <- row.names(d_psms)
    d_psms$scan_id <- str_replace(d_psms$scan_id_sub, "\\.[1-9]$", "")
    
    rmpc <- function(x) x[ !names(x) %in% "peaks"]
    d_scans <- do.call(rbind.data.frame, lapply(json$scans, rmpc))
    d_scans$scan_id <- row.names(d_scans)
    
    d_table <- merge(d_scans, d_psms, by='scan_id', all=T)
    
    if(!is.null(file))
        d_table$scan_id = paste0(file, "_", d_table$scan_id)
    
    return( d_table )
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
cat(" reading files")
files <- list.files(path = ui_path, pattern = paste0(ui_pattern, ".*json$"), full.names = T)
for ( file in files ){
    
    df <- c()
    d_json <- c()
    d_table <- c()
    
    cat("\n                             :", basename(file))
    
    attach_file <- NULL
    if(length(files) > 1)
        attach_file <- str_replace(basename(file), "\\.json$", "")
    
    d_list <- rjson::fromJSON(file=file)
    d_list_names <- setdiff(names(d_list), c('file', 'scans', 'jaevalen'))
    for( d_list_name in d_list_names){
        d_list[d_list_name] <- NULL    
    }
    
    d_table <- jsonlite::fromJSON(file)
    d_table <- jvln_tibble(d_table)
    cat("\n")
    cat(" compute discriminant score  :")
    d_table <- jvln_dscore(d_table)
    cat(" apply ...")
    for( scan_id in names(d_list$jaevalen$psms)){
        t_ds <- d_table[d_table$scan_id == scan_id,]
        for( i in 1:nrow(t_ds)){
            d_list$jaevalen$psms[[scan_id]][[i]]$psm_dscore <- t_ds$psm_dscore[i]
            d_list$jaevalen$psms[[scan_id]][[i]]$hit_rank <- t_ds$hit_rank[i]
        }
    }
    cat(" done\n")
    
    json_txt <- rjson::toJSON(d_list)    
    cat(" write out file              :")
    write(json_txt, file)
    cat(" done\n\n")
    
}

