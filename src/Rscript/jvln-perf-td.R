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

 jvlnds.R --path=<path> --fdr=0.05


DESCRIPTION

 calculates the q-value [psm_qvalue], [psm_eprob] and [psm_prob] for each
 peptide-spectrum-match.

 [psm_eprob] the posterior error probability proposed by Kall, L. et. al J.
 Proteome Res. 7, 40-44 (2008).
 
 [psm_prob] the peptide PSM probability calculated using the dscore
 value according to Keller, A., et. al. Anal. Chem. 74, 5383-5392
 (2002).
 
 [psm_qvalue] the specific FDR associated with a give score value computed
 from a population of suspected positives (target) and known negatives
 (decoy)


COMMAND LINE

 --path    : directory path to json file/files [default: ./]
 --pattern : file pattern to restrict analysis to [default: NULL]

 --metric  : metric to use in error estimation [default: psm_dscore]
 --fdr     : target fdr value [default: 0.05]


EXAMPLE

 Rscript jvlnds.R --path=./

"

###############################################################################
# USER INPUT
ui_path                      <- "./"
ui_pattern                   <- NULL
ui_fdr                       <- 0.05
ui_metric                    <- "psm_dscore"
ui_model                     <- "target-decoy"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--path", arg) ) ui_path <- arg_value
    if( grepl("--pattern", arg) ) ui_pattern <- arg_value
    if( grepl("--fdr", arg) ) ui_fdr <- arg_value
    if( grepl("--metric", arg) ) ui_metric <- arg_value
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
cat(" fdr                         :", ui_fdr, "\n")
cat(" metric                      :", ui_metric, "\n")
cat(" model                       :", ui_model, "\n")

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
    d_psms$scan_id <- row.names(d_psms)
    d_psms$scan_id <- str_replace(d_psms$scan_id, "\\.[1-9]$", "")
    
    rmpc <- function(x) x[ !names(x) %in% "peaks"]
    d_scans <- do.call(rbind.data.frame, lapply(json$scans, rmpc))
    d_scans$scan_id <- row.names(d_scans)
    
    d_table <- merge(d_scans, d_psms, by='scan_id', all=T)
    
    if(!is.null(file))
        d_table$scan_id = paste0(file, "_", d_table$scan_id)
    
    return( d_table )
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
    
    d_table <- jsonlite::fromJSON(file)
    d_table <- jvln_tibble(d_table)
    d_table <- d_table[!is.na(d_table[,ui_metric]),]
    
    df <- d_table
    cat("\n")
    
    
    #
    # compute the discriminant score
    # estimate the probabilities for
    #
    #                  p(F|+)p(+)
    # p(+|F) = _________________________
    #           p(F|+)p(+) + p(F|-)p(-)
    #
    
    cat(" estimate psm performance    :")
    dscore_bins <- 2^8
    dscore_min <- min(df[,ui_metric])
    dscore_max <- max(df[,ui_metric])
    dscore_bwd <- (dscore_max - dscore_min) / (dscore_bins)
    
    df$ui_metric_bin <- round((df[,ui_metric] - dscore_min) / dscore_bwd)
    df$psm_metric <- df[,ui_metric]
    df$correct <- df$hit_rank == 1
    
    #
    # calculate the components of the psm probabilities
    #
    df_model_f <- df[df$correct == FALSE,]
    df_model_t <- df[df$correct == TRUE,]
    
    d_fn <- nrow(df_model_f)
    d_tn <- nrow(df_model_t)
    
    p_pos <- d_tn / (d_fn + d_tn)
    p_neg <- d_fn / (d_fn + d_tn)
    
    #
    # estimate the densities of each kernel
    #
    d_f <- density(df_model_f[,ui_metric])
    d_t <- density(df_model_t[,ui_metric])
    
    #
    # smooth the density estimations
    #
    d_min <- min(d_t$x, d_f$x)
    d_max <- max(d_t$x, d_f$x)
    
    d_t <- ksmooth(d_t$x, d_t$y, bandwidth = dscore_bwd, range.x = c(d_min,d_max), n.points = dscore_bins)
    d_f <- ksmooth(d_f$x, d_f$y, bandwidth = dscore_bwd, range.x = c(d_min,d_max), n.points = dscore_bins)
    
    d_t$y[is.na(d_t$y)] <- 0
    d_f$y[is.na(d_f$y)] <- 0
    
    #
    # create the function to return the density x estimation
    #
    pdist_t <- approxfun(d_t$x, d_t$y)
    pdist_f <- approxfun(d_f$x, d_f$y)
    
    #
    # create the continous distribution functions
    #
    fdist_t <- Vectorize(function(x){
        max(0,min(1, integrate(pdist_t, x, d_max, stop.on.error = F)$value))    
    })
    
    fdist_f <- Vectorize(function(x){
        max(0,min(1, integrate(pdist_f, x, d_max, stop.on.error = F)$value)) 
    })
    
    p_S_pos <- fdist_t(df[,ui_metric])
    p_S_neg <- fdist_f(df[,ui_metric])
    
    #
    # vector of posterior-error probabilities
    #
    psm_eprob <- pdist_f(df[,ui_metric])  / (pdist_f(df[,ui_metric]) + pdist_t(df[,ui_metric]))
    
    #
    # vector of FDR q-values
    #
    psm_qvalue <- p_S_neg / (p_S_pos + p_S_neg)
    
    #
    # vector of probabilites for correct peptide assignment 
    #
    psm_prob <- (p_S_pos * p_pos) / (p_S_pos * p_pos + p_S_neg * p_neg)
    
    
    df$psm_eprob <- psm_eprob
    df$psm_prob <- psm_prob
    df$psm_qvalue <- psm_qvalue
    
    df_rank1 <- df[df$correct == T,]
    row.names(df_rank1) <- df_rank1$scan_id
    df_rank1 <- df_rank1[,c('peptide','psm_metric','psm_eprob','psm_prob','psm_qvalue')]
    
    df_g <- melt(df_rank1, 
                      id.vars=c('peptide', 'psm_metric'),
                      variable.name = "estimate", value.nsm="value")
    
    g_estimates <- ggplot(df_g, aes(psm_metric, value, color=estimate)) + 
        geom_line() + ggtitle("Jaevalen Discriminat Estimates", subtitle = ui_path) + xlab(ui_metric)
    
    ds_test_cut <- df_rank1[order(-df_rank1$psm_qvalue),]
    ds_test_cut <- ds_test_cut[ds_test_cut$psm_qvalue <= ui_fdr,]
    
    df_rank1 <- df_rank1[,c('psm_eprob','psm_prob','psm_qvalue')]
    colnames(df_rank1)[colnames(df_rank1) == 'psm_metric'] <- ui_metric
    ds_results <- setNames(split(df_rank1, seq(nrow(df_rank1))), rownames(df_rank1))
    
    ds_summary <- list(file = basename(file),
                       modeling_method = ui_model,
                       fdr_target = 0.05,
                       fdr_observed = ds_test_cut$psm_qvalue[1],
                       cutoff_metric = ui_metric,
                       cutoff_value = ds_test_cut$psm_metric[1],
                       peptide_pass_n = nrow(ds_test_cut), 
                       peptide_pass_r = round(nrow(ds_test_cut)/nrow(df_rank1) * 100, 2))
    
	cat(" done\n");
    cat(" write out file              :")
    d_list$performance <- list()
    d_list$performance$summary <- ds_summary
    d_list$performance$psms <- ds_results
    
    json_txt <- rjson::toJSON(d_list)
    
    write(json_txt, file)
    cat(" done\n\n")
    
}

