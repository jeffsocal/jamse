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
 Proteome Res. 7, 40–44 (2008).
 
 [psm_prob] the peptide PSM probability calculated using the dscore
 value according to Keller, A., et. al. Anal. Chem. 74, 5383–5392
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
ui_model                     <- "kernel"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--path", arg) ) ui_path <- arg_value
    if( grepl("--pattern", arg) ) ui_pattern <- arg_value
    if( grepl("--fdr", arg) ) ui_fdr <- arg_value
    if( grepl("--metric", arg) ) ui_metric <- arg_value
    if( grepl("--model", arg) ) ui_model <- arg_value
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
# FUNCTION: update the data frame with the modeled discriminant score
#
jvln_dscore <- function(df){
    library(e1071)
    model <- readRDS("/var/jvlnse/bin/dscore_svm_model.rds")
    fea <- names(model$coefficients)
    df <- df[!is.na(df$psm_score),]
    df$abs_pre_mma <- abs(df$pre_mma)
    df$pre_mz_log2 <- log2(df$pre_mz)
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
files <- list.files(path = ui_path, pattern = ui_pattern, full.names = T)
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
    d_table <- jvln_dscore(d_table)
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
	
	cat(" compute discriminant score  :")
	dscore_bins <- 2^8
	dscore_min <- min(df[,ui_metric])
	dscore_max <- max(df[,ui_metric])
	dscore_bwd <- (dscore_max - dscore_min) / (dscore_bins)
	
	df$ui_metric_bin <- round((df[,ui_metric] - dscore_min) / dscore_bwd)
	df$psm_metric <- df[,ui_metric]
	
	df_test <- df[df$hit_rank == 1,]
	df_test$correct <- TRUE
	df_null <- df[df$hit_rank == 2,]
	
	all_incorrect <- FALSE
	for(i in dscore_bins:1 ){
	    w <- which(df_test$ui_metric_bin == i)
	    # subtract by raw values
	    df_null_frac <- dim(df_null[df_null$ui_metric_bin == i,])[1]
	    n <- min(length(w), df_null_frac)
	    if(n == 0) next()
	    # provide a method to remove true values left of null
	    if(n == length(w)) all_incorrect <- TRUE
	    if(all_incorrect == TRUE) n <- length(w)
	    if(n == 0) next()
	    w <- sample(w, n)
	    df_test$correct[w] <- FALSE
	}
	cat(" done\n")
	
	p2 <- ggplot(df_test, aes(psm_metric)) + 
	    geom_histogram(binwidth = dscore_bwd*4, position='identity',
	                   aes(fill=as.factor(correct), alpha=1/2)) + xlab(ui_metric)
	cat(" estimate probabilities      :")
	
	false_gaussian_mean <- mean(df_null[,ui_metric])
	false_gaussian_sd <- sd(df_null[,ui_metric])
	false_gaussian <- data.frame(psm_metric = rnorm(50000, mean=false_gaussian_mean, sd=false_gaussian_sd),
	                             correct = TRUE)
	
	true_gaussian_mean <- mean(df_test[df_test$correct == TRUE,][,ui_metric])
	true_gaussian_sd <- sd(df_test[df_test$correct == TRUE,][,ui_metric])
	true_gaussian <- data.frame(psm_metric = rnorm(50000, mean=true_gaussian_mean, sd=true_gaussian_sd),
	                            correct = TRUE)
	
	# gaussian model method
	# calculate the components of the psm probabilities
	df_model_f <- df_test[df_test$correct == FALSE,]
	df_model_t <- df_test[df_test$correct == TRUE,]
	
	d_fn <- nrow(df_model_f)
	d_tn <- nrow(df_model_t)
	
	p_pos <- d_tn / (d_fn + d_tn)
	p_neg <- d_fn / (d_fn + d_tn)
	
	df_test_org <- df_test
	g_distributions <- ggplot(df_test_org[df_test_org$correct == T,], aes(psm_metric)) + 
	    # DIST CORRECT PSMS
	    geom_density(position='identity', fill='dodgerblue', color=NA, alpha=1/3) + 
	    geom_density(data = true_gaussian, position='identity', linetype = 2, color='blue') +
	    # DIST FALSE PSMS
	    geom_density(data = df_null, position='identity', fill='red', color=NA, alpha=1/3) + 
	    geom_density(data = false_gaussian, position='identity', linetype = 2, color='red') +
	    scale_fill_brewer(palette = 'Set1') +
	    ggtitle("Jaevalen Discriminat Densities", subtitle = ui_path) +
	    xlab(ui_metric)
	
	df_test <- df_test_org
	cat("", ui_model)
	if(ui_model == 'gaussian'){
	    p_S_pos <- pnorm(df_test[,ui_metric], true_gaussian_mean, true_gaussian_sd, lower.tail = F)
	    p_S_neg <- pnorm(df_test[,ui_metric], false_gaussian_mean, false_gaussian_sd, lower.tail = F)
	    
	    d_S_pos <- dnorm(df_test[,ui_metric], true_gaussian_mean, true_gaussian_sd)
	    d_S_neg <- dnorm(df_test[,ui_metric], false_gaussian_mean, false_gaussian_sd) 
	    
	    psm_eprob <- d_S_neg  / (d_S_pos + d_S_neg)
	} else {
	    # d_f <- density(rnorm(d_fn, mean = mean(df_model_t[,ui_metric]), sd=sd(df_model_t[,ui_metric])))
	    d_f <- density(df_model_f[,ui_metric])
	    d_t <- density(df_model_t[,ui_metric])
	    
	    # smooth the density estimations
	    d_min <- min(d_t$x, d_f$x)
	    d_max <- max(d_t$x, d_f$x)
	    
	    d_t <- ksmooth(d_t$x, d_t$y, bandwidth = dscore_bwd, range.x = c(d_min,d_max), n.points = dscore_bins)
	    d_f <- ksmooth(d_f$x, d_f$y, bandwidth = dscore_bwd, range.x = c(d_min,d_max), n.points = dscore_bins)
	    
	    d_t$y[is.na(d_t$y)] <- 0
	    d_f$y[is.na(d_f$y)] <- 0
	    
	    # create the function to return the density x estimation
	    pdist_t <- approxfun(d_t$x, d_t$y)
	    pdist_f <- approxfun(d_f$x, d_f$y)
	    
	    # create the continous distribution functions
	    fdist_t <- Vectorize(function(x){
	        max(0,min(1, integrate(pdist_t, x, d_max, stop.on.error = F)$value))    
	    })
	    
	    fdist_f <- Vectorize(function(x){
	        max(0,min(1, integrate(pdist_f, x, d_max, stop.on.error = F)$value)) 
	    })
	    
	    p_S_pos <- fdist_t(df_test[,ui_metric])
	    p_S_neg <- fdist_f(df_test[,ui_metric])
	    
	    psm_eprob <- pdist_f(df_test[,ui_metric])  / (pdist_f(df_test[,ui_metric]) + pdist_t(df_test[,ui_metric]))
	}
	
	psm_qvalue <- p_S_neg / (p_S_pos + p_S_neg)
	psm_prob <- (p_S_pos * p_pos) / (p_S_pos * p_pos + p_S_neg * p_neg)
	
	df_test$psm_eprob <- psm_eprob
	df_test$psm_prob <- psm_prob
	df_test$psm_qvalue <- psm_qvalue
	
	row.names(df_test) <- df_test$scan_id
	df_test <- df_test[,c('peptide','psm_metric','psm_eprob','psm_prob','psm_qvalue')]
	
	df_test_g <- melt(df_test, 
	                  id.vars=c('peptide'),
	                  variable.name = "estimate", value.nsm="value")
	
	g_estimates <- ggplot(df_test_g, aes(psm_metric, value, color=estimate)) + 
	    geom_line() + ggtitle("Jaevalen Discriminat Estimates", subtitle = ui_path) + xlab(ui_metric)
	
	ds_test_cut <- df_test[order(-df_test$psm_qvalue),]
	ds_test_cut <- ds_test_cut[ds_test_cut$psm_qvalue <= 0.05,]
	
	colnames(df_test)[colnames(df_test) == 'psm_metric'] <- ui_metric
	ds_results <- setNames(split(df_test, seq(nrow(df_test))), rownames(df_test))
	
	ds_summary <- list(files = basename(files),
	                   modeling_method = ui_model,
	                   fdr_target = 0.05,
	                   fdr_observed = ds_test_cut$psm_qvalue[1],
	                   cutoff_metric = ui_metric,
	                   cutoff_value = ds_test_cut$psm_metric[1],
	                   peptide_pass_n = nrow(ds_test_cut), 
	                   peptide_pass_r = round(nrow(ds_test_cut)/nrow(df_test) * 100, 2))
	
	cat("\n")
	
    cat(" write out file              :")
    d_list$performance <- list()
    d_list$performance$summary <- ds_summary
    d_list$performance$psms <- ds_results
    
    json_txt <- rjson::toJSON(d_list)
    
    write(json_txt, file)
    cat(" done\n\n")

}

