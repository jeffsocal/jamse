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

	calculates the q-value [fdr_qvalue], [fdr_eprob] and [fdr_prob] for each
	peptide-spectrum-match.

	[fdr_eprob] the posterior error probability proposed by Kall, L. et. al J.
	Proteome Res. 7, 40-44 (2008).

	[fdr_prob] the peptide PSM probability calculated using the dscore
	value according to Keller, A., et. al. Anal. Chem. 74, 5383-5392
	(2002).

	[fdr_qvalue] the specific FDR associated with a give score value computed
	from a population of suspected positives (target) and known negatives
	(decoy)


COMMAND LINE

	--file    : path to and file name

	--metric  : metric to use in error estimation [default: fdr_dscore]
	--fdr     : target fdr value [default: 0.05]


EXAMPLE

 	Rscript jvlnds-csv.R --file=/var/jvlnse/dscore/mz5da74e2882b45.csv

"

###############################################################################
# USER INPUT
ui_file                      <- "./"
ui_fdr                       <- 0.05
ui_metric                    <- "psm_dscore"
ui_model                     <- "kernel"

for (arg in commandArgs()){
    arg_value <- as.character(sub("--[a-z]*\\=", "", arg))
    if( grepl("--file", arg) ) ui_file <- arg_value
    if( grepl("--fdr", arg) ) ui_fdr <- arg_value
    if( grepl("--metric", arg) ) ui_metric <- arg_value
    if( grepl("--model", arg) ) ui_model <- arg_value
    if( grepl("--help", arg) ) stop(help_text)
}

###############################################################################
# INPUT VALIDATION
message <- NULL
if(!file.exists(ui_file)) message <- paste0(message, "  file does not exist\n")
if(!grepl(".csv$", ui_file)) message <- paste0(message, "  csv file (--csv) not a supported format\n")
if(!is.null(message)) stop("ERROR\n", message)

cat("jvln discriminant score\n")
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
# read data using jsonlite package
# save data using the rjson package
#
df <- c()
d_table <- c()

cat(" read file                   :", basename(ui_file))
d_table <- read.csv(ui_file)
cat("\n")

df <- d_table


#
# compute the discriminant score
# estimate the probabilities for
#
#                  p(F|+)p(+)
# p(+|F) = _________________________
#           p(F|+)p(+) + p(F|-)p(-)
#

cat(" estimate probabilities      :")
dscore_bins <- 2^8
dscore_min <- min(df[,ui_metric])
dscore_max <- max(df[,ui_metric])
dscore_bwd <- (dscore_max - dscore_min) / (dscore_bins)

df$ui_metric_bin <- round((df[,ui_metric] - dscore_min) / dscore_bwd)
df$fdr_metric <- df[,ui_metric]
df$correct <- df$hit_rank == 1

df_test <- df[df$hit_rank == 1,]
df_null <- df[df$hit_rank %in% 3:4,]

all_incorrect <- FALSE
for(i in dscore_bins:1 ){
    w <- which(df_test$ui_metric_bin == i)
    # subtract by raw values
    df_null_frac <- nrow(df_null[df_null$ui_metric_bin == i,])
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

p2 <- ggplot(df_test, aes(fdr_metric)) +
    geom_histogram(binwidth = dscore_bwd*4, position='identity',
                   aes(fill=as.factor(correct), alpha=1/2)) + xlab(ui_metric)

false_gaussian_mean <- mean(df_null[,ui_metric])
false_gaussian_sd <- sd(df_null[,ui_metric])
false_gaussian <- data.frame(fdr_metric = rnorm(50000, mean=false_gaussian_mean, sd=false_gaussian_sd),
                             correct = TRUE)

true_gaussian_mean <- mean(df_test[df_test$correct == TRUE,][,ui_metric])
true_gaussian_sd <- sd(df_test[df_test$correct == TRUE,][,ui_metric])
true_gaussian <- data.frame(fdr_metric = rnorm(50000, mean=true_gaussian_mean, sd=true_gaussian_sd),
                            correct = TRUE)

#
# calculate the components of the psm probabilities
#
df_model_f <- df_test[df_test$correct == FALSE,]
df_model_t <- df_test[df_test$correct == TRUE,]

d_fn <- nrow(df_model_f)
d_tn <- nrow(df_model_t)

p_pos <- d_tn / (d_fn + d_tn)
p_neg <- d_fn / (d_fn + d_tn)

df_test_org <- df_test
g_distributions <- ggplot(df_test_org[df_test_org$correct == T,], aes(fdr_metric)) +
    # DIST CORRECT PSMS
    geom_density(position='identity', fill='dodgerblue', color=NA, alpha=1/3) +
    geom_density(data = true_gaussian, position='identity', linetype = 2, color='blue') +
    # DIST FALSE PSMS
    geom_density(data = df_null, position='identity', fill='red', color=NA, alpha=1/3) +
    geom_density(data = false_gaussian, position='identity', linetype = 2, color='red') +
    scale_fill_brewer(palette = 'Set1') +
    ggtitle("Jaevalen Discriminat Densities", subtitle = basename(ui_file)) +
    xlab(ui_metric)



df_test <- df_test_org
if(ui_model == 'gaussian'){
    p_S_pos <- pnorm(df_test[,ui_metric], true_gaussian_mean, true_gaussian_sd, lower.tail = F)
    p_S_neg <- pnorm(df_test[,ui_metric], false_gaussian_mean, false_gaussian_sd, lower.tail = F)

    d_S_pos <- dnorm(df_test[,ui_metric], true_gaussian_mean, true_gaussian_sd)
    d_S_neg <- dnorm(df_test[,ui_metric], false_gaussian_mean, false_gaussian_sd)

    fdr_eprob <- d_S_neg  / (d_S_pos + d_S_neg)
} else {
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

    p_S_pos <- fdist_t(df_test[,ui_metric])
    p_S_neg <- fdist_f(df_test[,ui_metric])

    #
    # vector of posterior-error probabilities
    #
    fdr_eprob <- pdist_f(df_test[,ui_metric])  / (pdist_f(df_test[,ui_metric]) + pdist_t(df_test[,ui_metric]))
}

#
# vector of FDR q-values
#
fdr_qvalue <- p_S_neg / (p_S_pos + p_S_neg)

#
# vector of probabilites for correct peptide assignment
#
fdr_prob <- (p_S_pos * p_pos) / (p_S_pos * p_pos + p_S_neg * p_neg)

df_test$fdr_eprob <- fdr_eprob
df_test$fdr_prob <- fdr_prob
df_test$fdr_qvalue <- fdr_qvalue

df_test <- df_test[,c('pk','fdr_metric','fdr_eprob','fdr_prob','fdr_qvalue')]
df_test_g <- melt(df_test,
                  id.vars=c('pk', 'fdr_metric'),
                  variable.name = "estimate", value.nsm="value")

g_estimates <- ggplot(df_test_g, aes(fdr_metric, value, color=estimate)) +
    geom_line() + ggtitle("Jaevalen Discriminat Estimates", subtitle = basename(ui_file)) + xlab(ui_metric)

ds_test_cut <- df_test[order(-df_test$fdr_qvalue),]
ds_test_cut <- ds_test_cut[ds_test_cut$fdr_qvalue <= ui_fdr,]

ds_summary <- list(file = sub(".csv", "", basename(ui_file)),
                   modeling_method = ui_model,
                   fdr_target = 0.05,
                   fdr_observed = ds_test_cut$fdr_qvalue[1],
                   cutoff_metric = ui_metric,
                   cutoff_value = ds_test_cut$fdr_metric[1],
                   peptide_pass_n = nrow(ds_test_cut),
                   peptide_pass_r = round(nrow(ds_test_cut)/nrow(df_test) * 100, 2))


cat(" write out file              :")
write.csv(df_test, ui_file, row.names = FALSE)
json_txt <- rjson::toJSON(ds_summary)
write(json_txt, sub(".csv", "_summary.json", ui_file))
cat(" done\n\n")
