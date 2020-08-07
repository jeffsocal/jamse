######################################################################
# AUTH  Jeff Jones -- SoCal Bioinformatics Inc.
# DATE  2018.07.23
# DESC  create the svm model locally
######################################################################

rm(list=ls())

library(e1071)
library(ROCR)

df_jm <- read.csv("./yeast_discr-data_20191003.csv")
df_jm$psm_score <- df_jm$score
df_jm$pre_mz_log2 <- log2(df_jm$pre_mz)

cols <- c('binary',
          'psm_score',
          'a_bins',
          'int_ntp',
          'psm_o',
          'psm_r',
          'abs_pre_mma',
          'pre_mz_log2'
)

############################################################
# ln model
fm <- as.formula(paste(cols[1],"~",paste(c(cols[-1],0), collapse = " + ")))
model <- lm(data=df_jm, formula = fm)

df_jm$pred_model <- (predict(model, df_jm) + 0.4 ) / 2.2

# ggplot(df_jm, aes(ppv, fill=as.factor(hit_rank))) + geom_density(alpha=1/2)

df_jm[df_jm$psm_z == 0,]$psm_z <- 1
df_jm$pre_mz_log2 = ceiling(log2(df_jm$pre_mz * df_jm$psm_z))
df_jm[df_jm$pre_mz_log2 >= 13, ]$pre_mz_log2  <- 12
df_jm$avg_length = ceiling((df_jm$pre_mz * df_jm$psm_z) / 228)

pred <- prediction(df_jm$score, df_jm$binary)
perf <- performance(pred, 'auc')
cat('score     - auc', unlist(perf@y.values), "\n")


pred <- prediction(df_jm$pred_model, df_jm$binary)
perf <- performance(pred, 'auc')
cat('new model - auc', unlist(perf@y.values), "\n")

print(summary(model))

saveRDS(model, file="./jvln_discr-model_20191003.rds")