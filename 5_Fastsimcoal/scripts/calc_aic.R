df <- read.delim("all_bestlhoods_raw_6k_clean_bestlhoods_per_scenario.txt",header=T,sep='\t')
df$AIC <- 2*df$params-2*(df$MaxEstLhood/log10(exp(1)))
write.table(df,"all_bestlhoods_raw_6k_clean_bestlhoods_per_scenario_AIC.txt",sep='\t',quote=FALSE,row.names=F)
