require(dplyr)
files = list()
index = 1
j = 10
for(i in seq(1,21, by = 10)) {
	
		files[[index]] = readRDS(paste0("biogrid_p_", i, "_", j, ".rds"))
j = j+10
index = index+1

}

saveRDS(do.call(rbind, files), "Complete_10k_Biogrid.rds")