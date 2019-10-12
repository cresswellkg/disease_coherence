require(dplyr)
files = list()
index = 1
j = 10
for(i in seq(1,71, by = 10)) {
	
		files[[index]] = readRDS(paste0("string_filt_p_", i, "_", j, ".rds"))
j = j+10
index = index+1

}

saveRDS(do.call(rbind, files), "Complete_10k_String_Filt.rds")