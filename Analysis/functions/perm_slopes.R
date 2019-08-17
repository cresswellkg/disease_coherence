slope_perm = function(Conn_Frame, Disease_Name, Norm_Reg = NULL, database, permutations = 5000, Disease = TRUE) {
  
  require(dplyr)
  source("./functions/GeneSample.R")
  
  
  if (Disease) {
  Gene_Frame = Conn_Frame %>% filter(Disease == Disease_Name)
  } else {
  Gene_Frame = Conn_Frame %>% filter(Category == Disease_Name)  
  }
  
  
  if (!is.null(Norm_Reg)) {
    norm_coef_reg = predict(Norm_Reg, newdata = data.frame(num_genes = nrow(Gene_Frame)))
    
    Gene_Frame = Gene_Frame %>% mutate(Internal = Internal/norm_coef_reg) 

  }
  
  #Calculate slope for original frame
  
  reg = lm(External ~ Internal -1,  data = Gene_Frame)
  slope = reg$coefficients[1]
  
  #Get slope from comparison

  real_diff = slope
  
  #Perform Permutation test
  
  Comb_Frame = bind_rows(Gene_Frame, Rand_Frame)
  
  diff = c()
  i = 1
  
 while(i<permutations) {
    
  rand_over = GeneSample(nrow(Gene_Frame), database, directed = FALSE)
  
  rand_over = bind_rows(rand_over$internal_connectivity, rand_over$external_connectivity)
  rand_over = cbind.data.frame(colnames(rand_over), t(rand_over[1,]), t(rand_over[2,]))
  colnames(rand_over) = c("Genes", "Internal", "External")
  
  rand_over = rand_over %>%  mutate(Pred_Edges = predict(pred_edges, data.frame(num_genes = n())), Internal = sqrt(Internal/Pred_Edges), External = sqrt(External))
  
  #Randomly sample
  
    
  new_reg = lm(External ~ Internal-1, data = rand_over)
  new_slope = new_reg$coefficients[1]
  
  if (is.na(new_slope)) {
    next
  } else {
    diff[i] = new_slope
    i = i+1
    
  }
  
  }
  
  p_val = (sum(diff<real_diff) + 1)/(permutations+1)
  return(p_val)
}


