# read the signature table
signature_table = read.csv(
  paste0(intermediate_dir, "SuperGeneSets_5_2023.csv"), 
  header = T, fill = T, check.names = F)
signature_table[signature_table==""] = NA

# convert a table to signature list
signature_list_master = as.list(signature_table)
signature_list_master = lapply(signature_list_master, function(x) x[!is.na(x)])
# refine names for each vector in list
names(signature_list_master) = sapply(names(signature_list_master), function(x){
  x = stringr::str_replace_all(x, " ", "_")
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
})
signature_list_master = lapply(signature_list_master, function(x) trimws(x, which = "both"))
signature_list_master$Complex_V = sapply(signature_list_master$Complex_V, function(x){
  if(grepl("\\(", x)){
    a = stringr::str_split_fixed(x, "\\(", 2)[1,2]
    a = stringr::str_split_fixed(a, "\\)", 2)[1,1]
  }else{
    a = x
  }
  return(a)
})
signature_list_master$Complex_V = as.character(signature_list_master$Complex_V)