check_dir <- function(dir_path){
  # Check if the directory exists
  if (!file.exists(dir_path)) {
    # If it doesn't exist, create it
    dir.create(dir_path)
    cat("Directory created successfully.\n")
  } else {
    cat("Directory already exists.\n")
  }
}
# create label tags based on cutoffs
create_label_tag = function(break_cutoffs){
  result = c()
  for(x in seq(1, length(break_cutoffs))){
    if(x == 1 & break_cutoffs[1] == -Inf){
      result[x] = paste0("<=", break_cutoffs[x+1])
      next
    }
    if(x == length(break_cutoffs)-1 & break_cutoffs[length(break_cutoffs)] == Inf){
      result[x] = paste0(">", break_cutoffs[x])
      break
    }
    result[x] = paste0(break_cutoffs[x], " to ", break_cutoffs[x+1])
  }
  return(result)
}