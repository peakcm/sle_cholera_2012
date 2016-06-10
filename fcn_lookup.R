fcn_lookup <- function(query_1, query_2 = NA, reference, value_column, transformation = "none"){
  out <- 0
  reference <- data.frame(reference)
  if (is.na(query_2)==0){
    row_1 <- which(reference[,1] %in% query_1) 
    row_2 <- which(reference[,2] %in% query_2)
    if (length(row_1) > 0 & length(row_2) > 0){
      out <- value_column[intersect(row_1, row_2)]
    }
  } else {
    out <- value_column[which(reference == query_1)]
  }
  if (length(out)==1){
    if (transformation == "as.character"){
      return(as.character(out))
    } else {return(out)}
  } else {
    return(0)
  }
}