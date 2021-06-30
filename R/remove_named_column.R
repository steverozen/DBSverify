# Function to delete columns by name
remove_named_column <- function(df, col) {
  for (cname in col) {
    ii <- which(colnames(df) == cname)
    if (length(ii) > 0) {
      df <- df[ , -ii]
    }
  }
  return(df)
}
