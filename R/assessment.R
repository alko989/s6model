
fitWL <- function(df, colnames = c("Weight", "Length")) {
  df <- setNames(df[colnames], c("Weight", "Length"))
  fit <- lm(log(Weight) ~ log(Length), data = df) 
  a <- exp(fit$coefficients[1])
  b <- fit$coefficients[2]
  list(a = a, b = b)
}

l2w <- function(l, a, b) {
  a * l ^ b
}

addWeight <- function(df, a, b, lengthcol = "Length") {
  if(is.null(df$Weight))
    df$Weight <- l2w(l = df[lengthcol], a = a, b = b)
  else
    warning("Column `Weight` exists, the original data frame `df` is returned")
  df
}