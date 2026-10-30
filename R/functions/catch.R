# Catch function
# Author: Domagoj Kifer
# Date: May 10, 2020
# Amended HO


# Catch error: https://stackoverflow.com/a/24569739
ch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
#  list(v=value, w=warn, e=err)
  list( w=warn, e=err)

}