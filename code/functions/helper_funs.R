################################################################################
###############---------------- HELPER FUNCTIONS ----------------###############
################################################################################

### Helper functions used in subsequent functions and code

## Population standard deviation -----------------------------------------------
# @ x: A numerical vector of values

sdP <- function(x){
  sqrt(mean(x^2) - (mean(x))^2) 
}

## Standard error --------------------------------------------------------------
# @ x: A numerical vector of values

se <- function(x) sqrt(var(x) / length(x))

## Mode ------------------------------------------------------------------------
# @ x: A numerical vector of values

getmode <- function(x) {
  uniqX <- unique(x)
  uniqX[which.max(tabulate(match(x, uniqX)))]
  
}

## Task data standardization ---------------------------------------------------
#@ x: A numerical vector of values
#@ ref: A single numerical values to be used as the reference
#@ levels: Possible values the data to be standardized can take

standardization <- function(x, ref = 0, levels = 1:4){
  (x - ref) / sdP(levels)
}

## Shuffle subject IDs ---------------------------------------------------------
#@ x: Character string

id_shuffle <- function(x){
  new_id <- unlist(strsplit(x, ""))
  new_id <- new_id[1:10]
  set.seed(1)
  new_id <- sample(new_id, length(new_id), FALSE)
  new_id <- paste(new_id, collapse = "")
  new_id
}


## Pass through piping ---------------------------------------------------------
#@ data: data
#@ fun: function to apply

pass_through <- function(data, fun){
  fun(data)
  data
}




