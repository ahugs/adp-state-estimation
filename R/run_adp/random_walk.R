transition_func <- function(y_t, y_tm1){
  if((y_t == y_tm1 + 1)|(y_t == y_tm1 - 1)|(y_t == y_tm1)){
    return(log(1/3))
  }
  else {
    return(log(0))
  }
}
