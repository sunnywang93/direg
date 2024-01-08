

log_transform <- function(a, b, c) {

  if(a > b & b > c) {
    (log(a - b) / (2*log(2))) - (log(b - c) / (2*log(2)))
  } else {
    1
  }

}

H_replace <- function(theta_2delta, theta_delta, sigma) {

  if(theta_2delta > 2*sigma**2 & theta_delta > 2*sigma**2) {
    (log(theta_2delta - 2*sigma**2) / (2*log(2))) -
      (log(theta_delta - 2*sigma**2) / (2*log(2)))
  } else {
    1
  }

}



