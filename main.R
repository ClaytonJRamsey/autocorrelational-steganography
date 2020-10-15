################ Introduction ########################

# Autocorrelational Steganography

# A technique to hide information in noise data by encoding that data in the autocorrelations over a particular interval.

# Basic algorithm to encode:
# 1. Data is converted to raw binary
# 2. Each bit is used in turn to generate a burst of auto-regressive time series data, where the parameter is based on
#    the bit value. The number of samples and parameters are specified.
# 3. The data is concatenated to form a large time series object.
# 4. The time series object is converted to a sound file.

# Basic algorithm to decode:
# 1. Read the sound file as a time series object.
# 2. Break it into appropriate length samples
# 3. Feed each sample into a function that will extract the PACF at lag 1 and then decide if it reprenta 1 or 0
# 4. Convert the string of binary into the data again.

# NOTE TO SELF: How robust is the autocorrelation of the data if subjected to file format changes such as
# conversion to a compressed format and back?

############## Packages  #######################
if(!require(gtools)){
  install.packages("gtools")
  library(gtools)
}

############## Functions #######################

# A function to generate autoregressive data with a certain parameter.
# As this is an AR(1) the PACF will have value p at lag 1.
gen_ar_data <- function(n, p){
  data_out <- double(n)
  data_out[1] <- runif(1)
  for (i in 2:n){
    data_out[i] <- p*data_out[i-1] + runif(1)
  }
  return(data_out)
}

# A function to convert an ASCII string into a vector of binary values.
ascii_to_bin <- function(m){
  b <- integer()
  codes <- asc(m)[,1] # convert message into dec numeric ASCII codes
  for(i in codes){
    b <- append(b, as.integer(intToBits(i))[1:8])
  }
  return(b)
}

# a function to convert a vector of binary values into an ASCII string
bin_to_ascii <- function(b){
  bin_starts <- seq(1, length(b)-7, 8)
  j <- 1
  m <- integer()
  for(i in bin_starts){
    m[j] <- b[i]+2*b[i+1]+4*b[i+2]+8*b[i+3]+16*b[i+4]+32*b[i+5]+64*b[i+6]+128*b[i+7]
    j <- j+1
  }
  return(paste(chr(m), collapse = ""))
}

# a function to convert a vector of binary values into a time series object in which the samples
# have one of two different autocorrelations for a certain number of samples.
ts_generation <- function(x, n, p0, p1){
  # x is the vector of binary values
  # n is the number of samples of ts data for each bit
  # p0 is the parameter associated with a zero
  # p1 is the parameter associated with a one
  
  t <- ts()
  for(i in x){
    t <- append(t, gen_ar_data(n, ifelse(i, p1, p0)))
  }
  return(t[-1])
}

# a function to take a time series object and return 0 or 1 based on the PACF
ts_to_bit <- function(x, p0, p1){
  p <- pacf(x)
  p <- p$acf[1]
  dist_0 <- abs(p - p0)
  dist_1 <- abs(p - p1)
  bit_out <- ifelse(dist_1 < dist_0, 1, 0)
  return(bit_out)
}

# a wrapper function for the above function that converts the ts into binary piece by pieve
ts_to_bin_converter <- function(x, n, p0, p1){
  # x is the ts
  # n is the size of each chunk that represents a single binary digit
  bin_out <- numeric()
  indices <- seq(1, length(x), n)
  for(i in indices){
    x_current <- x[i:(i+99)]
    bit_current <- ts_to_bit(x_current, p0, p1)
    bin_out <- append(bin_out, bit_current)
  }
  return(bin_out)
}

################################ Diagnostics #################################

# test the ts_to_bit() function to check the error rate for combinations of 
# p0 and p1 that are separated by different amounts.
n <- 100 # samples per bit
k <- 100 # no. of bits
b <- sample(c(0, 1), k, replace = TRUE)
p1seq <- p0seq <- seq(0.1, 0.9, 0.1)
param0 <- param1 <- error_rate <- numeric()

for(p0 in p0seq){
  for(p1 in p1seq){
    param0 <- append(param0, p0)
    param1 <- append(param1, p1)
    t <- ts_generation(b, n, p0, p1)
    b_new <- ts_to_bin_converter(t, n, p0, p1)
    error_r <- sum(b != b_new)/k
    error_rate <- append(error_rate, error_r)
  }
}
e_r_results <- data.frame(param0, param1, error_rate)

################################ Test Message ################################

# Conversion Parameters:
n <- 1000
p0 <- 0.5
p1 <- 0.3

test_message <- "The crow flies at midnight."
binary_version <- ascii_to_bin(test_message)
ts_data <- ts_generation(binary_version, n, p0, p1)

binary_decode_version <- ts_to_bin_converter(ts_data, n, p0, p1)
decode_message <- bin_to_ascii(binary_decode_version)
cat(decode_message)
