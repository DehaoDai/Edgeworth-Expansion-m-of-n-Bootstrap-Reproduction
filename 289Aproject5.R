# setwd("~/Documents/R/math289")

set.seed(2952) # Set seed for reproducibility

# tdistribution <- curve(dt(x, df=5), from=-4, to=4, 
#       main = 't Distribution (df = 5)', #add title
#       ylab = 'Density', #change y-axis label
#       lwd = 2, #increase line width to 2
#       col = 'steelblue') #change line color to steelblue
# 
# curve(df(x, df1=10, df2=10), from=-4, to=4, 
#       main = 'F Distribution (df1 = 10, df2 = 10)', #add title
#       ylab = 'Density', #change y-axis label
#       lwd = 2, #increase line width to 2
#       col = 'red') #change line color to red


# For the F(10,10) distribution, the median is 1
f_theta <- df(x = 1, df1 = 10, df2 = 10) # the density at the median
sigma <- 1 / (2 * f_theta) # asymptotic standard deviation of the sample median
B <- 1000 # number of bootstrap repeats
num_repeat <- 1000
ns <- c(100, 200, 400)
betas <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

# sample(x = dat, size = round(length(dat)^(2/5)), replace = TRUE)

# estimated density at the median
# density(dat, kernel = "epanechnikov", n = 1, 
#         from = median(dat), to = median(dat))$y

library(tidyverse)

create_tibble <- function() {
  tibble(
    lower = numeric(num_repeat),
    upper = numeric(num_repeat),
    length = numeric(num_repeat),
    contain = logical(num_repeat)
  )
}

generate_f <- function(n) {
  rf(n = n, df1 = 10, df2 = 10)
}

generate_t <- function(n) {
  rt(n = n, df = 5)
}

n_out_of_n <- function(n, theta, generate) {
  result_percentile <- create_tibble()
  result_root <- create_tibble()
  for (i in 1:num_repeat) {
    dat <- generate(n)
    theta_hat <- median(dat)
    bootstrap_dist <- numeric(B)
    for (b in 1:B) {
      bootstrap_sample <- sample(x = dat, size = n, replace = TRUE)
      bootstrap_dist[b] <- median(bootstrap_sample)
    }
    result_root$lower[i] <- 2 * theta_hat - quantile(bootstrap_dist, probs = 0.975)
    result_root$upper[i] <- 2 * theta_hat - quantile(bootstrap_dist, probs = 0.025)
    result_root$length[i] <- result_root$upper[i] - result_root$lower[i]
    result_root$contain[i] <- (result_root$lower[i] <= theta) && 
      (result_root$upper[i] >= theta)
    
    result_percentile$lower[i] <- quantile(bootstrap_dist, probs = 0.025)
    result_percentile$upper[i] <- quantile(bootstrap_dist, probs = 0.975)
    result_percentile$length[i] <- result_percentile$upper[i] - result_percentile$lower[i]
    result_percentile$contain[i] <- (result_percentile$lower[i] <= theta) && 
      (result_percentile$upper[i] >= theta)
  }
  result <- list("result_percentile" = result_percentile,
                 "result_root" = result_root)
  return(result)
}

f_n_out_of_n <- function(n) {
  result <- n_out_of_n(n = n, theta = 1, generate = generate_f)
  write.csv(result$result_percentile, 
            str_c("f_percentile", format(n, scientific = FALSE), ".csv"), 
            row.names = FALSE)
  write.csv(result$result_root, 
            str_c("f_root", format(n, scientific = FALSE), ".csv"), 
            row.names = FALSE)
}

for(i in ns) {
  f_n_out_of_n(i)
}


t_n_out_of_n <- function(n) {
  result <- n_out_of_n(n = n, theta = 0, generate = generate_t)
  write.csv(result$result_percentile, 
            str_c("t_percentile", format(n, scientific = FALSE), ".csv"), 
            row.names = FALSE)
  write.csv(result$result_root, 
            str_c("t_root", format(n, scientific = FALSE), ".csv"), 
            row.names = FALSE)
}

for(i in ns) {
  t_n_out_of_n(i)
}

# Read results

for (n in ns) {
  result <- read.csv(str_c("f_percentile", format(n, scientific = FALSE), ".csv"))
  print(str_c("F distribution, percentile interval, n = ", n, 
        " Coverage probability = ", mean(result$contain),
        " Mean length = ", round(mean(result$length), 2)))
}

for (n in ns) {
  result <- read.csv(str_c("f_root", format(n, scientific = FALSE), ".csv"))
  print(str_c("F distribution, root interval, n = ", n, 
              " Coverage probability = ", mean(result$contain),
              " Mean length = ", round(mean(result$length), 2)))
}

for (n in ns) {
  result <- read.csv(str_c("t_percentile", format(n, scientific = FALSE), ".csv"))
  print(str_c("t distribution, percentile interval, n = ", n, 
              " Coverage probability = ", mean(result$contain),
              " Mean length = ", round(mean(result$length), 2)))
}

for (n in ns) {
  result <- read.csv(str_c("t_root", format(n, scientific = FALSE), ".csv"))
  print(str_c("t distribution, root interval, n = ", n, 
              " Coverage probability = ", mean(result$contain),
              " Mean length = ", round(mean(result$length), 2)))
}

m_out_of_n <- function(n, beta, theta, generate, replace) {
  result <- create_tibble()
  for (i in 1:num_repeat) {
    dat <- generate(n)
    theta_hat <- median(dat)
    bootstrap_dist <- numeric(B)
    m <- round(n^beta)
    for (b in 1:B) {
      bootstrap_sample <- sample(x = dat, size = m, replace = replace)
      bootstrap_dist[b] <- median(bootstrap_sample)
    }
    # MODIFICATION
    root_dist <- sqrt(m / n) * (bootstrap_dist - theta_hat)
    result$lower[i] <- theta_hat - quantile(root_dist, probs = 0.975)
    result$upper[i] <- theta_hat - quantile(root_dist, probs = 0.025)
    result$length[i] <- result$upper[i] - result$lower[i]
    result$contain[i] <- (result$lower[i] <= theta) && (result$upper[i] >= theta)
  }
  return(result)
}

f_m_out_of_n <- function(n, beta, replace) {
  result <- m_out_of_n(n = n, beta = beta, theta = 1,
                       generate = generate_f, replace = replace)
  write.csv(result, 
            str_c("f_n", format(n, scientific = FALSE),
                  "m", format(beta, scientific = FALSE),
                  ifelse(replace, "boot", "sub"),
                  ".csv"), 
            row.names = FALSE)
}

for (n in ns) {
  for (beta in betas) {
    f_m_out_of_n(n, beta, TRUE)
    f_m_out_of_n(n, beta, FALSE)
  }
}

t_m_out_of_n <- function(n, beta, replace) {
  result <- m_out_of_n(n = n, beta = beta, theta = 0,
                       generate = generate_t, replace = replace)
  write.csv(result, 
            str_c("t_n", format(n, scientific = FALSE),
                  "m", format(beta, scientific = FALSE),
                  ifelse(replace, "boot", "sub"),
                  ".csv"), 
            row.names = FALSE)
}

for (n in ns) {
  for (beta in betas) {
    t_m_out_of_n(n, beta, TRUE)
    t_m_out_of_n(n, beta, FALSE)
  }
}

for (n in ns) {
  for (beta in betas) {
    result <- read.csv(str_c("f_n", format(n, scientific = FALSE),
                             "m", format(beta, scientific = FALSE),
                             "boot", ".csv"))
    print(str_c("F distribution, bootstrap, n = ", n, " beta = ", beta, 
          " Coverage probability = ", mean(result$contain), 
          " Mean length = ", round(mean(result$length), 2)))
  }
}

for (n in ns) {
  for (beta in betas) {
    result <- read.csv(str_c("f_n", format(n, scientific = FALSE),
                             "m", format(beta, scientific = FALSE),
                             "sub", ".csv"))
    print(str_c("F distribution, subsampling, n = ", n, " beta = ", beta, 
                " Coverage probability = ", mean(result$contain), 
                " Mean length = ", round(mean(result$length), 2)))
  }
}

for (n in ns) {
  for (beta in betas) {
    result <- read.csv(str_c("t_n", format(n, scientific = FALSE),
                             "m", format(beta, scientific = FALSE),
                             "boot", ".csv"))
    print(str_c("t distribution, bootstrap, n = ", n, " beta = ", beta, 
                " Coverage probability = ", mean(result$contain), 
                " Mean length = ", round(mean(result$length), 2)))
  }
}

for (n in ns) {
  for (beta in betas) {
    result <- read.csv(str_c("t_n", format(n, scientific = FALSE),
                             "m", format(beta, scientific = FALSE),
                             "sub", ".csv"))
    print(str_c("t distribution, subsampling, n = ", n, " beta = ", beta, 
                " Coverage probability = ", mean(result$contain), 
                " Mean length = ", round(mean(result$length), 2)))
  }
}
