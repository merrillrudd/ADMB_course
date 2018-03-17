y <- c(1.4, 4.7, 5.1, 8.3, 9.0, 14.5, 14.0, 13.4, 19.2, 18.0)

n <- length(y)
mu <- mean(y)
RSS <- sum((y-mu)^2)
sigma <- sqrt(RSS/n)

RSS
# Residual sum of squares 320.624

0.5*n*log(2*pi) + n*log(sigma) + RSS/(2*sigma^2)
# Likelihood 31.52781


0.5*n*log(2*pi) + 0.5*n*log(RSS) - 0.5*n*log(n) + n/2;
# Concentrated 31.52781

-sum(dnorm(y, mu, sigma, log=TRUE))
# dnorm 31.52781
