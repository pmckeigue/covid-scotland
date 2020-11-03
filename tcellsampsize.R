cell.freqs <- function(freq.diff, p, q) {
    ## get cell freqs a, b, c, d given marginal freqs, freq.diff
    a <- p * q + freq.diff
    b <- p * (1 - q) - freq.diff
    c <- (1 - p) * q - freq.diff
    d <- (1 - p) * (1 - q) + freq.diff
    freqs <- c(a, b, c, d)
    stopifnot(all(freqs >= 0))
    names(freqs) <- c("a", "b", "c", "d")
    return(freqs)
}

logoddsratio.err <- function(freq.diff, p, q, theta) {
    ## returns objective function: squared difference of log odds ratio from target
    ## p prob for row 1
    ## q prob ror col 1
    ## freq.diff difference of cell freqs from value at beta=0
    freqs <- cell.freqs(freq.diff, p, q)
    oddsratio.calc <- freqs[1] * freqs[4] / (freqs[2] * freqs[3])
    return <- (log(oddsratio.calc) - theta)^2
}

range.freqs.diff <- function(p, q) {
    f <- cell.freqs(0, p, q)
    min.diff <- max(-f[c(1, 4)])
    max.diff <- min(f[c(2, 3)])
    return(c(min.diff, max.diff))
}

alpha <- 0.01 # threshold p-value = Type 1 error probability
beta <- 0.1   # Type 2 error probability = 1 - power
p = 0.60  # proportion of population with prior T-cell reactivity (row freqs)
q = 0.55  # secondary attack rate: proportion of exposed individuals who become infected (col freqs)
theta <- log(1/3)  # log odds ratio for association of prior T-cell reactivity with infection

interval.diff <- range.freqs.diff(p, q)
seq.diff <- seq(interval.diff[1], interval.diff[2], by=0.1)

## optimize freq.diff to minimize logoddsratio.err
## freq.diff must be the first argument of logoddsratio.err
freq.diff <- optimize(logoddsratio.err, interval=interval.diff, p=p, q=q, theta=theta)$minimum

freqs <- cell.freqs(freq.diff, p, q)

## column proportions in 2 x 2 table in which col 1 contains exposed uninfected, col 2 contains unexposed uninfected  
u.colfreqs <- c(
    c(freqs[2], freqs[4]) / (freqs[2] + freqs[4]),
    p, 1 - p)

print(matrix(u.colfreqs, ncol=2))

u.freqs <- 0.5 * u.colfreqs

## ratio of odds of T-cell reactivity among exposed uninfected to odds of T-cell reactivity among unexposed uninfected

u.oddsratio <- u.freqs[1] * u.freqs[4] / (u.freqs[2] * u.freqs[3])
u.beta <- log(u.oddsratio)

info.u.beta <- sum(1 / u.freqs)^-1


N <- (qnorm(1 - alpha/2) + qnorm(1 - beta))^2 / (info.u.beta * u.beta^2)

table.2x2 <- round(t(matrix(freqs, ncol=2)), 2)
table.2x2 <- cbind(table.2x2, rowSums(table.2x2))
colnames(table.2x2) <- c("Infected", "Uninfected", "All")
table.2x2 <- rbind(table.2x2, colSums(table.2x2))

rownames(table.2x2) <- c("Prior T-cell crossreactivity",
                         "No prior T-cell cross-reactivity", "All")

cat("\nProportion of population with prior T-cell reactivity (row freqs)", p, "\n")

cat("Secondary attack rate: proportion of exposed individuals who become infected (col props)", q, "\n")

print(table.2x2)

cat("\nProportion of exposed uninfected with prior T cell reactivity:",
    round(freqs[2] / (freqs[2] + freqs[4]), 3), "\n")

cat("\nTo detect odds ratio", exp(theta),
    "for association of prior T cell reactivity with infection, required total sample size for Type I error", alpha,
    "and Type II error", beta, "is", ceiling(N), "\n")
