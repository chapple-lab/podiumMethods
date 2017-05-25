pval <-
function(X, classlabel, teststat) {

    n1 <- rowSums(!is.na(X[,classlabel == 0]))
    n2 <- rowSums(!is.na(X[,classlabel == 1]))
    A <- apply(X[,classlabel == 0], 1, sd, na.rm=TRUE)^2/n1 ## sd(t(X[,classlabel == 0]), na.rm = TRUE)^2/n1
    B <- apply(X[,classlabel == 1], 1, sd, na.rm=TRUE)^2/n2 ## sd(t(X[,classlabel == 1]), na.rm = TRUE)^2/n2
    df <- (A+B)^2/(A^2/(n1-1)+B^2/(n2-1))

    pvalue <- 2 * (1 - pt(abs(teststat), df))
    invisible(pvalue)
}
