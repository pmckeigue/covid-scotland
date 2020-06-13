
## calculate propensity score for proton pump using hospitalized not severe
## select columns containing BNF subpara codes

subpara.colnames <- grep("^subpara\\.", colnames(cc.hosp), value=TRUE)
subpara.colnames <- grep("Proton pump", subpara.colnames, invert=TRUE, value=TRUE) 
train.cols <- match(c("y.protonpump", "AGE", "sex", subpara.colnames), colnames(cc.hosp))

train.data <- subset(cc.hosp, select=train.cols)

colnames(train.data) <-
    gsub("(subpara\\.[0-9]+)(.+)", "\\1", colnames(train.data))                

## select subparas with frequency > 0.1

subparanames <- grep("^subpara\\.", colnames(train.data), value=TRUE)
subparanames.freqs <- apply(subset(train.data, select=subparanames), 2,
                           function(x) sum(as.integer(x))) / nrow(train.data)

keep <- subparanames.freqs > 0.05
subparanames.kept <- subparanames[keep]

lower.formula <- as.formula("y.protonpump ~ AGE + sex")
upper.formula <- as.formula(paste("y.protonpump ~ AGE + sex +",
                            paste(subparanames.kept, collapse=" + ")))

full.model <- glm(formula=upper.formula, data=train.data, family="binomial")

cat("Stepwise regression to construct propensity score for proton pump ...")
propensity.stepwise.full <- step(full.model,
                       scope=list(lower=lower.formula, upper=upper.formula),
                      direction="both", method="approximate", trace=-1)
cat("done\n")

## get coeffs, dropping age and sex
coeffs.full <- summary(propensity.stepwise.full)$coefficients[, 1][-(1:3)]
rm(propensity.stepwise.full) 

names(coeffs.full) <- gsub("(.+)1$", "\\1", names(coeffs.full))
names(coeffs.full) <- subpara.colnames[match(names(coeffs.full),
                                             gsub("(subpara\\.[0-9]+)(.+)", "\\1",
                                                  subpara.colnames))]

newdata.x <- subset(cc.severe, select=match(names(coeffs.full), colnames(cc.severe)))
newdata.x <- as.matrix(newdata.x)
## convert from character to numeric
newdata.x <- matrix(as.numeric(newdata.x), nrow=nrow(newdata.x))
propensity <- newdata.x %*% matrix(coeffs.full, ncol=1)

coeffs.propensity <- data.frame(subpara=gsub("subpara\\.", "", names(coeffs.full)),
                                coefficient=as.numeric(coeffs.full))

rm(newdata.x)
rm(full.model)
rm(train.data)
