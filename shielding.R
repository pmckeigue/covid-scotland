## tables for shielding report

paste.vectortomatrix <- function(x, y) {
    matrix(paste(x, y), nrow=nrow(x), dimnames=dimnames(x))
}

## table by age group
shieldgroups.byage <- with(cc.all,
                           as.matrix(as.data.frame.matrix(table(shield.group, agegr3))))
shieldgroups.byage.rowpct <- round(100 * shieldgroups.byage /
                                   rowSums(shieldgroups.byage))
shieldgroups.byage.rowpct <- paste0("(", shieldgroups.byage.rowpct, "\\%)")
table.shieldgroups.byage <- paste.vectortomatrix(shieldgroups.byage, shieldgroups.byage.rowpct)


## table by case status: all and severe
table.shieldgroups <- univariate.tabulate("shield.group", outcome="CASE", data=cc.all)

table.shieldgroups.severe <- univariate.tabulate("shield.group", outcome="CASE", data=cc.severe)

## table by date of letter
table.date.letter <- as.data.frame.matrix(
    with(cc.all, table(shield.group, Date.Sent, exclude=NULL)))

table.date.letter <- data.frame(table.date.letter[, 1:4],
                                rowSums(table.date.letter[, 5:ncol(table.date.letter)]))
colnames(table.date.letter) <- format(as.Date(levels(as.factor(cc.all$Date.Sent))[1:5]),
                                      "%d %b %y")
colnames(table.date.letter)[5] <- paste(colnames(table.date.letter)[5], "or later")

table.date.letter <- as.matrix(table.date.letter[-1, ])
date.letter.rowpct <- round(100 * table.date.letter / rowSums(table.date.letter))
date.letter.rowpct <- paste0("(", date.letter.rowpct, "\\%)")

table.date.letter <- paste.vectortomatrix(table.date.letter, date.letter.rowpct)

## calculate estimated proportion of all cases infected after arrival of each batch of letters

severe.infected.date <- table(cc.severe[CASE==1, SPECIMENDATE - 7])
dateby.props <- cumsum(severe.infected.date) / sum(severe.infected.date)
dateby.props <- dateby.props[match(levels(as.factor(cc.all$Date.Sent)), names(dateby.props))]


## table by infection after letter
table.after.letter <- as.data.frame.matrix(
    with(cc.all[CASE==1], table(shield.group, after.letter, exclude=NULL)))
colnames(table.after.letter) <- c("Infection before letter", "Infection after letter")
table.after.letter <- as.matrix(table.after.letter[-1, ])
after.letter.rowpct <- round(100 * table.after.letter / rowSums(table.after.letter))
after.letter.rowpct <- paste0("(", after.letter.rowpct, "\\%)")
table.after.letter <- paste.vectortomatrix(table.after.letter, after.letter.rowpct)



## regression of case status on shielding group and after.letter

table.shielded.allcases <-
    tabulate.freqs.regressions(varnames=c("shield.group","after.letter"),
                               outcome="CASE", data=cc.all)
table.shielded.allcases["after.letter", 1:4] <- NA

## regression of severe case status on shielding group and after.letter

table.shielded.severecases <-
    tabulate.freqs.regressions(varnames=c("shield.group","after.letter"),
                               outcome="CASE", data=cc.severe)
table.shielded.severecases["after.letter", 1:4] <- NA
