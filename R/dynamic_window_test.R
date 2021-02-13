# library(slider)
# install.packages("slider")
#
# n <- dim(PROBES)[1]
# epimutation <- sample(c(0,1), replace=TRUE, size=n)
#
# df <- data.frame(PROBES, MUTATIONS = epimutation)
#
# df <- subset(df, CHR=="11")
#
# start <- sort(df$START, decreasing = FALSE)
#
# left <- slide_index(start, start, ~length(.x), .after=5000, .complete = TRUE)
# rigth <- slide_index(start, start, ~length(.x), .before = 10000, .complete = TRUE)
#
#
# tot <- unlist(left) + unlist(rigth)
#
# tt <- unlist(rigth)
# hist(tt)
# set.seed(123)
# table(tt)
#
# tt <- unlist(left)
# hist(tt)
# set.seed(123)
# table(tt)
#
#
#
#
# df <- data.frame(
#   y = rnorm(100),
#   x = rnorm(100),
#   i = as.Date("2019-08-15") + c(0, 2, 4, 6:102) # <- irregular
# )
#
# # 20 day rolling regression. Current day + 19 days back.
# # Additionally, set `.complete = TRUE` to not compute partial results.
# regr <- slide_index(df, df$i, ~lm(y ~ x, .x), .before = 19, .complete = TRUE)
#
#
# i <- c(2017, 2017, 2018, 2019, 2020, 2020)
# slide_index(i, i, ~.x)
#
#
#
