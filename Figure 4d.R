# This code reads the data created by "PO_example_simulations.ipynb"

add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

LCnos <- c(1,10,54,45,7,77,103,20,13,5,28,3)
LCs   <- paste("LC", LCnos, sep="")

LCcols <- ifelse(LCnos < 20, "gray40", "deeppink3")
names(LCcols) <- LCs

n <- length(LCnos)
L <- 9.5

ymax <- c(2,2,2,2)
ymin <- c(-4,-4,-4,-6)

dy <- 2

plot(0,0,
     type="n",
     xlim=c(0, n*L),
     ylim=c(0, 50),
     bty="n", yaxt="n", xaxt="n",
     xlab="Evolutionary time",
     ylab="Cell density")
axis(side=1, at=seq(0, n*L, L), labels=c(LCnos, NA))

for(q in 1:4){
  w <- 5 - q
  size_cutoff <- 10^(ymin[w])

  fy <- function(y){
    (10 + dy) * (q - 1) + 10 * (y - ymin[w]) / (ymax[w] - ymin[w])
  }

  invasions <- lapply(1:n, function(i){
    read.table(
      paste(c("Data/Figure 4/example_invasion_", w, "_", i, ".csv"), collapse=""),
      header=TRUE, stringsAsFactors=FALSE, sep=","
    )
  })

  t <- seq(ymin[w], ymax[w], 2)
  axis(side=2, at=sapply(t, fy), labels=10^t, las=2)

  for(i in 1:n){
    d <- subset(invasions[[i]], t >= 0.01)

    rect(L*(i-1), fy(ymin[w]), L*i, fy(ymax[w]),
         border=NA, col=add.alpha("gray", 0.15))

    for(l in LCs){
      if(l %in% colnames(d)){

        if(min(d[, l]) < size_cutoff){
          j  <- min(which(d[, l] < size_cutoff))
          d1 <- d[d[, l] >= size_cutoff, , drop=FALSE]

          points(L*(i-1) + log10(d1$t) + 2,
                 fy(log10(d1[, l])),
                 type="l", lwd=2, col=LCcols[l])

          points(L*(i-1) + log10(d$t[j]) + 2,
                 fy(log10(d[j, l])),
                 pch=4, col=LCcols[l], lwd=2)

        } else {
          d1 <- d[d[, l] >= size_cutoff, , drop=FALSE]

          points(L*(i-1) + log10(d1$t) + 2,
                 fy(log10(d1[, l])),
                 type="l", lwd=2, col=LCcols[l])
        }

        if(l == LCs[i]){
          points(L*(i-1),
                 fy(log10(d1[1, l])),
                 pch=21, col=LCcols[l], bg=LCcols[l], lwd=2)
        }
      }
    }
  }
}
