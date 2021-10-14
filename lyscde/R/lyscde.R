
diffexp <- function(data1, data2, plotting=FALSE) {
  res = filter_rows(data1, data2)
  res = diff_base(res[["data1"]], res[["data2"]])
  res = normalize_statistic(res, plotting=plotting)
  return(res)
}

diff_base <- function(data1, data2){
  tt = list()
  gmean = list()

  data = log2(1+cbind(data1, data2))
  statistic = apply(data,1,function(x) {t.test(x[1:ncol(data1)],x[(ncol(data1)+1):ncol(data)])$statistic})
  cmean = (log2(1+rowMeans(data1))+log2(1+rowMeans(data2)))/2

  return(data.frame(statistic, mean=cmean))
}

normalize_statistic <-function(data, plotting=FALSE){
  x = data[,"mean"]
  y = result[,"statistic"]
  names(x) = rownames(data)
  names(y) = rownames(data)
  oo <- order(x)
  x = x[oo]
  y = y[oo]

  ww = which(!is.na(y))
  x = x[ww]
  y = y[ww]

  lw1 <- loess(y ~ x)

  if(plotting){
    plot(x, y, pch=20, cex=0.5, cex.axis=2, cex.lab=2, xlab="average expression", ylab="diff statistic")
    lines(x, lw1$fitted, col="red", lwd=3)
  }

  normt = y-lw1$fitted
  rr = rollmean(abs(normt), k=100)

  wpos = which(normt > 0)
  wneg = which(normt < 0)

  mpos = loess(normt[wpos] ~ wpos, span=0.2)$fitted
  mneg = loess(normt[wneg] ~ wneg, span=0.2)$fitted

  offset = 1
  nnpos = normt[wpos][offset:length(mpos)]/mpos[offset:length(mpos)]
  nnneg = -normt[wneg][offset:length(mneg)]/mneg[offset:length(mneg)]

  xpos = wpos[offset:length(mpos)]
  xneg = wneg[offset:length(mneg)]

  if(plotting){
    plot(normt, pch=20, cex=0.2, ylim=c(-10,10))
    abline(h=0, col="black", lwd=4)
    lines(wpos, mpos, col="green", lwd=3)
    lines(wneg, mneg, col="red", lwd=3)
  }

  if(plotting){
    plot(1:length(nnpos), nnpos, pch=20, cex=0.1, col="blue", ylim=c(-7, 7))
    points(1:length(xneg), nnneg, pch=20, cex=0.2, col="red")

    ww1 = which(1-pnorm(nnpos) < 0.05)
    points(ww1, nnpos[ww1], col="magenta", pch=20, cex=0.5)

    ww2 = which(pnorm(nnneg) < 0.05)
    points(ww2, nnneg[ww2], col="magenta", pch=20, cex=0.5)
    abline(h=0, lwd=4)
  }

  pneg = pnorm(nnneg)
  ppos = 1-pnorm(nnpos)

  rrr = c(nnpos, nnneg)
  rnn = names(rrr)
  pvals = c(ppos, pneg)
  fdr = p.adjust(pvals, method="fdr")

  diffexp = data.frame(rnn, rrr, pvals, fdr)
  colnames(diffexp) = c("gene", "t", "pval", "fdr")
  diffexp = diffexp[order(diffexp[,"t"]),]

  return(diffexp)
}

filter_rows <- function(data1, data2, min_detect=0.005){
  inter = intersect(rownames(data1), rownames(data2))
  data = cbind(data1[inter,], data2[inter,])
  pass = list()

  pass = which(rowSums(data > 0)/ncol(data) >= min_detect)
  print(length(pass))
  result = list()
  result[["data1"]] = data.frame(data1[inter,][pass,])
  result[["data2"]] = data.frame(data2[inter,][pass,])
  return(result)
}

plotcounts <- function(gene, data1, data2, l1="data1", l2 = "data2"){
  print("hello")
}
