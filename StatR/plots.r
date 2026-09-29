#' BoxPlot
#'
#' Create a boxplot or comparative boxplots. It is an adaptation of the \link{boxplot} function to use type 6 quartiles.
#'
#' @param x for specifying data from which the boxplots are to be produced. Either a numeric vector, or a single list containing such vectors.
#' @param ... For the formula method, named arguments to be passed to the default method.
#'
#' @examples
#' library(IntroBioStats)
#' BoxPlot(rnorm(20))
#' BoxPlot(rnorm(20),rnorm(20,1.2,2.5))
#' ## y is the response, x is the group variable
#' y=c(rnorm(20),rnorm(10,1.2,2.5))
#' x=c(rep("Sample 1",20),rep("Sample 2",10))
#' ppnorm(y~x)
#' @export
BoxPlot <- function(x, ...)  UseMethod("BoxPlot")
#' @export
BoxPlot.default <- function (x, ..., range = 1.5, width = NULL, varwidth = FALSE,
                             notch = FALSE, outline = TRUE, names, plot = TRUE, border = par("fg"),
                             col = NULL, log = "", pars = list(boxwex = 0.8, staplewex = 0.5,
                                                               outwex = 0.5), horizontal = FALSE, add = FALSE, at = NULL)
{
  args <- list(x, ...)
  namedargs <- if (!is.null(attributes(args)$names))
    attributes(args)$names != ""
  else rep_len(FALSE, length(args))
  groups <- if (is.list(x))
    x
  else args[!namedargs]
  if (0L == (n <- length(groups)))
    stop("invalid first argument")
  if (length(class(groups)))
    groups <- unclass(groups)
  if (!missing(names))
    attr(groups, "names") <- names
  else {
    if (is.null(attr(groups, "names")))
      attr(groups, "names") <- 1L:n
    names <- attr(groups, "names")
  }
  cls <- sapply(groups, function(x) class(x)[1L])
  cl <- if (all(cls == cls[1L]))
    cls[1L]
  else NULL
  for (i in 1L:n) groups[i] <- list(BoxPlot.stats(unclass(groups[[i]]),
                                                  range))
  stats <- matrix(0, nrow = 5L, ncol = n)
  conf <- matrix(0, nrow = 2L, ncol = n)
  ng <- out <- group <- numeric(0L)
  ct <- 1
  for (i in groups) {
    stats[, ct] <- i$stats
    conf[, ct] <- i$conf
    ng <- c(ng, i$n)
    if ((lo <- length(i$out))) {
      out <- c(out, i$out)
      group <- c(group, rep.int(ct, lo))
    }
    ct <- ct + 1
  }
  if (length(cl) && cl != "numeric")
    oldClass(stats) <- cl
  z <- list(stats = stats, n = ng, conf = conf, out = out,
            group = group, names = names)
  if (plot) {
    if (is.null(pars$boxfill) && is.null(args$boxfill))
      pars$boxfill <- col
    do.call("bxp", c(list(z, notch = notch, width = width,
                          varwidth = varwidth, log = log, border = border,
                          pars = pars, outline = outline, horizontal = horizontal,
                          add = add, at = at), args[namedargs]))
    invisible(z)
  }
  else z
}
#' @export
BoxPlot.stats <- function (x,...,coef = 1.5, do.conf = TRUE, do.out = TRUE)
{
  if (coef < 0)
    stop("'coef' must not be negative")
  nna <- !is.na(x)
  n <- sum(nna)
  stats <- stats::quantile(x, na.rm = TRUE,type=6)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0)
    do.out <- FALSE
  else {
    out <- if (!is.na(iqr)) {
      x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef *
                                            iqr)
    }
    else !is.finite(x)
    if (any(out[nna], na.rm = TRUE))
      stats[c(1, 5)] <- range(x[!out], na.rm = TRUE)
  }
  conf <- if (do.conf)
    stats[3L] + c(-1.58, 1.58) * iqr/sqrt(n)
  list(stats = stats, n = n, conf = conf, out = if (do.out) x[out &
                                                                nna] else numeric())
}
#' @export
BoxPlot.formula <-
  function(formula, data = NULL, ..., subset, na.action = NULL,
           drop = FALSE, sep = ".", lex.order = FALSE)
  {
    if(missing(formula) || (length(formula) != 3L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m$... <- m$drop <- m$sep <- m$lex.order <- NULL

    m$na.action <- na.action # force use of default for this method
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")

    BoxPlot(split(mf[[response]], mf[-response],drop = drop, sep = sep),...)
  }

#' Normal probability plot
#'
#' Create a normal probablity plot or overlay normal probability plots
#'
#' @param x for specifying data from which the normal probablity plots are to be produced. Either a numeric vector, or a single list containing such vectors.
#' @param ... For the formula method, named arguments to be passed to the default method.
#'
#'#'
#' @examples
#' library(IntroBioStats)
#' ppnorm(rnorm(20))
#' ppnorm(rnorm(20),rnorm(20,1.2,2.5))
#' ## y is the response, x is the group variable
#' y=c(rnorm(20),rnorm(10,1.2,2.5))
#' x=c(rep("Sample 1",20),rep("Sample 2",10))
#' ppnorm(y~x)
#'
#' @export
ppnorm <- function(x, ...)  UseMethod("ppnorm")
#' @export
ppnorm.default <-function(x,...,names,main,xlab,ylab,french)
{

  args <- list(x, ...)
  namedargs <- if (!is.null(attributes(args)$names))
    attributes(args)$names != ""
  else rep_len(FALSE, length(args))
  groups <- if (is.list(x)) x else args[!namedargs]
  if (0L == (n <- length(groups)))
    stop("invalid first argument")
  if (length(class(groups)))
    groups <- unclass(groups)
  if (!missing(names))
    attr(groups, "names") <- names
  else {
    if (is.null(attr(groups, "names")))
    {
      na <-numeric()
      for (i in 1L:n) { na[i]<-paste("sample",i,sep=" ") }
      attr(groups, "names") <- na
    }
    names <- attr(groups, "names")
  }


  cls <- sapply(groups, function(x) class(x)[1L])
  cl <- if (all(cls == cls[1L]))
    cls[1L]
  else NULL
  for (i in 1L:n) groups[i] <- list(ppnorm.points(unclass(groups[[i]])))


  if (missing(french)) french <- FALSE
  if (!french)
  {
    if (missing(main)) main <- c("Normal Probability Plot")
    if (missing(xlab)) xlab <- c("Sample Quantiles")
    if (missing(ylab)) ylab <- c("Percentage")
  } else
  {
    if (missing(main)) main <- c("Diagramme \u{E0} \u{E9}chelle fonctionnelle normale")
    if (missing(xlab)) xlab <- c("Quantiles empiriques")
    if (missing(ylab)) ylab <- c("Pourcentage")

  }





  xlim <- ylim <- c(numeric(0),numeric(0))

  for (i in groups) {
    xlim <- c(min(i$x,xlim),  max(i$x,xlim))
    ylim <- c(min(i$y,ylim),  max(i$y,ylim))

  }
  ylim[1] <- min(-3.1,ylim[1])
  ylim[2] <- max(3.1,ylim[2])
  xrange=xlim[2]-xlim[1]
  xlim[1]=xlim[1]-xrange/5
  xlim[2]=xlim[2]+xrange/5


  ct<-1
  for (i in groups) {
    do.call("ppnorm.plot",list(i$x,...,lty=ct,pch=ct,xlim=xlim,ylim=ylim,main=main,ylab=ylab,xlab=xlab,french=french))
    ct<-ct+1
    par(new=TRUE)
  }
  par(new=FALSE)
  if (n>1) legend("bottomright",legend=names,pch=1:n,lty=1:n,bg='white',lwd=1.5)

}

#' @export
ppnorm.points<- function(x,...)
{
  nna <- !is.na(x)
  n <- sum(nna)
  stats <-c(mean(x[nna]),sd(x[nna]))
  yPoints<-qnorm(ppoints(n))
  xPoints<-sort(x)

  list(x=xPoints, y=yPoints)
}


#' @export
ppnorm.plot <- function(x,...,xlim,ylim,lty,pch,main,ylab,xlab,french)
{
  args <- list(x, ...)
  namedargs <- if (!is.null(attributes(args)$names))
    attributes(args)$names != ""
  else rep(FALSE, length(args))


  xPoints <- sort(x[!is.na(x)])
  n <- length(xPoints)
  yPoints <- qnorm(ppoints(n))

  xrange=max(x)-min(x)
  if (missing(xlim)) xlim <- c(min(x)-xrange/5,max(x)+xrange/5)
  if (missing(ylim)) ylim <- c(min(min(yPoints)[1],-3.1),max(max(yPoints)[1],3.1))
  if (missing(pch)) pch <- 1
  if (missing(lty)) lty <- 1
  if (missing(french)) french <- FALSE
  if (!french)
  {
    if (missing(main)) main <- c("Normal Probability Plot")
    if (missing(xlab)) xlab <- c("Sample Quantiles")
    if (missing(ylab)) ylab <- c("Percentage")
  } else
  {
    if (missing(main)) main <- c("Diagramme \u{E0} \u{E9}chelle fonctionnelle normale")
    if (missing(xlab)) xlab <- c("Quantiles empiriques")
    if (missing(ylab)) ylab <- c("Pourcentage")

  }





  do.call("plot",c(list(xPoints,yPoints,ylab=ylab,xlab=xlab,yaxt="n",main=main,
                        ylim=ylim,xlim=xlim,pch=pch),args[namedargs]))
  abline(-mean(x)/sd(x),1/sd(x),lty=lty,lwd=2)
  p <- c(.001,.01,.05,.1,.25,.5,.75, .9,.95,.99,.999)
  q <- qnorm(p)
  abline(h=q,lty=5)
  #   axis(2, at=q,labels=p*100, las=2)


  if (!french)
  {
    axis(2, at=q,labels=p*100, las=2)
  }
  else {
    axis(2, at=q,labels=formatC(p*100, decimal.mark=",", big.mark=" ", digits=1, format = "f"), las=2)
  }



}

#' @export
ppnorm.formula <-
  function(formula, data = NULL, ..., subset, na.action = NULL,
           drop = FALSE, sep = ".", lex.order = FALSE)
  {
    if(missing(formula) || (length(formula) != 3L))
      stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m$... <- m$drop <- m$sep <- m$lex.order <- NULL


    m$na.action <- na.action # force use of default for this method
    ## need stats:: for non-standard evaluation
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")


    ppnorm(split(mf[[response]], mf[-response],drop = drop, sep = sep),...)
  }











