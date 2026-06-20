#' @title Minimum and maximum metric values plot.
#' @name plotMetricsMinMax
#' @aliases plotMetricsMinMax
#' @description
#' It plots the minimum, maximum and standard deviation
#' values of the metrics in a \code{\link{SummarizedExperiment}} object.
#'
#' @inheritParams stability
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsMinMax(ontMetrics)
#'
plotMetricsMinMax <- function(data) {
  data <- assayAsDF(data)
  # Prepare data for plotting
  # Data matrix without descritive column
  matrix = data.matrix(data[,-1])
  maxs = matrixStats::colMaxs(matrix)
  mins = matrixStats::colMins(matrix)
  means = colMeans(matrix)
  sd = matrixStats::colSds(matrix)

  dataStats = matrix(NA, nrow=5, ncol = length(data[,-1]), byrow=TRUE,
                     dimnames = list(c("Metric", "Min","Max","Mean","Sd"),
                                     c(colnames(data[,-1]))))
  dataStats["Metric",] = colnames(data[,-1])
  dataStats["Min",] = mins
  dataStats["Max",] = maxs
  dataStats["Mean",] = means
  dataStats["Sd",] = sd
  dataStats.df = as.data.frame(dataStats)
  dataStats.df.t = as.data.frame(t(dataStats.df))
  # Factor to numeric conversion
  dataStats.df.t$Min = as.numeric(as.character(dataStats.df.t$Min))
  dataStats.df.t$Max = as.numeric(as.character(dataStats.df.t$Max))
  dataStats.df.t$Mean = as.numeric(as.character(dataStats.df.t$Mean))
  dataStats.df.t$Sd = as.numeric(as.character(dataStats.df.t$Sd))

  ## Plotting
  p <- ggplot(dataStats.df.t, aes(x=dataStats.df.t$Metric)) +
    geom_linerange(aes(ymin=dataStats.df.t$Min,ymax=dataStats.df.t$Max),
                   linetype=2,color="#4E84C4") +
    geom_point(aes(y=dataStats.df.t$Min),size=3,color="#00AFBB") +
    geom_point(aes(y=dataStats.df.t$Max),size=3,color="#FC4E07") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          #axis.text.y = element_blank(),
          text = element_text(size=15),
          axis.line = element_line(colour = "black",
                                   size = 1, linetype = "solid")
    ) +
    geom_errorbar(aes(ymin=(dataStats.df.t$Max-dataStats.df.t$Sd),
                      ymax=(dataStats.df.t$Max+dataStats.df.t$Sd)), width=.2,
                  position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=(dataStats.df.t$Min-dataStats.df.t$Sd),
                      ymax=(dataStats.df.t$Min+dataStats.df.t$Sd)), width=.2,
                  position=position_dodge(.9)) +
    scale_y_continuous(breaks=seq(round(min(dataStats.df.t$Min-dataStats.df.t$Sd)), # 10 ticks across min - max range
                                  round(max(dataStats.df.t$Max+dataStats.df.t$Sd)),
                                  round((max(dataStats.df.t$Max)-min(dataStats.df.t$Min)))/10),
                       labels=function(x) sprintf("%.2f", x)) + # Two decimals
    labs(x = "Metrics", y = "Metric value", title = "Min/max/sd values across metrics") +
    guides(fill=TRUE)

  print(p)
}


#' @title Metric values as a boxplot.
#' @name plotMetricsBoxplot
#' @aliases plotMetricsBoxplot
#' @description
#' It plots the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object as a boxplot.
#'
#' @inheritParams stability
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsBoxplot(ontMetrics)
#'
plotMetricsBoxplot <- function(data) {
  data <- assayAsDF(data)
  num_metrics_plot=20
  data.metrics = data[,-1] # Removing Description column

  metrics_length = length(colnames(data.metrics))
  num_iterations = round(metrics_length/num_metrics_plot)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
    i = 1
    rangeStart = (iteration*num_metrics_plot)+1
    rangeEnd = rangeStart+num_metrics_plot-1
    if (rangeEnd > metrics_length) {
      rangeEnd = metrics_length
    }
    suppressMessages({
      data.melt = melt(data.metrics[,rangeStart:rangeEnd])
    })
    # Melting 1 variable (e.g: data.metrics[,11:11])
    # won't create $variable column in data.melt.
    if (rangeStart == rangeEnd) {
      metricName = data[rangeStart, "Description"]
      data.melt$variable = rep(metricName, length(data.melt$value))
    }
    p <- ggplot(data.melt, aes(x=data.melt$variable, y=data.melt$value)) +
      geom_boxplot(
        #aes(fill=data.melt$variable), # Colors
        outlier.colour = "black",
        outlier.alpha = 0.7,
        outlier.shape = 21,
        show.legend = FALSE
      ) +
      #scale_y_continuous(limits = quantile(data.melt$value, c(0.1, 0.9))) +
      scale_color_grey() +
      theme_bw() +
      theme(
        text = element_text(size=20),
        axis.text.x = element_text(angle = 90)
      ) +
      labs(x = "Metrics", y="Metric value", fill="Metrics")
    # compute lower and upper whiskers
    #ylim1 = boxplot.stats(data.melt$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    #p1 = p + coord_cartesian(ylim = ylim1*1.05)
    print(p)
  }
}

#' @title Metric values clustering.
#' @name plotMetricsCluster
#' @aliases plotMetricsCluster
#' @description
#' It clusters the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object a an hclust dendogram from \code{\link{stats}}. By default distance is measured in 'euclidean'
#' and hclust method is 'ward.D20.
#'
#' @inheritParams stability
#' @param scale Boolean. If true input data is scaled. Default: FALSE.
#' @param k Integer. If not NULL a 'cutree' cut on the cluster is done. Default: NULL
#'
#' @return An hclust object.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsCluster(ontMetrics, scale=TRUE)
#'
plotMetricsCluster <- function(data, scale=FALSE, k=NULL) {
  data <- assayAsDF(data)
  data.metrics = data[,-1] # Removing Description column
  if (isTRUE(scale)) {
    data.metrics = base::scale(data.metrics)
  }
  d <- dist(t(data.metrics), method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward.D2")

  dend <- as.dendrogram(fit)

  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
  plot(dend, xlab = "", sub="", ylab = "Euclidean distance",
       main = paste0("Metrics dendrogram 'ward.D2'"), nodePar = nodePar)

  if (!is.null(k)) {
    dendextend::rect.dendrogram(dend , k=k, border="purple")
  }
  return(fit)
}

#' @title Metric values as violin plot.
#' @name plotMetricsViolin
#' @aliases plotMetricsViolin
#' @description
#' It plots the value of the metrics in a \code{\link{SummarizedExperiment}}
#' object as a violin plot.
#'
#' @inheritParams stability
#' @param nplots Positive integer. Number of metrics per violin plot. Default: 20.
#'
#' @return Nothing.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' plotMetricsViolin(ontMetrics)
#'
plotMetricsViolin <- function(data, nplots=20) {
  data <- assayAsDF(data)
  data.metrics = data[,-1] # Removing Description column
  nplots=20

  metrics_length = length(colnames(data.metrics))
  num_iterations = round(metrics_length/nplots)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
      i = 1
      rangeStart = (iteration*nplots)+1
      rangeEnd = rangeStart+nplots-1
      if (rangeEnd > metrics_length) {
        rangeEnd = metrics_length
      }
    suppressMessages({
      data.melt = melt(data.metrics[,rangeStart:rangeEnd])
    })
    # Melting 1 variable (11:11), won't create $variable column in data.melt.
    if (rangeStart == rangeEnd) {
      metricName = data[rangeStart, "Description"]
      data.melt$variable = rep(metricName, length(data.melt$value))
    }
    p <- ggplot(data.melt, aes(x=data.melt$variable, y=data.melt$value)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      #scale_y_continuous(limits = quantile(data.melt$value, c(0.1, 0.9))) +
      scale_color_grey() +
      theme_bw() +
      theme(
        text = element_text(size=20),
        axis.text.x = element_text(angle = 90)
      ) +
      labs(x = "Metrics", y="Metric value", fill="Metrics")
    # compute lower and upper whiskers
    #ylim1 = boxplot.stats(data.melt$value)$stats[c(1, 5)]
    # scale y limits based on ylim1
    #p1 = p + coord_cartesian(ylim = ylim1*1.05)
    print(p)
  }
}

