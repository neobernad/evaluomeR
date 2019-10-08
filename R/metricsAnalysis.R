#' @title Global statistic analysis
#' @name metricsAnalysis
#' @aliases metricsAnalysis
#' @description
#' Global statistic analysis.
#'
#' @param margins See \code{\link{par}}.
#'
#' @return A dataframe.
#'
#' @examples
#' # Using example data from our package
#' data("ontMetrics")
#' cor = metricsAnalysis(ontMetrics, getImages = TRUE, margins = c(1,0,5,11))
#'
metricsAnalysis <- function(data) {

  data <- as.data.frame(assay(data))
  cat("Plotting Min/Max/Sd plot...\n")
  plotMinMax(data)
  cat("Plotting boxplots...\n")
  suppressWarnings(plotBloxplot(data))
  cat("Plotting clustering...\n")
  plotCluster(data)


}

## Private functions ##
plotMinMax <- function(data) {
  # Prepare data for plotting
  # Data matrix without descritive column
  matrix = data.matrix(data[,-1])
  maxs = colMaxs(matrix)
  mins = colMins(matrix)
  means = colMeans(matrix)
  sd = colSds(matrix)

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
          axis.line = element_line(colour = "black",
                                   size = 1, linetype = "solid")
    ) +
    geom_errorbar(aes(ymin=(dataStats.df.t$Max-dataStats.df.t$Sd),
                      ymax=(dataStats.df.t$Max+dataStats.df.t$Sd)), width=.2,
                  position=position_dodge(.9)) +
    geom_errorbar(aes(ymin=(dataStats.df.t$Min-dataStats.df.t$Sd),
                      ymax=(dataStats.df.t$Min+dataStats.df.t$Sd)), width=.2,
                  position=position_dodge(.9)) +
    scale_y_continuous(breaks=seq(min(dataStats.df.t$Min-dataStats.df.t$Sd), # 15 ticks across min - max range
                  max(dataStats.df.t$Max+dataStats.df.t$Sd),
                  round((max(dataStats.df.t$Max)-min(dataStats.df.t$Min)))/15),
                  labels=function(x) sprintf("%.2f", x)) + # Two decimals
    labs(x = "Metrics", y = "Metric value", title = "Min/max/sd values across metrics") +
    guides(fill=TRUE)

    print(p)
}

plotBloxplot <- function(data) {
  num_metrics_plot=10
  data.metrics = data[,-1] # Removing Description column

  metrics_length = length(colnames(data.metrics))
  num_iterations = round(metrics_length/num_metrics_plot)
  if (num_iterations > 0) {
    num_iterations = num_iterations - 1
  }
  for (iteration in 0:num_iterations) {
    i = 1
    rangeStart = (iteration*num_metrics_plot)+1
    rangeEnd = rangeStart+num_metrics_plot
    if (rangeEnd > metrics_length) {
      rangeEnd = metrics_length
    }
    suppressMessages({
      data.melt = melt(data.metrics[,rangeStart:rangeEnd])
    })
    p <- ggplot(data.melt, aes(x=data.melt$variable, y=data.melt$value)) +
      geom_boxplot(
        aes(fill=data.melt$variable),
        outlier.colour = "black",
        outlier.alpha = 0.7,
        outlier.shape = 21
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(x = "Metrics", y="Metric value")
    print(p)
  }
}

plotCluster <- function(data) {
  data.metrics = data[,-1] # Removing Description column
  d <- dist(t(data.metrics), method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward.D2")
  theme_set(theme_bw())
  p <- ggdendrogram(fit, rotate = FALSE, size = 2) + # display dendogram
    labs(title="Metrics dendrogram")
  print(p)
}

