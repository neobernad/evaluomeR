#' @title Calculation of Pearson correlation coefficient.
#' @name correlations
#' @aliases correlations
#' @description
#' Calculation of Pearson correlation coefficient between every pair of
#' metrics available in order to quantify their interrelationship degree.
#' The score is in the range [-1,1]. Perfect correlations: -1 (inverse),
#' and 1 (direct).
#'
#' @inheritParams stability
#' @param margins See \code{\link{par}}.
#'
#' @return The Pearson correlation matrix.
#'
#' @examples
#' # Using example data from our package
#' metrics = loadSample("ont-metrics")
#' cor = correlations(data=metrics, getImages = TRUE)
#' # Changing figure margings
#' cor = correlations(data=metrics, getImages = TRUE, margins = c(1,0,5,11))
#'
correlations <- function(data, margins=c(0,10,9,11), getImages=TRUE,
                         label=NULL, path=NULL) {

  if (!is.null(label)) {
    isString(label)
  }

  cur.env <- new.env()

  MatCorr <- cor(data[,2:length(data)])

  if (getImages == TRUE) {
    if (!is.null(path)) {
      path <- checkDirectory(path)
    }
    runMetricsCorrelationIMG(data, MatCorr, margins, label, path)
  }

  return(MatCorr)
}

runMetricsCorrelationIMG <- function(data, correlations, margins, label, path) {
  datos.bruto=NULL; names.metr=NULL; names.index=NULL; k.min=NULL; k.max=NULL;
  datos.bruto = data
  MatCorr = correlations
  names.metr=names(datos.bruto)[-c(1)]  #metric names

  ##########################################################
  ancho=NULL; alto=NULL; ajuste=NULL; margenes=NULL;
  escala=NULL; escalax=NULL; escalal=NULL; escalat=NULL; escalap=NULL;
  ancho=5 #dimension de gr?fica  (pulgadas)
  alto=4 #dimension de gr?fica  (pulgadas)
  escala=0.6 #reescalamiento de texto
  escalax=escala #reescalamiento de ejes
  escalal=0.8 #reescalamiento de etiquetas de ejes
  ajuste=0.5 #ajuste 0=izq, 0.5=centrado(defecto), 1=derecha
  escalat=0.5
  escalap=0.4
  #margenes=c(1,0,5,11)
  margenes=margins
  ## para graficos de correlaciones
  col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                             "cyan", "#007FFF", "blue","#00007F"))
  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                             "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
  col3 <- colorRampPalette(c("red", "white", "blue"))
  col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F",
                             "cyan", "#007FFF", "blue","#00007F"))
  wb <- c("white","black")
  #par(ask = TRUE)
  ##########################################################

  listaFigurasCor=NULL
  contadorFiguras=1
  #Pattern: Correlations_X_metrics
  figurename="Correlations_"

  ##########################################################

  ## circle + colorful number
  corrplot(MatCorr, is.corr = TRUE,
           method = "number", #"circle", "square", "ellipse", "number", "shade", "color", "pie"
           order = "alphabet",  #"original", "AOE", "FPC", "hclust", "alphabet"
           diag = FALSE, type = "lower", #"full", "lower", "upper"
           tl.pos = "n", #pos text labels: lt", "ld", "td", "d", "n"
           tl.cex = 0.8, tl.col = "black", tl.offset = 0.4, tl.srt = 90, #ratio, col, pos & rotation
           col = col1(200), outline=FALSE, #col of symbol & outline of symbol
           cl.pos = "n", #pos col labels: "r", "b", "n"
           number.cex = 0.7,  #ratio to write corr in table
           number.font = 1, #font to write corr in table
           number.digits = 2,#NULL, #number of digits to write corr in table
           win.asp = 1, #ratio of symbol in table
           mar=margenes, #margenes del grafico
           bg = "white", title = "", fg="black",
           addgrid.col = NULL, #col of grid: NULL, NA, "red", ...
           addCoef.col = NULL, #col of coefficients: NULL, "red", ...
           addCoefasPercent = FALSE, #coefficients in %
           add = FALSE)

  corrplot(MatCorr, is.corr = TRUE,
           method = "circle", #"circle", "square", "ellipse", "number", "shade", "color", "pie"
           order = "alphabet",  #"original", "AOE", "FPC", "hclust", "alphabet"
           diag = TRUE, type = "upper", #"full", "lower", "upper"
           tl.pos = "t", #pos text labels: lt", "ld", "td", "d", "n"
           tl.cex = 0.8, tl.col = "black", #tl.offset = 0.4, tl.srt = 90, #ratio, col, pos & rotation
           col = col1(200), outline=FALSE, #col of symbol & outline of symbol
           cl.pos = "r", #pos col labels: "r", "b", "n"
           cl.cex = 0.8, cl.ratio = 0.15, cl.align.text = "r", cl.offset = -0.25,
           win.asp = 1, #ratio of symbol in table
           bg = "white", title = "", fg="black",
           #cl.pos = "n", #pos col labels: "r", "b", "n"
           #mar=margenes, #margenes del grafico
           add = TRUE)

  if (!is.null(path)) {

    if (!is.null(label)) {
      name <- paste(path, label, " - ", figurename, length(names.metr), "_metrics.png", sep="")
    } else {
      name <- paste(path, label, figurename, length(names.metr), "_metrics.png", sep="")
    }

    png(name)
    corrplot(MatCorr, is.corr = TRUE,
             method = "number", #"circle", "square", "ellipse", "number", "shade", "color", "pie"
             order = "alphabet",  #"original", "AOE", "FPC", "hclust", "alphabet"
             diag = FALSE, type = "lower", #"full", "lower", "upper"
             tl.pos = "n", #pos text labels: lt", "ld", "td", "d", "n"
             tl.cex = 0.8, tl.col = "black", tl.offset = 0.4, tl.srt = 90, #ratio, col, pos & rotation
             col = col1(200), outline=FALSE, #col of symbol & outline of symbol
             cl.pos = "n", #pos col labels: "r", "b", "n"
             number.cex = 0.7,  #ratio to write corr in table
             number.font = 1, #font to write corr in table
             number.digits = 2,#NULL, #number of digits to write corr in table
             win.asp = 1, #ratio of symbol in table
             mar=margenes, #margenes del grafico
             bg = "white", title = "", fg="black",
             addgrid.col = NULL, #col of grid: NULL, NA, "red", ...
             addCoef.col = NULL, #col of coefficients: NULL, "red", ...
             addCoefasPercent = FALSE, #coefficients in %
             add = FALSE)

    corrplot(MatCorr, is.corr = TRUE,
             method = "circle", #"circle", "square", "ellipse", "number", "shade", "color", "pie"
             order = "alphabet",  #"original", "AOE", "FPC", "hclust", "alphabet"
             diag = TRUE, type = "upper", #"full", "lower", "upper"
             tl.pos = "t", #pos text labels: lt", "ld", "td", "d", "n"
             tl.cex = 0.8, tl.col = "black", #tl.offset = 0.4, tl.srt = 90, #ratio, col, pos & rotation
             col = col1(200), outline=FALSE, #col of symbol & outline of symbol
             cl.pos = "r", #pos col labels: "r", "b", "n"
             cl.cex = 0.8, cl.ratio = 0.15, cl.align.text = "r", cl.offset = -0.25,
             win.asp = 1, #ratio of symbol in table
             bg = "white", title = "", fg="black",
             #cl.pos = "n", #pos col labels: "r", "b", "n"
             #mar=margenes, #margenes del grafico
             add = TRUE)
    par(new=FALSE)
    dev.off()
  }

}
