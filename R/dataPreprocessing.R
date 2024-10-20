#' @title Remove Non-Numeric Columns
#' @description
#' This function removes non-numeric columns from the dataset and returns only numeric columns.
#'
#' @param dataset A data frame from which non-numeric columns are to be removed.
#'
#' @return A data frame containing only numeric columns. Also prints non-numeric columns found.

removeNonNumericColumns <- function(dataset) {
  message("Removing non-numeric columns...")
  # Get only numeric columns
  numeric_cols = dataset[sapply(dataset, is.numeric)]
  # Get non-numeric columns
  non_numeric_cols <- names(dataset)[!sapply(dataset, is.numeric) & names(dataset) != "Description"]
  if (length(non_numeric_cols) > 0) {
    message("\tNon-numeric columns found:")
    message(paste0("\t\t", paste(non_numeric_cols, collapse = ", ")))
  } else {
    # 'Description' is non-numeric but we need it to know individual labels
    message("\tColumns are numeric")
  }

  # Get only numeric columns
  numeric_dataset <- dataset[sapply(dataset, is.numeric)]
  return(numeric_dataset)
}

#' @title Clean Dataset
#' @description
#' Cleans the dataset by removing non-numeric columns and those with perfect correlations.
#'
#' @param dataset A data frame containing the dataset to be cleaned.
#' @param correlation_threshold A numeric value between -1 and 1 that specifies the correlation threshold
#'                              for removing perfectly correlated columns.
#'
#' @return A list containing the cleaned dataset (\code{dataset}) with the description column
#' and the correlation matrix R (\code{R}).
#'
#' @export
cleanDataset <- function(dataset, correlation_threshold = 1) {
  if (correlation_threshold >= 1 || correlation_threshold < -1) {
    message("Error: 'correlation_threshold' must be in range [-1,1].")
    return(list(dataset = NULL,
                R = NULL))
  }
  if (ncol(dataset) < 2) {
    message("Error: Not enough columns for cleaning up (ncol < 2).")
    return(list(dataset = NULL,
                R = NULL))
  }


  message("Preprocessing dataset")
  # Get description column
  description_col = dataset[,"Description", drop=FALSE]

  numeric_dataset = removeNonNumericColumns(dataset)

  r_remove_correlations = removeCorrelations(numeric_dataset, correlation_threshold)
  R_cleaned = r_remove_correlations$R

  # Extract the remaining column names
  remaining_colnames = r_remove_correlations$remaining_colnames
  numeric_dataset <- numeric_dataset[, remaining_colnames, drop = FALSE]
  dataset_w_description <- cbind(description_col, numeric_dataset)


  return(list(dataset = dataset_w_description,
              R = R_cleaned))
}

#' @title Remove Correlations
#' @description
#' Identifies and removes  correlated columns from the numeric dataset, according to `correlation_threshold`.
#'
#' @param numeric_dataset A data frame consisting of only numeric columns.
#' @param correlation_threshold A numeric threshold for identifying perfect correlations.
#'
#' @return A list containing the cleaned dataset and the remaining column names, along with the correlation matrix R.
#'
removeCorrelations <- function(numeric_dataset, correlation_threshold = 1) {
  message("Removing correlations...")

  # Keep track of the initial dataset
  dataset_to_check <- numeric_dataset

  # Initialize R as the correlation matrix
  R <- cor(dataset_to_check)

  repeat {
    # Find perfect correlations
    perfect_corr_indices <- which((R >= correlation_threshold | R <= correlation_threshold*-1)
                                  & row(R) != col(R), arr.ind = TRUE)

    if (nrow(perfect_corr_indices) == 0) {
      break  # No more perfect correlations found
    }

    # Get the first pair of perfectly correlated columns
    i = perfect_corr_indices[1, 1]  # Row index (not using this, but just for the record)
    j <- perfect_corr_indices[1, 2]  # Column index of the column to remove

    message("\tRemoving correlated column: ", colnames(dataset_to_check)[j])

    # Remove the j-th column from the dataset
    dataset_to_check <- dataset_to_check[, -j, drop = FALSE]

    # Recalculate the correlation matrix with the updated dataset
    R <- cor(dataset_to_check)
  }

  return(list(cleaned_dataset = dataset_to_check,
              remaining_colnames = colnames(dataset_to_check),
              R=R))
}


#' @title PCA Suitability
#' @description
#' Performs Bartlett's test and KMO test to determine the suitability of PCA on the given correlation matrix R.
#'
#' @param R A correlation matrix of the dataset.
#' @param sig_level A numeric significance level for Bartlett's test, default is 0.05.
#'
#' @return A list indicating if PCA is suitable, along with the results of Bartlett's and KMO tests.
#'
#' @export
PCASuitability <- function(R, sig_level = 0.05) {
  message("Checking PCA suitability...")

  bartlett_test = psych::cortest.bartlett(R, n = 100)
  is_suitable = FALSE


  if (is.na(bartlett_test$p.value)) {
    message("\tPCA is not suitable. Bartlett's test produced NA for p-value.")
    return (list(pca_suitable=is_suitable,
                 bartlett.test=NULL,
                 kmo=NULL))
  }

  # Check p-value for Bartlett's test according to sig_level
  is_suitable = bartlett_test$p.value < sig_level

  if (!is_suitable) {
    message("\tPCA is not suitable. Bartlett's test p-value: ", bartlett_test$p.value)
    return (list(pca_suitable=is_suitable,
                 bartlett.test=bartlett_test,
                 kmo=NULL))
  }
  # Else: Suitable, perform KMO test

  kmo_result = psych::KMO(R)
  kmo_value = kmo_result$MSA

  # Interpret KMO value
  kmo_interpretation <- ""
  if (kmo_value < 0.5) {
    kmo_interpretation = "unacceptable"
  } else if (kmo_value < 0.6) {
    kmo_interpretation = "miserable"
  } else if (kmo_value < 0.7) {
    kmo_interpretation = "mediocre"
  } else if (kmo_value < 0.8) {
    kmo_interpretation = "middling"
  } else if (kmo_value < 0.9) {
    kmo_interpretation = "meritorious"
  } else {
    kmo_interpretation = "marvelous"
  }

  message("\tPCA is suitable.")
  message("\tBartlett's test p-value: ", bartlett_test$p.value)
  message("\tKMO value: ", kmo_value, " - ", kmo_interpretation)

  return (list(pca_suitable=is_suitable,
               bartlett.test=bartlett_test,
               kmo=kmo_result))
}

#' @title Perform PCA
#' @description
#' Executes PCA on the provided dataset and summarizes the results.
#'
#' @param dataset A data frame to be analyzed via PCA.
#' @param ncp An integer specifying the number of principal components to retain, default is 5.
#' @param scale A boolean indicating whether to scale the data, default is TRUE.
#'
#' @return A list containing the PCA results, summary of eigenvalues, contributions, coordinates, and more.
#'
#' @export
performPCA <- function(dataset, ncp = 5, scale = TRUE, visualize = FALSE) {

  description_col = dataset[[1]]
  if ("Description" %in% colnames(dataset)) {
    description_col = dataset[,"Description", drop=FALSE]
    dataset <- dataset[, !colnames(dataset) %in% "Description"]
  }

  pca_result <- FactoMineR::PCA(dataset, scale.unit = scale, ncp = ncp, graph = FALSE)
  pca_dimdesc = FactoMineR::dimdesc(pca_result, axes=c(1,2))

  summary(pca_result)

  # Mimic what summary(pca_result) would do
  eigenvalues <- pca_result$eig
  var_contributions <- pca_result$var$contrib
  var_coordinates <- pca_result$var$coord
  var_cos2 <- pca_result$var$cos2
  ind_contributions <- pca_result$ind$contrib
  ind_coordinates <- pca_result$ind$coord
  ind_cos2 <- pca_result$ind$cos2

  pca_summary <- list(
    eigenvalues = eigenvalues,
    var_contributions = var_contributions,
    var_coordinates = var_coordinates,
    var_cos2 = var_cos2,
    ind_contributions = ind_contributions,
    ind_coordinates = ind_coordinates,
    ind_cos2 = ind_cos2
  )

  dataset_ncp = as.data.frame(pca_result$ind$coord[, 1:nFactors])
  dataset_ncp <- cbind(description_col, dataset_ncp)

  return (list(dataset_ncp = dataset_ncp,
               pca = pca_result,
               summary = pca_summary,
               dimdesc = pca_dimdesc
  ))

}


#' @title Scree Plot for PCA
#' @description
#' Creates a scree plot to visualize the explained variance by principal components from a PCA result.
#'
#' @param pca_result The result of a PCA analysis (object returned from `performPCA`).
#' @param title An optional title for the scree plot. If NULL, a default title will be used.
#'
#' @return A ggplot object representing the scree plot.
#'
#' @export
plotPCA_fviz_screeplot <- function(pca_result, title=NULL) {
  if (is.null(title)) {
    title = "Scree Plot: Explained Variance by Principal Components"
  }

  factoextra::fviz_screeplot(pca_result,
                             addlabels = TRUE,
                             ylim = c(0, 100),
                             barcolor = "steelblue",
                             barfill = "lightblue",
                             ggtheme = ggplot2::theme_minimal(),
                             xlab = "Principal Components",
                             ylab = "Percentage of Variance (%)",
                             title = title,
                             labelsize = 4,
                             font.x = 14,
                             font.y = 14,
                             font.title = 16)

}

#' @title PCA Biplot
#' @description
#' Creates a biplot to visualize both the individuals and variables in a PCA result.
#'
#' @param pca_result The result of a PCA analysis (object returned from `performPCA`).
#' @param title An optional title for the biplot. If NULL, a default title will be used.
#'
#' @return A ggplot object representing the PCA biplot.
#'
#' @export
plotPCA_fviz_biplot <- function(pca_result, title=NULL) {
  if (is.null(title)) {
    title = "PCA Biplot"
  }

  factoextra::fviz_pca_biplot(pca_result,
                              repel = TRUE,
                              label = "var",
                              palette = "Dark2",
                              col.var = "steelblue",
                              col.ind = "orange",
                              pointshape = 21,
                              pointsize = 3,
                              arrowsize = 1,
                              fill.ind = "lightyellow",
                              fill.var = "lightgreen",
                              ggtheme = ggplot2::theme_minimal(),
                              title = "PCA Biplot",
                              font.x = 14,
                              font.y = 14,
                              font.title = 16)



}

#' @title Determine Number of Factors
#' @description
#' This function determines the optimal number of factors to extract from a dataset
#' using the eigenvalue-based scree test and parallel analysis.
#'
#' @param dataset A data frame from which factors are to be extracted.
#'                The function will ignore the "Description" column if it exists.
#'
#' @return An integer representing the number of factors to retain based on the
#'         scree test results.
#'
#' @export
determineNumberOfFactors <- function(dataset) {
  if ("Description" %in% colnames(dataset)) {
    dataset <- dataset[, !colnames(dataset) %in% "Description"]
  }

  eigenvalues = nFactors::eigenComputes(dataset)
  aparallel = nFactors::eigenBootParallel(dataset, quantile=0.95)$quantile
  results = nFactors::nScree(eigenvalues, aparallel = aparallel)
  num_components = results$Components[4]$nkaiser
  return (num_components)

}
