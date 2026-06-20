library(evaluomeR)

data("rnaMetrics")

dataFrame <- stability(data=rnaMetrics, k=4, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,4), bs=20, seed=13606, all_metrics = FALSE, getImages = FALSE)
result_stab <- as.data.frame(assay(dataFrame))
# Metric    Mean_stability_k_2  Mean_stability_k_3  Mean_stability_k_4
# [1,] "RIN"     "0.825833333333333" "0.778412698412698" "0.69625"
# [2,] "DegFact" "0.955595238095238" "0.977777777777778" "0.820833333333333"
stopifnot(isTRUE(all.equal(
  as.numeric(result_stab[result_stab[,"Metric"]=="RIN",    "Mean_stability_k_2"]),
  0.825833333333333, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_stab[result_stab[,"Metric"]=="RIN",    "Mean_stability_k_3"]),
  0.778412698412698, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_stab[result_stab[,"Metric"]=="RIN",    "Mean_stability_k_4"]),
  0.69625, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_stab[result_stab[,"Metric"]=="DegFact","Mean_stability_k_2"]),
  0.955595238095238, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_stab[result_stab[,"Metric"]=="DegFact","Mean_stability_k_3"]),
  0.977777777777778, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_stab[result_stab[,"Metric"]=="DegFact","Mean_stability_k_4"]),
  0.820833333333333, tolerance=1e-8)))

dataFrame <- stabilitySet(data=rnaMetrics, k.set=c(2,3,4), bs=20, seed=13606, all_metrics = FALSE, getImages = FALSE)

dataFrame <- quality(data=rnaMetrics, cbi="kmeans", k=3, seed=13606, all_metrics = FALSE, getImages = FALSE)
result_qual3 <- as.data.frame(assay(dataFrame))
# Metric    Cluster_1_SilScore  Cluster_2_SilScore  Cluster_3_SilScore
# [1,] "RIN"     "0.420502645502646" "0.724044583696066" "0.68338517747747"
# [2,] "DegFact" "0.876516605981734" "0.643613928123002" "0.521618857725795"
# Avg_Silhouette_Width Cluster_1_Size Cluster_2_Size Cluster_3_Size
# [1,] "0.627829396038413"  "4"            "4"            "8"
# [2,] "0.737191191352892"  "8"            "5"            "3"
stopifnot(isTRUE(all.equal(
  as.numeric(result_qual3[result_qual3[,"Metric"]=="RIN",    "Avg_Silhouette_Width"]),
  0.627829396038413, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_qual3[result_qual3[,"Metric"]=="DegFact","Avg_Silhouette_Width"]),
  0.737191191352892, tolerance=1e-8)))

dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,4), seed = 20, all_metrics = FALSE, getImages = FALSE)
result_qr2 <- as.data.frame(assay(getDataQualityRange(dataFrame, 2)))
result_qr4 <- as.data.frame(assay(getDataQualityRange(dataFrame, 4)))
# k=2
# [1,] "RIN"     "0.583166775069983" "0.619872562681118" "0.608402004052639"
# [2,] "DegFact" "0.664573423022171" "0.675315791048653" "0.666587617027136"
stopifnot(isTRUE(all.equal(
  as.numeric(result_qr2[result_qr2[,"Metric"]=="RIN",    "Avg_Silhouette_Width"]),
  0.608402004052639, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_qr2[result_qr2[,"Metric"]=="DegFact","Avg_Silhouette_Width"]),
  0.666587617027136, tolerance=1e-8)))
# k=4
stopifnot(isTRUE(all.equal(
  as.numeric(result_qr4[result_qr4[,"Metric"]=="RIN",    "Avg_Silhouette_Width"]),
  0.463905611516569, tolerance=1e-8)))
stopifnot(isTRUE(all.equal(
  as.numeric(result_qr4[result_qr4[,"Metric"]=="DegFact","Avg_Silhouette_Width"]),
  0.634170498361632, tolerance=1e-8)))

dataFrame1 <- qualitySet(data=rnaMetrics, k.set=c(2,3,4), seed=13606, all_metrics = FALSE, getImages = FALSE)


dataFrame <- metricsCorrelations(data=rnaMetrics, getImages = FALSE, margins = c(4,4,11,10))
assay(dataFrame, 1)


dataFrame <- stability(data=rnaMetrics, cbi="kmeans",    k=2, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="clara",     k=2, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="clara_pam", k=2, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="hclust",    k=2, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="pamk",      k=2, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame <- stability(data=rnaMetrics, cbi="pamk_pam",  k=2, bs=100, seed=13606, all_metrics = FALSE, getImages = FALSE)
#dataFrame <- stability(data=rnaMetrics, cbi="rskc", k=2, bs=100, all_metrics = TRUE, L1 = 2, alpha=0, getImages = FALSE)

# Supported CBIs:
evaluomeRSupportedCBI()

dataFrame <- qualityRange(data=rnaMetrics, k.range=c(2,10), seed=13606, all_metrics = FALSE, getImages = FALSE)
dataFrame

#dataFrame <- stabilityRange(data=rnaMetrics, k.range=c(2,8), bs=20, getImages = FALSE)
#assay(dataFrame)
