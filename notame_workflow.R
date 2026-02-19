# 1. Generate a peak abundance matrix from aligned metabolomics data (acquired
#    from the software handling the peak detection and alignment) and save it
#    in an Excel format or as a tab-separated text or CSV file.

# 2. Install renv, initialize, activate and restore the project environment.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Activate and restore environment from renv.lock
renv::activate()
renv::restore()

# 3. Load notame and other necessary packages, set paths, generate a folder for
#    figures.

library(notame)
library(notameViz)
library(notameStats)
library(dplyr)
library(tidyr)
library(pcaMethods)
library(Rtsne)
library(missForest)
library(igraph)
library(lmerTest)
library(PK)
library(MuMIn)
library(MUVR2)

# Create the path for output data
dir.create(file.path("data", "figures"))

# Set path for all data
ppath <- file.path("data")

# 4. Choose from the following alternatives based on your data.
# a) Load the Excel data containing all modes into R environment and create the
#    SummarizedExperiment data containers.
se <- import_from_excel(
  file = system.file("extdata", "toy_notame_set.xlsx", package = "notame"),
  sheet = 1,
  split_by = "Mode"
)

# b) if the modes (in this example, four modes named hilic_neg.xlsx, etc.) are
#    in separate spreadsheets, load them separately while ensuring they have the
#    same run sequence (the same number and order of samples in the columns). It
#    is recommended to merge them as well at this point into one
#    SummarizedExperiment object.

modesList <- c("HILIC_neg", "HILIC_pos", "RP_neg", "RP_pos")
modes <- list()
for (mode in modesList) {
  # Read single mode, set it in the list
  modes[[mode]] <- import_from_excel(
    file = file.path(ppath, mode, ".xlsx"),
    name = mode
  )
}
se <- merge_notame_sets(object = modes)

# c) Explore and test the package with toy_notame_set.
data(toy_notame_set)
# And then assign it to data (remove the comment #)
se <- toy_notame_set

# 5. Classify the data as (metabolite) abundances. Create any necessary missing
#    columns for the pheno and feature data, clean the object, and split the
#    object by mode using fix_object.

names(assays(se)) <- "abundances"
modes <- fix_object(object = se, split_data = TRUE, assay.type = "abundances")

# 6. Ensure that the objects contain all the data in the correct form.

names(modes)
sapply(modes, class)

# 7. Take several processor cores into use. Leave two cores for other purposes.

library(BiocParallel)
register(SnowParam(workers = snowWorkers() - 2))

# 8. Perform the contaminant flagging, quality check, and drift correction on
#    each mode separately within a loop and create visualisations of the
#    various phases of the process in the figures folder of the working
#    directory.

# Initialize empty list for processed objects
processed <- list()
for (i in seq_along(modes)) {
  name <- names(modes)[i]
  mode <- modes[[i]]
  mode <- mark_nas(mode, value = 0) # Set all zero abundances to NA
  mode <- flag_detection(mode, group = "Group") # Flag features with low detection rate
  mode <- flag_contaminants(
    mode,
    blank_col = "QC",
    blank_label = "Blank",
    blank_type = "mean",
    sample_type = "max",
    flag_thresh = 5
  ) # Flag contaminant features
  mode <- mode[, mode$QC != "Blank"] # Remove blanks
  save_QC_plots(
    mode,
    prefix = file.path(ppath, "figures", paste0(name, "_ORIG")),
    perplexity = 5,
    group = "Group",
    time = "Time",
    id = "Subject_ID",
    color = "Group"
  ) # Visualize data before drift correction
  corrected <- correct_drift(mode) # Correct drift
  save_QC_plots(
    corrected,
    prefix = file.path(ppath, "figures", paste0(name, "_DRIFT")),
    perplexity = 5,
    group = "Group",
    time = "Time",
    id = "Subject_ID",
    color = "Group"
  ) # Visualize data after drift correction
  corrected <- corrected %>%
    assess_quality() %>%
    flag_quality() # Flag low-quality features
  save_QC_plots(
    corrected,
    prefix = file.path(ppath, "figures", paste0(name, "_CLEAN")),
    perplexity = 5,
    group = "Group",
    time = "Time",
    id = "Subject_ID",
    color = "Group"
  ) # Visualize data after removal of low-quality features
  processed[[i]] <- corrected # Save iterated results
}

# Release the cores
register(SerialParam())

# 9. Merge the processed data and create visualisations of the whole data.

merged <- merge_notame_sets(object = processed)
save_QC_plots(
  merged,
  prefix = file.path(ppath, "figures", paste0(name, "_FULL")),
  group = "Group",
  time = "Time",
  id = "Subject_ID",
  color = "Group"
)

# 10. Set a seed number for reproducibility, such as 10, for the random forest
#     imputation.

set.seed(10)

# 11. Perform the first round of imputation with the random forest method on
#     the good-quality features. Please note that RF imputation may generate
#     values that change the biological interpretation of the data and it
#     should be applied with caution (see https://zenodo.org/records/18621244).

imputed <- impute_rf(object = merged, all_features = FALSE)

# 12. Perform the second round of imputation by imputing 1 for the remaining
#     features with more than 50% missing values.

imputed <- impute_simple(object = imputed, value = 1, na_limit = 0.5)

# 13. Perform the third round of imputation with the random forest method for
#     the remaining features.

imputed <- impute_rf(object = imputed, all_features = TRUE)

# 14. Perform batch correction if the samples were analysed in more than one
#     batch. By default, the QC samples are named as “QC”.

library(batchCorr)
batch_corrected <- batchCorr::normalizeBatches(
  peakTableCorr = imputed,
  batches = "Batch",
  sampleGroup = "QC",
  refGroup = "QC",
  population = "all",
  assay.type = "abundances",
  name = "normalized"
)

# 15. Remove the QC sample information and save visualisations once more
# without QC samples.
merged_no_qc <- drop_qcs(imputed) # change object if batch correction was applied
save_QC_plots(
  merged_no_qc,
  prefix = file.path(ppath, "figures", paste0(name, "FULL_NO_QC")),
  group = "Group",
  time = "Time",
  id = "Subject_ID",
  color = "Group"
)

# 16. Optional step: Cluster the molecular features. Use the same time unit as
#     the retention time information in the data.
clustered <- cluster_features(
  object = merged_no_qc,
  rt_window = 1 / 60,
  corr_thresh = 0.9,
  d_thresh = 0.8
)

# 17. Calculate the summary statistics for each molecular feature.

summary_statistics <- summary_statistics(object = imputed, grouping_cols = NULL)

# 18. Calculate the area under the curve (AUC) for each feature and study
#     subject in studies where samples from several time points are included.
#     The resulting object contains the AUC values for each subject and it can
#     be used for statistics, such as t-test (step 21).

aucs <- perform_auc(
  object = imputed,
  time = "Time",
  subject = "Subject_ID",
  group = "Group"
)

# 19. Calculate fold changes between the study groups. If there are more than
#     two study groups, the function will by default calculate the fold changes
#     between all possible pairs. group_col(object) is the group given in step
#     5 and can be replaced with the name of another column containing sample
#     grouping.

fc <- fold_change(object = imputed, group = "Group")

# 20. Calculate Cohen’s d between the study groups. This can be done for two
#     specified groups at a time. If id (for subject ID) and time are given,
#     Cohen’s d will be computed for the change in time.

cohens_d <- cohens_d(
  object = imputed,
  group = "Group",
  id = "Subject_ID",
  time = "Time"
)

# 21. Perform Welch’s t-test between two or more study groups. In case more
#     than two study groups exist, all possible pairs will be tested
#     separately. This function and several other statistical functions use a
#     formula interface, where Feature is replaced by the ID of the molecular
#     feature in each iteration. If the group variances are known to be equal,
#     a Student’s t-test can be performed using the same function by setting
#     the attribute var.equal = TRUE.

t_test_results <- perform_t_test(
  object = imputed,
  formula_char = "Feature ~ Group"
)

# 22. Perform Mann–Whitney U test between two study groups.

mann_whitney_results <- perform_non_parametric(
  object = imputed,
  formula_char = "Feature ~ Group"
)

# 23. Perform paired parametric t-test between two study groups.

paired_t_results <- perform_t_test(
  object = imputed,
  formula_char = "Feature ~ Group",
  is_paired = TRUE,
  id = "Subject_ID"
)

# 24. Perform paired non-parametric t-test (Wilcoxon signed-rank test) between
#     two study groups.

wilcoxon_results <- perform_non_parametric(
  object = imputed,
  formula_char = "Feature ~ Group",
  is_paired = TRUE,
  id = "Subject_ID"
)

# 25. Perform Welch’s ANOVA to compare the averages of two or more study groups.

oneway_anova_results <- perform_oneway_anova(
  object = imputed,
  formula_char = "Feature ~ Group"
)

# 26. Perform Kruskal–Wallis test to compare the averages of two or more study
#     groups.

kruskal_wallis_results <- perform_kruskal_wallis(
  object = imputed,
  formula_char = "Feature ~ Group"
)

# 27. Perform a linear model to study whether the molecular features predict
#     the difference in selected fixed effects, such as group and time, as well
#     as their interaction. An equivalent formula would be "Feature ~ Group *
#     Time".

lm_results <- perform_lm(
  object = imputed,
  formula_char = "Feature ~ Group + Time + Group:Time"
)

# 28. Perform a linear mixed model to study whether the molecular features
#     predict the difference in selected fixed effects, such as group and time,
#     and the contribution of random effects, such as subject ID. The
#     confidence interval for the parameters is by default calculated with the
#     Wald method, which is suitable for large datasets with normally
#     distributed effects. In other cases, profile or boot (bootstrapped)
#     confidence interval should be considered.

lmer_results <- perform_lmer(
  object = imputed,
  formula_char = "Feature ~ Group + Time + (1 | Subject_ID)",
  ci_method = "Wald"
)

# 29. Perform MUVR to select relevant molecular features for the target
#     variable y to predict. The number of iterations (nRep) is recommended to
#     be at least 30.
rf_model <- muvr_analysis(
  object = imputed,
  y = "Group",
  nRep = 30,
  method = "RF"
)

# 30. Perform PLS to identify a stable set of metabolites that best explain the
#     outcome y, while controlling overfitting.

pls_opt <- mixomics_pls_optimize(
  object = object,
  y = "Group",
  ncomp = 3,
  nrepeat = 5
)

# 31. Combine the statistics results (needed for manual annotation of
#     metabolites) into the main preprocessed object. In this example, results
#     from the Mann–Whitney U-test and linear model were chosen. Export the
#     combined data object into an Excel table.

with_results <- join_rowData(imputed, cohens_d)
with_results <- join_rowData(with_results, fc)
with_results <- join_rowData(with_results, mann_whitney_results)
with_results <- join_rowData(with_results, lm_results)

write_to_excel(with_results, file = file.path(ppath, "imputed_statistics.xlsx"))

# 32. Perform manual annotation of metabolites and add manual metabolite ID
#     (column Curated_ID), MSI ID level (column ID_level), and manually chosen
#     representative ion columns (column Representative_ion, where the
#     representative ions have been marked with x). Remove the extra rows, so
#     that Feature_ID etc. become the column names (top rows). Save the Excel
#     file under a new name, such as curated_data.xlsx.

# 33. Join the manual annotation results from the Excel file (via openxlsx
#     package) as extra data columns to the with_results object. Ensure that
#     there is at least one redundant column common for both datasets, such as
#     the unique Feature ID generated by notame.

extra_data <- openxlsx::read.xlsx(file.path("data", "curated_data.xlsx"))
colnames(extra_data)
extra_data <- extra_data[, c(
  "Feature_ID",
  "Curated_ID",
  "ID_level",
  "Representative_ion"
)]
with_results <- join_rowData(object = with_results, extra_data)

# 34. Generate an object containing only the representative ions of each
#     annotated metabolite.

annotated <- with_results[
  !is.na(rowData(with_results)$Representative_ion) &
    rowData(with_results)$Representative_ion == "x",
]

# 35. Calculate Spearman correlation coefficients and p-values between the
#     annotated metabolites.

correlations <- perform_correlation_tests(
  object = annotated,
  x = rownames(annotated),
  method = "spearman"
)

# 36. Calculate Kendall’s τ correlation coefficients and p-values between the
#     annotated metabolites and two sample information variables: time and
#     injection order.

correlations_variables <- perform_correlation_tests(
  object = annotated,
  x = rownames(annotated),
  y = c("Time", "Injection_order"),
  method = "kendall"
)

# 37. Plot a PCA including the first two principal components, unit variance
#     scaling, density curves, Group defining the colour, and Time defining the
#     shape.

plot_pca(
  object = imputed,
  pcs = c(1, 2),
  scale = "uv",
  density = TRUE,
  color = "Group",
  shape = "Time"
)

# 38. Plot an arrow t-SNE including the first two principal components, unit
#     variance scaling, density curves, the colour defining grouping variable
#     1, the time point for constructing the arrows, and subject ID. If the
#     sample size is large and the image becomes crowded, it is possible to
#     show the groups in separate faceted plots by adding + facet_wrap(~Group)
#     at the end of the code.

plot_tsne_arrows(
  object = imputed,
  scale = "uv",
  perplexity = 10,
  color = "Group",
  time = "Time",
  subject = "Subject_ID",
  arrow_style = arrow(
    angle = 30,
    length = unit(0.15, "inches"),
    ends = "last",
    type = "open"
  ),
  line_width = 0.75
) +
  facet_wrap(~Group)

# 39. Generate an effect heatmap with hierarchical clustering from the
#     correlation results (object correlations generated in Step 35) between
#     the annotated metabolites.

plot_effect_heatmap(
  correlations,
  x = "X",
  y = "Y",
  effect = "Correlation_coefficient",
  p = "Correlation_P_FDR",
  point_size_range = c(2, 8),
  discretize_effect = TRUE,
  breaks = 7,
  lower_tri = TRUE
)

# 40. Generate a volcano plot from a pairwise comparison by using the Cohen’s d
#     and FDR-corrected p-values from the pairwise statistical test,
#     log2-transform the x-axis, centre the zero effect point, and label those
#     metabolites with manually curated identifications and FDR < 0.05.

volcano_plot(
  object = with_results,
  x = "B_vs_A_FC",
  p = "GroupB.Time2.p.value",
  color = "B_vs_A_2_minus_1_Cohen_d",
  log2_x = TRUE,
  center_x_axis = TRUE,
  label = "Curated_ID",
  label_limit = 0.05,
  title = NULL
) +
  labs(
    x = "Fold change 2 vs 1",
    y = "p-value",
    color = "Cohen's d"
  )

# 41. Save boxplots of all features by group, with “Curated ID” column used to
#     label the metabolites, into a single PDF file.

save_group_boxplots(
  object = annotated,
  file_path = file.path(ppath, "figures", "group_boxplots.pdf"),
  format = "pdf",
  x = "Group",
  title = "Curated_ID",
  color = "Group"
)

# 42. Save beeswarm plots of unique identified features by group into separate
#     PNG files.

save_beeswarm_plots(
  object = annotated,
  file_path = file.path(ppath, "figures", "beeswarm_plots"),
  format = "png",
  x = "Group",
  color = "Group"
)
