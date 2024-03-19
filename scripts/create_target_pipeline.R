# The ultimate goal of the script is to process, analyze, and visualize this intricate data to identify genes that are expressed differently in males and females. This kind of analysis is crucial for understanding biological differences between sexes at a molecular level, which can have important implications for medical research, including the development of targeted therapies and personalized medicine.
#
# To achieve this, the script performs several key tasks:
#
# Data Preparation: It starts by organizing and preparing the sequencing data, making it ready for analysis.
# This involves filtering the data to focus on specific types of cells or tissues and formatting it into a
# structure suitable for detailed analysis.
#
# Comprehensive Analysis: The script then employs advanced statistical methods to sift through the data,
# looking for genes that show significant differences in expression between male and female samples.
# This step is complex and requires careful handling of the data to ensure accurate results.
#
# Visualizing the Results: After the analysis, the script creates various graphical representations
# of the findings. These visualizations make it easier to understand and interpret the results,
# highlighting the key differences in gene expression between sexes.
#
# Summarizing Insights: Finally, the script compiles the results into a comprehensive format,
# such as a detailed chart or graph, which can be used for further scientific research or reporting.
#
# In summary, this script is a powerful tool in the field of genomics, enabling researchers to delve
# deep into the genetic differences between males and females at a cellular level.
# The insights gained from this analysis have the potential to contribute significantly to our
# understanding of human biology and disease.


library(tidyverse)
library(targets)
library(glue)
library(tictoc)
library(patchwork)
# # Get input from other workflow
result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"
root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap"
store=glue("{result_directory}/_targets__pseudobulk_non_immune_split3")

pseudobulk_df_tissue = tar_read(pseudobulk_df_tissue, store=store)
pseudobulk_df_tissue = pseudobulk_df_tissue |> filter(tissue_harmonised == "blood")
pseudobulk_df_tissue = pseudobulk_df_tissue |> group_by(cell_type_harmonised) |> tar_group()
pseudobulk_df_tissue |> qs::qsave("pseudobulk_df_tissue_blood.qs")

root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc"


store=glue("_targets_tidyomics")


tar_script({

  result_directory = "./"

  #-----------------------#
  # Input
  #-----------------------#
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(qs)
  library(crew)
  library(crew.cluster)

  #-----------------------#
  # Future SLURM
  #-----------------------#

  small_slurm =
    crew_controller_slurm(
      name = "small_slurm",
      slurm_memory_gigabytes_per_cpu = 20,
      slurm_cpus_per_task = 2,
      workers = 200,
      verbose = T

    )

  small_slurm_memory =
    crew_controller_slurm(
      name = "small_slurm_memory",
      slurm_memory_gigabytes_per_cpu = 50,
      slurm_cpus_per_task = 1,
      workers = 10,
      verbose = T

    )

  slurm_3 =
    crew_controller_slurm(
      name = "slurm_3",
      slurm_memory_gigabytes_per_cpu = 20,
      slurm_cpus_per_task = 3,
      workers = 400,
      verbose = T

    )

  big_slurm =
    crew_controller_slurm(
      name = "big_slurm",
      slurm_memory_gigabytes_per_cpu = 4,
      slurm_cpus_per_task = 5,
      workers = 300,
      verbose = T

    )

  #-----------------------#
  # Packages
  #-----------------------#
  tar_option_set(
    packages = c(
      "CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
      "Seurat", "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment", "edgeR", "jascap", "crew", "magrittr", "digest", "glmmSeq"
    ),

    garbage_collection = TRUE,
    #trust_object_timestamps = TRUE,
    memory = "transient",
    storage = "worker",
    retrieval = "worker",
    error = "continue",
    format = "qs",
    controller = crew_controller_group(small_slurm, small_slurm_memory, big_slurm, slurm_3),
    resources = tar_resources(crew = tar_resources_crew("small_slurm"))
    # ,
    # debug = "estimates_sex_tissue_14f58586", # Set the target you want to debug.
    # cue = tar_cue(mode = "never") # Force skip non-debugging outdated targets.
  )

  # Define a function to process and split metadata
  split_metadata = function(metadata_DB_path){
    # Read metadata from an RDS file and add a column for age in days
    metadata = readRDS(metadata_DB_path) |> mutate(age_days = age_days_original)

    # Start processing the metadata
    get_metadata() |>
      # Filter samples based on a condition
      filter(
        sample_ %in% (
          !!metadata |>
            distinct(sample_) |>
            tidyr::extract(sample_, "sample_", "([a-zA-Z0-9]+)_.+") |>
            pull(sample_)
        )) |>

      # Join with additional metadata about cell types
      left_join(
        read_csv("~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv") |>
          replace_na(list(lineage_1 = "other_non_immune")) |>
          mutate(is_immune = as.character(lineage_1 == "immune")),
        by = join_by(cell_type),
        copy = TRUE
      ) |>

      # Filter and format the dataset
      filter(!is.na(lineage_1)) |>
      filter(!cell_type_harmonised %in% c("platelet", "immune_unclassified")) |>
      as_tibble() |>

      # Filter out rare cell types (less than 100 occurrences)
      add_count(cell_type_harmonised, name = "count_cell_type_harmonised") |>
      filter(count_cell_type_harmonised > 100) |>
      select(-count_cell_type_harmonised) |>

      # Format the assay field by replacing spaces and dashes
      mutate(assay = assay |> str_replace_all(" ", "_") |> str_replace_all("-", "_")  |> str_remove_all("'")) |>

      # Simplify ethnicity information
      mutate(
        ethnicity = case_when(
          ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
          ethnicity |> str_detect("African") ~ "African",
          TRUE ~ ethnicity
        )
      ) |>

      # Filter out unknown developmental stages
      filter(development_stage != "unknown") |>

      # Simplify ethnicity categorization
      mutate(ethnicity_simplified = case_when(
        ethnicity %in% c("European", "Chinese", "African", "Hispanic or Latin American") ~ ethnicity,
        TRUE ~ "Other"
      )) |>
      mutate(
        ethnicity_simplified =
          ethnicity_simplified |>
          fct_relevel(c("European", "Chinese", "African", "Hispanic or Latin American", "Other"))
      ) |>

      # Simplify assay information
      mutate(assay_simplified = if_else(assay |> str_detect("10x"), "10x", assay)) |>
      mutate(assay_simplified = factor(assay_simplified)) |>

      # Set a baseline for disease status
      mutate(disease = if_else(disease == "normal", "aaa_normal", disease)) |>

      # Select specific columns for further analysis
      select(
        cell_, cell_type_harmonised, sample_,  file_id, file_id_db,
        age_days, development_stage, sex, tissue_harmonised, disease,
        ethnicity_simplified, assay_simplified, is_immune
      ) |>

      # Filter out certain cell types and further categorize cell types
      filter(cell_type_harmonised != "animal_cell") |>
      mutate(cell_type_harmonised = case_when(
        cell_type_harmonised |> str_detect("fibro") ~ "stromal_cell",
        cell_type_harmonised |> str_detect("chondro") ~ "connective_tissue_cell",
        cell_type_harmonised |> str_detect("adipoc") ~ "fat_cell",
        cell_type_harmonised |> str_detect("hematopoietic") ~ "hematopoietic_cell",
        cell_type_harmonised |> str_detect("epithe") ~ "epithelial_cell",
        cell_type_harmonised |> str_detect("hepato") ~ "hepatic_cell",
        cell_type_harmonised |> str_detect("keratino") ~ "keratinocyte",
        cell_type_harmonised == "neuron" ~ "neural_cell",
        cell_type_harmonised |> str_detect("progenitor") ~ "hematopoietic_cell",
        cell_type_harmonised |> str_detect("stem") ~ "hematopoietic_cell",
        cell_type_harmonised |> str_detect("tendon") ~ "connective_tissue_cell",
        TRUE ~ cell_type_harmonised
      )) |>

      # Split datasets with many samples into chunks of maximum 100 samples
      nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_)) |>
      with_groups(
        c(cell_type_harmonised, tissue_harmonised, file_id, is_immune),
        ~ .x |>
          arrange(file_id) |>
          mutate(sample_chunk = row_number() |>
                   divide_by(100) |>
                   ceiling()
          )
      ) |>
      unnest(data) |>

      # Nest the data again for final group-wise processing
      nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_chunk))

  }

  # Define a function 'get_sce' that takes metadata as input
  get_sce = function(tissue_cell_type_metadata) {

    # Process each item in the metadata
    tissue_cell_type_metadata |>
      mutate(data = pmap(
        # Map over data, cell_type_harmonised, and tissue_harmonised columns
        list(data, cell_type_harmonised, tissue_harmonised),
        ~ ..1 |> # Unpack the first argument (data)
          # Fetch single cell experiment data from a cache directory
          get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated") |>
          # Mutate to create a unique identifier for each sample
          mutate(sample_se =
                   # Use the glue function to combine sample, disease, cell type, and tissue into a unique ID
                   glue("{sample_}___{disease}___{..2}___{..3}") |>
                   # Replace spaces and slashes in the ID for standardization
                   str_replace_all(" ", "_") |>
                   str_replace_all("/", "__")
          )
      ))

  }

  # Define a function 'get_pseudobulk' that takes a data frame 'sce_df' as input
  get_pseudobulk = function(sce_df) {

    # Process the input data frame
    sce_df |>
      # Use 'mutate' to apply a function to each element of the 'data' column
      mutate(data = map(
        data,
        # For each element in 'data', apply a function
        ~ .x |> # '.x' refers to the current element of 'data'
          # Aggregate single-cell data into pseudobulk data
          # The function 'aggregate_cells' is from the tidySingleCellExperiment package
          # It aggregates data based on the 'sample_se' column which is a sample identifier
          tidySingleCellExperiment::aggregate_cells(.sample = sample_se)
      ))

  }

  # In this function, the relationship between two specified columns (.col1 and .col2) in a dataset is analyzed to detect a potential confounder. The function nests the data based on each column in turn and calculates the number of unique values in the other column for each group. This can help in understanding how these two variables are related and whether they might confound each other in a statistical analysis. For example, if .col1 represents "ethnicity" and .col2 represents "assay type", this function would help determine how many unique assay types there are for each ethnicity and vice versa.
  # Define a function 'nest_detect_complete_confounder' that takes a dataset and two column names as input
  nest_detect_complete_confounder = function(.data, .col1, .col2){
    # Capture the column names as quosures (to allow for non-standard evaluation)
    .col1 = enquo(.col1)
    .col2 = enquo(.col2)

    # Begin processing the input data
    .data |>
      # Nest data except for the two specified columns
      nest(se_data = -c(!!.col1, !!.col2)) |>

      # Analyze the first column against the second
      # First, nest the data except for the first column
      nest(data = -!!.col1) |>
      # Count the unique values in the second column for each group of the first column
      mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
      # Unnest the data to revert to the original format
      unnest(data) |>

      # Analyze the second column against the first
      # Nest the data except for the second column
      nest(data = -!!.col2) |>
      # Count the unique values in the first column for each group of the second column
      mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
      # Unnest the data to revert to the original format
      unnest(data)
  }

  # This function processes a dataset by first determining the number of unique values in .col2 for each unique value in .col1 and vice versa. Then, it joins this information back to the original dataset. This approach can be useful for understanding the interaction between two variables and checking if one variable is a complete confounder of the other. For instance, in a dataset where .col1 is "ethnicity" and .col2 is "assay type", this function helps to assess how many assay types are used per ethnicity and how many ethnicities are involved per assay type, which could influence the analysis results.
  # Define a function 'left_join_detect_complete_confounder' that takes a dataset and two column names as input
  left_join_detect_complete_confounder = function(.data, .col1, .col2){
    # Capture the column names as quosures (to allow for non-standard evaluation)
    .col1 = enquo(.col1)
    .col2 = enquo(.col2)

    # Calculate counts for the unique values in the two columns
    counts =
      .data |>
      # Select distinct pairs of values from the two columns
      distinct(!!.col1, !!.col2) |> # Quicker than nesting

      # Nest data except for the first column
      nest(data = -!!.col1) |>
      # Count the unique values in the second column for each group of the first column
      mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
      # Unnest the data
      unnest(data) |>

      # Nest data except for the second column
      nest(data = -!!.col2) |>
      # Count the unique values in the first column for each group of the second column
      mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
      # Unnest the data
      unnest(data)

    # Join the original dataset with the calculated counts
    .data |>
      left_join(counts) # Join on the two specified columns
  }

  # In this function, the left_join_detect_complete_confounder function is used to join the original data with counts of unique values in the two specified columns (.col1 and .col2). It then filters out rows where the sum of these counts (n1 + n2) is not greater than 2. This filtering step aims to remove data points that could be complete confounders, thus potentially improving the quality and reliability of subsequent analyses. The function finally removes the count columns from the dataset, leaving only the relevant data for analysis.
  # Define a function 'drop_samples_complete_confounder' that takes a dataset and two column names as input
  drop_samples_complete_confounder = function(.data, .col1, .col2){
    # Capture the column names as quosures (to allow for non-standard evaluation)
    .col1 = enquo(.col1)
    .col2 = enquo(.col2)

    # Process the input dataset
    .data |>
      # Join the dataset with counts of unique values for the two specified columns
      # This function is defined to detect complete confounders
      left_join_detect_complete_confounder(!!.col1, !!.col2) |>
      # Filter out rows where the sum of unique counts for both columns is not greater than 2
      # This step aims to remove rows where there might be complete confounding
      filter(n1 + n2 > 2) |>
      # Remove the count columns (n1 and n2) from the dataset
      select(-n1, -n2)
  }

  # Define a function to process a dataset 'se' (likely a SingleCellExperiment object)
  samples_NOT_complete_confounders_for_ethnicity_assay = function(se){
    # Process the dataset
    se =
      se |>
      # Nest data except for assay and ethnicity
      nest(se_data = -c(assay_simplified, ethnicity_simplified)) |>

      # Count distinct ethnicities for each assay type
      nest(data = -assay_simplified) |>
      mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |>
      unnest(data) |>

      # Count distinct assays for each ethnicity
      nest(data = -ethnicity_simplified) |>
      mutate(n2 = map_int(data, ~ .x |> distinct(assay_simplified) |> nrow())) |>
      unnest(data)

    # Select an assay as a dummy value based on the highest sum of n1 and n2
    dummy_assay = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(assay_simplified)

    # Replace assay type in cases with less than 3 distinct values in either ethnicity or assay
    se |>
      mutate(assay_simplified = if_else(n1 + n2 < 3, dummy_assay, assay_simplified)) |>
      # Clean up the dataset by removing n1 and n2 columns
      select(-n1, -n2) |>
      # Unnest the summarized experiment data
      unnest_summarized_experiment(se_data)
  }

  # Define a function to process a dataset 'se' considering ethnicity and disease
  samples_NOT_complete_confounders_for_ethnicity_disease = function(se){
    # Process the dataset
    se =
      se |>
      # Nest data except for disease and ethnicity
      nest(se_data = -c(disease, ethnicity_simplified)) |>

      # Count distinct ethnicities for each disease
      nest(data = -disease) |>
      mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |>
      unnest(data) |>

      # Count distinct diseases for each ethnicity
      nest(data = -ethnicity_simplified) |>
      mutate(n2 = map_int(data, ~ .x |> distinct(disease) |> nrow())) |>
      unnest(data)

    # Select an ethnicity as a dummy value based on the highest sum of n1 and n2
    dummy_ethnicity = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(ethnicity_simplified)

    # Replace ethnicity in cases with less than 3 distinct values in either ethnicity or disease
    se |>
      mutate(ethnicity_simplified = if_else(n1 + n2 < 3, dummy_ethnicity, ethnicity_simplified)) |>
      # Clean up the dataset by removing n1 and n2 columns
      select(-n1, -n2) |>
      # Unnest the summarized experiment data
      unnest_summarized_experiment(se_data)
  }

  # Define a function 'aggregate' for processing and aggregating single-cell experiment data
  aggregate = function(se_df){
    # Print a message indicating the start of the aggregation process
    print("Start aggregate")
    # Perform garbage collection to free up memory
    gc()

    # Process the input data frame
    se_df |>
      # Mutate to apply a function to each element of the 'data' column
      mutate(data = pmap(
        list(data, cell_type_harmonised, tissue_harmonised, file_id),
        ~ {
          # Add relevant columns to the SingleCellExperiment object
          se =
            ..1 |>
            mutate(cell_type_harmonised = ..2, tissue_harmonised = ..3, file_id = ..4) |>
            select(-any_of(c("file_id_db", ".cell", "original_cell_id")))

          # Identify samples with a high number of genes (more than 5000)
          sample_with_many_genes =
            se |>
            assay("counts") |>
            apply(2, function(x) (x>0) |> which() |> length()) |>
            enframe() |>
            mutate(value = as.character(value), name = as.character(name)) |>
            filter(value > 5000) |>
            pull(name)
          # Filter the dataset to keep only samples with many genes
          se = se[,sample_with_many_genes, drop=FALSE]

          # Further filter the dataset to keep samples with more than 10 cells
          se = se |> filter(.aggregated_cells > 10)
        },
        .progress=TRUE
      )) |>

      # Nest the data for output preparation
      nest(data = -c(tissue_harmonised, cell_type_harmonised)) |>

      # Bind the nested data
      mutate(data = map(
        data,
        ~ {
          # Combine data from different SingleCellExperiment objects
          se = do.call(cbind, .x |> pull(data))

          # Filter out genes that are very rare (all zero counts or total counts less than 15)
          all_zeros = assay(se, "counts") |> rowSums() |> equals(0)
          se = se[!all_zeros,]
          lower_than_total_counts = assay(se, "counts") |> rowSums() < 15
          se = se[!lower_than_total_counts,]

          # Convert the counts matrix to a sparse matrix to optimize memory usage
          se@assays@data$counts = as(se@assays@data$counts, "sparseMatrix")

          # Return the processed SingleCellExperiment object
          se
        }
      ))

  }

  # Define a function 'map_quantile_scale_abundance' for scaling and normalizing abundance in single-cell experiment data
  map_quantile_scale_abundance = function(se_df){
    # Print a message indicating the start of the scaling process
    print("Start scale abundance")
    # Perform garbage collection to free up memory
    gc()

    # Process the input data frame
    se_df |>
      # Apply quantile normalization to each element of the 'data' column
      mutate(data = map(data, quantile_normalise_abundance, method = "preprocesscore_normalize_quantiles_use_target")) |>

      # Convert the scaled counts to a sparse matrix format
      mutate(data = map(data, ~ {
        # Convert the 'counts_scaled' matrix to a sparse matrix
        .x@assays@data$counts_scaled = as(.x@assays@data$counts_scaled, "sparseMatrix")

        # Clean up the environment attribute to avoid memory leaks
        attr(.x, "internals")$tt_columns$.abundance_scaled |> attr(".Environment") = NULL # new_environment()

        # Return the processed SingleCellExperiment object
        .x
      }))
  }

  # Define a function 'map_keep_abundant' for filtering abundant genes in single-cell experiment data
  map_keep_abundant = function(se_df){
    # Print a message indicating the start of the process
    print("Start keep abundant")
    # Perform garbage collection to free up memory
    gc()

    # Process the input data frame
    se_df |>
      mutate(data = map(
        data,
        ~ {
          # Check if the number of columns (samples) is greater than 1000
          # This is likely a check for a specific type of data or condition, like blood
          if(ncol(.x) > 1000)
            # For datasets with many samples, apply a specific filtering function
            .x |>
            keep_abundant(
              .abundance = counts_scaled,
              factor_of_interest = c(sex, ethnicity_simplified),
              minimum_counts = 500,
              minimum_proportion = 0.9
            )

          else {
            # For datasets with fewer samples, apply a different strategy

            # Remove duplicate columns
            .x = .x[, !colnames(.x) |> duplicated()]

            # Determine the abundant genes based on nested conditions
            abundant_genes =
              .x |>
              nest(data = -cell_type_harmonised) |>
              mutate(abundant_genes = map(
                data,
                ~ {
                  # Apply the keep_abundant function to each nested dataset
                  se = .x
                  se |>
                    keep_abundant(
                      .abundance = counts_scaled,
                      factor_of_interest = c(sex, ethnicity_simplified),
                      minimum_counts = 50
                    ) |>
                    rownames()
                },
                .progress = TRUE
              )) |>
              pull(abundant_genes) |>
              unlist() |>
              unique()

            # Filter the original dataset to keep only the abundant genes
            .x |> filter(.feature %in% abundant_genes)
          }

        }
      ))
  }

  # This function performs several steps to calculate dispersion:
  # Filters out samples with only one instance in the categories of sex or ethnicity to avoid issues with statistical indeterminacy.
  # Drops samples that are complete confounders based on sex and ethnicity.
  # Calculates dispersion estimates differently based on the size of the dataset:
  # For smaller datasets (<1000 samples), it uses all samples to calculate dispersion.
  # For larger datasets, it samples a subset of the data (up to 2000 samples) and then calculates dispersion on this subset to reduce computational load.
  # Dispersion is calculated using functions from RNA-seq analysis packages (like estimateDisp and estimateTrendedDisp), which handle the complexity of dispersion estimation in gene expression data. This approach is critical in differential expression analysis, as it helps to account for the variability inherent in such data.
  # # Define a function 'se_add_dispersion' to add dispersion estimates to single-cell experiment data
  se_add_dispersion = function(se_df){
    # Print a message indicating the start of the process
    print("Start add dispersion")
    # Perform garbage collection to free up memory
    gc()

    # Process the input data frame
    se_df |>
      mutate(data = map(
        data,
        ~ {
          # Assign the current SingleCellExperiment object to 'se'
          se = .x

          # Check if there are enough columns (samples) to proceed
          if(ncol(se) > 0) {
            # Filter out ethnicities with only one sample
            ethnicity_to_keep = se |> pivot_sample() |> count(ethnicity_simplified) |> filter(n > 1) |> pull(ethnicity_simplified)
            se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)

            # Filter out sex categories with only one sample
            sex_to_keep = se |> pivot_sample() |> count(sex) |> filter(n > 1) |> pull(sex)
            se = se |> filter(sex %in% sex_to_keep)
          }

          # Drop complete confounders based on sex and ethnicity
          if(ncol(se) > 0) {
            se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
          }

          # Handle cases where the SingleCellExperiment object is empty
          if(ncol(se) == 0) {
            rowData(se)$dispersion = rep(NA, nrow(se))
          }
          else if(ncol(se) < 1000) {
            # Calculate dispersion for datasets with fewer than 1000 columns
            # Create a model design matrix based on sex and ethnicity
            factors =
              c("sex", "ethnicity_simplified") |>
              enframe(value = "factor") |>
              mutate(n = map_int(
                factor, ~ se |> select(.x) |> distinct() |> nrow()
              )) |>
              filter(n > 1) |>
              pull(factor) |>
              str_c(collapse = " + ")

            my_design = glue("~ {factors}") |> as.formula() |> model.matrix(data = colData(se) |> droplevels())
            rowData(se)$dispersion = se |> assay("counts_scaled") |> estimateDisp(design = my_design) %$% tagwise.dispersion
          }
          else {
            # For larger datasets, sample a subset of columns
            sampled_samples = sample(seq_len(ncol(se)), size = min(ncol(se), 2000))

            # Create a model design matrix for the sampled subset
            factors =
              c("sex", "ethnicity_simplified") |>
              enframe(value = "factor") |>
              mutate(n = map_int(
                factor, ~ se[,sampled_samples, drop=FALSE] |> select(.x) |> distinct() |> nrow()
              )) |>
              filter(n > 1) |>
              pull(factor) |>
              str_c(collapse = " + ")

            my_design = model.matrix(~ sex + ethnicity_simplified, data = se[,sampled_samples, drop=FALSE] |> colData() |> droplevels())

            # Estimate trended dispersion based on the sampled subset
            rowData(se)$dispersion =
              assay(se[,sampled_samples, drop=FALSE], "counts_scaled") |>
              estimateTrendedDisp(design = my_design, subset=1000, rowsum.filter=10)
          }

          # Return the processed SingleCellExperiment object
          se
        }
      ))

  }

  # In this function, the number of gene chunks is calculated based on the number of distinct samples (sample_se) in each SingleCellExperiment object. The calculation involves dividing the number of distinct samples by 50, multiplying the result by 10, and then rounding up to the nearest whole number. There's also a safeguard to ensure a minimum of one chunk and a maximum limit of 10,000 chunks. This method of chunking is typically used to balance computational efficiency with the need to process large datasets in manageable segments.
  # Define a function 'add_number_of_gene_chunks' to calculate and add the number of gene chunks in single-cell experiment data
  add_number_of_gene_chunks = function(se_df){
    # Print a message indicating the start of the chunk calculation process
    print("Start add number of chunks")
    # Perform garbage collection to free up memory
    gc()

    # Process the input data frame
    se_df |>
      mutate(number_of_chunks = map_int(
        data,
        ~ {
          # Process the current SingleCellExperiment object
          se = .x
          # Count the distinct sample identifiers (sample_se)
          distinct_samples = se |>
            distinct(sample_se) |>
            nrow()

          # Calculate the number of chunks
          # Ensure there's at least one chunk even if there are no samples
          # Divide the number of samples by 50 and scale it up by a factor of 10
          # The ceiling function ensures the result is rounded up to the nearest whole number
          number_of_chunks = max(1, distinct_samples) |>
            divide_by(50) |>
            multiply_by(10) |>
            ceiling()

          # Limit the number of chunks to a maximum of 10000
          # This is likely a safeguard to prevent excessive computational load
          min(number_of_chunks, 10000)
        }
      ))
  }

  map_analyse_sex_dataset = function(se_df, max_rows_for_matrix_multiplication = NULL, cores = 1){

    se_df |>
      mutate(data = map(
        data,
        ~ {

          # Filter
          se =
            .x |>
            filter(sex != "unknown") |>


            # Eliminate complete confounders
            samples_NOT_complete_confounders_for_ethnicity_assay() |>
            samples_NOT_complete_confounders_for_ethnicity_disease()

          rm(.x)
          gc()

          # Filter disease
          se =
            se |>
            filter(disease %in% (
              se |>
                distinct(disease, sex) |>
                count(disease) |>
                filter(n>1) |>
                pull(disease)
            ))

          # Filter file_id
          se =
            se |>
            filter(file_id %in% (
              se |>
                distinct(file_id, sex) |>
                count(file_id) |>
                filter(n>1) |>
                pull(file_id)
            ))


          # # Filter tissue that has two sexes
          # if(ncol(se)>0)
          # 	se =
          # 	se |>
          # 	filter(tissue_harmonised %in% (
          # 		se |>
          # 			distinct(tissue_harmonised, sex) |>
          # 			count(tissue_harmonised) |>
          # 			filter(n>1) |>
          # 			pull(tissue_harmonised)
          # 	))
          #

          # Return prematurely
          if(ncol(se) == 0) return(se)
          if(se |> distinct(sex, ethnicity_simplified) |> count(ethnicity_simplified) |> pull(n) |> max() == 1) return(se)

          # Build the formula
          factors =
            c("age_days", "sex", "ethnicity_simplified", "assay_simplified",  ".aggregated_cells", "disease") |>
            enframe(value = "factor") |>
            mutate(n = map_int(
              factor, ~ se |> select(.x) |> distinct() |> nrow()
            )) |>
            filter(n>1) |>
            pull(factor) |>
            str_c(collapse = " + ")

          random_effects =

            # SOME DO NOT WORK LIKE THIS
            # c("age_days", "sex", "ethnicity_simplified") |>

            c("age_days", "sex") |>
            enframe(value = "factor") |>
            mutate(n = map_int(
              factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
            ))   |>
            filter(n>1) |>
            pull(factor) |>
            str_c(collapse = " + ")

          # The default
          my_formula = glue("~ {factors}")
          method = "edgeR_quasi_likelihood"

          if(
            se |> distinct(file_id) |> nrow() > 1 &
            length(random_effects) > 0
          ) {
            my_formula = glue("{my_formula} + (1 + {random_effects} | file_id)")
            method = "glmmseq_lme4"
          }

          # Add the interaction
          my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")

          se =
            se |>

            # Scale continuous variables
            mutate(age_days = scale(age_days) |> as.numeric()) |>

            # otherwise I get error for some reason
            mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised", "cell_type_harmonised")), as.character)) |>
            mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European"))

          # Drop random effect grouping with no enough data
          combinations_to_keep =
            se |>
            distinct(sex, ethnicity_simplified) |>
            add_count(ethnicity_simplified) |>
            filter(n>1)

          se =
            se |>
            right_join(combinations_to_keep)

          # Filter file_id to keep
          file_id_to_keep =
            se |>
            distinct(sample_, file_id) |>
            count(  file_id) |>
            filter(n > 3) |>
            pull(file_id)

          se = se |>  filter(file_id %in% file_id_to_keep)

          # Use fast method but does not have dispersion
          if(ncol(se) > 2000) method = "glmmseq_glmmTMB"

          se = se |>

            # Test
            test_differential_abundance(
              as.formula(my_formula),
              .abundance = counts_scaled,
              method = method,
              cores = min(nrow(se), cores),
              max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
              .dispersion = dispersion
            )

          attr(se, "internals")$glmmseq_glmmTMB = NULL

          se
        }))

  }


  #-----------------------#
  # Pipeline
  #-----------------------#
  list(

    # Grouped file
    tar_target(
      pseudobulk_df_blood,
      qs::qread(glue("{result_directory}/pseudobulk_df_tissue_blood.qs")),
      deployment = "main"
      #resources = tar_resources(crew = tar_resources_crew("small_slurm_memory"))
    ),

    # Group samples
    tarchetypes::tar_group_by(
      pseudobulk_df_blood_cell_type,
      pseudobulk_df_blood,
      cell_type_harmonised,
      deployment = "main"
    ),

    # Add dispersion
    tar_target(
      sce_df_split_by_gene,
      pseudobulk_df_blood_cell_type |>
        aggregate() |>
        map_quantile_scale_abundance() |>
        map_keep_abundant() |>
        se_add_dispersion() |>
        add_number_of_gene_chunks() |>
        map_split_se_by_gene(data, number_of_chunks),
      pattern = map(pseudobulk_df_blood_cell_type),
      iteration = "group",
      resources = tar_resources(crew = tar_resources_crew("small_slurm_memory"))
      # ,
      #  deployment = "main"
      #resources = small_slurm_resource
    ),

    # Parallelise rows
    tarchetypes::tar_group_by(
      sce_df_split_by_gene_grouped,
      sce_df_split_by_gene,
      se_md5,
      deployment = "main",
      memory = "persistent"
    ) ,


    # Analyse
    tar_target(
      estimates_sex_tissue,
      sce_df_split_by_gene_grouped |> map_analyse_sex_dataset(max_rows_for_matrix_multiplication = 10000, cores = 6),
      pattern = map(sce_df_split_by_gene_grouped),
      iteration = "group",
      resources = tar_resources(crew = tar_resources_crew("slurm_3")),
      cue = tar_cue(mode = "never")
      # ,
      # retrieval = "main"
    )
  )


}, ask = FALSE, script = glue("{store}.R"))



tar_make(
  script = glue("{store}.R"),
  store = store
  # ,
  # callr_function = NULL
)



# Loading required libraries
library(furrr)  # For parallel processing
library(tidybulk)  # For bulk RNA-seq data analysis

# Setting up parallel processing plan with 18 workers
plan(multisession, workers = 18)
# Setting the maximum size for global objects passed to future expressions
options(future.globals.maxSize = 200000 * 1024^2)

# Retrieving and preprocessing blood data
blood =
  tar_meta(  store=store  ) |>  # Accessing metadata from a target store
  filter(name |> str_detect("estimates_sex_tissue")) |>  # Filtering for specific dataset names
  filter(!is.na(data)) |>  # Removing entries with NA data
  filter(name != "estimates_sex_tissue") |>  # Excluding a specific dataset
  mutate(se = future_map(
    name,
    ~{  # Applying function to each element in 'name'
      #print(.x)  # Uncomment to print the current name
      .x |>
        tar_read_raw(store=store ) |>  # Reading raw data from target store
        mutate(data = map(data, pivot_transcript))  # Pivoting transcript data
    },
    .progress = TRUE  # Showing progress
  ))

# Preparing data for differential expression analysis
de =
  blood |>
  select(se) |>  # Selecting the 'se' column
  unnest(se) |>  # Expanding list columns
  unnest(data) |>  # Further expanding nested data
  filter(cell_type_harmonised != "thymocyte")  # Filtering out specific cell types

# Saving the prepared data to an RDS file
de |> saveRDS("de_blood.rds")

# Sourcing a script for cell type color mapping
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/f0b6bf9f59847c8b9f0a638262a6b8dd697affb7/color_cell_types.R")
# Getting colors for cell types
cell_type_color =
  de |>
  pull(cell_type_harmonised) |>  # Extracting cell type column
  unique() |>  # Getting unique cell types
  get_cell_type_color()  # Applying function to get colors

# Adjusting p-values for multiple testing
de =
  de |>
  with_groups(cell_type_harmonised, ~ .x |> mutate(
    P_sex_adjusted = p.adjust(P_sex, method = "BH"),
    P_age_days_adjusted = p.adjust(P_age_days, method = "BH"),
    P_age_days.sex_adjusted = p.adjust(P_age_days.sex, method = "BH")
  ))

# Determining the most differentially expressed genes by sex in each cell type
df_sex_cell_type_most_de =
  de |>
  filter(!is.na(P_sex_adjusted)) |>  # Filtering out NA adjusted p-values
  add_count(cell_type_harmonised, name = "gene_number") |>  # Counting genes per cell type
  mutate(de_only_in_sex = P_sex_adjusted < 0.05 | P_age_days.sex_adjusted < 0.05 & P_age_days_adjusted > 0.05) |>  # Determining differential expression by sex
  dplyr::count(
    cell_type_harmonised, gene_number,
    de_only_in_sex
  ) |>  # Counting differential expressions
  mutate(proportion_of_significant = n/gene_number) |>  # Calculating proportion of significant genes
  filter(de_only_in_sex)  # Filtering for differential expression by sex only

# Creating a bar plot for the proportion of significant genes by cell type
plot_sex_cell_type_most_de =
  df_sex_cell_type_most_de |>
  ggplot(aes(fct_reorder(cell_type_harmonised,proportion_of_significant, .desc = TRUE), proportion_of_significant, fill = cell_type_harmonised)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_type_color) +
  ylab("Proportion of significant genes") +
  xlab("Cell types") +
  guides(fill = "none") +
  scale_x_discrete(labels = function(x) x |> str_replace("terminal effector", "eff") |> str_remove("_cell") |> str_remove("cyte") ) +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Loading ggupset for creating UpSet plots
library(ggupset)
# Creating an UpSet plot to visualize intersections of cell types and features
plot_sex_cell_type_upset =
  de |>
  filter(!is.na(P_sex_adjusted)) |>  # Filtering out NA adjusted p-values
  filter( P_sex_adjusted < 0.05 | P_age_days.sex_adjusted < 0.05 & P_age_days_adjusted > 0.05) |>  # Filtering for significant genes
  filter(cell_type_harmonised %in% (df_sex_cell_type_most_de |> arrange(desc(proportion_of_significant)) |> head(9) |> pull(cell_type_harmonised))) |>  # Filtering top cell types
  select(cell_type_harmonised, .feature) |>  # Selecting relevant columns
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("terminal effector", "eff")) |>  # Modifying cell type names
  nest(cell_types = cell_type_harmonised) |>  # Nesting data by cell types
  mutate(cell_types = map(cell_types, pull, cell_type_harmonised)) |>  # Transforming nested data
  ggplot(aes(x=cell_types)) +
  geom_bar() +
  scale_x_upset(n_intersections = 20) +
  ylab("Significant gene count") +
  theme_multipanel +
  theme(axis.title.x = element_blank()) +
  theme_combmatrix(
    combmatrix.panel.point.size = 0.5,
    combmatrix.panel.line.size = 0.3, combmatrix.label.height = unit(15, "mm")
  )

# Reading chromosome data for genes
gene_chr = read_csv("~/PostDoc/immuneHealthyBodyMap/symbol_chr.csv")

# Calculating the proportion of interaction-only significant genes
proportion_of_interaction_only =
  de |>
  filter(cell_type_harmonised != "thymocyte") |>
  filter(!is.na(P_sex_adjusted)) |>  # Filtering out NA adjusted p-values
  filter(!.feature %in% gene_chr$ID) |>  # Excluding genes listed in gene_chr
  add_count(cell_type_harmonised, name = "gene_number") |>  # Counting genes per cell type
  mutate(de_only_in_sex = P_sex_adjusted < 0.05 & P_age_days.sex_adjusted > 0.05 & P_age_days_adjusted > 0.05) |>  # Defining differential expression by sex
  mutate(de_only_in_interaction = P_sex_adjusted > 0.05 & P_age_days.sex_adjusted < 0.05 & P_age_days_adjusted > 0.05) |>  # Defining differential expression by interaction
  dplyr::count(
    cell_type_harmonised, gene_number,
    de_only_in_sex, de_only_in_interaction
  ) |>  # Counting differential expressions
  complete(nesting(cell_type_harmonised, gene_number), de_only_in_sex, de_only_in_interaction, fill = list(n=0)) |>
  filter(de_only_in_sex + de_only_in_interaction == 1) |>  # Filtering for exclusive categories
  mutate(proportion_of_significant = n/gene_number) |>  # Calculating proportion of significant genes
  mutate(label = if_else(de_only_in_sex, "de_only_in_sex", "de_only_in_interaction")) |>  # Labeling categories
  ggplot(aes(label, proportion_of_significant, fill=label)) +
  geom_boxplot(outlier.shape = NA, lwd=0.3, fatten=0.3) +
  geom_jitter(width = 0.2, size=0.2) +
  scale_fill_brewer(palette="Set2") +
  xlab("Proportion of significant genes") +
  guides(fill="none") +
  coord_flip() +
  theme_multipanel +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(angle=90, hjust = 0.5))

# Calculating the contribution of age-interaction in the number of sex-dependent genes
de |>
  filter(!is.na(P_sex_adjusted)) |>
  filter(!.feature %in% gene_chr$ID) |>
  add_count(cell_type_harmonised, name = "gene_number") |>
  mutate(de_only_in_sex = P_sex_adjusted < 0.05 &	P_age_days.sex_adjusted > 0.05 & P_age_days_adjusted > 0.05) |>
  mutate(de_only_in_interaction = P_sex_adjusted > 0.05 &	P_age_days.sex_adjusted < 0.05 & P_age_days_adjusted > 0.05) |>
  dplyr::count(
    cell_type_harmonised, gene_number,
    de_only_in_sex, de_only_in_interaction
  ) |>
  filter(de_only_in_sex + de_only_in_interaction == 1) |>
  mutate(proportion_of_significant = n/gene_number) |>
  mutate(label = if_else(de_only_in_sex, "de_only_in_sex", "de_only_in_interaction")) |>
  select(cell_type_harmonised, proportion_of_significant, label) |>
  pivot_wider(names_from = label , values_from = proportion_of_significant) |>
  mutate(contribution = de_only_in_interaction / (de_only_in_sex + de_only_in_interaction) ) |>
  pull(contribution) |>
  mean(na.rm=TRUE)

p = ((
  plot_spacer() | plot_sex_cell_type_most_de
) + plot_layout(  width = c(2, 1) )) / ((
  plot_sex_cell_type_upset | proportion_of_interaction_only | plot_spacer()
) + plot_layout(  width = c(1, 1, 1) ) ) +
  plot_layout(  height = c(30, 20) )

ggsave(
  "plot_bio_application.pdf",
  plot = p,
  units = c("mm"),
  width = 70 *2,
  height = 55 *2 ,
  limitsize = FALSE
)
t
