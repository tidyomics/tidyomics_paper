library(tidyverse)
library(targets)
library(glue)
library(tictoc)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")
library(patchwork)
# # Get input from other workflow
# result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap/pseudobulk_0.2.3.5_non_immune"
# root_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/immuneHealthyBodyMap"
# store=glue("{result_directory}/_targets__pseudobulk_non_immune_split3")
#
# pseudobulk_df_tissue = tar_read(pseudobulk_df_tissue, store=store)
# pseudobulk_df_tissue = pseudobulk_df_tissue |> filter(tissue_harmonised == "blood")
# pseudobulk_df_tissue = pseudobulk_df_tissue |> group_by(cell_type_harmonised) |> tar_group()
# pseudobulk_df_tissue |> qs::qsave("pseudobulk_df_tissue_blood.qs")

result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/tidyomics"
store=glue("{result_directory}/_targets_tidyomics")


tar_script({

  result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/tidyomics"

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
      #,
      #script_lines = "module load R/4.2.1",
      #host = "spartan.hpc.unimelb.edu.au"
    )

  small_slurm_memory =
    crew_controller_slurm(
      name = "small_slurm_memory",
      slurm_memory_gigabytes_per_cpu = 50,
      slurm_cpus_per_task = 1,
      workers = 10,
      verbose = T
      #,
      #script_lines = "module load R/4.2.1",
      #host = "spartan.hpc.unimelb.edu.au"
    )

  slurm_3 =
    crew_controller_slurm(
      name = "slurm_3",
      slurm_memory_gigabytes_per_cpu = 20,
      slurm_cpus_per_task = 3,
      workers = 400,
      verbose = T
      #,
      #script_lines = "module load R/4.2.1",
      #host = "spartan.hpc.unimelb.edu.au"
    )

  big_slurm =
    crew_controller_slurm(
      name = "big_slurm",
      slurm_memory_gigabytes_per_cpu = 4,
      slurm_cpus_per_task = 5,
      workers = 300,
      verbose = T
      #,
      #script_lines = "module load R/4.2.1",
      #host = "spartan.hpc.unimelb.edu.au"
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

  split_metadata = function(metadata_DB_path){
    metadata = readRDS(metadata_DB_path) |> mutate(age_days = age_days_original)

    get_metadata() |>
      filter(
        sample_ %in% (
          !!metadata |>
            distinct(sample_) |>
            tidyr::extract(sample_, "sample_", "([a-zA-Z0-9]+)_.+") |>
            pull(sample_)
        )) |>

      # Attach lineage
      left_join(
        read_csv("~/PostDoc/immuneHealthyBodyMap/metadata_cell_type.csv") |>
          replace_na(list(lineage_1 = "other_non_immune")) |>
          mutate(is_immune = as.character(lineage_1 == "immune")),
        by = join_by(cell_type),
        copy = TRUE
      ) |>
      #filter(is_immune == "FALSE" & !is.na(lineage_1)) |>
      filter(!is.na(lineage_1)) |>
      filter(!cell_type_harmonised %in% c("platelet", "immune_unclassified")) |>
      as_tibble() |>

      # Filter rare cell types
      add_count(cell_type_harmonised, name = "count_cell_type_harmonised") |>
      filter(count_cell_type_harmonised>100) |>
      select(-count_cell_type_harmonised) |>

      # Format covatriates
      mutate(assay = assay |> str_replace_all(" ", "_") |> str_replace_all("-", "_")  |> str_remove_all("'")) |>
      mutate(
        ethnicity = case_when(
          ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
          ethnicity |> str_detect("African") ~ "African",
          TRUE ~ ethnicity
        )
      ) |>

      # Mutate days
      filter(development_stage!="unknown") |>

      # Establish the baseline for simplified ethnicity. European as it is the most represented
      # This is so I have a tight intercept term for data simulation
      mutate(ethnicity_simplified = case_when(
        ethnicity %in% c("European", "Chinese", "African", "Hispanic or Latin American") ~ ethnicity,
        TRUE ~ "Other"
      )) |>
      mutate(
        ethnicity_simplified =
          ethnicity_simplified |>
          fct_relevel(c("European", "Chinese", "African", "Hispanic or Latin American", "Other")
          )) |>

      # Establish the baseline for simplified assay
      # Summarise assays to get more stable data simulations
      # 10x as baseline
      mutate(assay_simplified = if_else(assay |> str_detect("10x"), "10x", assay)) |>
      mutate(assay_simplified = factor(assay_simplified)) |>

      # Establish the baseline for disease
      mutate(disease = if_else(disease == "normal", "aaa_normal", disease)) |>

      # Select few columns to make things lighter
      select(
        cell_, cell_type_harmonised, sample_,  file_id, file_id_db,
        age_days, development_stage, sex, tissue_harmonised, disease,
        ethnicity_simplified, assay_simplified, is_immune
      ) |>

      # Cell type for non immune are not summarised ernought I'm loosing a lot of samples
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

      # Split very big datasets with a lot of samples, with maximum 100 samples
      nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_)) |>
      with_groups(
        c(cell_type_harmonised, tissue_harmonised, file_id, is_immune),
        ~ .x |>
          arrange(file_id) |>
          mutate(sample_chunk = row_number() |>

                   # MAX SAMPLES
                   divide_by(100) |>
                   ceiling()
          )
      ) |>
      unnest(data) |>

      # group
      nest(data = -c(cell_type_harmonised, tissue_harmonised, file_id, is_immune, sample_chunk))

  }

  get_sce = 	function(tissue_cell_type_metadata) {

    tissue_cell_type_metadata |>
      mutate(data = pmap(
        list(data, cell_type_harmonised, tissue_harmonised),
        ~ ..1 |>
          get_single_cell_experiment(cache_directory = "/vast/projects/cellxgene_curated") |>
          mutate(sample_se =

                   # I need to fix Curated CellAtlas with disease sample, duplication for
                   # file_id=="cc3ff54f-7587-49ea-b197-1515b6d98c4c", cell_type_harmonised=="stromal_cell"
                   # for lung
                   glue("{sample_}___{disease}___{..2}___{..3}") |>
                   str_replace_all(" ", "_") |>
                   str_replace_all("/", "__")
          )
      ))

  }

  get_pseudobulk = 	function(sce_df) {

    sce_df |>
      mutate(data = map(
        data,
        ~ .x |>
          tidySingleCellExperiment::aggregate_cells( .sample = sample_se	)
      ))

  }

  nest_detect_complete_confounder = function(.data, .col1, .col2){

    .col1 = enquo(.col1)
    .col2 = enquo(.col2)

    .data |>

      nest(se_data = -c(!!.col1, !!.col2)) |>

      # How many ethnicity per assay
      nest(data = -!!.col1) |>
      mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
      unnest(data) |>

      # How many assay per ethnicity
      nest(data = - !!.col2) |>
      mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
      unnest(data)
  }

  left_join_detect_complete_confounder = function(.data, .col1, .col2){

    .col1 = enquo(.col1)
    .col2 = enquo(.col2)

    counts =
      .data |>


      # nest(se_data = -c(!!.col1, !!.col2)) |>
      distinct(!!.col1, !!.col2) |> # Quicker

      # How many ethnicity per assay
      nest(data = -!!.col1) |>
      mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
      unnest(data) |>

      # How many assay per ethnicity
      nest(data = - !!.col2) |>
      mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
      unnest(data)

    .data |>
      left_join(counts) #, by = join_by(!!.col1, !!.col2))
  }

  drop_samples_complete_confounder = function(.data, .col1, .col2){

    .col1 = enquo(.col1)
    .col2 = enquo(.col2)

    .data |>
      left_join_detect_complete_confounder(!!.col1, !!.col2) |>
      filter(n1 + n2 > 2) |>
      select(-n1, -n2)

  }

  samples_NOT_complete_confounders_for_ethnicity_assay = function(se){



    se =
      se |>
      # distinct(assay_simplified, ethnicity_simplified, .sample) |>
      #
      nest(se_data = -c(assay_simplified, ethnicity_simplified)) |>

      # How many ethnicity per assay
      nest(data = -assay_simplified) |>
      mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |>
      unnest(data) |>

      # How many assay per ethnicity
      nest(data = - ethnicity_simplified) |>
      mutate(n2 = map_int(data, ~ .x |> distinct(assay_simplified) |> nrow())) |>
      unnest(data)

    # Replace ethnicity
    dummy_assay = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(assay_simplified)

    se |>
      mutate(assay_simplified = if_else(n1 + n2 < 3, dummy_assay, assay_simplified)) 	|>

      # # Filter
      # filter(!(n1==1 & n2==1)) |>
      select(-n1, -n2) |>

      unnest_summarized_experiment(se_data)
    # |>
    # 	pull(.sample) |>
    # 	unique()
  }

  samples_NOT_complete_confounders_for_ethnicity_disease = function(se){



    se =
      se |>
      #distinct(disease, ethnicity_simplified, .sample) |>

      nest(se_data = -c(disease, ethnicity_simplified)) |>

      # How many ethnicity per assay
      nest(data = -disease) |>
      mutate(n1 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |>
      unnest(data) |>

      # How many assay per ethnicity
      nest(data = - ethnicity_simplified) |>
      mutate(n2 = map_int(data, ~ .x |> distinct(disease) |> nrow())) |>
      unnest(data)


    # Replace ethnicity
    dummy_ethnicity = se |> arrange(desc(n1 + n2)) |> slice(1) |> pull(ethnicity_simplified)

    se |>
      mutate(ethnicity_simplified = if_else(n1 + n2 < 3, dummy_ethnicity, ethnicity_simplified)) 	|>

      # # Filter
      # filter(!(n1==1 & n2==1)) |>
      select(-n1, -n2) |>

      unnest_summarized_experiment(se_data)
    # |>
    # 	pull(.sample) |>
    # 	unique()
  }

  aggregate = function(se_df){

    print("Start aggregate")
    gc()

    se_df |>

      # Add columns and filter
      mutate(data = pmap(
        list(data, cell_type_harmonised, tissue_harmonised, file_id),
        ~ {
          # Add columns
          se =
            ..1 |>
            mutate(cell_type_harmonised = ..2, tissue_harmonised = ..3, file_id = ..4) |>
            select(-any_of(c("file_id_db", ".cell", "original_cell_id")))


          # Identify samples with many genes
          sample_with_many_genes =
            se |>
            assay("counts") |>
            apply(2, function(x) (x>0) |> which() |> length()) |>
            enframe() |>
            mutate(value = as.character(value), name = as.character(name)) |>
            filter(value > 5000) |>
            pull(name)
          se = se[,sample_with_many_genes, drop=FALSE]

          # Filter samples with too few cells
          se = se |> filter(.aggregated_cells > 10)

        },
        .progress=TRUE
      )) |>

      # nest for output
      nest(data = -c(tissue_harmonised, cell_type_harmonised)) |>

      # bind
      mutate(data = map(
        data,
        ~ {

          se = do.call(cbind,  .x |> pull(data))

          # Filter very rare gene-transcripts
          all_zeros = assay(se, "counts") |> rowSums() |> equals(0)
          se = se[!all_zeros,]
          lower_than_total_counts = assay(se, "counts") |> rowSums() < 15
          se = se[!lower_than_total_counts,]

          # Make it sparse
          se@assays@data$counts = as(se@assays@data$counts, "sparseMatrix")

          se
        }
      ))


  }

  map_quantile_scale_abundance = function(se_df){

    print("Start scale abundance")
    gc()

    se_df |>

      # Quantile tranformation
      mutate(data = map(data,	quantile_normalise_abundance, method = "preprocesscore_normalize_quantiles_use_target")) |>

      # convert to sparse again
      mutate(data = map(data, ~ {

        .x@assays@data$counts_scaled = as(.x@assays@data$counts_scaled, "sparseMatrix")

        attr(.x, "internals")$tt_columns$.abundance_scaled |> attr(".Environment") = NULL # new_environment()

        .x
      }))

  }

  map_keep_abundant = function(se_df){

    print("Start keep abundant")
    gc()


    se_df |>
      mutate(data = map(
        data,
        ~ {

          # For blood
          if(ncol(.x) > 1000)
            .x |>
            keep_abundant(
              .abundance = counts_scaled,
              factor_of_interest = c(sex, ethnicity_simplified),
              minimum_counts = 500,
              minimum_proportion = 0.9
            )

          else {
            # Select abundant genes within tissues and unite

            # DEDUPLLICATE
            # I have to investigate e80ca11e47c3f5db861058eabac363c5___aaa_normal___hematopoietic_cell___blood
            .x = .x[,!colnames(.x) |> duplicated()]

            abundant_genes =
              .x |>
              nest(data = -cell_type_harmonised) |>

              # Filter if only one sex
              mutate(abundant_genes = map(
                data,
                ~ {

                  se = .x

                  # # Avoid indeterminability
                  # if(ncol(se) > 0) {
                  #
                  # 	ethnicity_to_keep = se |> pivot_sample() |> count(ethnicity_simplified) |> filter(n >1) |> pull(ethnicity_simplified)
                  # 	se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)
                  #
                  # 	sex_to_keep = se |> pivot_sample() |> count(sex) |> filter(n >1) |> pull(sex)
                  # 	se = se |> filter(sex %in% sex_to_keep)
                  #
                  # 	se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
                  #
                  # }
                  #
                  # factors =
                  # 	c( "sex", "ethnicity_simplified") |>
                  # 	enframe(value = "factor") |>
                  # 	mutate(n = map_int(
                  # 		factor, ~ se |> select(.x) |> distinct() |> nrow()
                  # 	)) |>
                  # 	filter(n>1) |>
                  # 	pull(factor) |>
                  # 	map(sym) |>
                  # 	unlist()

                  se |>
                    keep_abundant(.abundance = counts_scaled, factor_of_interest = c(sex, ethnicity_simplified), minimum_counts = 50) |>
                    rownames()

                },
                .progress = TRUE
              )) |>
              pull(abundant_genes) |>
              unlist() |>
              unique()

            .x |> filter(.feature %in% abundant_genes)
          }

        }
      ))
    # se_df |>
    # 	mutate(data = map(
    # 		data,
    # 		~ {
    #
    # 			# Select abundant genes within tissues and unite
    # 			abundant_genes =
    # 				.x |>
    # 				nest(data = -cell_type_harmonised) |>
    #
    # 				# Filter if only one sex
    # 				mutate(abundant_genes = map(
    # 					data,
    # 					~ {
    #
    # 						se = .x
    #
    # 						# # Avoid indeterminability
    # 						# if(ncol(se) > 0) {
    # 						#
    # 						# 	ethnicity_to_keep = se |> pivot_sample() |> count(ethnicity_simplified) |> filter(n >1) |> pull(ethnicity_simplified)
    # 						# 	se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)
    # 						#
    # 						# 	sex_to_keep = se |> pivot_sample() |> count(sex) |> filter(n >1) |> pull(sex)
    # 						# 	se = se |> filter(sex %in% sex_to_keep)
    # 						#
    # 						# 	se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
    # 						#
    # 						# }
    # 						#
    # 						# factors =
    # 						# 	c( "sex", "ethnicity_simplified") |>
    # 						# 	enframe(value = "factor") |>
    # 						# 	mutate(n = map_int(
    # 						# 		factor, ~ se |> select(.x) |> distinct() |> nrow()
    # 						# 	)) |>
    # 						# 	filter(n>1) |>
    # 						# 	pull(factor) |>
    # 						# 	map(sym) |>
    # 						# 	unlist()
    #
    # 						se |>
    # 						keep_abundant(.abundance = counts_scaled, factor_of_interest = c(sex, ethnicity_simplified), minimum_counts = 50) |>
    # 						rownames()
    #
    # 						},
    # 					.progress = TRUE
    # 				)) |>
    # 				pull(abundant_genes) |>
    # 				unlist() |>
    # 				unique()
    #
    # 			.x |> filter(.feature %in% abundant_genes)
    #
    # 		}
    # 	))
  }

  se_add_dispersion = function(se_df){

    print("Start add dispersion")
    gc()

    se_df |>
      mutate(data = map(
        data,
        ~ {

          # Because I have nested map
          se = .x

          # Avoid indeterminability
          if(ncol(se) > 0) {

            ethnicity_to_keep = se |> pivot_sample() |> count(ethnicity_simplified) |> filter(n >1) |> pull(ethnicity_simplified)
            se = se |> filter(ethnicity_simplified %in% ethnicity_to_keep)

            sex_to_keep = se |> pivot_sample() |> count(sex) |> filter(n >1) |> pull(sex)
            se = se |> filter(sex %in% sex_to_keep)
          }
          if(ncol(se) > 0) {
            se = se |> drop_samples_complete_confounder(sex, ethnicity_simplified)
          }

          # If SE empty add dummy dispersion
          if(ncol(se) == 0) rowData(se)$dispersion = rep(NA, nrow(se))
          else if(ncol(se)<1000) {

            factors =
              c( "sex", "ethnicity_simplified") |>
              enframe(value = "factor") |>
              mutate(n = map_int(
                factor, ~ se |> select(.x) |> distinct() |> nrow()
              )) |>
              filter(n>1) |>
              pull(factor) |>
              str_c(collapse = " + ")

            my_design = glue("~ {factors}") |> as.formula() |> model.matrix( data = colData(se) |> droplevels())
            rowData(se)$dispersion =  se |> assay("counts_scaled") |> estimateDisp(design = my_design) %$% tagwise.dispersion
          }
          else {

            sampled_samples = sample(seq_len(ncol(se)), size = min(ncol(se), 2000))

            factors =
              c( "sex", "ethnicity_simplified") |>
              enframe(value = "factor") |>
              mutate(n = map_int(
                factor, ~ se[,sampled_samples, drop=FALSE]  |> select(.x) |> distinct() |> nrow()
              )) |>
              filter(n>1) |>
              pull(factor) |>
              str_c(collapse = " + ")

            my_design = model.matrix(~ sex + ethnicity_simplified, data = se[,sampled_samples, drop=FALSE] |> colData() |> droplevels())

            rowData(se)$dispersion =
              assay(se[,sampled_samples, drop=FALSE] , "counts_scaled") |>
              estimateTrendedDisp(design = my_design, subset=1000, rowsum.filter=10)
          }

          se
        }
      ))

  }

  add_number_of_gene_chunks =	function(se_df){

    print("Start add number of chunks")
    gc()

    se_df |>
      mutate(number_of_chunks = map_int(
        data,
        ~ .x |>
          distinct(sample_se) |>
          nrow() |>

          # Avoid 0 because it will cause error
          # And because ifI have 0 samples I still have one chunk of genes
          max(1) |>
          divide_by(50) |>
          multiply_by(10) |>
          ceiling() |>

          # Otherwise it takes more than 20 minutes
          min(10000)
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

    # tar_target(
    #   sce_df_split_by_gene_grouped_list,
    #   sce_df_split_by_gene_grouped |> list(),
    #   pattern = map(sce_df_split_by_gene_grouped),
    #   iteration = "list",
    #   deployment = "main",
    #   memory = "persistent",
    #   storage = "main",
    #   retrieval = "main",
    # ) ,
    #
    #
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



library(furrr)
plan(multisession, workers = 18)
options(future.globals.maxSize = 200000 * 1024^2)
library(tidybulk)

blood =
  tar_meta(  store=store  ) |>
  filter(name |> str_detect("estimates_sex_tissue")) |>
  filter(!is.na(data)) |>
  filter(name != "estimates_sex_tissue") |>
  mutate(se = future_map(
    name,
    ~{
      #print(.x)
      .x |>
        tar_read_raw(store=store ) |>
        mutate(data = map(data, pivot_transcript))
    },
    .progress = TRUE
  ))

de =
  blood |>
  select(se) |>
  unnest(se) |>
  unnest(data) |>
  filter(cell_type_harmonised != "thymocyte")

de |> saveRDS("de_blood.rds")

source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/f0b6bf9f59847c8b9f0a638262a6b8dd697affb7/color_cell_types.R")
cell_type_color =
  de |>
  pull(cell_type_harmonised) |>
  unique() |>
  get_cell_type_color()

de =
  de |>
  with_groups(cell_type_harmonised, ~ .x |> mutate(
    P_sex_adjusted = p.adjust(P_sex, method = "BH"),
    P_age_days_adjusted = p.adjust(P_age_days, method = "BH"),
    P_age_days.sex_adjusted = p.adjust(P_age_days.sex, method = "BH")
  ))

df_sex_cell_type_most_de =
  de |>
  filter(!is.na(P_sex_adjusted)) |>
  add_count(cell_type_harmonised, name = "gene_number") |>
  mutate(de_only_in_sex = P_sex_adjusted < 0.05 |	P_age_days.sex_adjusted < 0.05 & P_age_days_adjusted > 0.05) |>
  dplyr::count(
    cell_type_harmonised, gene_number,
    de_only_in_sex
  ) |>
  mutate(proportion_of_significant = n/gene_number) |>
  filter(de_only_in_sex)

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

library(ggupset)
plot_sex_cell_type_upset =
  de |>
  filter(!is.na(P_sex_adjusted)) |>
  filter( P_sex_adjusted < 0.05 |	P_age_days.sex_adjusted < 0.05 & P_age_days_adjusted > 0.05) |>
  filter(cell_type_harmonised %in% (df_sex_cell_type_most_de |> arrange(desc(proportion_of_significant)) |> head(9) |> pull(cell_type_harmonised))) |>
  select(cell_type_harmonised, .feature) |>
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("terminal effector", "eff")) |>
  nest(cell_types = cell_type_harmonised) |>
  mutate(cell_types = map(cell_types, pull, cell_type_harmonised)) |>
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

gene_chr = read_csv("~/PostDoc/immuneHealthyBodyMap/symbol_chr.csv")

proportion_of_interaction_only =
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
  ggplot(aes(label, proportion_of_significant, fill=label)) +
  geom_boxplot(outlier.shape = NA, lwd=0.3, fatten=0.3) +
  geom_jitter(width = 0.2, size=0.2) +
  scale_fill_brewer(palette="Set2") +
  xlab("Proportion of significant genes") +
  guides(fill="none") +
  coord_flip() +
  theme_multipanel +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(angle=90, hjust = 0.5))

# Contribution of age-interaction in number of sex-deendent genes
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
  mutate(label = if_else(de_only_in_sex, "de_only_in_sex", "de_only_in_interaction")) |> select(cell_type_harmonised, proportion_of_significant, label) |> pivot_wider(names_from = label , values_from = proportion_of_significant) |> mutate(contribution = de_only_in_interaction / (de_only_in_sex + de_only_in_interaction) ) |> pull(contribution) |> mean(na.rm=TRUE)

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
