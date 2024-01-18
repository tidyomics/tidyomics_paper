library(targets)
library(tidyverse)
library(targets)
library(glue)
library(tidySummarizedExperiment)
library(magrittr)
library(tidybulk)
library(tidyseurat)
library(Seurat)
# Get input


result_directory = "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab_projects/people/mangiola.s/PostDoc/tidyomics"
store=glue("{result_directory}/_targets_tidyomics")
pseudobulk = tar_read(pseudobulk_df_blood_cell_type, store=store)


aggregate_no_filter = function(se_df){
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
				# Make it sparse
				se@assays@data$counts = as(se@assays@data$counts, "sparseMatrix")
				se
			}
		))
}
pseudobulk =
	pseudobulk  |>
	group_split(cell_type_harmonised) |>
	map(	aggregate_no_filter ,	.progress = TRUE	)
pseudobulk = do.call(cbind, pseudobulk |> bind_rows() |> pull(data))

# I have one sample that is duplicated. Need to investigate,
# but among 28K samples it is not a big deal
pseudobulk = pseudobulk[,!pseudobulk |> colnames() |> duplicated()]

# Filter very rare gene-transcripts
all_zeros = assay(pseudobulk, "counts") |> rowSums() |> equals(0)
pseudobulk = pseudobulk[!all_zeros,]
lower_than_total_counts = assay(pseudobulk, "counts") |> rowSums() < 500
pseudobulk = pseudobulk[!lower_than_total_counts,]

# Normalise quantile as the data property is very heterogeneous
pseudobulk = pseudobulk |> quantile_normalise_abundance(method = "preprocesscore_normalize_quantiles_use_target")
attr(pseudobulk, "internals")$tt_columns$.abundance_scaled |> attr(".Environment") = NULL # new_environment()

# Keep abundant genes
pseudobulk =
	pseudobulk |>
	keep_abundant(
		.abundance = counts_scaled,
		factor_of_interest = c(sex, cell_type_harmonised),
		minimum_counts = 500,
		minimum_proportion = 0.9
	)
# pseudobulk |> qs::qsave("pseudobulk_cell_type_blood.qs")

pseudobulk = qs::qread("pseudobulk_cell_type_blood.qs")

# Use Seurat for integration a plotting
seu = pseudobulk |> as("SingleCellExperiment") |>  as.Seurat(data = NULL)
seu <- seu |> SplitObject( split.by = "file_id")
seu  <- seu[map_int(seu, ncol)>50] |>
	map(~ .x |>
				NormalizeData() |>
				FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
	)
features <- SelectIntegrationFeatures(object.list = seu, nfeatures = 100)
seu |> qs::qsave("pseudobulk_cell_type_blood_features.qs")
anchors <-
	seu |>
	FindIntegrationAnchors( anchor.features = features, dims = 1:10, k.filter = 10) #, reference = 1139, dims = 1:5,   k.filter = 5, reduction = "rpca")

seu <- IntegrateData(anchorset = anchors, k.weight = 10)
DefaultAssay(seu) <- "integrated"

# Run the standard workflow for visualization and clustering
seu <-
	seu |>
	ScaleData( verbose = FALSE) |>
	RunPCA( npcs = 30, verbose = FALSE) |>
	RunUMAP( reduction = "pca", dims = 1:10)

seu |> qs::qsave("pseudobulk_cell_type_blood_seurat.qs")

seu = qs::qread("pseudobulk_cell_type_blood_seurat.qs")

# Load cell type colors
source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/a175f7d0fe95ce663a440ecab0023ca4933e5ab8/color_cell_types.R")
cell_type_color =
	seu |>
	pull(cell_type_harmonised) |>
	unique() |>
	get_cell_type_color()

plot_umap =
	seu |>
	DimPlot(group.by = "cell_type_harmonised", cols = cell_type_color, label = TRUE) + NoLegend() +
	theme(
		axis.title.x = element_blank(), axis.title.y = element_blank(),
		axis.text.x = element_blank(), axis.text.y = element_blank(),
		axis.ticks.x = element_blank(), 	axis.ticks.y = element_blank()
	)

plot_umap |> qs::qsave("plot_umap.qs")

plot_umap = qs::qread("plot_umap.qs")

ggsave(plot_umap, file = "plot_umap.pdf", width = 150, height = 120, units = "mm", device = "pdf")

# Files for Michael

# Save pseudobulk dataset
library(SummarizedExperiment)
library(tidySummarizedExperiment)
pseudobulk_with_stats = qs::qread("pseudobulk_cell_type_blood.qs")
rowData(pseudobulk_with_stats)  = rowData(pseudobulk_with_stats) |> as_tibble() |> left_join(
  de |>
    select(.feature, cell_type_harmonised, sexmale, age_days, age_days.sexmale, P_sex_adjusted ,  P_age_days.sex_adjusted ,   P_age_days_adjusted ) |>
    pivot_wider(names_from = cell_type_harmonised, values_from = -c(.feature, cell_type_harmonised), names_sep = "___"),
  join_by(feature == .feature)
) |> as("DataFrame")


pseudobulk_only_cd4_naive = pseudobulk_with_stats[,colData(pseudobulk_with_stats)$cell_type_harmonised == "cd4 naive"]
rowData(pseudobulk_only_cd4_naive) = rowData(pseudobulk_only_cd4_naive) |> as_tibble() |> select(feature, .abundant, contains("cd4.naive")) |> as("DataFrame")

pseudobulk_only_cd4_naive |>
  saveRDS("pseudobulk_for_michael_cd4_naive.rds", compress = "xz")


