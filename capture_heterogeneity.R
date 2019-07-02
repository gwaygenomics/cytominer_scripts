#!/usr/bin/env Rscript

'capture_heterogeneity

Usage:
  capture_heterogeneity.R -s <sqlite_file> -n <num_comp> -o <file>

Options:
  -h --help                                       Show this screen.
  -s <sqlite_filte> --sqlite_file=<sqlite_file>   Location of the single cell data
  -n <num_comp> --num_components=<num_comp>       How many components to randomly project [default: 3000]
  -o <file> --output=<file>                       per-well aggregated data.' -> doc

suppressWarnings(suppressMessages(library(docopt)))

suppressWarnings(suppressMessages(library(dplyr)))

suppressWarnings(suppressMessages(library(magrittr)))

suppressWarnings(suppressMessages(library(stringr)))

opts <- docopt(doc)

sqlite_file <- opts[["sqlite_file"]]
num_components <- opts[["num_comp"]]
output_file <- opts[["file"]]

db <- dplyr::src_sqlite(path = sqlite_file)

image_df <- dplyr::tbl(src = db, "image") %>%
  dplyr::select(TableNumber, ImageNumber, Image_Metadata_Plate, Image_Metadata_Well)

load_compartment <- function(db, image, compartment) {
  object_df <- dplyr::tbl(src = db, compartment)

  object_df %<>% dplyr::inner_join(image, by = c("TableNumber", "ImageNumber"))

  return(object_df)
}

# Merge all compartments of sc profiles together
merge_cols <- c("TableNumber",
                "ImageNumber",
                "ObjectNumber",
                "Image_Metadata_Plate",
                "Image_Metadata_Well")

sc_profiles_df <-
  load_compartment(db = db,
                   image = image_df,
                   compartment = "Cells") %>%
  inner_join(
    load_compartment(db = db,
                     image = image_df,
                     compartment = "Cytoplasm"),
    by = merge_cols
  ) %>%
  inner_join(
    load_compartment(db = db,
                     image = image_df,
                     compartment = "Nuclei"),
    by = merge_cols
  ) %>% dplyr::collect()

# Prepare data for input to covariance
futile.logger::flog.info(str_c("Calculating covariance"))

cp_features <-
    colnames(sc_profiles_df) %>%
    stringr::str_subset("^Nuclei_|^Cells_|^Cytoplasm_")

# Apply covariance function
sc_profile_cov_df <- cytominer::covariance(population = sc_profiles_df,
                                           variables = cp_features)

# Apply sparse random projection function
futile.logger::flog.info(str_c("Performing random projection"))
cp_cov_features <- colnames(sc_profile_cov_df)

random_projection_df <-
  cytominer::sparse_random_projection(population = sc_profile_cov_df,
                                      variables = cp_cov_features,
                                      n_components = num_components)

futile.logger::flog.info(paste0("Writing aggregated to ", output_file))

random_projection_df %>% readr::write_csv(output_file)
