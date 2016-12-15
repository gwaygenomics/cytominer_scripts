#!/usr/bin/env Rscript

'annotate

Usage: 
  annotate.R -b <id> -p <id> [-d -j <file>]

Options:
  -b <id>, --batch_id=<id>                 Batch ID.
  -p <id>, --plate_id=<id>                 Plate ID.
  -d, --format_broad_cmap                  Add columns to make compatible with Broad CMap naming conventions.
  -j <file>, --external_metadata=<file>    External metadata to join with.' -> doc

suppressWarnings(suppressMessages(library(docopt)))

suppressWarnings(suppressMessages(library(dplyr)))

suppressWarnings(suppressMessages(library(magrittr)))

opts <- docopt(doc)

batch_id <- opts[["batch_id"]]

external_metadata <- opts[["external_metadata"]]

format_broad_cmap <- opts[["format_broad_cmap"]]

plate_id <- opts[["plate_id"]]

metadata_dir <- paste("../..", "metadata", batch_id, sep = "/")

backend_dir <- paste("../..", "backend", batch_id, plate_id, sep = "/")

# read profiles and rename column names 

profiles <- suppressMessages(readr::read_csv(paste(backend_dir, paste0(plate_id, ".csv"), sep = "/")))

profiles %<>% setNames(names(profiles) %>% stringr::str_replace_all("^Image_Metadata", "Metadata"))

# read and join metadata map

metadata_map <- suppressMessages(readr::read_csv(paste(metadata_dir, "barcode_platemap.csv", sep = "/")))

testthat::expect_true("Assay_Plate_Barcode" %in% colnames(metadata_map))

metadata_map %<>% setNames(names(metadata_map) %>% stringr::str_replace_all("^", "Metadata_"))

profiles %<>% mutate(Metadata_Assay_Plate_Barcode = Metadata_Plate)

profiles %<>% inner_join(metadata_map, by = c("Metadata_Assay_Plate_Barcode"))

# read and join platemap

platemap_name <- profiles %>% select(Metadata_Plate_Map_Name) %>% distinct() %>% extract2("Metadata_Plate_Map_Name")

testthat::expect_equal(length(platemap_name), 1)

platemap <- suppressMessages(readr::read_tsv(paste(metadata_dir, "platemap", paste0(platemap_name, ".txt"), sep = "/")))

testthat::expect_true("well_position" %in% colnames(platemap))

cnames <- colnames(platemap)

cnames %<>% stringr::str_replace_all("^", "Metadata_")

names(platemap) <- cnames

profiles %<>% mutate(Metadata_well_position = Metadata_Well)

profiles %<>% inner_join(platemap, by = c("Metadata_well_position"))

# format_broad_cmap
if (format_broad_cmap) {
    profiles %<>% 
        mutate(Metadata_broad_sample_type = ifelse(is.na(Metadata_broad_sample), "control", "trt"),
               Metadata_broad_sample = ifelse(Metadata_broad_sample_type =="control", "DMSO", Metadata_broad_sample),
               Metadata_mg_per_ml = ifelse(Metadata_broad_sample_type =="control", 0, Metadata_mg_per_ml),
               Metadata_mmoles_per_liter = ifelse(Metadata_broad_sample_type =="control", 0, Metadata_mmoles_per_liter),
               Metadata_pert_id = stringr::str_extract(Metadata_broad_sample, "(BRD-[A-Z0-9]+)"),
               Metadata_pert_mfc_id = Metadata_broad_sample,
               Metadata_pert_well = Metadata_Well,
               Metadata_pert_vehicle = Metadata_solvent,
               Metadata_pert_idose = Metadata_mmoles_per_liter,
               Metadata_pert_id_vendor = "")
}

# external_metadata
if(!is.null(external_metadata)) {
    external_metadata_df <- suppressMessages(readr::read_csv(external_metadata))

    # Check whether the columns have "Metadata" prefix; if not, assume that all columsn need the suffix
    if (length(grep("Metadata_", colnames(external_metadata_df))) == 0) {
        external_metadata_df %<>% setNames(names(external_metadata_df) %>% stringr::str_replace_all("^", "Metadata_"))

    }
    
    profiles %<>% 
        left_join(
            external_metadata_df %>% 
            distinct()
            )
}


# format_broad_cmap: columns that may be added after joining with external metadata
if (format_broad_cmap) {

    if ("Metadata_pert_iname" %in% colnames(profiles)) {
        profiles %<>% 
            mutate(Metadata_pert_mfc_desc = Metadata_pert_iname)

    }

}


# save 

profiles_augmented <- paste(backend_dir, paste0(plate_id, "_augmented.csv"), sep = "/")

metadata_cols <- stringr::str_subset(names(profiles), "^Metadata_")

feature_cols <- stringr::str_subset(names(profiles), "^Cells_|^Cytoplasm_|^Nuclei_")

all_cols <- c(metadata_cols, feature_cols)

profiles[all_cols] %>% readr::write_csv(profiles_augmented)