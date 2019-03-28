# NOTE: This script is used to generate Zonation setups for public release.
# See XXXXXX for a more detailed explanation.
#
# NOTE: you will need the latest version for this to work
# zonator > 0.5.0
# devtools::install_github("cbig/zonator")
library(dplyr)
library(raster)
library(readr)
library(readxl)
library(zonator)
library(rgdal)

source("R/00_lib/utils.R")

build_zonation_params = list()

# Generate variants for all taxa ------------------------------------------

# Define names for variants. "[ID]" is a placeholder for running id, "[TX]" is
# for taxon codes.

build_zonation_params$variant_templates = c("[TX]_caz",
                                            "[TX]_abf",
                                            "[TX]_caz_wgt",
                                            "[TX]_abf_wgt",
                                            "[TX]_caz_wgt_con",
                                            "[TX]_abf_wgt_con",
                                            "[TX]_caz_wgt_con_hm3",
                                            "[TX]_abf_wgt_con_hm3")

build_zonation_params$data_dir <- paste0(path.expand('~'), "/workdir/species_data/zenodo/")
build_zonation_params$zsetup_root <- "zsetup"
build_zonation_params$feature_data_file <- paste0(build_zonation_params$data_dir, "species_params/biodiv_features.csv")

build_zonation_params$agg_weights_file <- paste0(build_zonation_params$data_dir, "species_params/aggregate_weights.csv")
build_zonation_params$ppa_raster_file <- paste0(build_zonation_params$data_dir, "group_layers/prefecture.tif")
build_zonation_params$ppa_cmd_string <- paste(c("LSM", build_zonation_params$ppa_raster_file, 0, -1, 0), collapse = " ")

build_zonation_params$condition_raster_file <- paste0(build_zonation_params$data_dir, "group_layers/hii_rescaled_simple.tif")
# Hierarchical masks
build_zonation_params$hm3_raster_file <- paste0(build_zonation_params$data_dir, "group_layers/PA_3_levels_simple.tif")

build_zonation_params$w_field <- "weight"
build_zonation_params$condition_config_file <- "condition_config.txt"
build_zonation_params$ppa_config_file <- "ppa_config.txt"
build_zonation_params$use_groups = rep(1, 100)
build_zonation_params$cell_removal_rule = rep(c(1, 2), 50)

# Make a list to hold information on the different taxa names, short codes
# and Excel sheet names (sheet names not actually needed).

# build_zonation_params$species_data_template <- list("amphibians" = list("code" = "amp", "sheet" = "amphibians"),
#                                                 "birds" = list("code" = "bir", "sheet" = "birds"),
#                                                 "freshwater_fish" = list("code" = "frf", "sheet" = "freshwater_fish"),
#                                                 "mammals" = list("code" = "mam", "sheet" = "mammals"),
#                                                 "plants" = list("code" = "pla", "sheet" = "plants"),
#                                                 "reptiles" = list("code" = "rep", "sheet" = "reptiles"))

build_zonation_params$species_group_names = c('amphibians', 'birds', 'freshwater_fish', 'mammals', 'plants', 'reptiles')
build_zonation_params$species_group_codes = c('amp', 'bir', 'frf', 'mam', 'pla', 'rep')
build_zonation_params$species_data_template = setNames(lapply(seq_along(build_zonation_params$species_group_names),
                                                              function(species_group_ind) setNames(list(build_zonation_params$species_group_codes[species_group_ind], 
                                                                                                        build_zonation_params$species_group_names[species_group_ind]),
                                                                                                   c('code', 'sheet'))),
                                                              build_zonation_params$species_group_names)

build_zonation_params$variant_cycle = rep(1:(length(build_zonation_params$species_data_template) + 1), 
                                          each = length(build_zonation_params$variant_templates)) 
build_zonation_params$z_control_template = "templates/run_Z_template.sh"
build_zonation_params$dat_template_file = "templates/template.dat"
build_zonation_params$z_control_file = "run_biodiversa_source.sh"


species_data <- as.data.frame(readr::read_csv(build_zonation_params$feature_data_file))

# Read in aggregate weights
agg_weights <- readr::read_csv(build_zonation_params$agg_weights_file) 

# Set up the project
zproject <- initiate_zproject(zsetup_root = build_zonation_params$zsetup_root,
                                    build_zonation_params$variant_templates,
                                    spp_data = build_zonation_params$species_data_template,
                                    data_dir = build_zonation_params$data_dir,
                                    prefix_spp_paths = "../..", 
                                    dat_template_file = build_zonation_params$dat_template_file)

# Set run configuration parameters ----------------------------------------

# In principle, the preprocessing of the same type of variant for each taxon
# (e.g. 01_amp_caz, 06_bir_caz, 11_frf_caz etc.) is the same. Thus  let's
# automate it a bit.

# Use a counting variable to identify variant IDs.


lapply(seq_along(build_zonation_params$species_data_template), 
       function(group_ind) write(build_zonation_params$ppa_cmd_string,
                                 file.path(build_zonation_params$zsetup_root, 
                                           build_zonation_params$species_data_template[[group_ind]]$sheet, 
                                           build_zonation_params$ppa_config_file)))

variants <- lapply(seq_along(zproject@variants), 
                   function(variant_id) get_variant(zproject, 
                                                    variant_id))

variants <- lapply(seq_along(zproject@variants), 
                   function(variant_id) set_dat_param(variants[[variant_id]], 
                                                      "use groups", 
                                                      build_zonation_params$use_groups[variant_id])) 

variants <- lapply(seq_along(zproject@variants), 
                   function(variant_id) set_dat_param(variants[[variant_id]], 
                                                      "groups file", 
                                                      file.path(variants[[variant_id]]@name, paste0(variants[[variant_id]]@name, "_groups.txt")))) 

variants <- lapply(seq_along(zproject@variants), 
                   function(variant_id) set_dat_param(variants[[variant_id]], 
                                                      "post-processing list file",
                                                      build_zonation_params$ppa_config_file))

variants <- lapply(seq_along(zproject@variants), 
                   function(variant_id) set_dat_param(variants[[variant_id]], 
                                                      "removal rule",
                                                      build_zonation_params$cell_removal_rule[variant_id] ))

# lapply(seq_along(zproject@variants), 
#        function(variant_id) save_zvariant(variants[[variant_id]], 
#                                           dir = file.path(build_zonation_params$zsetup_root, 
#                                                           build_zonation_params$species_data_template[[build_zonation_params$variant_cycle[group_ind]]]$sheet),
#                                           overwrite = TRUE, 
#                                           debug_msg = FALSE))

z_control_files = lapply(seq_along(zproject@variants), 
                         function(variant_id) create_sh_file(variants[[variant_id]], 
                                                             remove_bat = TRUE))

build_zonation_params$control_template = "templates/run_Z_template.sh"
build_zonation_params$dat_template_file = "templates/template.dat"
build_zonation_params$Z_sh_runfile = "Z_source.sh"

z_template <- readLines(build_zonation_params$z_control_template)

file_ID <- file(build_zonation_params$z_control_file)
writeLines(c(z_template, paste0("/bin/sh ", z_control_files)), file_ID)
close(file_ID)
Sys.chmod(build_zonation_params$z_control_file)




#   # rl_groups <- lapply(seq_along(variants), function(variant_id) lookup_rl_group(species_data, featurenames(variants[[variant_id]])))
# aa = lapply(seq_along(variants), function(variant_id) groups(variants[[variant_id]]) <- as.vector(rl_groups[[variant_id]]))
#   # 
# variant_id <- 1
# 
# for (taxon in spp_data) {
# 
#   variant <- get_variant(zproject, variant_id)
# 
#   # Generate groups
#   # Parse groups based on the Red-list status
#   groups <- lookup_rl_group(species_data,  featurenames(variant))
# 
#   # [XX]_[TX]_caz -------------------------------------------------------------
# 
#   message("Editing variant: ", variant@name)
#   # Set groups
#   groups(variant) <- as.vector(groups)
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set post-processing (LSM). First, let's create the file itself (zonator
#   # can't handle this yet). The file needs to be created only once per raxon
#   # since all the variants can use the same file.
#   ppa_file_name <- file.path(build_zonation_params$zsetup_root, taxon$sheet, build_zonation_params$ppa_config_file)
#   ppa_cmd_string <- paste(c("LSM", build_zonation_params$ppa_raster_file, 0, -1, 0), collapse = " ")
#   write(ppa_cmd_string, ppa_file_name)
#   # Need to define ppa_config.txt relative to the bat-file (same dir)-
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_abf -------------------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Manipulate the spp data
#   variant_spp_data <- sppdata(variant)
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with caz. NOTE: This must be done after setting
#   # the spp data
#   
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 2)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_caz_wgt ---------------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Set weights
#   wgts <- lookup_weight(species_data,  featurenames(variant), w_field = build_zonation_params$w_field)
#   # sppweights()<- not implemented yet.
#   variant_spp_data <- sppdata(variant)
#   variant_spp_data$weight <- wgts
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with previous
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 1)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_abf_wgt ---------------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Set weights
#   wgts <- lookup_weight(species_data,  featurenames(variant), build_zonation_params$w_field)
#   # sppweights()<- not implemented yet.
#   variant_spp_data <- sppdata(variant)
#   variant_spp_data$weight <- wgts
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with previous
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 2)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_caz_wgt_con -----------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Set weights
#   wgts <- lookup_weight(species_data,  featurenames(variant), w_field = build_zonation_params$w_field)
#   # sppweights()<- not implemented yet.
#   variant_spp_data <- sppdata(variant)
#   variant_spp_data$weight <- wgts
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with previous
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 1)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Set condition. First, let's create the file itself (zonator
#   # can't handle this yet). The file needs to be created only once per raxon
#   # since all the variants can use the same file.
#   condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon$sheet,
#                                    build_zonation_params$condition_config_file)
#   condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
#   write(condition_cmd_string, condition_file_name)
#   # Need to define ppa_config.txt relative to the bat-file (same dir)-
#   variant <- set_dat_param(variant, "use condition layer", 1)
#   variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
#   # We also need to (manually) set condition reference in the groups data
#   variant@groups$condition.group <- 1
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_abf_wgt_con -----------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Set weights
#   wgts <- lookup_weight(species_data,  featurenames(variant), w_field = build_zonation_params$w_field)
# 
#   # sppweights()<- not implemented yet.
#   variant_spp_data <- sppdata(variant)
#   variant_spp_data$weight <- wgts
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with previous
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 2)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Set condition. First, let's create the file itself (zonator
#   # can't handle this yet). The file needs to be created only once per raxon
#   # since all the variants can use the same file.
#   condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon$sheet,
#                                    build_zonation_params$condition_config_file)
#   condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
#   write(condition_cmd_string, condition_file_name)
#   # Need to define ppa_config.txt relative to the bat-file (same dir)-
#   variant <- set_dat_param(variant, "use condition layer", 1)
#   variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
#   # We also need to (manually) set condition reference in the groups data
#   variant@groups$condition.group <- 1
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_caz_wgt_con_hm3 -------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Set weights
#   wgts <- lookup_weight(species_data,  featurenames(variant), w_field = build_zonation_params$w_field)
#   # sppweights()<- not implemented yet.
#   variant_spp_data <- sppdata(variant)
#   variant_spp_data$weight <- wgts
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with previous
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 1)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Set condition. First, let's create the file itself (zonator
#   # can't handle this yet). The file needs to be created only once per raxon
#   # since all the variants can use the same file.
#   condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon$sheet,
#                                    build_zonation_params$condition_config_file)
#   condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
#   write(condition_cmd_string, condition_file_name)
#   # Need to define ppa_config.txt relative to the bat-file (same dir)-
#   variant <- set_dat_param(variant, "use condition layer", 1)
#   variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
#   # We also need to (manually) set condition reference in the groups data
#   variant@groups$condition.group <- 1
# 
#   # Set hierarchical mask hm2_raster_file
#   variant <- set_dat_param(variant, "use mask", 1)
#   variant <- set_dat_param(variant, "mask file", build_zonation_params$hm3_raster_file)
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   # [XX]_[TX]_abf_wgt_con_hm3 -------------------------------------------------
# 
#   variant_id <- variant_id + 1
#   variant <- get_variant(zproject, variant_id)
#   message("Editing variant: ", variant@name)
# 
#   # Set weights
#   wgts <- lookup_weight(species_data,  featurenames(variant), w_field = build_zonation_params$w_field)
#   # sppweights()<- not implemented yet.
#   variant_spp_data <- sppdata(variant)
#   variant_spp_data$weight <- wgts
#   sppdata(variant) <- variant_spp_data
# 
#   # Set groups, use the same as with previous
#   groups(variant) <- as.vector(groups)
# 
#   # Set cell removal rule
#   variant <- set_dat_param(variant, "removal rule", 2)
# 
#   # Set groups use and groups file
#   variant <- set_dat_param(variant, "use groups", 1)
#   # Note that groups file location is always relative to the bat file
#   groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
#   variant <- set_dat_param(variant, "groups file", groups_file)
# 
#   # Set PPA file
#   variant <- set_dat_param(variant, "post-processing list file",
#                            build_zonation_params$ppa_config_file)
# 
#   # Set condition. First, let's create the file itself (zonator
#   # can't handle this yet). The file needs to be created only once per raxon
#   # since all the variants can use the same file.
#   condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon$sheet,
#                                    build_zonation_params$condition_config_file)
#   condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
#   write(condition_cmd_string, condition_file_name)
#   # Need to define ppa_config.txt relative to the bat-file (same dir)-
#   variant <- set_dat_param(variant, "use condition layer", 1)
#   variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
#   # We also need to (manually) set condition reference in the groups data
#   variant@groups$condition.group <- 1
# 
#   # Set hierarchical mask hm2_raster_file
#   variant <- set_dat_param(variant, "use mask", 1)
#   variant <- set_dat_param(variant, "mask file", build_zonation_params$hm3_raster_file)
# 
#   # Save variant
#   save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon$sheet),
#                 overwrite = TRUE, debug_msg = FALSE)
#   create_sh_file(variant)
# 
#   variant_id <- variant_id + 1
#   
# }

# taxon <- "taxa_all"
# 
# # Generate variants for all species together -------------------------------
# 
# variant <- get_variant(zproject, variant_id)
# 
# # Generate groups
# # Parse groups based on the Red-list status
# groups <- lookup_rl_group(species_data,  featurenames(variant))
# # We need to manually fix few things in groups
# #groups["Salamandrella_keyserlingii_comp"] <- groups["Salamandrella_keyserlingii"]
# 
# groups[-is.na(groups)] <- 1
# # [XX]_caz
# 
# message("Editing variant: ", variant@name)
# # Set groups
# groups(variant) <- as.vector(groups)
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set post-processing (LSM). First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# ppa_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$ppa_config_file)
# ppa_cmd_string <- paste(c("LSM", build_zonation_params$ppa_raster_file, 0, -1, 0), collapse = " ")
# write(ppa_cmd_string, ppa_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_abf
# 
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Set z-value on taxon level. First, we need to know how many species are in
# # what taxon. HACK: this is not the safest way of doing the z-value allocation
# # but will do.
# 
# spp_rasters <- list.files("../Data.150928/",
#                           pattern = "[A-Z][a-z]+_.+\\.(tif|img)$",
#                           recursive = TRUE)
# # Split by path separator; subfolder is the taxon.
# spp_rasters <- as.data.frame(do.call("rbind",
#                                      strsplit(spp_rasters, .Platform$file.sep)))
# names(spp_rasters) <- c("taxon_name", "species")
# 
# # Manipulate the spp data
# variant_spp_data <- sppdata(variant)
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with caz
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 2)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_caz_wgt
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# dat_sk <- dat[dat$st.species == "Salamandrella_keyserlingii",]
# dat_sk$st.species <- "Salamandrella_keyserlingii_comp"
# dat <- rbind(species_data,  dat_sk)
# 
# # sppweights()<- not implemented yet.
# variant_spp_data <- sppdata(variant)
# 
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 1)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_abf_wgt
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 2)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_[TX]_caz_wgt_con
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 1)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Set condition. First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$condition_config_file)
# condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
# write(condition_cmd_string, condition_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "use condition layer", 1)
# variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
# # We also need to (manually) set condition reference in the groups data
# variant@groups$condition.group <- 1
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_[TX]_abf_wgt_con
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 2)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Set condition. First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$condition_config_file)
# condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
# write(condition_cmd_string, condition_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "use condition layer", 1)
# variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
# # We also need to (manually) set condition reference in the groups data
# variant@groups$condition.group <- 1
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_[TX]_caz_wgt_con_hm2
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 1)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Set condition. First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$condition_config_file)
# condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
# write(condition_cmd_string, condition_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "use condition layer", 1)
# variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
# # We also need to (manually) set condition reference in the groups data
# variant@groups$condition.group <- 1
# 
# # Set hierarchical mask hm2_raster_file
# variant <- set_dat_param(variant, "use mask", 1)
# variant <- set_dat_param(variant, "mask file", hm2_raster_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_[TX]_abf_wgt_con_hm2
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 2)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Set condition. First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$condition_config_file)
# condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
# write(condition_cmd_string, condition_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "use condition layer", 1)
# variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
# # We also need to (manually) set condition reference in the groups data
# variant@groups$condition.group <- 1
# 
# # Set hierarchical mask hm2_raster_file
# variant <- set_dat_param(variant, "use mask", 1)
# variant <- set_dat_param(variant, "mask file", hm2_raster_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_[TX]_caz_wgt_con_hm3
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 1)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Set condition. First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$condition_config_file)
# condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
# write(condition_cmd_string, condition_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "use condition layer", 1)
# variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
# # We also need to (manually) set condition reference in the groups data
# variant@groups$condition.group <- 1
# 
# # Set hierarchical mask hm2_raster_file
# variant <- set_dat_param(variant, "use mask", 1)
# variant <- set_dat_param(variant, "mask file", build_zonation_params$hm3_raster_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
# 
# # [XX]_[TX]_abf_wgt_con_hm3
# variant_id <- variant_id + 1
# variant <- get_variant(zproject, variant_id)
# message("Editing variant: ", variant@name)
# 
# # Manipulate sppdata
# variant_spp_data <- sppdata(variant)
# variant_spp_data$weight <- dat[[build_zonation_params$w_field]]
# sppdata(variant) <- variant_spp_data
# 
# # Set groups, use the same as with previous
# groups(variant) <- as.vector(groups)
# 
# # Set cell removal rule
# variant <- set_dat_param(variant, "removal rule", 2)
# 
# # Set groups use and groups file
# variant <- set_dat_param(variant, "use groups", 1)
# # Note that groups file location is always relative to the bat file
# groups_file <- file.path(variant@name, paste0(variant@name, "_groups.txt"))
# variant <- set_dat_param(variant, "groups file", groups_file)
# 
# # Set PPA file
# variant <- set_dat_param(variant, "post-processing list file",
#                          build_zonation_params$ppa_config_file)
# 
# # Set condition. First, let's create the file itself (zonator
# # can't handle this yet). The file needs to be created only once per raxon
# # since all the variants can use the same file.
# condition_file_name <- file.path(build_zonation_params$zsetup_root, taxon, build_zonation_params$condition_config_file)
# condition_cmd_string <- paste(c(1, build_zonation_params$condition_raster_file), collapse = " ")
# write(condition_cmd_string, condition_file_name)
# # Need to define ppa_config.txt relative to the bat-file (same dir)-
# variant <- set_dat_param(variant, "use condition layer", 1)
# variant <- set_dat_param(variant, "condition file", build_zonation_params$condition_config_file)
# # We also need to (manually) set condition reference in the groups data
# variant@groups$condition.group <- 1
# 
# # Set hierarchical mask hm2_raster_file
# variant <- set_dat_param(variant, "use mask", 1)
# variant <- set_dat_param(variant, "mask file", build_zonation_params$hm3_raster_file)
# 
# # Save variant
# save_zvariant(variant, dir = file.path(build_zonation_params$zsetup_root, taxon), overwrite = TRUE,
#               debug_msg = FALSE)
# create_sh_file(variant)
