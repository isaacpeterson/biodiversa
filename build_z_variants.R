# NOTE: This script is used to generate Zonation setups for public release.
# See XXXXXX for a more detailed explanation.
#
# devtools::install_github("cbig/zonator")
library(dplyr)
library(raster)
library(readr)
library(readxl)
library(rgdal)

build_variant_files <- function(variants_object, build_type, z_params){
  
  for (variant_ind in seq(nrow(variants_object$z_variants))){
    
    current_variant <- paste0('z_variant_', variant_ind)
    
    if (build_type == 'local_shell'){
      
      data_to_write = c(z_params$local_sh_preamble, 
                        paste("zig4 -r", 
                              paste0(current_variant, '/', current_variant, '.dat'),
                              paste0(current_variant, '/', current_variant, '.spp'), 
                              paste0(current_variant, '/out/', current_variant, '_out.txt'), 
                              z_params$base_z_control,
                              z_params$output_control_string))
      current_file = variants_object$z_variant_files[variant_ind]
    } else if (build_type == 'dat'){
      data_to_write = z_params$z_datfile_control_template
      current_file = paste0(z_params$workdir, current_variant, '/', current_variant, '.dat')
      
    } else if (build_type == 'species_list'){

      current_species_group = variants_object$links_variants[[which(names(variants_object$links_variants) == variants_object$z_variants$group_names[variant_ind])]]
      data_to_write = paste0(z_params$base_species_control, " ", z_params$species_data_rel_to_variant_path, current_species_group)
      current_file = paste0(z_params$workdir, current_variant, '/', current_variant, '.spp')
      
    }
    
    file_ID <- file(current_file)
    writeLines(data_to_write, file_ID)
    close(file_ID)
    
    Sys.chmod(current_file)
    
  }
  
}

params_object = list()

params_object$variants$run_type = c('caz', 'abg')
params_object$variants$use_weights = c(FALSE)
params_object$links_fraction = seq(10,100, by = 10)
params_object$variants$group_names = formatC(params_object$links_fraction, width = 3, format = "d", flag = "0")

params_object$z_params$output_control_string = list("--grid-output-formats=compressed-tif --image-output-formats=png")
params_object$z_params$base_species_control = '1 1 1 1 0.25'
params_object$z_params$base_z_control = '0.0 0 1.0 0'
params_object$z_params$local_sh_preamble = c("#!/bin/sh", 'cd "$(dirname $0)"')
params_object$z_params$species_data_rel_to_variant_path = '../species_data/5km/link_layers/'
params_object$z_params$workdir = 'z_variants/'
params_object$z_params$global_sh_preamble = c('#!/bin/bash -l',
                                            '# created: ',
                                            '# author: peterson',
                                            '#SBATCH -J variant_101',
                                            '#SBATCH --constraint="snb|hsw"',
                                            '#SBATCH -o variant101_output',
                                            '#SBATCH -e variant101_error',
                                            '#SBATCH -p serial',
                                            '#SBATCH -n 1',
                                            '#SBATCH -t 40:00:00',
                                            '#SBATCH --mem-per-cpu=130000',
                                            '#SBATCH --mail-type=END',
                                            '#SBATCH --mail-user=isaac.peterson@helsinki.fi',
                                            'module load geo-env')

params_object$z_params$global_z_control_file = "run_global_z.sh"
params_object$z_params$z_datfile_control_template = c('[Settings]',
                                                      'removal rule = 1', 
                                                      'warp factor = 1000', 
                                                      'edge removal = 1',
                                                      'add edge points = 0',
                                                      'initial removal percent = 0.0')

variants_object <- list()
variants_object$z_variants <- expand.grid(params_object$variants, stringsAsFactors=FALSE)
variants_object$z_variant_files <- paste0(params_object$z_params$workdir, "z_variant_", seq(nrow(variants_object$z_variants)), '.sh')
variants_object$links_files <- list.files(paste0('species_data/5km/link_layers/'), recursive = TRUE)

variants_object$links_variants <- setNames(lapply(seq_along(params_object$links_fraction), 
                                         function(i) sample(variants_object$links_files, 
                                                            round(length(variants_object$links_files)*params_object$links_fraction[[i]]/100))),
                                         params_object$variants$group_names)

file_ID <- file(params_object$z_params$global_z_control_file)
  writeLines(c(params_object$z_params$global_sh_preamble, paste("/bin/sh ", variants_object$z_variant_files)), file_ID)
close(file_ID)

Sys.chmod(params_object$z_params$global_z_control_file)

# build file variants
for (variant_ind in seq(nrow(variants_object$z_variants))){
  
  current_variant <- paste0('z_variant_', variant_ind)
  
  if (!dir.exists(paste0(params_object$z_params$workdir, current_variant, '/out/'))){
    dir.create(paste0(params_object$z_params$workdir, current_variant, '/out/'), recursive = TRUE)
  }
  
}

build_variant_files(variants_object, build_type = 'local_shell', params_object$z_params)
build_variant_files(variants_object, build_type = 'dat', params_object$z_params)
build_variant_files(variants_object, build_type = 'species_list', params_object$z_params)



# build_zonation_params$data_dir <- paste0('../../species_data/zenodo/')
# build_zonation_params$zsetup_root <- "zsetup"
# build_zonation_params$feature_data_file <- paste0(build_zonation_params$data_dir, "/species_params/biodiv_features.csv")
# 
# build_zonation_params$agg_weights_file <- paste0(build_zonation_params$data_dir, "/species_params/aggregate_weights.csv")
# build_zonation_params$ppa_raster_file <- paste0(build_zonation_params$data_dir, "/group_layers/prefecture.tif")
# build_zonation_params$ppa_cmd_string <- paste(c("LSM", build_zonation_params$ppa_raster_file, 0, -1, 0), collapse = " ")
# 
# build_zonation_params$condition_raster_file <- paste0(build_zonation_params$data_dir, "/group_layers/hii_rescaled_simple.tif")
# 
# build_zonation_params$hm3_raster_file <- paste0(build_zonation_params$data_dir, "/group_layers/PA_3_levels_simple.tif")
# 
# build_zonation_params$w_field <- "weight"
# build_zonation_params$condition_config_file <- "condition_config.txt"
# build_zonation_params$ppa_config_file <- "ppa_config.txt"
# build_zonation_params$use_groups = rep(1, 100)
# build_zonation_params$cell_removal_rule = rep(c(1, 2), 50)
# 
# 
# 
# 
# 
# build_zonation_params$species_group_codes = c('amp', 'bir', 'frf', 'mam', 'pla', 'rep')
# build_zonation_params$species_data_template = setNames(lapply(seq_along(build_zonation_params$group_names),
#                                                               function(species_group_ind) setNames(list(build_zonation_params$species_group_codes[species_group_ind], 
#                                                                                                         build_zonation_params$group_names[species_group_ind]),
#                                                                                                    c('code', 'sheet'))),
#                                                               build_zonation_params$group_names)
# 
# build_zonation_params$variant_cycle = rep(1:(length(build_zonation_params$species_data_template) + 1), 
#                                           each = length(build_zonation_params$variant_templates)) 
# 
# build_zonation_params$z_control_template = "templates/run_Z_template.sh"
# build_zonation_params$dat_template_file = "templates/template.dat"
# build_zonation_params$z_control_file = "run_biodiversa_source.sh"
# 
# 
# species_data <- as.data.frame(readr::read_csv(build_zonation_params$feature_data_file))
# 
# # Read in aggregate weights
# agg_weights <- readr::read_csv(build_zonation_params$agg_weights_file) 
# 
# # Set up the project
# zproject <- initiate_zproject(zsetup_root = build_zonation_params$zsetup_root,
#                                     build_zonation_params$variant_templates,
#                                     spp_data = build_zonation_params$species_data_template,
#                                     data_dir = build_zonation_params$data_dir,
#                                     prefix_spp_paths = '../', 
#                                     dat_template_file = build_zonation_params$dat_template_file)
# 
# # Set run configuration parameters ----------------------------------------
# 
# lapply(seq_along(build_zonation_params$species_data_template), 
#        function(group_ind) write(build_zonation_params$ppa_cmd_string,
#                                  file.path(build_zonation_params$zsetup_root, 
#                                            build_zonation_params$species_data_template[[group_ind]]$sheet, 
#                                            build_zonation_params$ppa_config_file)))
# 
# variants <- lapply(seq_along(zproject@variants), 
#                    function(variant_id) get_variant(zproject, 
#                                                     variant_id))
# 
# variants <- lapply(seq_along(zproject@variants), 
#                    function(variant_id) set_dat_param(variants[[variant_id]], 
#                                                       "use groups", 
#                                                       build_zonation_params$use_groups[variant_id])) 
# 
# variants <- lapply(seq_along(zproject@variants), 
#                    function(variant_id) set_dat_param(variants[[variant_id]], 
#                                                       "groups file", 
#                                                       file.path(variants[[variant_id]]@name, paste0(variants[[variant_id]]@name, "_groups.txt")))) 
# 
# variants <- lapply(seq_along(zproject@variants), 
#                    function(variant_id) set_dat_param(variants[[variant_id]], 
#                                                       "post-processing list file",
#                                                       build_zonation_params$ppa_config_file))
# 
# variants <- lapply(seq_along(zproject@variants), 
#                    function(variant_id) set_dat_param(variants[[variant_id]], 
#                                                       "removal rule",
#                                                       build_zonation_params$cell_removal_rule[variant_id] ))
# 
# # lapply(seq_along(zproject@variants), 
# #        function(variant_id) save_zvariant(variants[[variant_id]], 
# #                                           dir = file.path(build_zonation_params$zsetup_root, 
# #                                                           build_zonation_params$species_data_template[[build_zonation_params$variant_cycle[group_ind]]]$sheet),
# #                                           overwrite = TRUE, 
# #                                           debug_msg = FALSE))
# 
# z_control_files = lapply(seq_along(zproject@variants), 
#                          function(variant_id) create_sh_file(variants[[variant_id]], 
#                                                              remove_bat = TRUE))
# 
# build_zonation_params$control_template = "templates/run_Z_template.sh"
# build_zonation_params$dat_template_file = "templates/template.dat"
# build_zonation_params$Z_sh_runfile = "Z_source.sh"
# 
# z_template <- readLines(build_zonation_params$z_control_template)
# 
# file_ID <- file(build_zonation_params$z_control_file)
# writeLines(c(z_template, paste0("/bin/sh ", z_control_files)), file_ID)
# close(file_ID)
# Sys.chmod(build_zonation_params$z_control_file)
# 
# 
# 
