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

      data_to_write = '[Settings]'
      data_to_write = c(data_to_write, paste0("removal rule = ", variants_object$z_variants$removal_rule[variant_ind]))
      data_to_write = c(data_to_write, paste0("warp factor = ", variants_object$z_variants$warp_factor[variant_ind]))
      data_to_write = c(data_to_write, paste0("edge removal = ", variants_object$z_variants$edge_removal[variant_ind]))

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

params_object$variants$removal_rule = c(2)
params_object$variants$warp_factor = 1000
params_object$variants$edge_removal = 1

params_object$variants$use_weights = c(FALSE)
params_object$links_fraction = c(1, 5, 20, 50, 100, 10)
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
                                            '#SBATCH -J variant_block',
                                            '#SBATCH --constraint="snb|hsw"',
                                            '#SBATCH -o variant_block_output',
                                            '#SBATCH -e variant_block_error',
                                            '#SBATCH -p serial',
                                            '#SBATCH -n 1',
                                            '#SBATCH -t 40:00:00',
                                            '#SBATCH --mem-per-cpu=130000',
                                            '#SBATCH --mail-type=END',
                                            '#SBATCH --mail-user=isaac.peterson@helsinki.fi',
                                            'module load geo-env')

params_object$z_params$global_z_control_file = "run_global_z.sh"

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
