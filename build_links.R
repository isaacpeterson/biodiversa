library(Matrix)
library(raster)
library(foreign)
library(reshape2)


build_current_raster <- function(raster_to_write, species_loc_mapping, elements_to_replace, current_filename){
  values(raster_to_write)[species_loc_mapping] <- elements_to_replace
  writeRaster(raster_to_write, current_filename, format='GTiff', overwrite = TRUE)
  
}

build_link_layers <- function(pp_links, sites_by_species, study_region, species_loc_mapping, write_dir){
  
  if (!file.exists(write_dir)){
    dir.create(write_dir, recursive = TRUE)
  }
  
  links_to_use <- which(pp_links > 0)
  link_num = length(links_to_use)
  row_num = dim(pp_links)[1]
  
  for (link_ind in seq(link_num)){
    
    current_pred = ((links_to_use[link_ind] - 1) %% row_num) + 1
    current_prey = floor((links_to_use[link_ind] - 1) / row_num) + 1
    current_pp_link_per_site = input_data$sites_by_species[, current_pred] * input_data$sites_by_species[, current_prey]
    current_pred_folder = paste0(write_dir, 'pred_', current_pred, '/')
    
    if (!file.exists(current_pred_folder)){
      dir.create(current_pred_folder, recursive = TRUE)
    }
    
    build_current_raster(raster_to_write = study_region, 
                         species_loc_mapping, 
                         elements_to_replace = current_pp_link_per_site, 
                         current_filename = paste0(current_pred_folder, 'pred_', current_pred, '_prey_', current_prey, '.tif'))
    
    print(paste('link', link_ind, 'of', link_num, 'done'))
    
  }
  
}


build_species_layers <- function(site_by_species, study_region, species_loc_mapping, write_dir){
  if (!file.exists(write_dir)){
    dir.create(write_dir, recursive = TRUE)
  }
  
  species_num = dim(site_by_species)[2]
  
  for (species_ind in seq(species_num)){
    
    build_current_raster(raster_to_write = study_region, 
                         species_loc_mapping, 
                         elements_to_replace = site_by_species[, species_ind], 
                         current_filename = paste0(write_dir, 'species_', species_ind, '.tif'))
    
    print(paste('species', species_ind, 'of', species_num, 'done'))
  }
  
}


build_links_params = list()
build_links_params$build_link_layers = FALSE
build_links_params$build_species_layers = FALSE
build_links_params$workdir = '~/species_data/biodiversa/5km/'
build_links_params$site_by_species_filename = paste0(build_links_params$workdir, 'Site_By_Species_21Feb.rds')
build_links_params$pred_by_prey_filename = paste0(build_links_params$workdir, 'BARM_binary_without_eggs.rds')
build_links_params$grid_ref_filename = paste0(build_links_params$workdir, 'reference_grid_5km.img')
build_links_params$ref_grid_vals = paste0(build_links_params$workdir, 'reference_grid_5km.img.vat.dbf')



input_data = list()
input_data$study_region = raster(build_links_params$grid_ref_filename)
input_data$species_locs = read.dbf(build_links_params$ref_grid_vals)
input_data$pred_by_prey = readRDS(build_links_params$pred_by_prey_filename)
input_data$sites_by_species = readRDS(build_links_params$site_by_species_filename)

output_data_objects <- list()
output_data_objects$species_loc_mapping = match(seq(dim(input_data$sites_by_species)[1]), values(input_data$study_region))
output_data_objects$co_occurances <- t(input_data$sites_by_species) %*% input_data$sites_by_species
output_data_objects$co_occurances <- as.matrix(output_data_objects$co_occurances)
output_data_objects$co_occurances[output_data_objects$co_occurances > 1] <- 1
diag(output_data_objects$co_occurances) <- 0

output_data_objects$symmetrised_pp = input_data$pred_by_prey[ , match(rownames(input_data$pred_by_prey), colnames(input_data$pred_by_prey))]
output_data_objects$symmetrised_pp[is.na(output_data_objects$symmetrised_pp)] <- 0
diag(output_data_objects$symmetrised_pp) <- 0

sp_overlaps <- match(colnames(input_data$sites_by_species), rownames(output_data_objects$symmetrised_pp), nomatch = 0)
output_data_objects$symmetrised_pp <- output_data_objects$symmetrised_pp[sp_overlaps, sp_overlaps]
output_data_objects$pp_links <- output_data_objects$co_occurances * output_data_objects$symmetrised_pp

if (build_links_params$build_species_layers == TRUE){
  build_species_layers(input_data$sites_by_species, 
               input_data$study_region, 
               output_data_objects$species_loc_mapping, 
               write_dir = paste0(build_links_params$workdir, '/species_layers/'))
}


if (build_links_params$build_link_layers == TRUE){
  build_link_layers(output_data_objects$pp_links, 
                    input_data$sites_by_species, 
                    input_data$study_region, 
                    output_data_objects$species_loc_mapping, 
                    write_dir = paste0(build_links_params$workdir, '/link_layers/'))
}

link_filenames = list.files(paste0(build_links_params$workdir, '/link_layers/'))


