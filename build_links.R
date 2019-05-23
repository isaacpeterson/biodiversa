library(Matrix)
library(raster)
library(foreign)
library(reshape2)

aggregate_rasters <- function(write_dir, current_filenames){

  for (raster_ind in seq_along(current_filenames)){
    current_species_filename = paste0(write_dir, current_filenames[raster_ind])
    current_raster = raster(current_species_filename)
    if (raster_ind == 1){
      raster_stack = current_raster
    } else{
      raster_stack = stack(raster_stack, current_raster)
      raster_stack = sum(raster_stack)
    }
    
  }
  
  return(raster_stack) 
  
}

write_aggregated_rasters <- function(pp_links, raster_dir){

  if (!file.exists(paste0(raster_dir, 'aggregated/'))){
    dir.create(paste0(raster_dir, 'aggregated/'), recursive = TRUE)
  }
  
  for (predator_ind in seq_along(rownames(pp_links))){
    current_links_filenames = list.files(path = raster_dir, pattern = paste0('pred_', predator_ind, '_'))
    
    if (length(current_links_filenames) > 0){

      raster_to_write <- aggregate_rasters(raster_dir, current_links_filenames)
      writeRaster(raster_to_write, paste0(paste0(raster_dir, 'aggregated/'), 'pred_', predator_ind, '_aggregated.tif'), format='GTiff', overwrite = TRUE)
      print(paste('predator', predator_ind, 'done'))
    }
    
    
  }
}

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
build_links_params$workdir = 'species_data/5km/'
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

write_aggregated_rasters(output_data_objects$pp_links, raster_dir = paste0(build_links_params$workdir, '/link_layers/'))

pdf(file = '~/GitHub/biodiversa/pp_links.pdf')
par(pty="s")
image(raster(output_data_objects$symmetrised_pp), col = c("#000000", "#E6E6E6"))
graphics.off()




