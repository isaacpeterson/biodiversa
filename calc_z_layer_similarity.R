# devtools::install_github("cbig/zonator")
library(Matrix)
library(raster)
library(foreign)

assess_mean_sq <- function(layer_a, layer_b){
  
  apply(abs(layer_a - layer_b), 'sum')
  
}

z_rank_filenames = paste0('~/GitHub/biodiversa/tiff_ranks/', list.files('~/GitHub/biodiversa/tiff_ranks/', pattern = '.tif'))

raster_stack = stack((z_rank_filenames[1:(length(z_rank_filenames))]))
full_links_ranked = raster(z_rank_filenames[length(z_rank_filenames) - 1])
threshold_set = full_links_ranked
#threshold_set[which(values(threshold_set) < 0.90)] <- NA

threshold_set = which(!is.na(values(threshold_set)))

corrs = lapply(seq(dim(raster_stack)[3]), 
               function(i) cor(values(raster_stack)[threshold_set, i], values(full_links_ranked)[threshold_set], use = 'complete')
)

threshold_set = full_links_ranked
threshold_set[which(values(threshold_set) < 0.67)] <- NA
threshold_set = which(!is.na(values(threshold_set)))

diffs_full = lapply(seq(dim(raster_stack)[3]), 
                  function(i) sum(abs(values(raster_stack)[threshold_set, i] - values(full_links_ranked)[threshold_set])))

par(pty="s")
plot_vec = c(0.01, 0.05, seq(0.1, 0.7, by = 0.1), 1)
pdf('~/GitHub/biodiversa/diff_plot_overlay_upper_99.pdf', height = 5, width = 4.75)
#plot(plot_vec, diffs_full[1:10], ylab = "", xlab = "", type = 'l')
#plot(plot_vec, diffs_95[1:10], ylab = "", xlab = "", col = 'blue', type = 'l')
plot(plot_vec, diffs_99[1:10], ylab = "", xlab = "", col = 'red', type = 'l')

#lines(plot_vec, diffs_67[1:10], ylab = "", xlab = "", col = 'darkgreen')
graphics.off()

# temp_chl_s_nb <- full_links_ranked
# values(temp_chl_s_nb) <- 1:ncell(full_links_ranked)
#
# focal_cor <- focal(
#   x = temp_chl_s_nb,
#   w = matrix(1, 5, 5),
#   fun = function(x, y = raster_stack){
#     cor(values(y)[x, 1], values(y)[x, dim(raster_stack)[3]],
#         use = "na.or.complete")
#   })

# full_links_ranked = as.vector(full_links)
# full_links_ranked = full_links_ranked[!is.na(full_links_ranked)]
# full_links_ranked[full_links_ranked > 0.8]
# current_links_ranked = as.vector(subset(raster_stack, 1))
# current_links_ranked = current_links_ranked[!is.na(current_links_ranked)]

#shift_metric = abs(full_links_ranked[full_links_ranked > 0.8] - current_links_ranked[current_links_ranked > 0.8])
#shift_metric = shift_metric[!is.na(shift_metric)]shift_metric
