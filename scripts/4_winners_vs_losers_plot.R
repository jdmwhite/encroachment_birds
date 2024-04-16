#________________________________________________#
##### Load libraries #############################
#________________________________________________#
library(tidyverse)
library(grid)
library(gridExtra)
library(magick)

#________________________________________________#
##### Load processed colext dataset ##############
#________________________________________________#
all_spp_est <- read_csv('processed/all_spp_estimates.csv')

# select relevant columns
plot_est <- all_spp_est %>% dplyr::select(species:outcome)

# identify lambda vs z-score conflicts
plot_est %>% mutate(lambda_zscore_conflict = case_when(
                     lambda_category == 'decreasing' & z_score > 1.64 ~ 'conflict',
                     lambda_category == 'increasing' & z_score < -1.64 ~ 'conflict',
                     lambda_category == 'increasing' & z_score < 0 & z_score > -1.64 ~ 'possible_conflict',
                     lambda_category == 'decreasing' & z_score > 0 & z_score < 1.64 ~ 'possible_conflict',
                     TRUE ~ 'non-conflict')) -> plot_est
table(plot_est$lambda_zscore_conflict)

#________________________________________________#
##### Prepare plot annotations ###################
#________________________________________________#
text_negative <- textGrob("Negative", gp=gpar(fontsize=10))
text_wc_stable <- textGrob("Neutral", gp=gpar(fontsize=10))
text_positive <- textGrob("Positive", gp=gpar(fontsize=10))
text_increasing <- textGrob("Increasing", gp=gpar(fontsize=10), rot = 270)
text_occ_stable <- textGrob("Stable", gp=gpar(fontsize=10), rot = 270)
text_decreasing <- textGrob("Decreasing", gp=gpar(fontsize=10), rot = 270)
text_OR_franc_179 <- textGrob('Orange River\nfrancolin', gp=gpar(fontsize=10, fontface = 'italic'))
text_SDC_sub_760 <- textGrob('Southern \ndouble-collared sunbird', gp=gpar(fontsize=10, fontface = 'italic'))

# Load in grobs for png images
orf_png <- magick::image_read('bird_images/ORF.png') %>%
  magick::image_background('none') %>%
  grid::rasterGrob()

sdcs_png <- magick::image_read('bird_images/SDCS.png') %>%
  magick::image_background('none') %>%
  magick::image_flop() %>%
  grid::rasterGrob()

# separate out winners/losers and indifferent
plot_est %>% filter(lambda_zscore_conflict == 'non-conflict' & outcome %in% c('loser','winner')) -> wl
plot_est %>% filter(lambda_zscore_conflict == 'non-conflict'  & outcome == 'indifferent') -> indifferent

# change order of levels
plot_est$outcome <- factor(plot_est$outcome, levels = c('loser', 'indifferent', 'winner'))

#________________________________________________#
############## Make plot #########################
#________________________________________________#
plot_est %>% 
  filter(lambda_zscore_conflict == 'non-conflict') %>%
  ggplot(aes(x = occ_est, y = lambda_mean)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 1) +
  geom_errorbar(data = indifferent, aes(x = occ_est, ymin = lambda_0.05, ymax = lambda_0.95), col = 'gray90') +
  geom_errorbarh(data = indifferent, aes(xmin = occ_CI_0.05, xmax = occ_CI_0.95, y = lambda_mean), col = 'gray90') +
  geom_errorbar(data = wl, aes(x = occ_est, ymin = lambda_0.05, ymax = lambda_0.95), col = 'gray65') +
  geom_errorbarh(data = wl, aes(xmin = occ_CI_0.05, xmax = occ_CI_0.95, y = lambda_mean), col = 'gray65') +
  geom_point(aes(fill = outcome), shape = 21, alpha = 0.7, size = 2.5) +
  scale_fill_manual(values = c('red','grey', 'blue'), name = 'Class', labels = c('Losers','Indifferent', 'Winners')) +
  xlab('Woody cover coefficient (initial occupancy)') + ylab('Occupancy trend') +
  annotation_custom(text_positive,xmin=6.4,xmax=6.4,ymin=1.1,ymax=1.1) + 
  annotation_custom(text_wc_stable,xmin=0,xmax=0,ymin=1.1,ymax=1.1) + 
  annotation_custom(text_negative,xmin=-4.6,xmax=-4.6,ymin=1.1,ymax=1.1) + 
  annotation_custom(text_increasing,xmin=7.3,xmax=7.3,ymin=1.06,ymax=1.06) + 
  annotation_custom(text_occ_stable,xmin=7.3,xmax=7.3,ymin=1,ymax=1) + 
  annotation_custom(text_decreasing,xmin=7.3,xmax=7.3,ymin=0.77,ymax=0.77) + 
  annotation_custom(text_OR_franc_179,xmin=-3.5,xmax=-3.5,ymin=0.85,ymax=0.85) +
  annotation_custom(text_SDC_sub_760,xmin=3.6,xmax=3.6,ymin=1.075,ymax=1.075) +
  annotation_custom(orf_png, xmin = -3.8, xmax = -2.8, ymin = 0.86, ymax = 0.91) +
  annotation_custom(sdcs_png, xmin = 3, xmax = 4, ymin = 0.995+0.025, ymax = 1.045+0.025) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(color = 'black', size = 12),
        axis.text = element_text(color = 'black', size = 10),
        axis.ticks = element_line(color = 'black'),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 4),
        plot.margin = unit(c(2,2,2,2), "lines"))

# save output
ggsave('output/fig2_winners_vs_losers.png', dpi = 'retina', width = 8.61, height = 6.33)

# END
