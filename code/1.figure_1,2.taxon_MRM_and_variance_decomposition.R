# Clear console and environment
cat("\014")
rm(list = ls())

# Load required packages
library(readr)
library(dplyr)
library(tidyr)
library(betapart)
library(geosphere)
library(ecodist)
library(ggpubr)
library(tibble)
library(forcats)

# Set working directory
setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()  # Confirm working directory

# Load species checklist data
data <- read_csv("./data/Campus_checklist_status.csv")

# Function to compute beta diversity for a group of universities
get_beta_div <- function(df, status_name, all_univ) {
  # Create presence-absence matrix
  mat <- df %>%
    mutate(species = gsub(" ", "_", taxon_name)) %>%
    distinct(univ.links.uni.abbrev02, species) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = species, values_from = presence, values_fill = 0) %>%
    column_to_rownames(var = "univ.links.uni.abbrev02")
  
  # Add missing universities with zero-filled rows
  missing_rows <- setdiff(all_univ, rownames(mat))
  if (length(missing_rows) > 0) {
    zero_matrix <- matrix(0, nrow = length(missing_rows), ncol = ncol(mat))
    rownames(zero_matrix) <- missing_rows
    colnames(zero_matrix) <- colnames(mat)
    mat <- rbind(mat, zero_matrix)
  }

  # Compute beta diversity components
  beta_result <- beta.pair(mat, index.family = "sorensen")

  # Return result matrices
  return(list(
    turnover = as.matrix(beta_result$beta.sim),
    nestedness = as.matrix(beta_result$beta.sne),
    total = as.matrix(beta_result$beta.sor),
    status = status_name
  ))
}

# Subset data by species types
invasive.data <- data %>% filter(invasive == 1)
non_invasive.data <- data %>% filter(non_invasive == 1)
native.data <- data %>% filter(native == 1)
all.data <- data

# Get list of all university abbreviations
all_univ <- sort(unique(data$univ.links.uni.abbrev02))

# Compute beta diversity for each group
beta_invasive <- get_beta_div(invasive.data, "invasive", all_univ)
beta_noninvasive <- get_beta_div(non_invasive.data, "noninvasive", all_univ)
beta_native <- get_beta_div(native.data, "native", all_univ)
beta_all <- get_beta_div(all.data, "all", all_univ)

# Store beta diversity matrices
beta_list <- list(
  invasive_turnover = beta_invasive$turnover,
  invasive_nestedness = beta_invasive$nestedness,
  invasive_total = beta_invasive$total,

  noninvasive_turnover = beta_noninvasive$turnover,
  noninvasive_nestedness = beta_noninvasive$nestedness,
  noninvasive_total = beta_noninvasive$total,
  
  native_turnover = beta_native$turnover,
  native_nestedness = beta_native$nestedness,
  native_total = beta_native$total,
  
  all_turnover = beta_all$turnover,
  all_nestedness = beta_all$nestedness,
  all_total = beta_all$total
)

# When saving (optional)
# beta_list_dist <- lapply(beta_list, function(m) stats::as.dist(m))
# saveRDS(beta_list_dist, "./results/beta_list_taxon_dist.rds")

# When reading (optional)
# beta_list_dist <- readRDS("./results/beta_list_taxon_dist.rds")
# beta_list <- lapply(beta_list_dist, function(d) as.matrix(d))

# Load campus-level environmental and spatial data
campus_data <- read_csv("./data/Driving_factors.csv") %>%
  mutate(campus02 = as.character(campus02))

# Reorder to match matrix names
univ_order <- rownames(beta_invasive$turnover)
campus_data_ordered <- campus_data %>%
  filter(campus02 %in% univ_order) %>%
  arrange(factor(campus02, levels = univ_order))

# Extract coordinates and compute geographic distance matrix (km)
coords <- campus_data_ordered %>% dplyr::select(lng, lat)
geo_dist_matrix <- distm(coords, fun = distHaversine) / 1000
rownames(geo_dist_matrix) <- colnames(geo_dist_matrix) <- univ_order

# Extract and prepare environmental variables
env_vars <- campus_data_ordered[, c("est.time", "area_ha", "dem_mean", "map_mean", "mat_mean", "gdp_mean", "pd_mean", "ul_mean")]

# Compute Euclidean distance matrices
dist_est.time <- dist(env_vars$est.time) %>% as.matrix()
dist_area_ha <- dist(env_vars$area_ha) %>% as.matrix()
dist_dem_mean <- dist(env_vars$dem_mean) %>% as.matrix()
dist_map_mean <- dist(env_vars$map_mean) %>% as.matrix()
dist_mat_mean <- dist(env_vars$mat_mean) %>% as.matrix()
# dist_wet_area <- dist(env_vars$wet_area) %>% as.matrix()
# dist_wet_num <- dist(env_vars$wet_num) %>% as.matrix()
dist_gdp_mean <- dist(env_vars$gdp_mean) %>% as.matrix()
dist_pd_mean <- dist(env_vars$pd_mean) %>% as.matrix()
dist_ul_mean <- dist(env_vars$ul_mean) %>% as.matrix()

# Global standardization function
standardize_matrix <- function(mat) {
  mat_mean <- mean(mat)
  mat_sd <- sd(mat)
  return((mat - mat_mean) / mat_sd)
}

# Standardize all matrices
geo_dist_std <- standardize_matrix(geo_dist_matrix)
est.time_std <- standardize_matrix(dist_est.time)
area_ha_std <- standardize_matrix(dist_area_ha)
dem_mean_std <- standardize_matrix(dist_dem_mean)
map_mean_std <- standardize_matrix(dist_map_mean)
mat_mean_std <- standardize_matrix(dist_mat_mean)
# wet_area_std <- standardize_matrix(dist_wet_area)
# wet_num_std <- standardize_matrix(dist_wet_num)
gdp_mean_std <- standardize_matrix(dist_gdp_mean)
pd_mean_std <- standardize_matrix(dist_pd_mean)
ul_mean_std <- standardize_matrix(dist_ul_mean)

# Assign consistent names
matrix_list <- list(geo_dist_std, est.time_std, area_ha_std, dem_mean_std, map_mean_std, mat_mean_std, 
                     gdp_mean_std, pd_mean_std, ul_mean_std)
matrix_list <- lapply(matrix_list, function(m) {
  rownames(m) <- colnames(m) <- univ_order
  return(m)
})

# Prepare predictor dataframe (distance vectors)
predictors_df <- data.frame(
  geo1 = as.vector(as.dist(geo_dist_std)),
  est.time1 = as.vector(as.dist(est.time_std)),
  area_ha1 = as.vector(as.dist(area_ha_std)),
  dem_mean1 = as.vector(as.dist(dem_mean_std)),
  map_mean1 = as.vector(as.dist(map_mean_std)),
  mat_mean1 = as.vector(as.dist(mat_mean_std)),
  # wet_area1 = as.vector(as.dist(wet_area_std)),
  # wet_num1 = as.vector(as.dist(wet_num_std)),
  gdp_mean1 = as.vector(as.dist(gdp_mean_std)),
  pd_mean1 = as.vector(as.dist(pd_mean_std)),
  ul_mean1 = as.vector(as.dist(ul_mean_std))
)

# Perform MRM for all combinations of group and beta components
groups <- c("invasive", "noninvasive", "native", "all")
components <- c("total", "turnover", "nestedness")

results <- list()

for (group in groups) {
  for (component in components) {
    beta_name <- paste0(group, "_", component)
    beta_vec <- as.vector(as.dist(beta_list[[beta_name]]))
    
    mrm_model <- MRM(beta_vec ~ geo1 + est.time1 + area_ha1 + dem_mean1 + map_mean1 + 
                     mat_mean1 + gdp_mean1 + pd_mean1 + ul_mean1,
                     data = predictors_df,
                     nperm = 999)
    
    coef_table <- as.data.frame(mrm_model$coef)
    colnames(coef_table) <- c("Coefficient", "p_value")
    coef_table$Variable <- rownames(coef_table)
    
    results[[beta_name]] <- cbind(
      data.frame(
        Group = group,
        Component = component,
        R_squared = mrm_model$r.squared[1]
      ),
      coef_table
    )
  }
}

# Combine all results
final_results <- do.call(rbind, results)
# Add significance stars based on p-value thresholds
final_results$Significance <- cut(final_results$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", ""),
  right = FALSE
)

# Output updated result
print(final_results)

write.csv(final_results, "./results/MRM_results_taxon.csv", row.names = FALSE)

# Determine significance
final_results <- final_results %>%
  mutate(signif = if_else(p_value < 0.05, "significant", "non-significant"))

# Add beta_type and beta_index classifications to final_results
final_results <- final_results %>%
  mutate(
    beta_type = case_when(
      Group == "all" ~ "(A) All\nspecies",
      Group == "native" ~ "(B) Native\nspecies",
      Group == "noninvasive" ~ "(C) Non-invasive\nspecies",
      Group == "invasive" ~ "(D) Invasive\nspecies",
      TRUE ~ NA_character_
    ),
    beta_index = case_when(
      Component == "turnover" ~ "Turnover",
      Component == "nestedness" ~ "Nestedness",
      Component == "total" ~ "Beta diversity",
      TRUE ~ NA_character_
    ),
    signif = ifelse(p_value < 0.05, "1", "0"),
    beta_signif = ifelse(signif==1, beta_index, signif),
    variable = case_when(
      Variable == "geo1" ~ "Campus distance (km)",
      Variable == "est.time1" ~ "Campus age (years)",
      Variable == "area_ha1" ~ "Campus area (km²)",
      Variable == "dem_mean1" ~ "Elevation (m)",
      Variable == "mat_mean1" ~ "MAT (°C)",
      Variable == "map_mean1" ~ "MAP (mm)",
      # Variable == "wet_area1" ~ "City wetland area (km²)",
      # Variable == "wet_num1" ~ "City wetland number",
      Variable == "gdp_mean1" ~ "GDP per area\n(10,000 CNY/km²)",
      Variable == "pd_mean1" ~ "Population density\n(persons/km²)",
      Variable == "ul_mean1" ~ "Urbanization level (%)"
    )
  ) %>%
  filter(!is.na(variable)) %>%
  mutate(
    # Set custom order
    variable = factor(variable, levels = c(
      "Campus age (years)",
      "Campus area (km²)",
      "Campus distance (km)",
      "Elevation (m)",
      "MAT (°C)",
      "MAP (mm)",
      # "City wetland area (km²)",
      # "City wetland number",
      "GDP per area\n(10,000 CNY/km²)",
      "Population density\n(persons/km²)",
      "Urbanization level (%)"
    )),
    beta_index = factor(beta_index, levels = c("Beta diversity", "Turnover", "Nestedness")),
    beta_type = factor(beta_type, levels = c("(A) All\nspecies", "(B) Native\nspecies", "(C) Non-invasive\nspecies", "(D) Invasive\nspecies"))
  )

# Customize Graphic Style
beta_type_shape <- c("Beta diversity" = 21, "Turnover" = 22, "Nestedness" = 23)
beta_type_color <- c("Beta diversity" = "#1B9E77", "Turnover" = "#7570B3", "Nestedness" = "#D95F02", "0" = "gray")

# Obtain the legend (categorized by beta_type)
plot_legend_beta_type <- ggplot(final_results) + 
  facet_wrap(~ beta_type, scale = "free_x", ncol = 4) + 
  geom_point(aes(y = variable, x = Coefficient, shape = fct_rev(beta_index), fill = beta_index, color = beta_index),
             position = position_dodge(0.6), size = 3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  labs(y = NULL, x = "Coefficient") +  
  scale_shape_manual(name = NULL, values = beta_type_shape) + 
  scale_color_manual(name = NULL, values = beta_type_color) + 
  scale_fill_manual(name = NULL, values = beta_type_color) + 
  guides(shape = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(legend.position = "bottom", 
        legend.background = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        legend.text = element_text(size = 14))

legend_beta_type <- ggpubr::get_legend(plot_legend_beta_type)

# main plot(no legend)
beta_plot <- ggplot(final_results) + 
  facet_wrap(~ beta_type, scale = "free_x", ncol = 4) + 
  geom_point(aes(y = variable, x = Coefficient, shape = fct_rev(beta_index), fill = beta_signif, color = beta_signif),
             position = position_dodge(0.6), size = 3) + 
  geom_vline(xintercept = 0, linetype = 2, size = 0.3) +
  labs(y = NULL, x = "Coefficient") +  
  scale_shape_manual(name = NULL, values = beta_type_shape) + 
  scale_color_manual(name = NULL, values = beta_type_color) + 
  scale_fill_manual(name = NULL, values = beta_type_color) + 
  scale_x_continuous(
  limits = c(-0.11, 0.11),
  breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
  labels = c("-0.1", "-0.05", "0", "0.05", "0.1")) +
  theme_bw() +
  theme(legend.position = "n",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,0,1), "mm"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14))

# combine main plot and legend plot
final_plot <- cowplot::plot_grid(
  cowplot::plot_grid(NULL, legend_beta_type, rel_widths = c(0.2, 1)),
  beta_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

print(final_plot)

ggexport(final_plot, filename = "./results/plots.taxon.png",
         width = 3000,
         height = 2000,
         pointsize = 12,
         res = 300)

# ----------------- extraction R2 -----------------
get_R2 <- function(m){
  if(!is.null(m$R2)) return(as.numeric(m$R2))
  if(!is.null(m$r.squared)) return(as.numeric(m$r.squared[1]))
  stop("Cannot find R2 in MRM object.")
}

# ----------------- Define two variable blocks with names -----------------
blockA_vars <- c("geo1","dem_mean1","mat_mean1","map_mean1")  # Environmental
blockB_vars <- c("est.time1","area_ha1","gdp_mean1","pd_mean1","ul_mean1")  # Anthropogenic

blockA_name <- "Environmental"
blockB_name <- "Anthropogenic"

# ----------------- Variance decomposition function (for a beta matrix) -----------------
partition_two_blocks <- function(beta_mat, predictors_df, vars_A, vars_B, nperm_full = 999){
  y <- as.vector(as.dist(beta_mat))

  m_full <- MRM(y ~ ., data = predictors_df[, c(vars_A, vars_B)], nperm = nperm_full)
  m_A    <- MRM(y ~ ., data = predictors_df[, vars_A],             nperm = 0)
  m_B    <- MRM(y ~ ., data = predictors_df[, vars_B],             nperm = 0)

  R2_full <- get_R2(m_full)
  R2_A    <- get_R2(m_A)
  R2_B    <- get_R2(m_B)

  pure_A  <- R2_full - R2_B
  pure_B  <- R2_full - R2_A
  unexpl  <- 1 - R2_full
  shared  <- R2_A + R2_B - R2_full

tibble(
    R2_full     = R2_full,
    pure_blockA = pure_A,        
    pure_blockB = pure_B,       
    unexplained = unexpl,
    shared      = shared,
    pureA_pos   = pmax(0, pure_A),
    pureB_pos   = pmax(0, pure_B),
    unexpl_pos  = pmax(0, unexpl),
    shared_pos  = pmax(0, shared)
  )
}

# ----------------- Run batch operations for all group × component combinations -----------------
groups     <- c("invasive", "noninvasive", "native", "all")
components <- c("total", "turnover", "nestedness")

part_list <- list()

for(g in groups){
  for(comp in components){
    nm <- paste0(g, "_", comp)
    beta_mat <- beta_list[[nm]]
    res <- partition_two_blocks(beta_mat, predictors_df, blockA_vars, blockB_vars, nperm_full = 999)
    res$Group <- g
    res$Component <- comp
    part_list[[nm]] <- res
  }
}

partition_table <- bind_rows(part_list) %>%
  relocate(Group, Component)

partition_table0 <- partition_table %>%
  select(Group, Component, R2_full, pureA_pos, pureB_pos, unexpl_pos, shared_pos)

write.csv(partition_table0, "./results/MRM_taxon_Variance decomposition.csv", row.names = FALSE)

blockA_name <- "Environmental"
blockB_name <- "Anthropogenic"

# Use rename() to change the column name
plot_df <- partition_table %>%
  select(Group, Component, pureA_pos, pureB_pos, unexpl_pos, shared_pos) %>%
  dplyr::mutate(dplyr::across(c(pureA_pos, pureB_pos, unexpl_pos, shared_pos), ~ .x * 100)) %>%
  rename(
    PureB = pureB_pos, 
    PureA = pureA_pos   
  ) %>%
  tidyr::pivot_longer(
    cols = c(PureA, PureB, shared_pos, unexpl_pos),
    names_to = "Part",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Component = factor(
      Component,
      levels = c("total", "turnover", "nestedness"),
      labels = c("Beta diversity", "Turnover", "Nestedness")),
    Group = factor(
      Group,
      levels = c("all", "native", "noninvasive", "invasive"),
      labels = c("(A) All species", "(B) Native species", "(C) Non-invasive species", "(D) Invasive species")),
    Part = factor(
      Part,
      levels = c("unexpl_pos", "PureA", "PureB", "shared_pos"),
      labels = c("Unexplained", blockA_name, blockB_name, "Shared")
    ))

fill_colors <- c(
    "Unexplained"           = "#F2F2F2",  
    "Environmental"             = "#E41A1C", 
    "Anthropogenic"             = "#377EB8",  
    "Shared"                = "#FF7F00" 
  )
legend_breaks <- names(fill_colors)

# variance_decomposition plot
vd_plot <- ggplot(plot_df, aes(x = Component, y = Value, fill = Part)) +
  geom_col(width = 0.6) +
  facet_wrap(~ Group, nrow = 1) +
  scale_fill_manual(values = fill_colors,breaks = legend_breaks) +  
  scale_y_continuous(limits = c(0, 100.1),
                     breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = "Proportion of variance (%)", x = NULL, fill = NULL,
       title = NULL) +
  theme_bw()+
  theme(legend.position = "top",
        legend.direction = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1,1,0,1), "mm"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black",
                                 angle = 35, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14, colour = "black"),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 1, ncol = 4, byrow = TRUE)) 

ggexport(vd_plot, filename = "./results/plots.taxon.vd.png",
         width = 3000,
         height = 2400,
         pointsize = 12,
         res = 300)

