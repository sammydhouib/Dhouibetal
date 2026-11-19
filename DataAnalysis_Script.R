library(readxl)
library(tidyverse)
library(Hmisc)
library(EFAtools)
library(rmcorr)

# Set paths
zenk_path <- "/Users/sammydhouib/Desktop/Manuscript/Data/WTSP_Data.xlsx"                        # ZENK counts + MoveTime
mr_path   <- "/Users/sammydhouib/Desktop/Manuscript/Data/Migratory restlessness raw data.csv"   # repeated-measures data
figdir    <- "/Users/sammydhouib/Desktop/Manuscript/Figures"
dir.create(figdir, showWarnings = FALSE)

# Import and assemble

#Import hippocampus data 
HpRaw <- read_excel(
  path  = zenk_path,
  sheet = "HippocampusRaw")  %>% 
  mutate(
    filename_clean = filename %>%
      str_replace_all("\\s*\\(RGB\\)", "") %>%  
      str_replace_all("-", "_"),                 
    
    band = str_extract(filename_clean, "\\d+"),
    hemisphere = str_extract(filename_clean, "(L|R)(?=\\.tif$)"),
    m = str_match(filename_clean, "^Hp\\d+_\\d+_([aip])(d|dl|vm)hp"),
    caudrost = m[, 2],
    subregion = m[, 3]
  ) %>%
  select(band, caudrost, subregion, hemisphere, count)

HpMean <- HpRaw %>%
  group_by(band, caudrost, subregion) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop"
  )

HpMean <- HpMean %>%
  mutate(band = as.double(band))

#Import movement data
MoveTime <- read_excel(
  path  = zenk_path,
  sheet = "MoveTime")  

HpMean <- HpMean %>%
  left_join(MoveTime, by = "band")

HpMean_collapsed <- HpMean %>%
  group_by(band, caudrost, subregion) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

HpMean_collapsed <- HpMean_collapsed %>%
  mutate(region = paste(caudrost, subregion, sep = "_")) %>%
  select(-caudrost, -subregion) %>%
  relocate(region, .before = 2)

#export hp data to corrmatrixdf
corrmatrixdf <- HpMean_collapsed %>%
  select(band, region, mean_count) %>%
  group_by(band, region) %>%
  summarise(mean_count = mean(mean_count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = band,
    names_from = region,
    values_from = mean_count
  ) %>%
  arrange(band)

#Import Raw NCL data 
NCL_raw <- read_excel(
  path  = zenk_path,
  sheet = "NCLRaw")

NCLsorted <- NCL_raw %>%
  mutate(
    fname = filename %>%
      str_trim() %>%
      str_replace("NCLm\\s+_", "NCLm_") %>%
      str_remove("\\s*\\(RGB\\)$"),
    band = as.integer(str_extract(fname, "^NCL_(\\d+)") %>% str_remove("^NCL_")),
    
    subregion = str_match(fname, "_NCL([A-Za-z])_")[,2],
    
    hemisphere = str_match(fname, "_\\d+-\\d+([A-Za-z])\\.tif$")[,2],

    slide = as.integer(str_match(fname, "_(\\d+)-\\d+[A-Za-z]\\.tif$")[,2]),
    slice = as.integer(str_match(fname, "_\\d+-(\\d+)[A-Za-z]\\.tif$")[,2])
  ) %>%
  arrange(band, subregion, hemisphere, slide, slice)

NCLsorted <- NCLsorted %>%
  select(-filename, -fname)

NCLmeans <- NCLsorted %>%
  group_by(band, subregion) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop"
  )

NCLmeans <- NCLmeans %>%
  mutate(subregion = dplyr::recode(subregion,
                            "d" = "NCLd",
                            "m" = "NCLm",
                            "v" = "NCLv"
  ))

#export ncl data to corrmatrixdf
corrmatrixdf <- corrmatrixdf %>%
  left_join(
    NCLmeans %>%
      select(band, subregion, mean_count) %>%
      pivot_wider(
        names_from  = subregion,
        values_from = mean_count
      ),
    by = "band"
  )

#Import Cluster N data
clusternmean <- read_excel(
  path  = zenk_path,
  sheet = "ClusterN_Pooled")   # exact sheet name

clusternmean <- clusternmean %>%
  rename(ClN = mean_cell_count)

corrmatrixdf <- corrmatrixdf %>%
  left_join(clusternmean %>% select(band, ClN), by = "band") %>%
  arrange(band)

# Correlation matrix (Fig 5)

X  <- corrmatrixdf %>% select(where(is.numeric), -band) %>% as.matrix()
rc <- Hmisc::rcorr(X, type = "pearson")

var_order <- colnames(X)
n_var     <- length(var_order)

var_lab <- c("Anterior Dorsal Hippocampus","Anterior Dorsolateral Hippocampus",
             "Anterior Ventromedial Hippocampus",
             "Intermediate Dorsal Hippocampus","Intermediate Dorsolateral Hippocampus",
             "Intermediate Ventromedial Hippocampus",
             "Posterior Dorsal Hippocampus","Posterior Dorsolateral Hippocampus",
             "Posterior Ventromedial Hippocampus",
             "Dorsal NCL","Medial NCL","Ventral NCL","Cluster N")
names(var_lab) <- var_order

# 78 unique pairs, with BH-adjusted p values
pairs_df <- as.data.frame(rc$r) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
  left_join(
    as.data.frame(rc$P) %>%
      rownames_to_column("var1") %>%
      pivot_longer(-var1, names_to = "var2", values_to = "p"),
    by = c("var1", "var2")
  ) %>%
  mutate(i1 = match(var1, var_order), i2 = match(var2, var_order)) %>%
  filter(i1 > i2) %>%                                   # lower triangle, no diagonal
  mutate(q_BH    = p.adjust(p, method = "BH"),
         fdr_sig = q_BH < 0.05,
         var1    = var_lab[var1],
         var2    = var_lab[var2]) %>%
  select(var1, var2, r, p, q_BH, fdr_sig)

# ---- add the r = 1 diagonal, build plotting grid -------------
diag_df <- tibble(var1 = var_lab, var2 = var_lab, r = 1, fdr_sig = NA)

tiles <- bind_rows(pairs_df, diag_df) %>%
  mutate(i1 = match(var1, var_lab), i2 = match(var2, var_lab),
         x = pmin(i1, i2), y = -pmax(i1, i2),
         star = ifelse(!is.na(fdr_sig) & fdr_sig, "*", ""),
         cell = paste0(sprintf("%.2f", r), star))

# ---- figure -------------------------------------------------
pal_coral <- c("#FFFDF7", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#E36243")
figure5 <- ggplot(tiles, aes(x = x, y = y)) +
  geom_tile(aes(fill = r), colour = "white", linewidth = 0.9) +
  geom_text(aes(label = cell), size = 2.6) +
  scale_fill_gradientn(colours = pal_coral, limits = c(0, 1),
                       breaks = seq(0, 1, 0.1),
                       labels = sprintf("%.1f", seq(0, 1, 0.1)), name = NULL) +
  scale_x_continuous(breaks = 1:n_var, labels = var_lab,
                     limits = c(0.5, n_var + 0.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = -(1:n_var), labels = var_lab,
                     limits = c(-(n_var + 0.5), -0.5), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  labs(title = "",
       x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barheight = unit(78, "mm"), barwidth = unit(4, "mm"),
                                frame.colour = "grey70", frame.linewidth = 0.3,
                                ticks.colour = "white")) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, colour = "black"),
        axis.text.y = element_text(size = 8.5, colour = "black"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 10, 6, 10))

figure5

# Principal components (Tables 1-2)

pca_df <- corrmatrixdf %>%
  drop_na()

pca_res <- prcomp(pca_df[, colnames(corrmatrixdf) != "band"],
                  scale. = TRUE)

summary(pca_res)
pca_res$sdev^2
pca_res$rotation
eigenvalues <- pca_res$sdev^2
eigenvalues

pca_data <- pca_df %>%
  select(-band) %>%
  scale()

ekc <- EKC(pca_data)
ekc

scores <- as.data.frame(pca_res$x)

scores$band <- pca_df$band

#PC loadings reported in Table 2

round(pca_res$rotation[, "PC2"], 2)
round(pca_res$rotation[, "PC3"], 2)

#Multiple regression (Fig 6, S1 Table)

reg_data <- scores %>%
  dplyr::left_join(MoveTime %>% select(band, Movement, PerMove_4, PerMove_430, PerMove_5, PerMove_530), by = "band")

modelpc1 <- lm(PC1 ~ PerMove_4 + PerMove_430 + PerMove_5 + PerMove_530, data = reg_data)
modelpc2 <- lm(PC2 ~ PerMove_4 + PerMove_430 + PerMove_5 + PerMove_530, data = reg_data)

summary(modelpc1)
summary(modelpc2)

# assumption checks: residuals vs fitted, Q-Q, scale-location, residuals vs leverage
par(mfrow = c(2, 2))
plot(modelpc1)
plot(modelpc2)
par(mfrow = c(1, 1))

# collinearity among the four behavioural predictors (reported in Methods)
car::vif(modelpc1)

# full coefficient tables with 95% CIs (S1 Table)
cbind(summary(modelpc1)$coefficients, confint(modelpc1))
cbind(summary(modelpc2)$coefficients, confint(modelpc2))

plot_pc1 <- reg_data %>%
  select(PC1, PerMove_4, PerMove_430, PerMove_5, PerMove_530) %>%
  pivot_longer(
    cols = starts_with("PerMove"),
    names_to = "Time",
    values_to = "Movement"
      )

    # relabel time
    plot_pc1$Time <- dplyr::recode(
      plot_pc1$Time,
      "PerMove_4"   = "4:00–4:30 pm",
      "PerMove_430" = "4:30–5:00 pm",
      "PerMove_5"   = "5:00–5:30 pm",
      "PerMove_530" = "5:30–6:00 pm"
    )
    
    
    # keep panels ordered
    plot_pc1$Time <- factor(
      plot_pc1$Time,
      levels = c(
        "4:00–4:30 pm",
        "4:30–5:00 pm",
        "5:00–5:30 pm",
        "5:30–6:00 pm"
      )
    )
    
    # color: last panel blue, others orange
    plot_pc1$line_color <- ifelse(
      plot_pc1$Time == "5:30–6:00 pm",
      "blue",
      "orange"
    )
    
    # extract stats from PC1 multiple regression model
    # marginal (simple) correlations, one per time period — these match the fitted lines
    coef_df1 <- plot_pc1 %>%
      group_by(Time) %>%
      summarise(
        r = cor(Movement, PC1),
        p = cor.test(Movement, PC1)$p.value,
        .groups = "drop"
      ) %>%
      mutate(
        label = paste0("r = ", round(r, 2),
                       "\np = ", signif(p, 2))
      )

    # plot PC1
    figure6 <- ggplot(plot_pc1, aes(x = Movement, y = PC1)) +
      geom_point(size = 2) +
      geom_smooth(aes(color = line_color), method = "lm", se = FALSE) +
      geom_text(
        data = coef_df1,
        aes(x = Inf, y = Inf, label = label),
        hjust = 1.1,
        vjust = 1.2,
        size = 3,
        inherit.aes = FALSE
      ) +
      scale_color_identity() +
      facet_wrap(~ Time, scales = "free_x") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        strip.text = element_text()
      ) +
      labs(
        x = "Movement or Migratory Restlessness (%)",
        y = "PC1"
      )

# Repeated-measures correlation (Fig 7, S2 Table)

MRtime <- read.csv(
  mr_path)

MRtime$Date <- as.factor(MRtime$Date)
MRtime$Time <- as.factor(MRtime$Time)
MRtime$Bird.ID <- as.factor(MRtime$Bird.ID)
MRtime$Treatment <- as.factor(MRtime$Treatment)

MRtime <- MRtime %>%
  filter(Treatment %in% c("Control", "Warm", "Cold"))

MRtime$Treatment <- droplevels(MRtime$Treatment)

sub <- subset(MRtime, Time == "17:30")
sub$Bird.ID <- droplevels(factor(sub$Bird.ID))

# repeated measures correlation
rmc <- rmcorr(participant = Bird.ID,
              measure1 = Percentage.MR,
              measure2 = PerMove,
              dataset = sub)

# common within-bird slope, anchored at the grand mean for plotting
slope     <- unname(coef(lm(PerMove ~ Bird.ID + Percentage.MR,
                            data = sub))["Percentage.MR"])
intercept <- mean(sub$PerMove) - slope * mean(sub$Percentage.MR)

coef_df <- data.frame(
  label = paste0("r = ", round(rmc$r, 2),
                 "\np = ", signif(rmc$p, 1)))

pred <- data.frame(
  Percentage.MR = seq(min(sub$Percentage.MR), max(sub$Percentage.MR),
                      length.out = 100))
pred$PerMove <- intercept + slope * pred$Percentage.MR

figure7 <- ggplot(sub, aes(x = Percentage.MR, y = PerMove)) +
  geom_point(size = 2) +
  geom_line(data = pred, aes(color = "lightblue"), linewidth = 1) +
  geom_text(
    data = coef_df,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.2, size = 3,
    inherit.aes = FALSE
  ) +
  scale_color_identity() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Migratory Restlessness from 5:30 - 6:00 pm (%)",
       y = "Movement from 5:00 - 5:30 pm (%)")

#Save figures

ggsave(file.path(figdir, "Fig6.tif"), plot = figure6, device = "tiff", width = 6.5, 
       height = 5.2, units = "in", dpi = 300, bg = "white", compression = "lzw")

ggsave(file.path(figdir, "Fig5.tif"), plot = figure5, device = "tiff",
       width = 7.5, height = 7.5, units = "in", dpi = 300,
       bg = "white", compression = "lzw")

ggsave(file.path(figdir, "Fig7.tif"), plot = figure7, device = "tiff", width = 5.5, height = 4.4,
       units = "in", dpi = 300, bg = "white", compression = "lzw")

