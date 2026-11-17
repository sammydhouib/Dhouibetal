library(readxl)
library(janitor)
library(ggplot2)
library(rstatix)
library(stringr)
library(readr)
library(tidyverse)
library(dplyr)
library(export)
library(afex)  
library(Hmisc)
library(scales)

#Import Hp data
HpRaw <- read_excel(
  path  = "/Users/sammydhouib/Desktop/Honors Thesis/Data/WTSP_Data.xlsx",
  sheet = "HippocampusRaw")  %>% # exact sheet name 
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
  path  = "/Users/sammydhouib/Desktop/Honors Thesis/Data/WTSP_Data.xlsx",
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
HpMean_collapsed

Hp_subregionsonly <- HpMean_collapsed %>%
  filter(subregion %in% c("dl", "vm", "d")) %>%
  select(band, caudrost, subregion, mean_count) %>%
  pivot_wider(
    names_from = subregion,
    values_from = mean_count
  ) %>%
  arrange(band)

HpMean_caudrost <- HpMean_collapsed %>%
  group_by(band, caudrost) %>%
  summarise(
    mean_count = mean(mean_count, na.rm = TRUE),
    .groups = "drop"
  )

HpMean_collapsed <- HpMean_collapsed %>%
  mutate(region = paste(caudrost, subregion, sep = "_")) %>%
  select(-caudrost, -subregion) %>%
  relocate(region, .before = 2)
HpMean_collapsed

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
  path  = "/Users/sammydhouib/Desktop/Honors Thesis/Data/WTSP_Data.xlsx",
  sheet = "NCLRaw")
NCL_raw

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

NCLsorted

NCLmeans <- NCLsorted %>%
  group_by(band, subregion) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop"
  )

NCLmeans

NCL_collapsed <- NCLmeans %>%
  select(band, subregion, mean_count) %>%
  pivot_wider(
    names_from = subregion,
    values_from = mean_count
  )

NCLmeans <- NCLmeans %>%
  mutate(subregion = recode(subregion,
                            "d" = "NCLd",
                            "m" = "NCLm",
                            "v" = "NCLv"
  ))

NCL_MR_joined <- NCL_collapsed %>%
  left_join(MoveTime, by = "band")

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
  path  = "/Users/sammydhouib/Desktop/Honors Thesis/Data/WTSP_Data.xlsx",
  sheet = "ClusterN_Pooled")   # exact sheet name
clusternmean

clusternmean <- clusternmean %>%
  rename(ClN = mean_cell_count)

clustern_joined <- MoveTime %>%
  left_join(clusternmean, by = "band") %>%
  arrange(band)

#export ClN to corrmatrixdf
corrmatrixdf <- corrmatrixdf %>%
  left_join(clusternmean %>% select(band, ClN), by = "band") %>%
  arrange(band)

#Correlation matrix
library(tidyverse)
library(Hmisc)

# ============================================================
# Correlation matrix: ZENK across hippocampus, NCL, Cluster N
# Requires: corrmatrixdf (one row per bird, numeric columns + band)
# ============================================================

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

# ---- 78 unique pairs, with BH-adjusted p values -------------
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
corr_plot <- ggplot(tiles, aes(x = x, y = y)) +
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

corr_plot

#Start pca with existing corrmatrix df
pca_df <- corrmatrixdf %>%
  drop_na()

pca_res <- prcomp(pca_df[, colnames(corrmatrixdf) != "band"],
                  scale. = TRUE)

summary(pca_res)
pca_res$sdev^2
pca_res$rotation
eigenvalues <- pca_res$sdev^2
eigenvalues
#parallel analysis
library(psych)
pca_data <- pca_df %>%
  select(-band) %>%
  scale()

fa.parallel(
  pca_data,
  fa = "pc",
  n.iter = 1000,
  quant = 0.5,
  show.legend = TRUE,
  main = "Parallel Analysis"
)

#empirical kaiser criterion
library(EFAtools)
ekc <- EKC(pca_data)
ekc

#Scree plot
plot(
  pca_res,
  type = "l",
  main = "Scree Plot"
)

#arrow plot
biplot(
  pca_res,
  scale = 0,
  cex = 0.7
)

#PC scores for each PC, for each bird 
scores <- as.data.frame(pca_res$x)

ggplot(scores, aes(PC1, PC3)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA: PC1 vs PC3",
    x = "PC1",
    y = "PC2"
  )

#add back band column that was removed during PCA to run correlation by band ID 
#(could do this with MR but band is easier on the eyes)
scores$band <- pca_df$band

#Sanity check that bands are added correctly
scores
pca_df

# Eigenvalues and kaiser criterion
eigenvalues <- pca_res$sdev^2
eigenvalues
kaiser_pcs <- which(eigenvalues > 1)
length(kaiser_pcs)

#add MR to the scores dataset and run correlation
reg_data <- scores %>%
  dplyr::left_join(MoveTime %>% select(band, Movement, PerMove_4, PerMove_430, PerMove_5, PerMove_530), by = "band")

reg_data_new <- reg_data %>%
  mutate(treatment = rowMeans(dplyr::select(., PerMove_4, PerMove_430), na.rm = TRUE)) %>%
  dplyr::select(-PerMove_4, -PerMove_430)

#Run models
modelpc1 <- lm(PC1 ~ PerMove_4 + PerMove_430 + PerMove_5 + PerMove_530, data = reg_data)
modelpc2 <- lm(PC2 ~ PerMove_4 + PerMove_430 + PerMove_5 + PerMove_530, data = reg_data)


ggplot(reg_data, aes(x = PerMove_5, y = PerMove_530)) +
  geom_point() +
  labs(
    x = "PerMove_5",
    y = "PerMove_530"
  ) +
  theme_classic()

cor.test(reg_data$PerMove_5, reg_data$PerMove_530, method = "spearman", exact = FALSE)

summary(modelpc1)
summary(modelpc2)


#assumption check of the MODEL
par(mfrow=c(2,2))
plot(modelpc1)
termplot(modelpc1, partial.resid = TRUE, se = TRUE)

#see which brain regions explain variance in PC1 and PC2
round(pca_res$rotation[, "PC1"], 2)
round(pca_res$rotation[, "PC2"], 2)
round(pca_res$rotation[, "PC3"], 2)
round(pca_res$rotation[, "PC4"], 2)


library(broom)

plot_data <- reg_data %>%
  select(PerMove_4, PerMove_430, PerMove_5, PerMove_530, PC1, PC2, PC3) %>%
  pivot_longer(
    cols = c(PerMove_4, PerMove_430, PerMove_5, PerMove_530),
    names_to = "Time",
    values_to = "Movement"
  ) %>%
  pivot_longer(
    cols = c(PC1, PC2, PC3),
    names_to = "Predictor",
    values_to = "Score"
  )

ggplot(plot_data, aes(x = Score, y = Movement)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  facet_grid(Time ~ Predictor, scales = "free_x") +
  labs(
    x = "PC score",
    y = "Movement (%)"
  ) +
  theme_classic()

#PC1
plot_pc1 <- reg_data %>%
  select(PC1, PerMove_4, PerMove_430, PerMove_5, PerMove_530) %>%
  pivot_longer(
    cols = starts_with("PerMove"),
    names_to = "Time",
    values_to = "Movement"
      )
    
    # relabel time
    plot_pc1$Time <- recode(
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
    coef_df1 <- tidy(modelpc1) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        partial_r = statistic / sqrt(statistic^2 + df.residual(modelpc1)),
        Time = recode(
          term,
          "PerMove_4"   = "4:00–4:30 pm",
          "PerMove_430" = "4:30–5:00 pm",
          "PerMove_5"   = "5:00–5:30 pm",
          "PerMove_530" = "5:30–6:00 pm"
        ),
        Time = factor(
          Time,
          levels = c("4:00–4:30 pm", "4:30–5:00 pm",
                     "5:00–5:30 pm", "5:30–6:00 pm")
        ),
        label = paste0(
          "r = ", round(partial_r, 2),
          "\np = ", signif(p.value, 2)
        )
      )
    
    # plot PC1
    ggplot(plot_pc1, aes(x = Movement, y = PC1)) +
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

#Follow-up
reg_data <- reg_data %>%
  mutate(low_move = factor(PerMove_530 <= 15, levels = c(TRUE, FALSE),
                           labels = c(" <= 15%", "> 15%")))

t.test(PC1 ~ low_move, data = reg_data)

ggplot(reg_data, aes(x = low_move, y = PC1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 2) +
  labs(x = "PerMove_530", y = "PC1") +
  theme_classic()

#rmANOVA
areas <- setdiff(colnames(corrmatrixdf), "band")

# long format: one row per bird x brain area
long15 <- corrmatrixdf %>%
  left_join(MoveTime %>% select(band, PerMove_530), by = "band") %>%
  mutate(migrate = factor(PerMove_530 > 15,
                          levels = c(FALSE, TRUE),
                          labels = c("No", "Yes"))) %>%
  pivot_longer(all_of(areas), names_to = "area", values_to = "zenk") %>%
  mutate(band = factor(band),
         area = factor(area, levels = areas))

# mixed ANOVA: area within, migrate between
fit15 <- aov_ez(id = "band", dv = "zenk", data = long15,
                between = "migrate", within = "area",
                anova_table = list(correction = "GG", es = "pes"))

fit15                                # GG-corrected table with partial eta-squared
summary(fit15)$sphericity.tests      # Mauchly's test

# profile plot
summ15 <- long15 %>%
  group_by(area, migrate) %>%
  summarise(mean = mean(zenk), se = sd(zenk)/sqrt(n()), .groups = "drop")

ggplot(summ15, aes(x = area, y = mean, colour = migrate, group = migrate)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.18) +
  labs(x = "Brain area", y = "Mean ZENK count (± SE)", colour = "Migrant") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#MR Histogram
ggplot(reg_data, aes(x = PerMove_530)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "black", fill = "grey70") +
  geom_vline(xintercept = 15, colour = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x = "MR 5:30-6:00 pm (%)", y = "Number of birds") +
  theme_classic()

#MR vs. Quiescence
library(emmeans)
library(readxl)
library(rmcorr)

MRtime <- read.csv(
  "/Users/sammydhouib/Desktop/Manuscript/MM:MR Linear Regression/Migratory restlessness raw data.csv")
Movement <- read_xlsx(
  "/Users/sammydhouib/Desktop/Manuscript/MM:MR Linear Regression/4-5pm observation data.xlsx",
  sheet = "Sheet1")

MRtime$Date <- as.factor(MRtime$Date)
MRtime$Time <- as.factor(MRtime$Time)
MRtime$Bird.ID <- as.factor(MRtime$Bird.ID)
MRtime$Treatment <- as.factor(MRtime$Treatment)

MRtime <- MRtime %>%
  filter(Treatment %in% c("Control", "Warm", "Cold"))

MRtime$Treatment <- droplevels(MRtime$Treatment)

str(MRtime)

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

rmcorrplot <- ggplot(sub, aes(x = Percentage.MR, y = PerMove)) +
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

rmcorrplot


