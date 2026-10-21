library (vegan)
library (dplyr)
library(ggplot2)

# ---- 1. loading data ----
paired <- readRDS ("data/clean/paired_hmo_metagenome.rds")
hmo_percent_unique <- readRDS ("data/clean/hmo_percent_unique.rds")

# ---- 2. set up variables ----
paired <- paired |> mutate(dpd = days)

hmo_cols <- setdiff(names(hmo_percent_unique),
                    c("time", "name", "enrollment_id", "secretor",
                      "diversity", "evenness", "days", "days_from_61", "sum"))

# ---- 3. data setup ----
dat <- paired |>
  filter(!is.na(seroconverted)) |> #removes 7 pairs with no response value, leaves 296
  mutate(seroconverted = factor(seroconverted)) # tells R that responded/non-responder are two categories

hmo_mat <- as.matrix(dat[, hmo_cols]) #pulls out only sugar columns
rownames(hmo_mat) <- dat$enrollment_id # labels each row with the infant's ID so you can trace rows later

# ---- snaity checks ----
dim(hmo_mat)                   # expect 296 rows
table(dat$seroconverted)       # how many in each group
rowSums(hmo_mat) |> summary()  # should be close to 100 if these are percentages
sum(hmo_mat == 0)              # how many sugars were absent (zeros)

# ---- 4. accounting for 'distances' before running permanova
set.seed(123)
dist_bc <- vegdist(hmo_mat, method = "bray")

# ---- 5. running permanova ----
perm <- adonis2(dist_bc ~ seroconverted, data = dat,
                permutations = 9999, by = "margin")
perm

# ---- 6. spread check ----
bd <- betadisper(dist_bc, dat$seroconverted)
permutest(bd, permutations = 9999)

# ---- 7. db-RDA ----
cap_country <- capscale(dist_bc ~ seroconverted + Condition(country), data = dat)
anova(cap_country, permutations = 9999)
RsquareAdj(cap_country)

#---- ordination plot----
pcoa <- cmdscale(dist_bc, k = 2)
pc_df <- data.frame(PC1 = pcoa[, 1], PC2 = pcoa[, 2],
                    response = dat$seroconverted, country = dat$country)

ggplot(pc_df, aes(PC1, PC2, color = response)) + geom_point(alpha = 0.7)
ggplot(pc_df, aes(PC1, PC2, color = country)) + geom_point(alpha = 0.7) +
  scale_color_manual(values = country_colors)

p_pcoa_response <- ggplot(pc_df, aes(PC1, PC2, color = response)) + geom_point(alpha = 0.7)
ggsave("output/pcoa_by_response.png", p_pcoa_response, width = 6, height = 5)




