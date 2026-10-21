library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom)

# ---- 1. load data ----
paired             <- readRDS("data/clean/paired_hmo_metagenome.rds")
hmo_percent_unique <- readRDS("data/clean/hmo_percent_unique.rds")

# ---- 2. set up variables ----
paired <- paired |> mutate(dpd = days)   # dpd = days post delivery

hmo_cols <- setdiff(names(hmo_percent_unique),
                    c("time", "name", "enrollment_id", "secretor",
                      "diversity", "evenness", "days", "days_from_61", "sum"))

hmo_cols
sapply(paired[hmo_cols], class)   # should all say "numeric"

country_colors <- c("Panama" = "#6BB8E0",   # light blue
                    "Peru"   = "#EE8DB3",   # light pink
                    "US"     = "#A98BD6")   # light purple

paired_long <- paired |>
  pivot_longer(all_of(hmo_cols), names_to = "hmo", values_to = "percent")

# ---- 3. one model by hand ----
fit <- lm(x2fl ~ 1 + country + dpd, data = paired)
summary(fit)

# ---- 4. plots ----
# one HMO by country
p_one <- ggplot(paired, aes(x = country, y = x2fl, fill = country)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  scale_fill_manual(values = country_colors) +
  labs(y = "2'FL (%)") +
  theme(legend.position = "none")
p_one

# one HMO by days post delivery
p_dpd <- ggplot(paired, aes(x = dpd, y = x2fl, color = country)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = country_colors)
p_dpd

# all HMOs by country
p_country <- ggplot(paired_long, aes(x = country, y = percent, fill = country)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ hmo, scales = "free_y") +
  scale_fill_manual(values = country_colors) +
  theme(legend.position = "none")
p_country

# ---- 5. regression for every HMO ----
fits <- map(hmo_cols, ~ lm(reformulate(c("country", "dpd"), response = .x), data = paired))
names(fits) <- hmo_cols

results <- map_dfr(fits, tidy, .id = "hmo") |>
  filter(term != "(Intercept)") |>
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

results |> arrange(p_adj)

# ---- 6. Coefficient plot ----
p_coef <- results |>
  filter(term != "dpd") |>
  mutate(lower = estimate - 1.96 * std.error,
         upper = estimate + 1.96 * std.error) |>
  ggplot(aes(x = estimate, y = hmo, color = term)) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Difference from Panama (percentage points)", y = NULL, color = "Country")
p_coef


# ---- 7. save figures ----
dir.create("output", showWarnings = FALSE)
ggsave("output/x2fl_by_country.png",  p_one,     width = 6,  height = 5)
ggsave("output/x2fl_by_dpd.png",      p_dpd,     width = 7,  height = 5)
ggsave("output/hmo_by_country.png",   p_country, width = 10, height = 8)
ggsave("output/hmo_coefficients.png", p_coef, width = 8, height = 7)