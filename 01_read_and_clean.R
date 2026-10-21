# reading and cleaning script

library(readxl)
library(readr)
library(janitor) #sheet cleaning
library(dplyr)
library(writexl) #export back to excel
library(stringr) # for str_extract
library(tidyr)# for dropping

# ---- 1. microbiome data ----
microbiome <- read_excel("data/raw/microbiome-ogu-table.xlsx")
glimpse(microbiome)
microbiome <- clean_names(microbiome)

# ---- 2. metegenome data ----
metagenome <- read_excel ("data/raw/metagenome-metadata.xlsx")
glimpse (metagenome)
metagenome <- clean_names(metagenome)

# ---- 3. hmo pecent data ----
hmo_percent <- read_tsv ("data/raw/hmo-percent.tsv.xls")
glimpse (hmo_percent)
hmo_percent <- clean_names(hmo_percent)



# --- 4. dropping 'NA' data in metagenome file ----

#metagenome file
dim(metagenome)
names(metagenome)[1:15]
sum(is.na(metagenome))
colSums(is.na(metagenome)) |> sort(decreasing = TRUE) |> head(30)

# checking if NAs all in the same rows
metagenome |>
  filter(is.na(enrollment_id)) |>
  summarise(n = n(),
            all_missing_country = all(is.na(country)),
            all_missing_visit   = all(is.na(visit)),
            all_missing_age     = all(is.na(age)))

# sample names:
metagenome |> filter(is.na(enrollment_id)) |> pull(sample_name) |> head(20)

metagenome |>
  filter(is.na(enrollment_id)) |>
  count(prefix = str_extract(sample_name, "^[A-Za-z]+"))

metagenome |>
  filter(is.na(enrollment_id), str_starts(sample_name, "P")) |>
  pull(sample_name)

# removed 93 non-participant controls (extraction blanks, PCR blanks, zymo mock communities and zymo PCR controls)
nrow(metagenome)   # 928
metagenome <- metagenome |> drop_na(enrollment_id, country, visit, age)
nrow(metagenome)   # expect 835

# ---- 5. Removing anything from hmo percents file if needed ----
dim(hmo_percent)
sum(is.na(hmo_percent))
colSums(is.na(hmo_percent)) |> sort(decreasing = TRUE) |> head(15)
hmo_percent |> filter(if_any(everything(), is.na)) |> View()

# no data in the 'time' column
hmo_percent <- hmo_percent |> select(-time)

# ---- 6. parsing 'name' into 'enrollment_id' to match metagenome file
head(hmo_percent$name, 10)
head(metagenome$enrollment_id, 10)

# parsing hmo names
hmo_percent <- hmo_percent |>
  mutate(enrollment_id = paste0("RB", str_pad(str_split_i(name, "\\.", 2), width = 3, pad = "0")))

# testing first few
hmo_percent |> select(name, enrollment_id) |> head(10)

# testig the IDs are clean before joining
n_distinct(hmo_percent$enrollment_id)       
sum(duplicated(hmo_percent$enrollment_id))  # 3 duplicate ids
table(nchar(hmo_percent$enrollment_id))    

#looking at the 3 duplicates:
hmo_percent |>
  filter(enrollment_id %in% enrollment_id[duplicated(enrollment_id)]) |>
  arrange(enrollment_id) |>
  select(name, enrollment_id)

#checking duplicate table values
dups <- hmo_percent |>
  filter(enrollment_id %in% c("RB051", "RB067", "RB072")) |>
  arrange(enrollment_id)

t(dups |> select(-name, -enrollment_id))

nrow(hmo_percent)                  # 345
nrow(distinct(hmo_percent))        # 343 if the two identical-name pairs are exact copies

# filtering the metagenome data to 2 months
table(metagenome$visit)                     # confirm the exact label
metagenome_2mo <- metagenome |> filter(visit == "2mo")
nrow(metagenome_2mo)
sum(duplicated(metagenome_2mo$enrollment_id))

#DEDUPLICATING HERE: keep sample with 61 days and first occurrences of the last two duplicate entries
hmo_percent <- hmo_percent |>
  mutate(days = as.numeric(str_split_i(name, "\\.", 3)),
         days_from_61 = abs(days - 61))

hmo_percent_unique <- hmo_percent |>
  group_by(enrollment_id) |>
  slice_min(days_from_61, n = 1, with_ties = FALSE) |>
  ungroup()

nrow(hmo_percent_unique)   # expect 342
sum(duplicated(hmo_percent_unique$enrollment_id))   # expect 0

#check cases
hmo_percent_unique |>
  filter(enrollment_id %in% c("RB051", "RB067", "RB072")) |>
  select(name, enrollment_id, days)

# ---- ALL PREP TO JOIN FILES (metadata and HMO)----
# duplicate check in the meta genome (number of mother infant pairs that have both layers):
length(intersect(hmo_percent_unique$enrollment_id, metagenome_2mo$enrollment_id))

#remaining 39 that don't match:
unmatched <- setdiff(hmo_percent_unique$enrollment_id, metagenome_2mo$enrollment_id)
length(unmatched)   # expect 39
head(unmatched)

# checking if those infant appear in the metagenome file at other visits:
sum(unmatched %in% metagenome$enrollment_id)

 #joining files
paired <- inner_join(hmo_percent_unique, metagenome_2mo, by = "enrollment_id")
nrow(paired)   # should be 303

#seroconversion check:
sum(is.na(paired$seroconverted))
table(paired$seroconverted, useNA = "ifany")
table(paired$country)

#save cleaned objects:
dir.create("data/clean", showWarnings = FALSE)
saveRDS(paired,               "data/clean/paired_hmo_metagenome.rds")
saveRDS(hmo_percent_unique,   "data/clean/hmo_percent_unique.rds")
saveRDS(metagenome_2mo,       "data/clean/metagenome_2mo_clean.rds")
saveRDS(microbiome, "data/clean/microbiome_clean.rds")


