library(tidyverse)
library(reshape2)
library(cowplot)
source("status_allocation_functions.r")
source("cmm2d-new-functions.r")
install.packages("cowplot")

# Data preparation ----

# Reading in the data. The file contains the Eurostat table `lfsa-egised`, downloaded from the Eurostat's website on May 2, 2024.

p <- read_delim(file = "estat_lfsa_egised.tsv", delim = "\t") 

# Reshaping the data to the long format

p <- p |>
  pivot_longer(cols = -1, names_to = "year") |>
  separate_wider_delim(cols = 1, names = c("freq", "age", "sex", "isco08", "isced11", "unit",
                                           "cnt"), delim = ",") |>
  mutate(value = str_remove_all(string = value, pattern = "[a-z]|\\s|\\:"))

# Filtering out the subset of interest

d <- p |>
  filter(age != "Y20-64",
         sex != "T",
         isco08 %in% paste0("OC", 1:9),
         isced11 %in% c("ED0-2", "ED3_4", "ED5-8"),
         !str_detect(string = cnt, pattern = "E[AU]")) |>
  mutate(across(.cols = c(year, value), .fns = as.double))

# Some further recoding

d <- d |>
  mutate(degree = factor(isced11, levels = c("ED5-8", "ED3_4", "ED0-2"), 
                         labels = c("High", "Medium", "Low")),
         status = case_when(
           isco08 %in% paste0("OC", 1:3) ~ "S1",
           isco08 %in% paste0("OC", 4:5) ~ "S2",
           isco08 %in% paste0("OC", 6:7) ~ "S3",
           isco08 %in% paste0("OC", 8:9) ~ "S4"),
         status = factor(status))

# The final tibble ----

t1 <- d |>
  count(year, cnt, sex, degree, status, wt = value, name = "freq") |>
  mutate(origin = paste(str_sub(string = degree, start = 1, end = 1), sex, sep = "-")) |>
  group_by(year, cnt) |>
  nest()

# Removing 'empty' country-year combinations (e.g., the LFS survey was not carried out in AT in 1983). Counts (frequencies) in status allocation tables corresponding to the empty combination sum up to 0

t1 <- t1 |>
  mutate(test = map_dbl(data, ~sum(.$freq))) |>
  filter(test != 0) |>
  select(-test) |>
  ungroup()

# One-dimensional analysis: education-based meritocracy ----

t1 <- t1 |>
  mutate(tb1ded = map(data, ~xtabs(freq ~ degree + status, data = .)),
         sa1ded = map(tb1ded, mar),
         cm1ded = map(sa1ded, merit1),
         dm1ded = map(sa1ded, merit2))

# One-dimensional analysis: gender-based meritocracy ----

t1 <- t1 |>
  mutate(tb1dsx = map(data, ~xtabs(freq ~ sex + status, data = .)),
         sa1dsx = map(tb1dsx, mar),
         cm1dsx = map(sa1dsx, merit1),
         dm1dsx = map(sa1dsx, merit2))

# One-dimensional analysis: gender-based meritocracy ----

t1 <- t1 |>
  mutate(tb1dsx_alt = map(tb1dsx, ~.[2:1, ]),
         sa1dsx_alt = map(tb1dsx_alt, mar),
         cm1dsx_alt = map(sa1dsx_alt, merit1),
         dm1dsx_alt = map(sa1dsx_alt, merit2))

# Two-dimensional analysis

t1 <- t1 |>
  mutate(sa2d = map(data, ~stall2d(d = ., f = "freq", dim1 = "degree", dim2 = "sex", status = "status")),
         cm2d = map(sa2d, cmm2d),
         dm2d = map(sa2d, dmm2d))

# Two-dimensional analysis --- alternate take: men over women

t1 <- t1 |>
  mutate(data = map(data, ~mutate(.data = ., sex1 = fct_relevel(sex, "M", "F"))),
         sa2d_alt = map(data, ~stall2d(d = ., f = "freq", dim1 = "degree", dim2 = "sex1", status = "status")),
         cm2d_alt = map(sa2d_alt, cmm2d),
         dm2d_alt = map(sa2d_alt, dmm2d))

save(t1, file = "js-zk-multdim-estat-results.RData")

# Comparison of gender-based models in terms of fit

t1 |>
  select(year, cnt, starts_with("cm1dsx"), starts_with("dm1dsx")) |>
  mutate(across(.cols = 3:4, .fns = ~map_dbl(., ~.[[3]])),
         across(.cols = 5:6, .fns = ~map_dbl(., ~.[[4]]))) |>
  pivot_longer(cols = 3:6, names_sep = "_", names_to = c("model", "ranking"),
               names_transform = list(ranking = ~replace(., is.na(.), "std"))) |>
  ggplot(mapping = aes(x = model, y = value, colour = ranking)) + 
  geom_boxplot()

# Comparison of the two 2-dim models in terms of fit

t1 <- t1 |>
  mutate(delta_cmm = map_dbl(.x = cm2d, .f = ~.$`Dissimilarity index`),
         delta_cmm_alt = map_dbl(.x = cm2d_alt, .f = ~.$`Dissimilarity index`),
         delta_dmm = map_dbl(.x = dm2d, .f = ~.$`Dissimilarity index`),
         delta_dmm_alt = map_dbl(.x = dm2d_alt, .f = ~.$`Dissimilarity index`)) |>
  relocate(starts_with("delta"), .before = "data")

t1 |>
  select(starts_with("delta")) |>
  pivot_longer(cols = starts_with("delta"),
               names_to = c("model", "ordering"),
               names_sep = 10,
               values_to = "delta",
               names_transform = list(
                 model = ~ifelse(str_detect(string = ., pattern = "cmm"), "Constant", "Differential"),
                 ordering = ~ifelse(. == "alt", "Men over women", "Women over men")
               )) |>
  ggplot(mapping = aes(x = ordering, y = delta, colour = model)) + 
  geom_boxplot() + 
  scale_colour_manual(values = c("#762a83", "#1b7837")) + 
  labs(x = NULL, y = "Dissimilarity index", colour = "Model") + 
  theme_minimal() + 
  theme(legend.position = "top")

# Distances between the various types of referent allocations and the observed allocation

t1 <- t1 |>
  mutate(distances = map(.x = sa2d,
                         .f = ~map2_dbl(.x = .[-1], 
                                        .y = .[1], 
                                        .f = ~sum(abs(.y - .x))/(2 * sum(.y))
                                        )
                         )
         )

t1 <- t1 |>
  mutate(distances = map(.x = distances, .f = ~tibble(
    AM = .x[1], AU = .x[2], AV = .x[3], AL = .x[4]
  )))

t1 |>
  select(cnt, year, distances) |>
  unnest(cols = distances) |>
  pivot_longer(cols = starts_with("A")) |>
  mutate(name = factor(name, levels = c("AM", "AU", "AV", "AL"),
                       labels = paste("Observed to", c("MM", "ML", "LM", "LL")))) |>
  ggplot(mapping = aes(x = name, y = value)) + 
  geom_boxplot() + 
  labs(x = NULL, y = "Dissimilarity index") + 
  theme_minimal() 

t1

# Assessing the fit of the 2D model ----



z <- t1$cm2d |>
  map(~.$`Mixing coefficient`$Estimate) |>
  map_df(~tibble(
    term = c("Education", "Gender"),
    estimate = .
  ))

z <- z |>
  mutate(cnt = rep(t1$cnt, each = 2),
         year = rep(t1$year, each = 2))

z |>
  ggplot(mapping = aes(x = estimate)) + 
  geom_density() + 
  facet_wrap(~term)

z |>
  pivot_wider(names_from = term, values_from = estimate) |>
  ggplot(mapping = aes(x = Education, y = Gender)) + 
  geom_point()

z1 <- z |>
  pivot_wider(names_from = term, values_from = estimate) |>
  mutate(year = year - 1992)

m1 <- lm(Education ~ cnt + year, data = z1)
summary(m1)
m2 <- update(m1, .~. + Gender)
summary(m2)
m3 <- update(m2, .~. + year:Gender)
summary(m3)
anova(m2, m3)


t1 <- t1 |>
  mutate(dist = map(.x = sa2d, .f = ~map2(.x = .[1], 
                                              .y = .[-1], 
                                              .f = ~norm(.x - .y, type = "f")/sum(.x))))

t1 <- t1 |>
  mutate(dist = map(.x = dist, .f = unlist),
         dist = map_df(.x = dist, 
                       .f = ~tibble(AM = .[1], AU = .[2], AV = .[3], AL = .[4]))) |>
  unnest(cols = dist)

t1 |>
  select(year, cnt, AM, AU, AV, AL) |>
  pivot_longer(starts_with("A")) |>
  mutate(name = factor(name, levels = c("AM", "AU", "AV", "AL"))) |>
  ggplot(mapping = aes(x = name, y = value)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = NULL, y = "Distance") +
  scale_x_discrete(labels = c("Observed to\n$\\langle M_xM_y\\rangle$",
                              "Observed to\n$\\langle M_xL_y\\rangle$",
                              "Observed to\n$\\langle L_xM_y\\rangle$",
                              "Observed to\n$\\langle L_xL_y\\rangle$"))

# permutations ---- 

d2 <- d |>
  count(year, cnt, sex, degree, status, wt = values, name = "freq") |>
  mutate(origin = paste(str_sub(string = degree, start = 1, end = 1), sex, sep = "-")) |>
  group_by(year, cnt) |>
  nest()

d2 <- d2 |>
  mutate(hml = map(data, ~xtabs(freq ~ degree + status, data = .)),
         mhl = map(hml, ~.[c(2, 1, 3), ]),
         mlh = map(hml, ~.[c(2, 3, 1), ]),
         hlm = map(hml, ~.[c(1, 3, 2), ]),
         lhm = map(hml, ~.[c(3, 1, 2), ]),
         lmh = map(hml, ~.[c(3, 2, 1), ]))


d2 <- d2 |>
  mutate(across(.cols = hml:lmh, .fns = list(
    cmm = ~map(., ~merit1(mar(.))),
    dmm = ~map(., ~merit2(mar(.)))
  )))

d2 <- d2 |>
  mutate(across(.cols = ends_with("cmm"), .fns = list(
    di = ~map_dbl(., ~.[[3]])
  )))

d2 <- d2 |>
  mutate(across(.cols = ends_with("dmm"), .fns = list(
    di = ~map_dbl(., ~.[[4]])
  )))

d2 <- d2 |>
  select(ends_with("_di")) |>
  pivot_longer(cols = ends_with("_di")) |>
  mutate(name = str_remove(string = name, pattern = "_di")) |>
  separate_wider_delim(cols = name, names = c("ordering", "model"), delim = "_") |>
  ungroup()

d2 |>
  mutate(model = factor(model, levels = c("cmm", "dmm"), labels = c("Constant", "Differential"))) |>
  ggplot(mapping = aes(x = ordering, y = value)) + 
  geom_boxplot() + 
  facet_wrap(~model) + 
  theme_ipsum() + 
  labs(x = "Ordering of the categories of education",
       y = "Dissimilarity index")

# lex ----

di <- t1 |>
  mutate(across(.cols = starts_with("cm2d"), .fns = list(
    di = ~map_dbl(., ~.[[3]])
  ))) |>
  select(ends_with("_di")) |>
  ungroup()

di <- di |>
  pivot_longer(cols = ends_with("_di")) |>
  mutate(primary = ifelse(str_detect(string = name, pattern = "lex"), "Sex", "Education"),
         sexrank = ifelse(str_detect(string = name, pattern = "alt"), "Men", "Women")) 

|>
  ggplot(mapping = aes(x = sexrank, y = value, colour = primary)) + 
  geom_boxplot()

di2 <- t1 |>
  mutate(across(.cols = starts_with("dm2d"), .fns = list(
    di = ~map_dbl(., ~.[[3]])
  ))) |>
  select(ends_with("_di")) |>
  ungroup()

di2 <- di2 |>
  pivot_longer(cols = ends_with("_di")) |>
  mutate(primary = ifelse(str_detect(string = name, pattern = "lex"), "Sex", "Education"),
         sexrank = ifelse(str_detect(string = name, pattern = "alt"), "Men", "Women")) 

|>
  ggplot(mapping = aes(x = sexrank, y = value, colour = primary)) + 
  geom_boxplot()

bind_rows(di, di2, .id = "model") |>
  select(-name) |>
  mutate(model = factor(model, levels = 1:2, labels = c("Constant", "Differential"))) |>
  ggplot(mapping = aes(x = sexrank, y = value, colour = primary)) + 
  geom_boxplot() + 
  facet_wrap(~model) + 
  theme_ipsum() +
  theme(legend.position = "top") + 
  labs(x = "Which is the higher ranking sex category?",
       y = "Dissimilarity index",
       colour = "Which is the primary dimension?")


# Regression models for the differential coefficients

t2 <- dmcoef |>
  group_by(model, dimension, term) |>
  nest()

t2 <- t2 |>
  mutate(m1 = map(data, ~lm(theta ~ cnt + I(year - 1992), data = .)))

screenreg(l = t2$m1[t2$dimension == "Education"], omit = "cnt", digits = 3)
screenreg(l = t2$m1[t2$dimension != "Education"], omit = "cnt", digits = 3)


dmcoef |>
  filter(term == "Male", model == "Two-dimensional") |>
  summarise(k = median(theta))


  