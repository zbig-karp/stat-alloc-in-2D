library(stat.alloc)
library(tidyverse)
library(knitr)
library(reshape2)
library(kableExtra)
library(gghalves)
library(Hmisc)
library(broom)
library(texreg)

# TABLE 1 IN THE TEXT ----

# A tibble with the original data ----

t1 <- matrix(c(4051, 1044, 299, 38, 7832, 1341, 845, 239, 1536, 3222, 1088, 131, 3031, 1738, 3386, 578, 1898, 8640, 5051, 589, 3257, 2503, 11080, 2010, 356, 1606, 3160, 601, 725, 693, 5526, 1430, 98, 325, 1846, 457, 452, 252, 3885, 1402), ncol = 4, byrow = TRUE)

dimnames(t1) <- list(
  origin = paste(rep(paste0("E", 1:5), each = 2), 
                 rep(c("F", "M"), times = 5), sep = "-"),
  status = paste0("S", 1:4)
)

t1 <- t1 |>
  melt(value.name = "freq", varnames = c("origin", "status")) |>
  separate_wider_delim(cols = origin, names = c("degree", "sex"), delim = "-", 
                       cols_remove = FALSE)

# Calculating the reference allocations (assuming women to be ranked over men) ----

t1_ref <- refall2d(d = t1, dim1 = "degree", dim2 = "sex", status = "status", f = "freq")

# Calculating the index of dissimilarity between each of the reference allocations and the observed one

t1_diss <- map_dbl(.x = t1_ref[-1], .f = ~sum(abs(t1_ref[[1]] - .))/(2 * sum(.)))

# Calculating the (Frobenius) distance between each of the reference allocations and the observed one

t1_dist <- map_dbl(.x = t1_ref[-1], .f = ~norm(t1_ref[[1]] - ., type = "F"))

# Calculating the reference allocations (assuming men to be ranked over women) ----

t1_ref_alt <- t1 |> 
  mutate(sex2 = fct_relevel(sex, "M")) |> 
  refall2d(dim1 = "degree", dim2 = "sex2", status = "status", f = "freq")

# Calculating the index of dissimilarity between each of the reference allocations and the observed one

t1_diss_alt <- map_dbl(.x = t1_ref_alt[-1], .f = ~sum(abs(t1_ref_alt[[1]] - .))/(2 * sum(.)))

# Calculating the distance between each of the reference allocations and the observed one

t1_dist_alt <- map_dbl(.x = t1_ref_alt[-1], .f = ~norm(t1_ref_alt[[1]] - ., type = "F"))

# A matrix of observed counts (Table 1 in the text)

latex(addmargins(t1_ref[[1]]), file = "", title = "", cgroup = c("Status", ""), n.cgroup = c(4, 1), rgroup = c("Origin", ""), n.rgroup = c(10, 1), big.mark = ",", booktabs = TRUE, insert.bottom = "\\medskip\\footnotesize{\\emph{Destination statuses} S1: Professional, technical and kindred workers; Managers and administrators, except farm; S2: Sales workers; Clerical and kindred workers; S3: Craft and kindred workers; Operatives, except transport; Transport equipment operatives; Service workers, except household; S4: Laborers, except farm; Farm workers; Private household workers. \\emph{Educational levels} E1: College, 4 years or more; E2: College, 1-3 years; E3: High school, 4 years; E4: High school, 1-3 years; E5: Elementary school, 8 years or less. \\emph{Genders}  M: Male; F: Female.}", table.env = FALSE)

# TABLE 2 IN THE TEXT ----

# Matrices representing the four reference allocations

kable(t1_ref[[2]], booktabs = TRUE, linesep = "", format.args = list(big.mark = ","), escape = FALSE)
kable(t1_ref[[3]], booktabs = TRUE, linesep = "", format.args = list(big.mark = ","), escape = FALSE)
kable(t1_ref[[4]], booktabs = TRUE, linesep = "", format.args = list(big.mark = ","), escape = FALSE)
kable(t1_ref[[5]], booktabs = TRUE, linesep = "", format.args = list(big.mark = ","), escape = FALSE)

# TABLE 3 IN THE TEXT ----

# Single-dimensional models ----

# Education-based status-allocation ----

t1_cmm_edu <- xtabs(freq ~ degree + status, data = t1) |>
  refall() |>
  cmm_mde()

# Gender-based status-allocation ----

t1_cmm_sex <- xtabs(freq ~ sex + status, data = t1) |>
  refall() |>
  cmm_mde()

# Two-dimensional status allocation model ----

# Reference allocations

t1_stall2d <- t1 |>
  refall2d(dim1 = "degree", dim2 = "sex", status = "status", f = "freq")

# Fitting the model
t1_cmm2d <- cmm2d_mde(dat = t1_stall2d)

# Fit comparison between 1D and 2D models ---

a1 <- t1_cmm_edu[[1]]$Estimate
a2 <- t1_cmm_sex[[1]]$Estimate
mpsaed <- a1 * t1_ref[[3]] + (1 - a1) * t1_ref[[5]]
mpsasx <- a2 * t1_ref[[4]] + (1 - a2) * t1_ref[[5]]
gf1ded <- sum(abs(prop.table(mpsaed) - prop.table(t1_ref[[1]])))/2
gf1dsx <- sum(abs(prop.table(mpsasx) - prop.table(t1_ref[[1]])))/2

# Estimates of model parameters (Table 3 in the text) ----

tibble(
  Model = c("Education", "Gender", "Education only", "Gender only", "Education & gender"),
  Education = c(t1_cmm_edu[[1]]$Estimate, NA, t1_cmm_edu[[1]]$Estimate, NA, t1_cmm2d[[1]]$Estimate[[1]]),
  Gender = c(NA, t1_cmm_sex[[1]]$Estimate, NA, t1_cmm_sex[[1]]$Estimate, t1_cmm2d[[1]]$Estimate[[2]]),
  `Dissimilarity index` = c(t1_cmm_edu[[3]], t1_cmm_sex[[3]], gf1ded, gf1dsx, t1_cmm2d[[3]])
) |>
  mutate(across(.cols = -1, .fns = ~round(., 2))) |>
  latex(file = "", rgroup = c("Single-dimensional", "Two-dimensional"), n.rgroup = c(2, 3), 
        rowlabel = "", rowname = "", table.env = FALSE, 
        cgroup = c("", "Mixing coefficient", ""), n.cgroup = c(1, 2, 1), 
        booktabs = TRUE, dcolumn = TRUE, insert.bottom = "\\footnotesize{\\raggedright NB: The mixing coefficients in the ``Education only'' and ``Gender only'' versions of the two-dimensional model come from relevant single-dimensional models, with the mixing coefficient on the other merit dimension set to 0, as explained in the text.}")

# Results from the Eurostat data ----

# Preparing the data for analysis ----

# Reading in the data. The file contains the Eurostat table `lfsa-egised`, downloaded from the Eurostat's website on May 2, 2024 (go to https://ec.europa.eu/eurostat/web/products-datasets/-/lfsa_egised).

p <- read_delim(file = "estat_lfsa_egised.tsv", delim = "\t")

# Reshaping the data to the long format

p <- p |>
  pivot_longer(cols = -1, names_to = "year") |>
  separate_wider_delim(cols = 1, names = c("freq", "age", "sex", "isco08", "isced11", "unit",
                                           "cnt"), delim = ",") |>
  mutate(value = str_remove_all(string = value, pattern = "[a-z]|\\s|\\:"),
         across(.cols = c(year, value), .fns = parse_number))

# Selecting the cases of interest

p <- p |>
  filter(age != "Y20-64",
         sex != "T",
         isco08 %in% paste0("OC", 1:9),
         isced11 %in% c("ED0-2", "ED3_4", "ED5-8"),
         !str_detect(string = cnt, pattern = "E[AU]"))

# Some further recoding

p <- p |>
  mutate(degree = factor(isced11, levels = c("ED5-8", "ED3_4", "ED0-2"), 
                         labels = c("High", "Medium", "Low")),
         status = case_when(
           isco08 %in% paste0("OC", 1:3) ~ "S1",
           isco08 %in% paste0("OC", 4:5) ~ "S2",
           isco08 %in% paste0("OC", 6:7) ~ "S3",
           isco08 %in% paste0("OC", 8:9) ~ "S4"),
         status = factor(status))

# The final tibble ----

d <- p |>
  count(year, cnt, sex, degree, status, wt = value, name = "freq") |>
  mutate(origin = paste(str_sub(string = degree, start = 1, end = 1), sex, sep = "-")) |>
  group_by(year, cnt) |>
  nest()

# Removing 'empty' country-year combinations (e.g., the LFS survey was not carried out in AT in 1983). Counts (frequencies) in status allocation tables corresponding to the empty combination sum up to 0

d <- d |>
  mutate(test = map_dbl(data, ~sum(.$freq))) |>
  filter(test != 0) |>
  select(-test) |>
  ungroup()

# FIGURE 1 IN THE TEXT ----

# Gender differences in educational attainment

d <- d |>
  mutate(uedu = map(data, ~{
    . |>
      count(sex, degree, wt = freq) |>
      group_by(sex) |>
      mutate(p = prop.table(n),
             cp = cumsum(p[n():1])[n():1]) |>
      summarise(uedu = (n() - sum(cp))/(n() - 1)) |>
      summarise(dU.edu = first(uedu) - last(uedu)) 
  }))

# Gender differences in status allocation

d <- d |>
  mutate(uocc = map(data, ~{
    . |>
      count(sex, status, wt = freq) |>
      group_by(sex) |>
      mutate(p = prop.table(n),
             cp = cumsum(p[n():1])[n():1]) |>
      summarise(uocc = (n() - sum(cp))/(n() - 1)) |>
      summarise(dU.occ = first(uocc) - last(uocc)) 
  }))

# Status differences within educational categories

d <- d |>
  mutate(data1 = map(data, ~mutate(., across(.cols = c(degree, status), .fns = fct_rev))),
         data1 = map(data1, ~arrange(., degree, sex, status)),
         data1 = map(data1, ~{
           . |>
             group_by(degree, sex) |>
             mutate(p = prop.table(freq),
                    cp = cumsum(p)) |>
             summarise(U = (n() - sum(cp))/(n() - 1)) |>
             summarise(dU = first(U) - last(U))
         }))

fig.dU <- d |>
  unnest(cols = c(uocc, uedu)) |>
  select(cnt, year, starts_with("dU.")) |>
  pivot_longer(cols = starts_with("dU."),
               names_to = "dim", 
               names_transform = list(
                 dim = ~str_remove(string = ., pattern = "dU\\.")
               )) |>
  group_by(year, dim) |>
  summarise(
    m = mean(value),
    s = sd(value)/sqrt(n()),
    n = n(),
    .groups = "drop"
  ) |>
  mutate(lwr = m + qt(p = 0.025, df = n - 1) * s,
         upr = m + qt(p = 0.975, df = n - 1) * s)

fig.dU |>
  ggplot(mapping = aes(x = year, y = m, colour = dim, fill = dim)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(mapping = aes(ymin = lwr, ymax = upr, colour = dim), position = position_dodge(width = 1)) +
  geom_point(pch = 21, size = 2, alpha = 0.75, , position = position_dodge(width = 1)) + 
  theme_bw() + 
  theme(legend.position = "top") + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2023, by = 3)) + 
  scale_colour_manual(values = c("#0C7BDC", "#FFC20A"), labels = c("Education", "Status")) + 
  scale_fill_manual(values = c("#0C7BDC", "#FFC20A"), labels = c("Education", "Status")) + 
  labs(x = "Year", y = "Average $\\Delta U$", colour = "Dimension", fill = "Dimension")

# FIGURE 2 IN THE TEXT ----

d <- d |>
  mutate(tb1ded = map(data, ~xtabs(freq ~ degree + status, data = .)),
         tb1ded_prm = map(tb1ded, 
                          ~list(hml = .,
                                hlm = .[c(1, 3, 2), ],
                                mhl = .[c(2, 1, 3), ],
                                mlh = .[c(2, 3, 1), ],
                                lhm = .[c(3, 1, 2), ],
                                lmh = .[3:1, ])),
         sa1ded_prm = map(tb1ded_prm, ~map(., refall)),
         cm1ded_prm = map(sa1ded_prm, ~map(., cmm_mde)))

d <- d |>
  mutate(cm1ded_prm = map(cm1ded_prm, ~map_df(., ~.[[3]])))

d <- d |>
  mutate(tb1dsx = map(data, ~xtabs(freq ~ sex + status, data = .)),
         tb1dsx_prm = map(tb1dsx, 
                          ~list(
                            f = .,
                            m = .[2:1, ]
                          )),
         sa1dsx_prm = map(tb1dsx_prm, ~map(., refall)),
         cm1dsx_prm = map(sa1dsx_prm, ~map(., cmm_mde)))

d <- d |>
  mutate(cm1dsx_prm = map(cm1dsx_prm, ~map_df(., ~.[[3]])))

cm1ded_prm <- d |>
  select(year, cnt, cm1ded_prm) |>
  unnest(cols = cm1ded_prm)

cm1dsx_prm <- d |>
  select(year, cnt, cm1dsx_prm) |>
  unnest(cols = cm1dsx_prm)

ed_prm <- cm1ded_prm |>
  pivot_longer(cols = 3:8) 

sx_prm <- cm1dsx_prm |>
  pivot_longer(cols = 3:4)

ed_prm <- ed_prm |>
  mutate(name = factor(name, 
                       levels = c("hml", "hlm", "mhl", "mlh", "lhm", "lmh"),
                       labels = c("$H>M>L$", "$H>L>M$", "$M>H>L$", "$M>L>H$", "$L>H>M$", "$L>M>H$"))) 

sx_prm <- sx_prm |>
  mutate(name = factor(name, 
                       levels = c("f", "m"),
                       labels = c("$F > M$", "$M > F$"))) 

bind_rows(ed_prm, sx_prm, .id = "dim") |>
  mutate(dim = factor(dim, levels = 1:2, labels = c("Education", "Gender"))) |>
  ggplot(mapping = aes(x = name, y = value)) + 
  geom_half_boxplot(fill = "#bdbdbd", alpha = 0.5) + 
  geom_half_violin(side = "r", fill = "#bdbdbd", alpha = 0.5) + 
  facet_grid(~dim, scales = "free_x", space = "free_x") + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.text = element_text(face = "bold", size = 12)) +
  labs(x = NULL, y = "Dissimilarity index")

# FIGURE 3 IN THE TEXT ----

d <- d |>
  mutate(sa2d = map(data, ~refall2d(d = ., f = "freq", dim1 = "degree", dim2 = "sex", status = "status")),
         dist = map(.x = sa2d, .f = ~map2(.x = .[1], 
                                          .y = .[-1], 
                                          .f = ~norm(.x - .y, type = "f")/sum(.x))))

d <- d |>
  mutate(dist = map(.x = dist, .f = unlist),
         dist = map_df(.x = dist, 
                       .f = ~tibble(AM = .[1], AU = .[2], AV = .[3], AL = .[4]))) |>
  unnest(cols = dist)

d |>
  select(year, cnt, AM, AU, AV, AL) |>
  pivot_longer(starts_with("A")) |>
  mutate(name = factor(name, levels = c("AM", "AU", "AV", "AL"))) |>
  ggplot(mapping = aes(x = name, y = value)) +
  geom_half_boxplot(fill = "#bdbdbd", alpha = 0.5) +
  geom_half_violin(side = "r", fill = "#bdbdbd", alpha = 0.5) +
  theme_minimal() +
  labs(x = NULL, y = "Distance") +
  scale_x_discrete(labels = c("Observed to\n$\\langle M_xM_y\\rangle$",
                              "Observed to\n$\\langle M_xL_y\\rangle$",
                              "Observed to\n$\\langle L_xM_y\\rangle$",
                              "Observed to\n$\\langle L_xL_y\\rangle$"))

# FIGURE 4 IN THE TEXT ----

d <- d |>
  mutate(sa1ded = map(tb1ded, refall),
         cm1ded = map(sa1ded, cmm_mde),
         sa1dsx = map(tb1dsx, refall),
         cm1dsx = map(sa1dsx, cmm_mde),
         tb1dsx_alt = map(tb1dsx, ~.[2:1, ]),
         sa1dsx_alt = map(tb1dsx_alt, refall),
         cm1dsx_alt = map(sa1dsx_alt, cmm_mde),
         sa2d = map(data, ~refall2d(d = ., f = "freq", dim1 = "degree", dim2 = "sex", status = "status")),
         cm2d = map(sa2d, cmm2d_mde),
         data = map(data, ~mutate(.data = ., sex1 = fct_relevel(sex, "M", "F"))),
         sa2d_alt = map(data, ~refall2d(d = ., f = "freq", dim1 = "degree", dim2 = "sex1", status = "status")),
         cm2d_alt = map(sa2d_alt, cmm2d_mde),
         alfa1_edu = map_dbl(cm1ded, ~.$`Mixing coefficient`$Estimate),
         alfa2_edu = map_dbl(cm2d, ~.$`Mixing coefficient`$Estimate[1]),
         alfa1_sex = map_dbl(cm1dsx, ~.$`Mixing coefficient`$Estimate),
         alfa2_sex = map_dbl(cm2d, ~.$`Mixing coefficient`$Estimate[2]))

# T-test of the difference between estimates of the mixing coefficients for gender in the 2D and 1D model

t1_sex_diff <- with(d, t.test(x = alfa2_sex, y = alfa1_sex, var.equal = TRUE)) |> tidy()

# Boxplot of the estimates mixing coefficients
d |>
  pivot_longer(cols = starts_with("alfa"), names_to = c("model", "dimension"),
               names_sep = "_",
               names_transform = list(model = ~paste0(parse_number(.), "D"))) |>
  ggplot(mapping = aes(x = dimension, y = value, colour = model, fill = model)) + 
  geom_half_boxplot(alpha = 0.5) +
  geom_half_violin(alpha = 0.5, side = "r", show.legend = FALSE) +
  scale_colour_manual(values = c("#0C7BDC", "#FFC20A"), 
                      labels = c("Single-dimensional", "Two-dimensional")) + 
  scale_fill_manual(values = c("#0C7BDC", "#FFC20A"), 
                    labels = c("Single-dimensional", "Two-dimensional")) + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  labs(x = NULL, y = "Estimated mixing coefficient", colour = "Model", fill = "Model") + 
  scale_x_discrete(labels = c("Education", "Gender"))

# FIGURE 5 IN THE TEXT ----

d <- d |>
  mutate(gf1ded = map2(.x = alfa1_edu, .y = sa2d, .f = ~.x * .y[[3]] + (1 - .x) * .y[[5]]),
         gf1ded = map2_dbl(.x = gf1ded, .y = sa2d, .f = ~sum(abs(.x - .y[[1]]))/(2 * sum(.x))),
         gf1dsx = map2(.x = alfa1_sex, .y = sa2d, .f = ~.x * .y[[4]] + (1 - .x) * .y[[5]]),
         gf1dsx = map2_dbl(.x = gf1dsx, .y = sa2d, .f = ~sum(abs(.x - .y[[1]]))/(2 * sum(.x))),
         gf2d = map_dbl(.x = cm2d, .f = ~.[[3]]))

t1_fit_diff <- with(d, t.test(x = gf1ded, y = gf2d, paired = TRUE, var.equal = TRUE)) |>
  glance()

d |>
  select(cnt, year, starts_with("gf")) |>
  pivot_longer(cols = starts_with("gf")) |>
  ggplot(mapping = aes(x = name, y = value)) + 
  geom_half_boxplot(fill = "#bdbdbd", alpha = 0.5) + 
  geom_half_violin(side = "r", fill = "#bdbdbd", alpha = 0.5) +
  theme_minimal() + 
  scale_x_discrete(labels = c("Single-dimensional\nEducation", "Single-dimensional\nGender",
                              "Two-dimensional")) + 
  labs(x = NULL, y = "Dissimilarity index")

# TABLE 4 IN THE TEXT ----

d <- d |>
  mutate(cnt = fct_relevel(cnt, "PL"))

list(
  m1 = lm(alfa1_edu ~ cnt + I(year - 1992), data = d),
  m2 = lm(alfa2_edu ~ cnt + I(year - 1992), data = d),
  m3 = lm(alfa1_sex ~ cnt + I(year - 1992), data = d),
  m4 = lm(alfa2_sex ~ cnt + I(year - 1992), data = d)
) |>
  texreg(omit = "cnt", digits = 3, 
         custom.coef.names = c("Intercept", "Year (ref. 1992)"),
         booktabs = TRUE, dcolumn = TRUE, use.packages = FALSE, table = FALSE,
         custom.model.names = c("One dim", "Two dims", "One dim", "Two dims"),
         custom.header = list("Education" = 1:2, "Gender" = 3:4),
         custom.gof.rows = list("Country fixed effects" = rep("Yes", 4))
  )

# FIGURE 6 IN THE TEXT ----

d |>
  select(cnt, year, data1, uocc) |>
  unnest(cols = c(data1, uocc)) |>
  mutate(degree = fct_rev(degree)) |>
  ggplot(mapping = aes(x = degree, y = dU)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_half_boxplot(fill = "#bdbdbd", alpha = 0.5) +
  geom_half_violin(side = "r", fill = "#bdbdbd", alpha = 0.5) +
  theme_bw() + 
  labs(x = "Educational attainment", y = "Status advantage of women over men")