## Fit a Topic Model to the Abernethy Forest data set

# install.packages("remotes")
# remotes::install_github("gavinsimpson/ggvegan")

# install.packages("servr")

## Packages
library("here")
library("topicmodels")
#library("stm")
library("readr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("patchwork")
# library("cowplot")
# theme_set(theme_bw())
library("analogue")
#library("LDAvis")
#library("ggvegan")
# library("DirichletReg") # using brms for the dirichlet
library("brms")
# library("splines") # not needed for brms

## Load data
aber <- read_rds(here("data", "abernethy-count-data.rds"))
## or:
# aber <- read_rds("https://bit.ly/abercount")

## take a subset of spp
take <- c("BETULA", "PINUS_SYLVESTRIS", "ULMUS", "QUERCUS", "ALNUS_GLUTINOSA",
          "CORYLUS_MYRICA", "SALIX", "JUNIPERUS_COMMUNIS", "CALLUNA_VULGARIS",
          "EMPETRUM", "GRAMINEAE", "CYPERACEAE", "SOLIDAGO_TYPE",
          "COMPOSITAE_LIGULIFLORAE", "ARTEMISIA",
          "CARYOPHYLLACEAE_UNDIFFERENTIATED",
          "SAGINA", "SILENE_CF_S_ACAULIS", "CHENOPODIACEAE", "EPILOBIUM_TYPE",
          "PAPILIONACEAE_UNDIFFERENTIATED", "ANTHYLLIS_VULNERARIA",
          "ASTRAGALUS_ALPINUS", "ONONIS_TYPE", "ROSACEAE_UNDIFFERENTIATED",
          "RUBIACEAE", "RANUNCULACEAE_UNDIFFERENTIATED", "THALICTRUM",
          "RUMEX_ACETOSA_TYPE", "OXYRIA_TYPE", "PARNASSIA_PALUSTRIS",
          "SAXIFRAGA_GRANULATA", "SAXIFRAGA_HIRCULUS_TYPE", "SAXIFRAGA_NIVALIS",
          "SAXIFRAGA_OPPOSITIFOLIA", "SAXIFRAGA_UNDIFFERENTIATED", "SEDUM",
          "URTICA", "VERONICA")
## Don't do this!
##take <- c(1,2,3,4,6,10,11,12,14,15,39,40,41,42,43,46,49,50,53,54,57,58,59,60,67,
##          69,70,72,74,75,83,85,86)
aber <- aber[, take]
## are any columns now all zeroes?
all_missing <- unname(vapply(aber, function(x) all(is.na(x)), logical(1)))
## drop those with all NA
aber <- aber[, !all_missing]
## change all the NA to 0s
aber <- tran(aber, method = "missing")
## check that all remaining values are numeric
stopifnot(all(vapply(aber, data.class, character(1), USE.NAMES = FALSE) == "numeric"))
## check all columns still have at least 1 positive count
cs <- colSums(aber) > 0L
aber <- aber[, cs]
names(aber) <- tolower(names(aber))
## aber ages
aber_age <- read_rds(here("data", "abernethy-sample-age.rds"))
## or:
# aber_age <- read_rds("https://bit.ly/aberage")

## Models to fit
k <- 2:10 # 2, 3, ... 10 associations / groups
## setting the same random seed for each model
reps <- length(k)
ctrl <- replicate(reps, list(seed = 42), simplify = FALSE)
## repeat the data n times to facilitate `mapply`
aberrep <- replicate(reps, aber, simplify = FALSE)
# fit the sequence of topic models
tms <- mapply(LDA, k = k, x = aberrep, control = ctrl)

## extract model fit in terms of BIC and plot
plot(k, sapply(tms, AIC, k = log(nrow(aber))))

## so 5 groups looks OK
k_take <- 6
## which is the 5 group model?
k_ind <- which(k == k_take)
## could also selected purely on BIC terms...
k_bic <- which.min(sapply(tms, AIC, k = log(nrow(aber))))
## but we'll take the model with 5 groups
aber_lda <- tms[[k_ind]]

## extract the posterior fitted distribution of the model
aber_posterior <- posterior(aber_lda)
## topic proportions
aber_topics <- aber_posterior$topics
## term proportions
aber_terms <- aber_posterior$terms

# Prep the data for modelling
# aber_topics needs to be a data frame with suitable names, currently a matrix
# with invalid names
topic_df <- tibble::as_tibble(aber_topics) |>
    setNames(paste0("assoc", seq_len(k_take))) |>
    bind_cols(aber_age) |>
    rename(age = Age) |>
    mutate(neg_age = -age,
        t = 1 - ((age - min(age)) / (max(age) - min(age))))

topic_df |>
    pivot_longer(cols = c(-age, -neg_age, -t),
        names_to = "association",
        names_pattern = "assoc([[:digit:]]{1})",
        values_to = "proportion") |>
    ggplot(aes(x = age, y = proportion, colour = association)) +
        geom_point() +
        scale_x_reverse()

# Fit the brms dirichlet model
bind <- function(...) cbind(...)
bind <- cbind

fml1 <- brmsformula(bind(assoc1, assoc2, assoc3, assoc4, assoc5, assoc6) ~
    s(t, k = 10))
fml2 <- brmsformula(bind(assoc1, assoc2, assoc3, assoc4, assoc5, assoc6) ~
    s(t, k = 8), phi ~ s(t, k = 15))

m <- brm(fml1, data = topic_df, family = dirichlet(refcat = "assoc3"),
    chains = 4, cores = 4, control = list(adapt_delta = 0.975), iter = 5000,
    backend = "cmdstanr", seed = 15)

write_rds(m, "~/aber-dirichlet.rds")

m2 <- brm(fml2, data = topic_df, family = dirichlet(refcat = "assoc3"),
    chains = 4, cores = 4, control = list(adapt_delta = 0.99), iter = 5000,
    backend = "cmdstanr", seed = 15)

write_rds(m2, "~/aber-dirichlet-phi.rds")

# m <- brm(bind(assoc1, assoc2, assoc3, assoc4, assoc5) ~ s(t, k = 10),
#     data = topic_df, family = dirichlet(),
#     chains = 4, cores = 4, control = list(adapt_delta = 0.9))

# load the model (saved - compiling takes a lot of RAM at moment)
m  <- read_rds(here("models/aber-dirichlet.rds"))
m2 <- read_rds(here("models/aber-dirichlet-phi.rds"))

cs  <- conditional_smooths(m)
cs2 <- conditional_smooths(m2)

plot(cs)

new_df <- with(topic_df, tibble(age = seq(min(age), max(age), length = 100))) |>
    mutate(neg_age = -age,
        t = 1 - ((age - min(age)) / (max(age) - min(age))))

# fitted values
fv1 <- fitted(m, newdata = new_df, scale = "response")
fv2 <- fitted(m2, newdata = new_df, scale = "response")

fv1_samps <- fitted(m, newdata = new_df, scale = "response", summary = FALSE)
fv2_samps <- fitted(m2, newdata = new_df, scale = "response", summary = FALSE)

write_rds(fv1, "~/aber-dirichlet-fitted-summary.rds")
write_rds(fv2, "~/aber-dirichlet-phi-fitted-summary.rds")
write_rds(fv1_samps, "~/aber-dirichlet-fitted-samples.rds")
write_rds(fv2_samps, "~/aber-dirichlet-phi-fitted-samples.rds")

fv1 <- read_rds(here("models/aber-dirichlet-fitted-summary.rds"))
fv2 <- read_rds(here("models/aber-dirichlet-phi-fitted-summary.rds"))
fv1_samps <- read_rds(here("models/aber-dirichlet-fitted-samples.rds"))
fv2_samps <- read_rds(here("models/aber-dirichlet-phi-fitted-samples.rds"))

# need to extract and wrangle the fvs into the required format
`extract_fv` <- function(x, data, k) {
    prob <- array_to_long_tibble(x) # posterior mean (default for fitted)
    lwr  <- array_to_long_tibble(x, var = "lower") # 2.5th quantile
    upr  <- array_to_long_tibble(x, var = "upper") # 97.5th quantile
    out <- bind_cols(prob, lwr[, "lower"], upr[, "upper"])

    # add on the data var
    out <- tibble::add_column(out,
        age = rep(data$age, each = k),
        neg_age = rep(data$neg_age, each = k),
        t = rep(data$t, each = k), .before = 1L)
    out
}

`array_to_long_tibble` <- function(x, var = c("proportion", "lower", "upper")) {
    var <- match.arg(var)
    take <- case_when(
        var == "proportion" ~ 1,
        var == "lower" ~ 3,
        var == "upper" ~ 4)
    out <- as_tibble(x[, take, ]) |>
        pivot_longer(cols = everything(), # pivot
            names_pattern = "assoc([[:digit:]]{1})",
            names_to = "association",
            values_to = var)
    out
}

fvs <- extract_fv(fv2, data = new_df, k = 6)

topic_df_long <- topic_df |>
    pivot_longer(cols = c(-age, -neg_age, -t),
        names_to = "association",
        names_pattern = "assoc([[:digit:]]{1})",
        values_to = "proportion")

fvs |>
    ggplot(aes(x = age, y = proportion, colour = association)) +
    geom_ribbon(aes(x = age, ymin = lower, ymax = upper,
        fill = association), inherit.aes = FALSE, alpha = 0.2) +
    geom_line() +
    geom_point(data = topic_df_long,
        mapping = aes(x = age, y = proportion, colour = association)) +
    scale_x_reverse() +
    labs(x = "Age (years BP)", y = "Relative abundance")

# derivatives
eps <- 1e-7

deriv_df_bkw <- with(topic_df, tibble(age = seq(min(age), max(age), length = 100))) |>
    mutate(neg_age = -age,
        t = 1 - ((age - min(age)) / (max(age) - min(age))),
        t = t - (eps / 2))

deriv_df_fwd <- with(topic_df, tibble(age = seq(min(age), max(age), length = 100))) |>
    mutate(neg_age = -age,
        t = 1 - ((age - min(age)) / (max(age) - min(age))),
        t = t + (eps / 2))

deriv_df <- deriv_df_fwd |>
    bind_rows(deriv_df_bkw)

fd1_samps <- fitted(m, newdata = deriv_df, scale = "response", summary = FALSE)
fd2_samps <- fitted(m2, newdata = deriv_df, scale = "response", summary = FALSE)

write_rds(fd1_samps, "~/aber-dirichlet-fitted-samples-derivs.rds")
write_rds(fd2_samps, "~/aber-dirichlet-phi-fitted-samples-derivs.rds")

fd1_samps <- read_rds(here("models/aber-dirichlet-fitted-samples-derivs.rds"))
fd2_samps <- read_rds(here("models/aber-dirichlet-phi-fitted-samples-derivs.rds"))

#fd1_fwd <- fd1_samps[, 1:100, 1]
#fd1_bkw <- fd1_samps[, 101:200, 1]
#fd1_dif <- (fd1_fwd - fd1_bkw) / eps

# work on the columns
#fd1_med <- apply(fd1_dif, 2, quantile, probs = 0.5)
#fd1_lwr <- apply(fd1_dif, 2, quantile, probs = 0.025)
#fd1_upr <- apply(fd1_dif, 2, quantile, probs = 0.975)

# deriv_data <- new_df |>
#     tibble::add_column(fd = fd1_med, lwr = fd1_lwr, upr = fd1_upr)

# arr1 <- array(1:27, dim = rep(3, 3))
# arr2 <- array(27:1, dim = rep(3, 3))


new_df <- with(topic_df, tibble(age = seq(min(age), max(age), length = 200))) |>
    mutate(neg_age = -age,
        t = 1 - ((age - min(age)) / (max(age) - min(age))))
fd_dif <- (fd2_samps[, 1:200, ] - fd2_samps[, 201:400, ]) / eps

# work on the columns
fd_med <- apply(fd_dif, 2:3, quantile, probs = 0.5) |>
    as_tibble() |>
        pivot_longer(cols = everything(), # pivot
            names_pattern = "assoc([[:digit:]]{1})",
            names_to = "association",
            values_to = "derivative")

fd_lwr <- apply(fd_dif, 2:3, quantile, probs = 0.025) |>
    as_tibble() |>
    pivot_longer(cols = everything(), # pivot
        names_pattern = "assoc([[:digit:]]{1})",
        names_to = "association",
        values_to = "lower")

fd_upr <- apply(fd_dif, 2:3, quantile, probs = 0.975) |>
    as_tibble() |>
        pivot_longer(cols = everything(), # pivot
            names_pattern = "assoc([[:digit:]]{1})",
            names_to = "association",
            values_to = "upper")

fd_df <- fd_med |>
    bind_cols(fd_lwr[, "lower"]) |>
    bind_cols(fd_upr[, "upper"]) |>
    tibble::add_column(age = rep(new_df$age, each = k_take),
        neg_age = rep(new_df$neg_age, each = k_take),
        t = rep(new_df$t, each = k_take), .before = 1L)

fd_df |>
    ggplot(aes(x = age, y = derivative, colour = association)) +
    geom_line() +
    scale_x_reverse()

fd_df |>
    ggplot(aes(x = age, y = abs(derivative), colour = association)) +
    geom_line() +
    scale_x_reverse()

fd_df |>
    group_by(age) |>
    summarise(rate_of_change = sum(abs(derivative))) |>
    ggplot(aes(x = age, y = rate_of_change)) +
    geom_line() +
    scale_x_reverse() +
    labs(x = "Age (years BP)", y = "Rate of compositional change")

library("RRatepol")

aber2 <- cbind.data.frame(sample_id = as.character(seq_len(nrow(aber))), aber)
rownames(aber2) <- NULL
aber_age2 <- aber_age |> rename(age = Age) |>
    tibble::add_column(aber_age,
        sample_id = as.character(seq_len(nrow(aber_age))), .before = 1L) |>
    as.data.frame()
rownames(aber_age2) <- NULL

aber_rp <- estimate_roc(
    data_source_community = aber2,
    data_source_age = aber_age2,
    smooth_method = "age.w",
    dissimilarity_coefficient = "chord",
    working_units = "levels",
    bin_size = 200
    )

# plot_roc(aber_rp)

aber_rp |>
    ggplot(aes(y = ROC, x = Age)) +
    geom_line(linewidth = 1.5) +
    scale_x_reverse()
