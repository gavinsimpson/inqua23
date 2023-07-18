# Script to analyse the small water C and N isotopes using a copula GAM

# Load packages
pkgs <- c("here", "tidyr", "GJRM", "mgcv", "gratia", "dplyr", "readr",
    "ggplot2")
vapply(pkgs, library, logical(1L),
    character.only = TRUE, logical.return = TRUE)

# load data
sw <- readr::read_csv(here("data", "small-water-isotopes.csv"),
    col_types = "ddddddd")

# create a tidy format
sw_tidy <- sw |>
    pivot_longer(c(-depth, -year),
        names_to = "variable", values_to = "values")

# plot the data
ok_pal <- palette.colors(3)[-1] |> unname()

iso_lab <- function(x) {
    c("d15n" = "delta^{15}*N", "d15n" = "delta^{13}*C")
}

pm_lab <- "\u2030"

iso_plt <- sw_tidy |>
    filter(variable %in% c("d13c", "d15n")) |>
    ggplot(aes(x = year, y = values, group = variable, colour = variable)) +
    geom_point(show.legend = FALSE) +
    facet_grid(variable ~ ., scales = "free_y",
      labeller = as_labeller(iso_lab, default = label_parsed)) +
    scale_colour_discrete(type = ok_pal) +
    labs(x = NULL, y = pm_lab)
iso_plt

sw |>
    ggplot(aes(x = d13c, y = d15n)) +
    geom_point()

# individual models

mu_13c    <- d13c ~ s(year, k = 10)
mu_15n    <- d15n ~ s(year, k = 10)
sigma_13c <- ~ 1 # + s(year, k = 4)
sigma_15n <- ~ 1 # +  s(year, k = 4)
theta     <- ~ s(year, k = 10)
nu        <- ~ 1
fml <- list(mu_13c, mu_15n, sigma_13c, sigma_15n, theta)
fml_t <- c(fml, nu)

m_gauss_cop <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "N", gamlssfit = TRUE)

conv.check(m_gauss_cop)
summary(m_gauss_cop)


tmp <- sw |>
    mutate(theta = m_gauss_cop$theta[, 1], tau = m_gauss_cop$tau[, 1])

tmp |>
    ggplot(aes(x = year, y = tau)) +
    geom_line()

# different copulas
m_cop_N     <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "N", gamlssfit = TRUE)
m_cop_F     <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "F", gamlssfit = TRUE)
m_cop_AMH   <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "AMH", gamlssfit = TRUE)
m_cop_FGM   <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "FGM", gamlssfit = TRUE)
m_cop_C0    <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "C0", gamlssfit = TRUE)
m_cop_C270  <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "C270", gamlssfit = TRUE)
m_cop_G0    <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "G0", gamlssfit = TRUE)
m_cop_G90   <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "G90", gamlssfit = TRUE)
m_cop_G270  <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "G270", gamlssfit = TRUE)
m_cop_J0    <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "J0", gamlssfit = TRUE)
m_cop_J90   <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "J90", gamlssfit = TRUE)
m_cop_J270  <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "J270", gamlssfit = TRUE)
m_cop_PL    <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "PL", gamlssfit = TRUE)
m_cop_HO    <- gjrm(fml, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "HO", gamlssfit = TRUE)
m_cop_T     <- gjrm(fml_t, data = sw, margins = c("N", "N"), Model = "B",
                    BivD = "T", gamlssfit = TRUE)

# J270 has lowest AIC but it is highly uncertain, so ignore it
# C90 fails to converge
AIC(m_cop_N, m_cop_F, m_cop_AMH, m_cop_FGM, m_cop_PL, m_cop_HO, m_cop_T,
    m_cop_C0, m_cop_C270, # m_cop_C90,
    m_cop_G0, m_cop_G90, m_cop_G180, m_cop_G270,
    m_cop_J0, m_cop_J90 # , m_cop_J270)
    ) |>
    data.frame() |>
        tibble::rownames_to_column(var = "model") |>
        arrange(AIC)

tmp <- sw |>
    mutate(tau_N = m_cop_N$tau[, 1],
        tau_F = m_cop_F$tau[, 1],
        tau_PL = m_cop_PL$tau,
        tau_G90 = m_cop_G90$tau[, 1],
        tau_J270 = m_cop_J270$tau[, 1])

tmp |>
    ggplot(aes(x = year, y = tau_J270)) +
    geom_line()

# going to need some example code to illustrate how the estimated copulas look
# given the different estimates of Tau over time
# Note that copula::rotCopula() can produced rotated (reflected) versions of
# copula.

## Packages
library("copula")
library("ggplot2")
library("dplyr")
library("tibble")

## Load Data

## Example of different copulas with Gaussian marginals
tau <- -0.7
th.n <- iTau(normalCopula(), tau = tau)
th.t <- iTau(tCopula(df = 3), tau = tau)
th.c <- iTau(rotCopula(claytonCopula(), flip = c(TRUE, FALSE)), tau = tau)
th.g <- iTau(rotCopula(gumbelCopula(), flip = c(TRUE, FALSE)), tau = tau)

## sample from objects
set.seed(271)
n <- 2000
N01m <- list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
X.n <- rMvdc(n, mvdc = mvdc(normalCopula(th.n),    c("norm", "norm"), N01m))
X.t <- rMvdc(n, mvdc = mvdc(tCopula(th.t, df = 3), c("norm", "norm"), N01m))
X.c <- rMvdc(n, mvdc = mvdc(rotCopula(claytonCopula(th.c),
    flip = c(TRUE, FALSE)), c("norm", "norm"), N01m))
X.g <- rMvdc(n, mvdc = mvdc(rotCopula(gumbelCopula(th.g),
    flip = c(TRUE, FALSE)), c("norm", "norm"), N01m))

## put into something I can plot with ggplot
cop_names <- c("Normal", "t", "Clayton", "Gumbel")
cops <- tibble(x1 = c(X.n[, 1], X.t[, 1], X.c[, 1], X.g[, 1]),
               x2 = c(X.n[, 2], X.t[, 2], X.c[, 2], X.g[, 2]),
               copula = factor(rep(cop_names, each = n), levels = cop_names))

## colours
pal <- c(Normal = "#0077bb", t = "#009988", Clayton = "#cc3311",
    Gumbel = "#ee3377")
ggplot(cops, aes(x = x1, y = x2, colour = copula)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~ copula, nrow = 1) +
    scale_colour_manual(values = pal, guide = FALSE) +
    theme_bw(base_size = 14, base_family = 'Fira Sans') +
    theme(strip.text = element_text(size = rel(1), face = "bold", hjust = 0),
          plot.title = element_text(face = "bold")) +
    labs(title = "Comparing copulas",
        caption = bquote("Kendall's" ~ tau == .(tau))) # "Kendall's τ = -0.7")

# plot d15n and d13 c
sw |>
ggplot(aes(x = d13c, y = d15n)) +
    geom_point()
