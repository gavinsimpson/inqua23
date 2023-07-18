# Analyse the d18O record from the Gulf of Taranto, Ionian Sea of Taricco et al
# (2016)
#
# Taricco, C., Alessio, S., Rubinetti, S., Vivaldo, G., Mancuso, S., 2016. A
# foraminiferal δ18O record covering the last 2,200 years. Scientific Data 3,
# 160042. https://doi.org/10.1038/sdata.2016.42

# load packages
library("pangaear")
library("dplyr")
library("ggplot2")
library("mgcv")
library("gratia")

# load the data from Pangaea
res <- pg_data(doi = "10.1594/PANGAEA.857573")
foram <- res[[1]]$data |>
    setNames(c("Depth", "Age_AD", "Age_kaBP", "d18O"))

max_agebp <- with(foram, max(-Age_kaBP))
min_agebp <- with(foram, min(-Age_kaBP))

foram <- foram |>
    mutate(t = (-Age_kaBP - min_agebp) / (max_agebp - min_agebp),
    neg_age = -Age_kaBP)

# some plotting constants
ylabel <- expression(delta^{18} * O ~ "[‰ VPDB]")
xlabel <- "Age [ka BP]"

# fit the model
m <- gam(d18O ~ s(neg_age, k = 100, bs = "ad"), data = foram, method = "REML")

# plot the data
ggplot(foram, aes(y = d18O, x = Age_kaBP)) +
    geom_path() +
    scale_x_reverse(sec.axis = sec_axis( ~ 1950 - (. * 1000),
        name = "Year [CE]")) +
    labs(y = ylabel, x = xlabel)

# fit the model
m <- gam(d18O ~ s(Age, k = 100, bs = "ad"), data = foram, method = "REML")

ds <- data_slice(m, neg_age = evenly(neg_age, 300))
fv <- fitted_values(m, data = ds) |>
    mutate(Age_kaBP = -neg_age)


fv |>
    ggplot(aes(y = fitted, x = Age_kaBP)) +
        geom_point(data = foram, aes(y = d18O, x = Age_kaBP), alpha = 0.5) +
        geom_ribbon(aes(x = Age_kaBP, ymin = lower, ymax = upper),
            alpha = 0.2) +
        geom_line(aes(x = Age_kaBP, y = fitted),
                  linewidth = 1) +
        scale_x_reverse(sec.axis = sec_axis( ~ 1950 - (. * 1000),
            name = "Age [AD]")) +
        #scale_y_reverse() +
        labs(y = ylabel, x = xlabel)

# compute first derivative of the smooth
fd <- derivatives(m, data = ds, type = "central")
fd |>
    mutate(data = -data) |>
    draw(add_change = TRUE, change_type = "sizer", lwd_change = 3) &
    labs(x = xlabel) &
    scale_x_reverse()