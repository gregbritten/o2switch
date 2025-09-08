rm(list=ls())

getwd()
setwd("~/Library/CloudStorage/GoogleDrive-julia.a.huggins@gmail.com/My Drive/WHOI/Modeling/R Model")

library(dplyr)
library(ggplot2)
library(scico)
library(tidyr)
library(ggpubr)
library(grid)
library(gtable)
library(gridExtra)
library(ggplotify)
library(cowplot)

# calculate pathway-specific ATP yield

# define set values

prod = 0.3           # ATPs produced (per proton pumped)
cost = 0             # ATPs cost of enzyme (per unit length)
Farad = 96.485       # Faradays constant (Kc/mol)
Rgas = 0.008314      # Gas constant (KJ/(mol*K))
Temp = 10            # Ambient Temp (deg. C)
pH = 8.3             # Ambient pH of solution (-log[H+])

# define concentrations (mol/L)
concentrations <- data.frame(
  row.names = c("O2", "NO3", "NO2", "N2O", "N2", "NH4", "CH2O", "HCO3", "H2O", "H"),
  conc.min = c(
    0        * 10^-6,         # O2
    10       * 10^-6,         # NO3
    10       * 10^-6,         # NO2
    10       * 10^-6,         # N2O
    5        * 10^-4,         # N2
    10       * 10^-6,         # NH4
    100      * 10^-6,         # CH2O
    2200     * 10^-6,         # HCO3
    1,                        # H2O
    10^-pH                    # H
  ),
  conc.max = c(
    5        * 10^-6,         # O2
    10       * 10^-6,         # NO3
    10       * 10^-6,         # NO2
    10       * 10^-6,         # N2O
    5        * 10^-4,         # N2
    10       * 10^-6,         # NH4
    100      * 10^-6,         # CH2O
    2200     * 10^-6,         # HCO3
    1,                        # H2O
    10^(-pH)                  # H
  )
)

# calculate the range of concentrations
concentrations$conc.range <- mapply(
  function(min, max) {
    if (min == max) {
      rep(min, 501)
    } else {
      seq(min, max, by = 0.01 * 10^-6)
    }
  },
  concentrations$conc.min,
  concentrations$conc.max,
  SIMPLIFY = FALSE
)

# make parameters table with set (known) values
parameters.set <- data.frame(
  row.names = c(
    "O2_LA",      # aerobic resp, low affinity
    "O2_HA",      # aerobic resp, high affinity
    "NO3_NO2",    # nitrate reduction to nitrite
    "NO2_N2O",    # nitrite reduction to N2O
    "N2O_N2",     # N2O reduction to N2
    "NO2_N2",     # denitrification from nitrite
    "NO3_N2",     # denitrification from nitrate
    "NO2_NH4",    # DNRA from nitrite
    "NO3_NH4"     # DNRA from nitrate
  ),
  num.e = c(      # number of electrons transferred per TEA
    4,            # O2_LA
    4,            # O2_HA
    2,            # NO3_NO2
    2,            # NO2_N2O
    1,            # N2O_N2
    3,            # NO2_N2
    5,            # NO3_N2
    6,            # NO2_NH4
    8             # NO3_NH4
  ),
  coeff.TEA = c(  # stoichiometric coefficient on TEA
    1,            # O2_LA
    1,            # O2_HA
    2,            # NO3_NO2
    2,            # NO2_N2O
    2,            # N2O_N2
    4,            # NO2_N2
    4,            # NO3_N2
    2,            # NO2_NH4
    1             # NO3_NH4
  ),
  coeff.TEAred = c(  # stoichiometric coefficient on product from TEA reduction
    2,            # O2_LA
    2,            # O2_HA
    2,            # NO3_NO2
    1,            # NO2_N2O
    2,            # N2O_N2
    2,            # NO2_N2
    2,            # NO3_N2
    2,            # NO2_NH4
    1             # NO3_NH4
  ),
  coeff.PED = c(  # stoichiometric coefficient on PED
    1,            # O2_LA
    1,            # O2_HA
    1,            # NO3_NO2
    1,            # NO2_N2O
    1,            # N2O_N2
    3,            # NO2_N2
    5,            # NO3_N2
    3,            # NO2_NH4
    2             # NO3_NH4
  ),
  coeff.PEDox = c(  # stoichiometric coefficient on product from PED oxidation
    1,            # O2_LA
    1,            # O2_HA
    1,            # NO3_NO2
    1,            # NO2_N2O
    1,            # N2O_N2
    3,            # NO2_N2
    5,            # NO3_N2
    3,            # NO2_NH4
    2             # NO3_NH4
  ),
  coeff.H.r = c(  # stoichiometric coefficient on H+ in the reactants
    0,            # O2_LA
    0,            # O2_HA
    0,            # NO3_NO2
    1,            # NO2_N2O
    0,            # N2O_N2
    1,            # NO2_N2
    0,            # NO3_N2
    1,            # NO2_NH4
    0             # NO3_NH4
  ),
  coeff.H.p = c(  # stoichiometric coefficient on H+ in the products
    1,            # O2_LA
    1,            # O2_HA
    1,            # NO3_NO2
    0,            # NO2_N2O
    1,            # N2O_N2
    0,            # NO2_N2
    1,            # NO3_N2
    0,            # NO2_NH4
    0             # NO3_NH4
  ),
  coeff.H2O.r = c(  # stoichiometric coefficient on H2O in the reactants
    2,            # O2_LA
    2,            # O2_HA
    0,            # NO3_NO2
    0,            # NO2_N2O
    0,            # N2O_N2
    0,            # NO2_N2
    0,            # NO3_N2
    2,            # NO2_NH4
    1             # NO3_NH4
  ),
  coeff.H2O.p = c(  # stoichiometric coefficient on H2O in the products
    0,            # O2_LA
    0,            # O2_HA
    0,            # NO3_NO2
    1,            # NO2_N2O
    0,            # N2O_N2
    2,            # NO2_N2
    2,            # NO3_N2
    0,            # NO2_NH4
    0             # NO3_NH4
  ),
  Ev = c(         # Electrode potential of reaction (standard state)
    1.17,         # O2_LA
    1.17,         # O2_HA
    0.74,         # NO3_NO2
    1.28,         # NO2_N2O
    1.66,         # N2O_N2
    1.41,         # NO2_N2
    1.14,         # NO3_N2
    0.79,         # NO2_NH4
    0.78          # NO3_NH4
  )
)

# export table of set parameters

parameters.set

# make parameters table with variable (unknown) values
parameters <- data.frame(
  row.names = c(
    "O2_LA",      # aerobic resp, low affinity
    "O2_HA",      # aerobic resp, high affinity
    "NO3_NO2",    # nitrate reduction to nitrite
    "NO2_N2O",    # nitrite reduction to N2O
    "N2O_N2",     # N2O reduction to N2
    "NO2_N2",     # denitrification from nitrite
    "NO3_N2",     # denitrification from nitrate
    "NO2_NH4",    # DNRA from nitrite
    "NO3_NH4"     # DNRA from nitrate
  ),
  num.e = c(      # number of electrons transferred per TEA
    4,            # O2_LA
    4,            # O2_HA
    2,            # NO3_NO2
    2,            # NO2_N2O
    1,            # N2O_N2
    3,            # NO2_N2
    5,            # NO3_N2
    6,            # NO2_NH4
    8             # NO3_NH4
  ),
  Ev = c(         # Electrode potential of reaction (standard state)
    1.17,         # O2_LA
    1.17,         # O2_HA
    0.74,         # NO3_NO2
    1.28,         # NO2_N2O
    1.66,         # N2O_N2
    1.41,         # NO2_N2
    1.14,         # NO3_N2
    0.79,         # NO2_NH4
    0.78          # NO3_NH4
  ),
  H.e = c(        # proton pumping yield per electron
    1,            # O2_LA
    1,            # O2_HA
    1,            # NO3_NO2
    1,            # NO2_N2O
    1,            # N2O_N2
    1,            # NO2_N2
    1,            # NO3_N2
    1,            # NO2_NH4
    1             # NO3_NH4
  ),
  Km = c(         # half-saturation constants (affinity)
    200  * 10^-9, # O2_LA
    10    * 10^-9, # O2_HA
    5    * 10^-6, # NO3_NO2
    5    * 10^-6, # NO2_N2O
    5    * 10^-6, # N2O_N2
    5    * 10^-6, # NO2_N2
    5    * 10^-6, # NO3_N2
    5    * 10^-6, # NO2_NH4
    5    * 10^-6  # NO3_NH4  
  ),
  Vmax = c(       # maximum reaction rate
    10  * 10^-9, # O2_LA
    7.5   * 10^-9, # O2_HA
    10  * 10^-9, # NO3_NO2      # rough guess that reaction rate = 10/#steps
    5  * 10^-9, # NO2_N2O
    10  * 10^-9, # N2O_N2
    3.3  * 10^-9, # NO2_N2
    2.5  * 10^-9, # NO3_N2
    10  * 10^-9, # NO2_NH4
    5  * 10^-9  # NO3_NH4  
  ),
  length = c(     # length of enzymes
    1,            # O2_LA
    1,            # O2_HA
    1,            # NO3_NO2
    1,            # NO2_N2O
    1,            # N2O_N2
    1,            # NO2_N2
    1,            # NO3_N2
    1,            # NO2_NH4
    1             # NO3_NH4  
  )
)

# join set parameters and variable parameters

parameters <- parameters %>%
  mutate(Reaction = rownames(.))   # Turn row names into a column
parameters.set <- parameters.set %>%
  mutate(Reaction = rownames(.))   # Turn row names into a column

parameters.full <- as.data.frame(
  full_join(parameters.set, parameters, "Reaction")
)

### parameters.full <- parameters.full %>%
###  row.names = parameters.full$Reaction
  

# calculate ATP generation from all reactions

PATP <- function(TEA,reaction){
  TEA_conc <- concentrations [[TEA, "conc.range"]]
  p <- as.list(parameters[reaction, ])
  y    <- p$num.e * p$H.e 
  c    <- cost * p$length
  Km   <- p$Km
  Vmax <- p$Vmax
  PATP <- y * Vmax * (TEA_conc/(Km + TEA_conc)) - c
  return(PATP)
}


### trying to find intersection point
TEA_conc <- concentrations [["O2", "conc.range"]]
TEA_conc

p1    <- as.list(parameters["O2_LA", ])
y1    <- p1$num.e * p1$H.e 
c1    <- cost * p1$length
Km1   <- p1$Km
Vmax1 <- p1$Vmax

p2    <- as.list(parameters["O2_HA", ])
y2    <- p2$num.e * p2$H.e 
c2    <- cost * p2$length
Km2   <- p2$Km
Vmax2 <- p2$Vmax

diff <- function(O2){
  PATP1 <- y1 * Vmax1 * (O2/(Km1 + O2)) - c1
  PATP2 <- y2 * Vmax2 * (O2/(Km2 + O2)) - c2
  return(abs(PATP1 - PATP2))
}

optimize(diff, lower=10^-10,upper=10)




PATP_O2_LA <- PATP("O2", "O2_LA")
PATP_O2_HA <- PATP("O2", "O2_HA")
PATP_NO3_NO2 <- PATP("NO3", "NO3_NO2")
PATP_NO2_N2O <- PATP("NO2", "NO2_N2O")
PATP_N2O_N2 <- PATP("N2O", "N2O_N2")
PATP_NO2_N2 <- PATP("NO2", "NO2_N2")
PATP_NO3_N2 <- PATP("NO3", "NO3_N2")
PATP_NO2_NH4 <- PATP("NO2", "NO2_NH4")
PATP_NO3_NH4 <- PATP("NO3", "NO3_NH4")





# calculate deltaG per reaction

dG <- function(TEA, TEAred, reaction){
  TEA_conc <- concentrations [[TEA, "conc.range"]]
  TEAred_conc <- concentrations [[TEAred, "conc.min"]]
  PED_conc <- concentrations [[CH2O, "conc.min"]]
  PEDox_conc <- concentrations [[HCO3, "conc.min"]]
  p <- as.list(parameters[reaction, ])
  num.e <- p$num.e
  Q <- 
  dG <-
  return(dG)
}

PATP_O2_LA <- PATP("O2", "O2_LA")
PATP_O2_HA <- PATP("O2", "O2_HA")
PATP_NO3_NO2 <- PATP("NO3", "NO3_NO2")
PATP_NO2_N2O <- PATP("NO2", "NO2_N2O")
PATP_N2O_N2 <- PATP("N2O", "N2O_N2")
PATP_NO2_N2 <- PATP("NO2", "NO2_N2")
PATP_NO3_N2 <- PATP("NO3", "NO3_N2")
PATP_NO2_NH4 <- PATP("NO2", "NO2_NH4")
PATP_NO3_NH4 <- PATP("NO3", "NO3_NH4")


##### plot ATP generation from all metabolisms


# specify colors for plotting

colors <- scico(n = 100, palette = "bamako")
reaction_colors <- c(
  "O2_LA"     = colors[90],
  "O2_HA"     = colors[80],
  "NO3_NO2"   = colors[60],
  "NO2_N2O"   = colors[55],
  "N2O_N2"    = colors[50],
  "NO2_N2"   = colors[45],
  "NO3_N2"   = colors[40],
  "NO2_NH4"   = colors[15],
  "NO3_NH4"   = colors[5]
)


# compile all metabolism info and format for plotting

O2.uM <- concentrations[["O2", "conc.range"]] * 10^6

plot_data <- data.frame(
  O2.uM = O2.uM,
  PATP_O2_LA = PATP_O2_LA,
  PATP_O2_HA = PATP_O2_HA,
  PATP_NO3_NO2 = PATP_NO3_NO2,
  PATP_NO2_N2O = PATP_NO2_N2O,
  PATP_N2O_N2 = PATP_N2O_N2,
  PATP_NO2_N2 = PATP_NO2_N2,
  PATP_NO3_N2 = PATP_NO3_N2,
  PATP_NO2_NH4 = PATP_NO2_NH4,
  PATP_NO3_NH4 = PATP_NO3_NH4
)

plot_data_long <- plot_data %>%
  pivot_longer(
    cols = starts_with("PATP"),
    names_to = "reaction",
    values_to = "PATP"
  ) %>%
  mutate(reaction = sub("^PATP_", "", reaction))%>%
  mutate(reaction = factor(reaction, levels = rownames(parameters)))


# specify a subset of the reactions to plot

plot_subset <- c(
  "O2_LA", 
  "O2_HA",
  "NO3_NO2",
#  "NO2_N2O",
#  "N2O_N2",
  "NO2_N2",
  "NO3_N2",
  "NO2_NH4",
  "NO3_NH4"
  )


# plot the subset of reactions 

plot_data_subset <- plot_data_long %>%
  filter(reaction %in% plot_subset)

label_data <- plot_data_subset %>%
  group_by(reaction) %>%
  filter(O2.uM == max(O2.uM))

plot <- ggplot(plot_data_subset, 
       aes(x = O2.uM, 
           y = PATP * 10^9, 
           color = reaction)) +
  geom_line(linewidth = 1.5, alpha = 0.7) +
  theme_classic() +
  scale_color_manual(values = reaction_colors) +
  geom_text(
    data = label_data,
    aes(label = reaction),
    hjust = 3,  # adjust position
    vjust = -.5,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    x = "O2 concentration (µM)",
    y = "production of ATP (nM/hr)",
    color = "Reaction"
  )

plot


# create a table of the parameter values for the subset reactions

param_subset <- parameters %>%
  as.data.frame() %>%
  mutate(Reaction = rownames(.)) %>%   # Turn row names into a column
  mutate(Km.uM = Km * 10^6) %>%
  mutate(Vmax.nM.hr = Vmax * 10^9) %>%  filter(Reaction %in% plot_subset) %>%      # Filter by row name
  select(Reaction, H.e, Km.uM, Vmax.nM.hr)

simple_theme <- ttheme_minimal(
  core = list(fg_params = list(fontface = "plain", cex = .7)),
  colhead = list(fg_params = list(fontface = "bold", cex = .7))
)

table_grob <- tableGrob(param_subset, rows = NULL, theme = simple_theme)
table_plot <- as.ggplot(table_grob)


# include the parameter values on the ggplot figure

legend <- get_legend(
  plot + theme(legend.position = "right")
)

legend_with_table <- plot_grid(
  legend,
  table_plot,
  ncol = 1,
  rel_heights = c(1, 1)  # adjust ratio to control table size
)

final_plot <- plot_grid(
  plot + theme(legend.position = "none"),  # remove legend from main plot
  legend_with_table,
  ncol = 2,
  rel_widths = c(1, 0.4)  # adjust width of the legend + table panel
)

final_plot



# include the substrate concentrations in the plot

conc_subset <- concentrations %>%
  mutate(conc.uM = conc.min * 10^6) %>%
  select (conc.uM) %>%
  t() %>%    # transpose rows to columns
  as.data.frame() %>%
  mutate(pH = -log10(H*10^-6))

conc_subset1 <- conc_subset %>%
  select (
    "NO3.uM" = "NO3", 
    "NO2.uM" = "NO2", 
    "N2O.uM" = "N2O", 
    "N2.uM" = "N2", 
    "NH4.uM" = "NH4"
    )
conc_subset2 <- conc_subset %>%
  select (
    "pH", 
    "CH2O.uM" = "CH2O",
    "HCO3.uM" = "HCO3"
    )

table_grob2 <- tableGrob(conc_subset1, rows = NULL, theme = simple_theme)
table_grob3 <- tableGrob(conc_subset2, rows = NULL, theme = simple_theme)
table_plot2 <- as.ggplot(table_grob2)
table_plot3 <- as.ggplot(table_grob3)

# boxed_table <- gtable_add_grob(
#   x = table_grob2,
#   grobs = rectGrob(gp = gpar(fill = NA, col = "black", lwd = 1)),
#   t = 1, l = 1, b = nrow(table_grob2), r = ncol(table_grob2)
# )

legend_with_tables <- plot_grid(
  legend,
  table_plot2,
  table_plot3,
  table_plot,
  ncol = 1,
  rel_heights = c(1, 0.2, 0.3, 1)  # Adjust these to control relative heights
)

final_plot <- plot_grid(
  plot + theme(legend.position = "none"),
  legend_with_tables,
  ncol = 2,
  rel_widths = c(1, 0.4)
)

final_plot


# TO DO

# Calculate and plot deltaG of metabolisms, calculate deltaG intersection values


str(parameters)



# Make plots of N metabolisms over gradient of NOx substrates

# Calculate intersection points, plot intersection values as functions of parameters

# Plot functions as per TEA or per PED



# 
# O2switch <- function(par,N,...){
#     #evaluate all terms except aerobic terms
#     #find largest non-aerobic term X
#     #find O2 where X - R(O2) = 0
# }

