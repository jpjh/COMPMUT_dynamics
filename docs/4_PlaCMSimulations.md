COMPMUT Dynamics 4: Simulations and analysis
================
jpjh
compiled Nov 2023, revised Aug 2024

## Simulate the plasmid compensatory mutations experiment

A general ODE model was developed to simulate competition between
plasmid-free, wild-type plasmid containing, and two different CMs. The
model is specified in `./model/COMPMOD_model.R`.

It describes:

- Z_f = no plasmid (i.e. SBW25)
- X_0 = plasmid-carrier (i.e. SBW25(pQBR57))
- X_1 = plasmid-carrier with CM type 1 (i.e. SBW25(pQBR57::plaCM))
- X_2 = plasmid-carrier with CM type 2 (i.e. SBW25::chrCM(pQBR57))
- X_0t = transconjugant from X_0
- X_1t = transconjugant from X_1 (i.e. wild-type recipients of
  pQBR57::plaCM)
- X_2t = transconjugant from X_2 (i.e. wild-type recipients of
  unameliorated plasmid from CM type 2)

### Identify suitable model parameters

See `3_ParameterEstimation.Rmd` for details.

``` r
(parameters <- read.table("./data/3_Parameters.txt", header=TRUE)) %>% kable(digits = 20)
```

| parameter  |        value |           sd |
|:-----------|-------------:|-------------:|
| alpha      | 5.412246e-01 | 5.535240e-02 |
| beta       | 8.180085e-01 | 4.290745e-02 |
| beta_chrCM | 9.690664e-01 | 6.987190e-02 |
| beta_plaCM | 9.449607e-01 | 4.131810e-02 |
| K          | 5.716250e+09 | 1.590994e+09 |
| gamma      | 4.573199e-12 | 1.013806e-12 |
| mu         | 4.125000e-02 | 0.000000e+00 |

### Analytical predictions

The analytical model is scaled to K, so the effective gamma is equal to
gamma \* K. The analytical model also describes the ’beta’s as the
overall growth rate of the plasmid-carrying populations, rather than the
relative effect as described in the parameters table. Add additional
variables describing this. Also, add a variable describing a plausible
minimum conjugation rate (gamma_min):

``` r
param_sd <- parameters %>% select(-value) %>% pivot_wider(names_from = "parameter", values_from = "sd")
param_values <- parameters %>% select(-sd) %>% 
  pivot_wider(names_from = "parameter", values_from = "value") %>%
  mutate(gamma_ = gamma * K,
         gamma_min = (gamma - (2*param_sd$gamma)) * (K - (2*param_sd$K)),
         beta_P = alpha * beta,
         beta_C = alpha * beta_chrCM,
         beta_Q = alpha * beta_plaCM)
```

Consider wild-type plasmid first.

Invasion if gamma \> μ(α-β)/(α-μ). Domination if gamma \> μ(α-β)/(β-μ).

``` r
(comparison_table <- data.frame(
  comparison = rep(c("no","chrCM","plaCM"),2),
  beta = rep(c(param_values$beta_P,param_values$beta_C,param_values$beta_Q),2),
  gamma = c(rep(param_values$gamma_,3), rep(param_values$gamma_min, 3)),
  alpha = param_values$alpha,
  mu = param_values$mu) %>%
   mutate(inv_threshold = (mu * (alpha - beta)/(alpha - mu)),
          inv_ratio = gamma / inv_threshold,
          invasion = ifelse(gamma > inv_threshold, "YES", "no"),
          dom_threshold = (mu * (alpha - beta)/(beta - mu)),
          dom_ratio = gamma / dom_threshold,
          domination = ifelse(gamma > dom_threshold, "YES", "no"))) %>% 
  select(-alpha, -mu) %>% kable()
```

| comparison | beta | gamma | inv_threshold | inv_ratio | invasion | dom_threshold | dom_ratio | domination |
|:---|---:|---:|---:|---:|:---|---:|---:|:---|
| no | 0.4427263 | 0.0261415 | 0.0081265 | 3.2168193 | YES | 0.0101203 | 2.583085 | YES |
| chrCM | 0.5244826 | 0.0261415 | 0.0013813 | 18.9254977 | YES | 0.0014291 | 18.291763 | YES |
| plaCM | 0.5114360 | 0.0261415 | 0.0024577 | 10.6366501 | YES | 0.0026134 | 10.002916 | YES |
| no | 0.4427263 | 0.0064512 | 0.0081265 | 0.7938434 | no | 0.0101203 | 0.637451 | no |
| chrCM | 0.5244826 | 0.0064512 | 0.0013813 | 4.6704148 | YES | 0.0014291 | 4.514022 | YES |
| plaCM | 0.5114360 | 0.0064512 | 0.0024577 | 2.6249015 | YES | 0.0026134 | 2.468509 | YES |

Predicts that the wild-type plasmid can invade and dominate a
plasmid-free population, which is in fact what we can see experimentally
e.g. in [Stevenson et al. 2017 doi:
10.1038/ismej.2017.42](https://www.nature.com/articles/ismej201742).
However, this is sensitive, and a 2.5x reduction in the product of gamma
and K can result in a mixed population, and a \>3.2x reduction can
result in the plasmid being lost. Such a pattern can plausibly occur,
e.g. with 2x SD reduction in both terms (see bottom lines of table).

Both CMs exceed the thresholds for invasion and domination, even with
the reduced estimate for conjugation rate. To examine the dynamics,
develop numerical simulations.

### Numerical simulations

The model `COMPMOD_model.R` has been developed as a Shiny app for
further investigation, accessible at
<https://jpjh.shinyapps.io/COMPMOD_shiny/>. Here, the model is run with
a range of parameters to investigate the competition between plaCM and
wild-type plasmid as described in the manuscript.

In this model, we set up the parameters so CM 2 is not present, and
fitness and conjugation rate for the CM 1 transconjugants are the same
as for the CM 1 donors.

First, run as a continuous model, i.e. mu = 0.04125 and dilFac = 1
(i.e. no dilution at each ‘transfer’, washout is continuous based on an
hourly turnover corresponding to 1:100 dilution every 24 hours). Run for
a long time and then can zoom in.

``` r
source("./model/COMPMOD_model.R")

parameters_mod0 <- c(alpha_Z_f = param_values$alpha,
                     alpha_X_0 = param_values$beta_P,
                     alpha_X_1 = param_values$beta_Q,
                     alpha_X_2 = param_values$beta_C,
                     alpha_X_0t = param_values$beta_P,
                     alpha_X_1t = param_values$beta_Q,
                     alpha_X_2t = param_values$beta_P,
                     gamma_X_0 =  param_values$gamma,
                     gamma_X_1 =  param_values$gamma,
                     gamma_X_2 =  param_values$gamma,
                     gamma_X_0t = param_values$gamma,
                     gamma_X_1t = param_values$gamma,
                     gamma_X_2t = param_values$gamma,
                     eta_Z_f = 0,
                     eta_X_0 = 0,
                     eta_X_1 = 0,
                     eta_X_2 = 0,
                     eta_X_0t = 0,
                     eta_X_1t = 0,
                     eta_X_2t = 0,
                     K = param_values$K,
                     mu = param_values$mu)

numTransfers <- 768
transferTime <- 24
dilFac <- 1

times <- seq(0, transferTime*numTransfers, by = 0.1) 

startfrac <- unname(parameters_mod0["K"] * 0.01)   # the population size at the start, 
# as estimated as ~ starting overnight cultures at saturation

initial_state <- c(Z_f = startfrac,   
                   X_0 = startfrac,  
                   X_1 = startfrac,  
                   X_2 = 0,
                   X_0t = 0,
                   X_1t = 0,
                   X_2t = 0)

state_mod01 <- initial_state * c(0.9,0.05,0.05,0,0,0,0)

event <- data.frame(var = rep(c("Z_f", "X_0", "X_1", "X_2", "X_0t", "X_1t", "X_2t"), numTransfers),
                    time = rep(seq(transferTime,transferTime*numTransfers,by=transferTime),each=7), 
                    value = c(dilFac), # set to 1 for the continuous model
                    method = c("mult"))

runModel <- function(state, func, times, parms, events, title="") {
  data.frame(ode(y = state, times = times,
                 func = func, parms = parms,
                 events = list(data = events))) %>%
    pivot_longer(cols=-time, names_to = "subpop", values_to = "density") %>%
    mutate(subpop = factor(subpop, levels=c('Z_f','X_0','X_1','X_2','X_0t','X_1t','X_2t'))) %>%
    group_by(time, subpop) %>% summarise(tot = sum(density), .groups='drop_last') %>%
    mutate(fraction = tot/sum(tot),
           model = title)
}

mod01 <- runModel(state = state_mod01,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod0,
                  title = "plaCM vs. wt 10:1 free")

state_mod02 <- initial_state * c(0.5,0.25,0.25,0,0,0,0)         

mod02 <- runModel(state = state_mod02,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod0,
                  title = "plaCM vs. wt 1:1 free")

state_mod03 <- initial_state * c(0.1,0.45,0.45,0,0,0,0)    

mod03 <- runModel(state = state_mod03,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod0,
                  title = "plaCM vs. wt 1:10 free")

state_mod04 <- initial_state * c(0,0.5,0.5,0,0,0,0)    

mod04 <- runModel(state = state_mod04,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod0,
                  title = "plaCM vs. wt")

mod0 <- rbind(mod01, mod02, mod03, mod04) %>%
  mutate(model = factor(model, levels=c("plaCM vs. wt 10:1 free",
                                        "plaCM vs. wt 1:1 free",
                                        "plaCM vs. wt 1:10 free",
                                        "plaCM vs. wt")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")),
         gamma = "gamma_plaCM = 1")

mod0 %>% filter(time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
  ggplot(aes(time, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#882255","#CC6677","#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "CM1 donor","CM1 transconjugant",
                             "CM2 donor","CM2 transconjugant",
                             "plasmid-free"),
                    name="subpopulation", drop=FALSE) +
  facet_grid(.~model) +
  ggtitle("as parameterised") +
  theme(legend.position="bottom")
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

These models poorly recapitulate the experimental observations.
Consistent with the analytical model, where plaCM has similar gamma to
the wild-type, plaCM wins.

Investigate modification of gamma_plaCM.

``` r
parameters_mod1_adj <- parameters_mod0
parameters_mod1_adj[c("gamma_X_1","gamma_X_1t")] <- 0.1*param_values$gamma

mod11 <- runModel(state = state_mod01,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 10:1 free")

mod12 <- runModel(state = state_mod02,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:1 free")

mod13 <- runModel(state = state_mod03,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:10 free")

mod14 <- runModel(state = state_mod04,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt")

parameters_mod1_adj[c("gamma_X_1","gamma_X_1t")] <- 0.01*param_values$gamma

mod21 <- runModel(state = state_mod01,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 10:1 free")

mod22 <- runModel(state = state_mod02,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:1 free")

mod23 <- runModel(state = state_mod03,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:10 free")

mod24 <- runModel(state = state_mod04,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt")

parameters_mod1_adj[c("gamma_X_1","gamma_X_1t")] <- 0.001*param_values$gamma

mod31 <- runModel(state = state_mod01,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 10:1 free")

mod32 <- runModel(state = state_mod02,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:1 free")

mod33 <- runModel(state = state_mod03,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:10 free")

mod34 <- runModel(state = state_mod04,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt")

parameters_mod1_adj[c("gamma_X_1","gamma_X_1t")] <- 0.05*param_values$gamma

mod41 <- runModel(state = state_mod01,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 10:1 free")

mod42 <- runModel(state = state_mod02,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:1 free")

mod43 <- runModel(state = state_mod03,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt 1:10 free")

mod44 <- runModel(state = state_mod04,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod1_adj,
                  title = "plaCM vs. wt")

mod1 <- rbind(mod11, mod12, mod13, mod14) %>%
  mutate(gamma = "gamma_plaCM = 0.1")

mod2 <- rbind(mod21, mod22, mod23, mod24) %>%
  mutate(gamma = "gamma_plaCM = 0.01")

mod3 <- rbind(mod31, mod32, mod33, mod34) %>%
  mutate(gamma = "gamma_plaCM = 0.001")

mod4 <- rbind(mod41, mod42, mod43, mod44) %>%
  mutate(gamma = "gamma_plaCM = 0.05")

mod <- rbind(mod0, mod1, mod2, mod3, mod4) %>%
  mutate(model = factor(model, levels=c("plaCM vs. wt 10:1 free",
                                        "plaCM vs. wt 1:1 free",
                                        "plaCM vs. wt 1:10 free",
                                        "plaCM vs. wt")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")))

(p_gamma_scan <- mod %>% filter(time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
    ggplot(aes(time, fraction, fill = subpop)) +
    geom_area() +
    scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                               "#882255","#CC6677","#DDDDDD"),
                      breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                      labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                               "CM1 donor","CM1 transconjugant",
                               "CM2 donor","CM2 transconjugant",
                               "plasmid-free"),
                      name="subpopulation", drop=FALSE) +
    facet_grid(gamma~model) +
    theme(legend.position="bottom"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
png("./figs/p11.png", width=7.2, height=3.6, units="in", res=300)
p_gamma_scan + 
  geom_vline(xintercept=192, linewidth=0.2, linetype="dotted") +
  theme_pub() + theme(legend.position="bottom")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("./figs/p12.png", width=3.6, height=3.6, units="in", res=300)
p_gamma_scan + 
  scale_x_continuous(limits=c(0,384),
                     breaks=c(0,192,384)) +
  geom_vline(xintercept=192, linewidth=0.2, linetype="dotted") +
  theme_pub()
```

    ## Warning: Removed 105280 rows containing non-finite outside the scale range
    ## (`stat_align()`).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Consistent with analysis, reducing gamma_plaCM by ~20x is sufficient to
recapitulate the general experimental results showing fluctuating
dynamics in the competition between wild-type and plaCM pQBR57.

Generate a suitable figure for the paper, using the 100x reduction.

``` r
library(patchwork)

p_gamma_fig <- mod %>% 
  filter(time %in% seq(0,transferTime*numTransfers,transferTime) &
           !(subpop %in% c("X_2", "X_2t")) &
           gamma == "gamma_plaCM = 0.01") %>%
  ggplot(aes(time/10000, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "plaCM donor","plaCM transconjugant",
                             "plasmid-free"),
                    name="Population") +
  labs(x="time (/10000)") + 
  facet_grid(.~model) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

p_gamma_fig_crop <- mod %>% 
  filter(time %in% seq(0,transferTime*numTransfers,transferTime) &
           time < 385 & 
           !(subpop %in% c("X_2", "X_2t")) & 
           gamma == "gamma_plaCM = 0.01") %>%
  ggplot(aes(time, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "plaCM donor","plaCM transconjugant",
                             "plasmid-free"),
                    name="Population") +
  labs(x="time") + 
  facet_grid(.~model) +
  theme(legend.position="bottom", strip.text.x=element_blank())

(doubleplot1 <- (p_gamma_fig_crop +
                   theme_pub() + theme(legend.position="right")) /
    (p_gamma_fig + theme_pub() + 
       theme(axis.text.x=element_text(angle=45, hjust=1), strip.text.x=element_blank())))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
tiff("./figs/FigS8.tiff", width=5.2, height=3.6, units="in", res=300)
doubleplot1
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Output the parameters used in this model.

``` r
parameters_mod1_adj %>% kable()
```

|            |            x |
|:-----------|-------------:|
| alpha_Z_f  | 5.412246e-01 |
| alpha_X_0  | 4.427263e-01 |
| alpha_X_1  | 5.114360e-01 |
| alpha_X_2  | 5.244826e-01 |
| alpha_X_0t | 4.427263e-01 |
| alpha_X_1t | 5.114360e-01 |
| alpha_X_2t | 4.427263e-01 |
| gamma_X_0  | 0.000000e+00 |
| gamma_X_1  | 0.000000e+00 |
| gamma_X_2  | 0.000000e+00 |
| gamma_X_0t | 0.000000e+00 |
| gamma_X_1t | 0.000000e+00 |
| gamma_X_2t | 0.000000e+00 |
| eta_Z_f    | 0.000000e+00 |
| eta_X_0    | 0.000000e+00 |
| eta_X_1    | 0.000000e+00 |
| eta_X_2    | 0.000000e+00 |
| eta_X_0t   | 0.000000e+00 |
| eta_X_1t   | 0.000000e+00 |
| eta_X_2t   | 0.000000e+00 |
| K          | 5.716250e+09 |
| mu         | 4.125000e-02 |

### Investigation of the transfer model.

Re-run these analyses, but rather than including mu as a continuous
turnover effect, model turnover directly through transfers.

Set mu to zero, and implement through ‘events’.

``` r
parameters_mod2 <- parameters_mod0
parameters_mod2["mu"] <- 0

event_batch <- event %>%
  mutate(value = c(0.01))

mod101 <- runModel(state = state_mod01,
                   func = COMPMOD_model,
                   times = times,
                   parms = parameters_mod2,
                   events = event_batch,
                   title = "plaCM vs. wt 10:1 free")

mod102 <- runModel(state = state_mod02,
                   func = COMPMOD_model,
                   times = times,
                   parms = parameters_mod2,
                   events = event_batch,
                   title = "plaCM vs. wt 1:1 free")

mod103 <- runModel(state = state_mod03,
                   func = COMPMOD_model,
                   times = times,
                   parms = parameters_mod2,
                   events = event_batch,
                   title = "plaCM vs. wt 1:10 free")

mod104 <- runModel(state = state_mod04,
                   func = COMPMOD_model,
                   times = times,
                   parms = parameters_mod2,
                   events = event_batch,
                   title = "plaCM vs. wt")

mod10 <- rbind(mod101, mod102, mod103, mod104) %>%
  mutate(model = factor(model, levels=c("plaCM vs. wt 10:1 free",
                                        "plaCM vs. wt 1:1 free",
                                        "plaCM vs. wt 1:10 free",
                                        "plaCM vs. wt")))

mod10 %>% filter(time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
  ggplot(aes(time, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#882255","#CC6677","#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "CM1 donor","CM1 transconjugant",
                             "CM2 donor","CM2 transconjugant",
                             "plasmid-free"),
                    name="subpopulation", drop=FALSE) +
  facet_grid(.~model) +
  theme(legend.position="bottom")
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

The plasmid is lost from the system if the mean measured parameters are
used in the transfer model. This is likely because the population is
maintained at lower capacity, which reduces the effective conjugation
rate, preventing maintenance of the plasmid through infectious transfer.

Adjust relevant parameters within realistic bounds to investigate. The
main values likely to result in maintenance are K, alpha, and gamma.

Increase alpha, gamma and K by 2.5x SD.

``` r
max_gamma <- param_values$gamma + 2.5*(param_sd$gamma)
max_K <- param_values$K + 2.5*(param_sd$K)
max_alpha_modifier <- (param_values$alpha + (2.5*param_sd$alpha))/param_values$alpha

parameters_mod3_adj <- c(alpha_Z_f = param_values$alpha * max_alpha_modifier,
                         alpha_X_0 = param_values$beta_P * max_alpha_modifier,
                         alpha_X_1 = param_values$beta_Q * max_alpha_modifier,
                         alpha_X_2 = param_values$beta_C * max_alpha_modifier,
                         alpha_X_0t = param_values$beta_P * max_alpha_modifier,
                         alpha_X_1t = param_values$beta_Q * max_alpha_modifier,
                         alpha_X_2t = param_values$beta_P * max_alpha_modifier,
                         gamma_X_0 =  max_gamma,
                         gamma_X_1 =  max_gamma,
                         gamma_X_2 =  max_gamma,
                         gamma_X_0t = max_gamma,
                         gamma_X_1t = max_gamma,
                         gamma_X_2t = max_gamma,
                         eta_Z_f = 0,
                         eta_X_0 = 0,
                         eta_X_1 = 0,
                         eta_X_2 = 0,
                         eta_X_0t = 0,
                         eta_X_1t = 0,
                         eta_X_2t = 0,
                         K = max_K,
                         mu = 0)

mod201 <- runModel(state = state_mod01,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 10:1 free")

mod202 <- runModel(state = state_mod02,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:1 free")

mod203 <- runModel(state = state_mod03,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:10 free")

mod204 <- runModel(state = state_mod04,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt")

parameters_mod3_adj[c("gamma_X_1","gamma_X_1t")] <- max_gamma*0.1

mod211 <- runModel(state = state_mod01,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 10:1 free")

mod212 <- runModel(state = state_mod02,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:1 free")

mod213 <- runModel(state = state_mod03,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:10 free")

mod214 <- runModel(state = state_mod04,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt")

parameters_mod3_adj[c("gamma_X_1","gamma_X_1t")] <- max_gamma*0.01

mod221 <- runModel(state = state_mod01,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 10:1 free")

mod222 <- runModel(state = state_mod02,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:1 free")

mod223 <- runModel(state = state_mod03,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:10 free")

mod224 <- runModel(state = state_mod04,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt")

parameters_mod3_adj[c("gamma_X_1","gamma_X_1t")] <- max_gamma*0.001

mod231 <- runModel(state = state_mod01,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 10:1 free")

mod232 <- runModel(state = state_mod02,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:1 free")

mod233 <- runModel(state = state_mod03,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:10 free")

mod234 <- runModel(state = state_mod04,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt")


mod20 <- rbind(mod201, mod202, mod203, mod204) %>%
  mutate(gamma = "gamma_plaCM = 1")

mod21 <- rbind(mod211, mod212, mod213, mod214) %>%
  mutate(gamma = "gamma_plaCM = 0.1")

mod22 <- rbind(mod221, mod222, mod223, mod224) %>%
  mutate(gamma = "gamma_plaCM = 0.01")

mod23 <- rbind(mod231, mod232, mod233, mod234) %>%
  mutate(gamma = "gamma_plaCM = 0.001")

mod_transfer_adj <- rbind(mod20, mod21, mod22, mod23) %>%
  mutate(model = factor(model, levels=c("plaCM vs. wt 10:1 free",
                                        "plaCM vs. wt 1:1 free",
                                        "plaCM vs. wt 1:10 free",
                                        "plaCM vs. wt")),
         gamma = factor(gamma,
                        levels=c("gamma_plaCM = 1", "gamma_plaCM = 0.1", "gamma_plaCM = 0.01", "gamma_plaCM = 0.001")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")))

(p_transfer_gamma_adj_scan <- mod_transfer_adj %>% 
    filter(time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
    ggplot(aes(time, fraction, fill = subpop)) +
    geom_area() +
    scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                               "#882255","#CC6677","#DDDDDD"),
                      breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                      labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                               "CM1 donor","CM1 transconjugant",
                               "CM2 donor","CM2 transconjugant",
                               "plasmid-free"),
                      name="subpopulation", drop=FALSE) +
    facet_grid(gamma~model) +
    theme(legend.position="bottom") +
    ggtitle("transfer model with upper limits for K and gamma"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
png("./figs/p21.png", width=7.2, height=3.6, units="in", res=300)
p_transfer_gamma_adj_scan + 
  geom_vline(xintercept=192, linewidth=0.2, linetype="dotted") +
  theme_pub() + theme(legend.position="bottom")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
png("./figs/p22.png", width=3.6, height=3.6, units="in", res=300)
p_transfer_gamma_adj_scan + 
  scale_x_continuous(limits=c(0,384),
                     breaks=c(0,192,384)) +
  geom_vline(xintercept=192, linewidth=0.2, linetype="dotted") +
  theme_pub()
```

    ## Warning: Removed 84224 rows containing non-finite outside the scale range
    ## (`stat_align()`).

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Here, we can clearly see the oscillatory effect when gamma_Q is reduced
to 0.1x gamma_P.

Make a figure for the manuscript describing these features. Use gamma_Q
that is 75x reduced compared with gamma_P.

``` r
parameters_mod3_adj[c("gamma_X_1","gamma_X_1t")] <- max_gamma*0.0133

mod241 <- runModel(state = state_mod01,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 10:1 free")

mod242 <- runModel(state = state_mod02,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:1 free")

mod243 <- runModel(state = state_mod03,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt 1:10 free")

mod244 <- runModel(state = state_mod04,
                   func = COMPMOD_model,
                   times = times, events = event_batch,
                   parms = parameters_mod3_adj,
                   title = "plaCM vs. wt")

mod_transfer_adj_fig <- rbind(mod241, mod242, mod243, mod244) %>%
  mutate(model = factor(model, levels=c("plaCM vs. wt 10:1 free",
                                        "plaCM vs. wt 1:1 free",
                                        "plaCM vs. wt 1:10 free",
                                        "plaCM vs. wt")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")))

library(patchwork)

p_transfer_gamma_adj_fig <- mod_transfer_adj_fig %>% 
  filter(time %in% seq(0,transferTime*numTransfers,transferTime) &
           !(subpop %in% c("X_2", "X_2t"))) %>%
  ggplot(aes(time/24, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "plaCM donor","plaCM transconjugant",
                             "plasmid-free"),
                    name="Population") +
  labs(x="time (d)") + 
  facet_grid(.~model) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

p_transfer_gamma_adj_fig_crop <- mod_transfer_adj_fig %>% 
  filter(time %in% seq(0,transferTime*numTransfers,transferTime) &
           time < 193 & 
           !(subpop %in% c("X_2", "X_2t"))) %>%
  ggplot(aes(time/24, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "plaCM donor","plaCM transconjugant",
                             "plasmid-free"),
                    name="Population") +
  labs(x="time (d)") + 
  facet_grid(.~model) +
  theme(legend.position="bottom", strip.text.x=element_blank())

(doubleplot <- (p_transfer_gamma_adj_fig_crop + theme_pub() +
                  theme(legend.position="right")) /
    (p_transfer_gamma_adj_fig +
       theme_pub() +
       theme(axis.text.x=element_text(angle=45, hjust=1), strip.text.x=element_blank())))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
png("./figs/p41.png", width=5.2, height=3.6, units="in", res=300)
doubleplot
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
tiff("./figs/FigS8_alt.tiff", width=5.2, height=3.6, units="in", res=300)
doubleplot
dev.off()
```

    ## quartz_off_screen 
    ##                 2

This plot shows the full dynamics in the bottom panels, revealing the
damped oscillatory behaviour, and the short-term dynamics equivalent to
the experimental results in the top panels.

Output the parameters used in this model.

``` r
parameters_mod3_adj %>% kable()
```

|            |            x |
|:-----------|-------------:|
| alpha_Z_f  | 6.796056e-01 |
| alpha_X_0  | 5.559232e-01 |
| alpha_X_1  | 6.422006e-01 |
| alpha_X_2  | 6.585830e-01 |
| alpha_X_0t | 5.559232e-01 |
| alpha_X_1t | 6.422006e-01 |
| alpha_X_2t | 5.559232e-01 |
| gamma_X_0  | 0.000000e+00 |
| gamma_X_1  | 0.000000e+00 |
| gamma_X_2  | 0.000000e+00 |
| gamma_X_0t | 0.000000e+00 |
| gamma_X_1t | 0.000000e+00 |
| gamma_X_2t | 0.000000e+00 |
| eta_Z_f    | 0.000000e+00 |
| eta_X_0    | 0.000000e+00 |
| eta_X_1    | 0.000000e+00 |
| eta_X_2    | 0.000000e+00 |
| eta_X_0t   | 0.000000e+00 |
| eta_X_1t   | 0.000000e+00 |
| eta_X_2t   | 0.000000e+00 |
| K          | 9.693735e+09 |
| mu         | 0.000000e+00 |

The continuous transfer model figure was selected as the example as this
was consistent with the plots provided in Fig S1.

## Simulate the experiments presented in Figures 3 and 5

To test that model outputs qualitatively conform to the experimental
results presented in Figures 3 and 5 (head-to-head competitions of plaCM
and chrCM), set up another set of simulations following these
experiments.

For these experiments, simulate using the continuous transfer model with
general parameters as above:

``` r
parameters_mod1_adj %>% kable()
```

|            |            x |
|:-----------|-------------:|
| alpha_Z_f  | 5.412246e-01 |
| alpha_X_0  | 4.427263e-01 |
| alpha_X_1  | 5.114360e-01 |
| alpha_X_2  | 5.244826e-01 |
| alpha_X_0t | 4.427263e-01 |
| alpha_X_1t | 5.114360e-01 |
| alpha_X_2t | 4.427263e-01 |
| gamma_X_0  | 0.000000e+00 |
| gamma_X_1  | 0.000000e+00 |
| gamma_X_2  | 0.000000e+00 |
| gamma_X_0t | 0.000000e+00 |
| gamma_X_1t | 0.000000e+00 |
| gamma_X_2t | 0.000000e+00 |
| eta_Z_f    | 0.000000e+00 |
| eta_X_0    | 0.000000e+00 |
| eta_X_1    | 0.000000e+00 |
| eta_X_2    | 0.000000e+00 |
| eta_X_0t   | 0.000000e+00 |
| eta_X_1t   | 0.000000e+00 |
| eta_X_2t   | 0.000000e+00 |
| K          | 5.716250e+09 |
| mu         | 4.125000e-02 |

Parameterise initially with gamma_plaCM = gamma_chrCM, with 50:50
CM:free and 50:50 plaCM:chrCM.

``` r
parameters_mod4_adj <- parameters_mod1_adj
parameters_mod4_adj[c("gamma_X_1","gamma_X_1t")] <- param_values$gamma

initial_state_mod4 <- c(Z_f = startfrac,   
                        X_0 = 0,  
                        X_1 = startfrac,  
                        X_2 = startfrac,
                        X_0t = 0,
                        X_1t = 0,
                        X_2t = 0)

state_mod4 <- initial_state_mod4 * c(0.5,0,0.25,0.25,0,0,0)

mod41 <- runModel(state = state_mod4,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod4_adj,
                  title = "plaCM vs. chrCM 1:1 free")

mod41 %>% filter(time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
  ggplot(aes(time, fraction, fill = subpop)) +
  geom_area() +
  scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                             "#882255","#CC6677","#DDDDDD"),
                    breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                    labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                             "CM1 donor","CM1 transconjugant",
                             "CM2 donor","CM2 transconjugant",
                             "plasmid-free"),
                    name="subpopulation", drop=FALSE) +
  facet_grid(.~model) +
  theme(legend.position="bottom")
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Reorder colours to resemble experiments and shorten timeframe. Note that
in the experimental data presented in Figure 3, the CM2 transconjugants
(pink above) would appear indistinguishable from plasmid-free (grey).
Use a darker shade of grey to indicate and distinguish.

``` r
(fig3a_modelled_equal <- mod41 %>% filter(time < 1000 & time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
   mutate(subpop = factor(subpop, levels=c("X_0","X_0t","X_2","X_1","X_1t","X_2t","Z_f"))) %>%
   ggplot(aes(time, fraction, fill = subpop)) +
   geom_area() +
   scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                              "#882255","#BBBBBB","#DDDDDD"),
                     breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                     labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                              "CM1 donor","CM1 transconjugant",
                              "CM2 donor","CM2 transconjugant",
                              "plasmid-free"),
                     name="subpopulation", drop=FALSE) +
   facet_grid(.~model) +
   theme(legend.position="right"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
png("./figs/p3a_modelled_equal.png", width=6.2, height=3.6, units="in", res=300)
fig3a_modelled_equal
dev.off()
```

    ## quartz_off_screen 
    ##                 2

This resembles the experimental result somewhat, resulting in ultimate
success of CM2 (parameterised as chrCM), but there are a lot more plaCM
transconjugants than expected given the experimental data.

However, this is assuming equal conjugation rates. Reduce conjugation
rate as discussed in the analyses above.

``` r
parameters_mod4_adj[c("gamma_X_1","gamma_X_1t")] <- param_values$gamma * 0.01

mod42 <- runModel(state = state_mod4,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod4_adj,
                  title = "plaCM vs. chrCM 1:1 free")

(fig3a_modelled <- mod42 %>% filter(time < 1000 & time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
    mutate(subpop = factor(subpop, levels=c("X_0","X_0t","X_2","X_1","X_1t","X_2t","Z_f"))) %>%
    ggplot(aes(time, fraction, fill = subpop)) +
    geom_area() +
    scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                               "#882255","#BBBBBB","#DDDDDD"),
                      breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                      labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                               "CM1 donor","CM1 transconjugant",
                               "CM2 donor","CM2 transconjugant",
                               "plasmid-free"),
                      name="subpopulation", drop=FALSE) +
    facet_grid(.~model) +
    theme(legend.position="right"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
png("./figs/p3a_modelled.png", width=6.2, height=3.6, units="in", res=300)
fig3a_modelled
dev.off()
```

    ## quartz_off_screen 
    ##                 2

This looks similar to Figure 3a: initial peak of plasmid-free, ultimate
invasion of CM, with chrCM dominating over plaCM when the recipients (or
wt-carrying strains) are outcompeted.

Using these conjugation rates, run across starting proportions to get a
figure analogous to Figure 5/6.

``` r
state_mod50 <- initial_state_mod4 * c(0.9,0,0.05,0.05,0,0,0)
state_mod51 <- initial_state_mod4 * c(0.5,0,0.25,0.25,0,0,0)
state_mod52 <- initial_state_mod4 * c(0.1,0,0.45,0.45,0,0,0)
state_mod53 <- initial_state_mod4 * c(0,0,0.5,0.5,0,0,0)

mod50 <- runModel(state = state_mod50,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod4_adj,
                  title = "10:1")

mod51 <- runModel(state = state_mod51,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod4_adj,
                  title = "1:1")

mod52 <- runModel(state = state_mod52,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod4_adj,
                  title = "1:10")

mod53 <- runModel(state = state_mod53,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod4_adj,
                  title = "none")

mod5X <- rbind(mod50, mod51, mod52, mod53) %>%
  mutate(model = factor(model, levels=c("10:1",
                                        "1:1",
                                        "1:10",
                                        "none")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")))

(fig5_modelled <- mod5X %>% 
    filter(time < 1000 & time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
    mutate(subpop = factor(subpop, levels=c("X_0","X_0t","X_2","X_2t","X_1","X_1t","Z_f"))) %>%
    ggplot(aes(time, fraction, fill = subpop)) +
    scale_x_continuous(breaks = c(0,300,600,900),
                       labels = c(0,300,600,900)) + 
    geom_area() +
    scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                               "#882255","#CC6677","#DDDDDD"),
                      breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                      labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                               "CM1 donor","CM1 transconjugant",
                               "CM2 donor","CM2 transconjugant",
                               "plasmid-free"),
                      name="subpopulation", drop=FALSE) +
    labs(x="time") + 
    facet_grid(.~model) +
    theme(legend.position="right"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
png("./figs/p5_modelled.png", width=8.2, height=4.6, units="in", res=300)
fig5_modelled
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Test the effect of plasmid ‘weaponisation’

To test whether costly plasmid conjugation from chrCM could be
responsible for the dynamics, repeat the previous analyses but without
conjugation from chrCM.

``` r
parameters_mod6_adj <- parameters_mod4_adj
parameters_mod6_adj["gamma_X_2"] <- 0

mod60 <- runModel(state = state_mod50,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod6_adj,
                  title = "10:1")

mod61 <- runModel(state = state_mod51,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod6_adj,
                  title = "1:1")

mod62 <- runModel(state = state_mod52,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod6_adj,
                  title = "1:10")

mod63 <- runModel(state = state_mod53,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod6_adj,
                  title = "none")

mod6X <- rbind(mod60, mod61, mod62, mod63) %>%
  mutate(model = factor(model, levels=c("10:1",
                                        "1:1",
                                        "1:10",
                                        "none")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")))

(p8 <- mod6X %>% 
    filter(time < 1000 & time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
    mutate(subpop = factor(subpop, levels=c("X_0","X_0t","X_2","X_2t","X_1","X_1t","Z_f"))) %>%
    ggplot(aes(time, fraction, fill = subpop)) +
    geom_area() +
    scale_x_continuous(breaks = c(0,300,600,900),
                       labels = c(0,300,600,900)) + 
    scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                               "#882255","#CC6677","#DDDDDD"),
                      breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                      labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                               "CM1 donor","CM1 transconjugant",
                               "CM2 donor","CM2 transconjugant",
                               "plasmid-free"),
                      name="subpopulation", drop=FALSE) +
    labs(x="time") + 
    facet_grid(.~model) +
    theme(legend.position="bottom"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
png("./figs/p8.png", width=8.2, height=4.6, units="in", res=300)
p8
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Both compensatory mutations are lost in this case, clearly illustrating
the importance of the plasmid ‘weaponisation’ on the outcomes.

Output the parameters for further exploration in the Shiny app.

``` r
parameters_mod6_adj
```

    ##    alpha_Z_f    alpha_X_0    alpha_X_1    alpha_X_2   alpha_X_0t   alpha_X_1t 
    ## 5.412246e-01 4.427263e-01 5.114360e-01 5.244826e-01 4.427263e-01 5.114360e-01 
    ##   alpha_X_2t    gamma_X_0    gamma_X_1    gamma_X_2   gamma_X_0t   gamma_X_1t 
    ## 4.427263e-01 4.573199e-12 4.573199e-14 0.000000e+00 4.573199e-12 4.573199e-14 
    ##   gamma_X_2t      eta_Z_f      eta_X_0      eta_X_1      eta_X_2     eta_X_0t 
    ## 4.573199e-12 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
    ##     eta_X_1t     eta_X_2t            K           mu 
    ## 0.000000e+00 0.000000e+00 5.716250e+09 4.125000e-02

Does the importance of weaponisation result hold if conjugation rates
are not reduced for plaCM?

``` r
parameters_mod7_adj <- parameters_mod6_adj

parameters_mod7_adj[c("gamma_X_1","gamma_X_1t")] <- param_values$gamma

mod70 <- runModel(state = state_mod50,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod7_adj,
                  title = "10:1")

mod71 <- runModel(state = state_mod51,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod7_adj,
                  title = "1:1")

mod72 <- runModel(state = state_mod52,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod7_adj,
                  title = "1:10")

mod73 <- runModel(state = state_mod53,
                  func = COMPMOD_model,
                  times = times,
                  events = event,
                  parms = parameters_mod7_adj,
                  title = "none")

mod7X <- rbind(mod70, mod71, mod72, mod73) %>%
  mutate(model = factor(model, levels=c("10:1",
                                        "1:1",
                                        "1:10",
                                        "none")),
         subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                          "X_2","X_2t","Z_f")))

(p8 <- mod7X %>% 
    filter(time < 1000 & time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
    mutate(subpop = factor(subpop, levels=c("X_0","X_0t","X_2","X_2t","X_1","X_1t","Z_f"))) %>%
    ggplot(aes(time, fraction, fill = subpop)) +
    scale_x_continuous(breaks = c(0,300,600,900),
                       labels = c(0,300,600,900)) + 
    geom_area() +
    scale_fill_manual(values=c("#ccb333","#e6d898","#44AA99","#92d3c8",
                               "#882255","#CC6677","#DDDDDD"),
                      breaks=c("X_0","X_0t","X_1","X_1t","X_2","X_2t","Z_f"),
                      labels=c("wild-type plasmid donor","wild-type plasmid transconjugant",
                               "CM1 donor","CM1 transconjugant",
                               "CM2 donor","CM2 transconjugant",
                               "plasmid-free"),
                      name="subpopulation", drop=FALSE) +
    labs(x="time") + 
    facet_grid(.~model) +
    theme(legend.position="bottom"))
```

![](4_PlaCMSimulations_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
parameters_mod7_adj
```

    ##    alpha_Z_f    alpha_X_0    alpha_X_1    alpha_X_2   alpha_X_0t   alpha_X_1t 
    ## 5.412246e-01 4.427263e-01 5.114360e-01 5.244826e-01 4.427263e-01 5.114360e-01 
    ##   alpha_X_2t    gamma_X_0    gamma_X_1    gamma_X_2   gamma_X_0t   gamma_X_1t 
    ## 4.427263e-01 4.573199e-12 4.573199e-12 0.000000e+00 4.573199e-12 4.573199e-12 
    ##   gamma_X_2t      eta_Z_f      eta_X_0      eta_X_1      eta_X_2     eta_X_0t 
    ## 4.573199e-12 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
    ##     eta_X_1t     eta_X_2t            K           mu 
    ## 0.000000e+00 0.000000e+00 5.716250e+09 4.125000e-02

Here, because plaCM can invade, the stage is set for chrCM to invade, in
a similar manner to the ‘none’ competition above. Make the supplementary
figure, using the reduced conjugation rates described above.

``` r
tiff("./figs/FigS7_alt.tiff", width=5.2, height=3.6, units="in", res=300)
(fig5_modelled + labs(tag = "A") + theme_pub() + theme(legend.position="right")) / 
  (p8 + labs(tag = "B") + theme_pub())
dev.off()
```

    ## quartz_off_screen 
    ##                 2

------------------------------------------------------------------------

**[Back to index.](../README.md)**
