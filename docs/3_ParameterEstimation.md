COMPMUT Dynamics 3: Estimation and gathering of relevant parameters
================
jpjh
compiled November 2023

## Background

#### Growth rates

Estimate growth rates from a plate-reader experiment, using the R
package `gcplyr` and values identified in a preliminary analysis.
Measurement cycles every 15 minutes.

``` r
library(gcplyr)
```

    ## ## 
    ## ## gcplyr (Version 1.6.0, Build Date: 2023-09-13)
    ## ## See http://github.com/mikeblazanin/gcplyr for additional documentation
    ## ## Please cite software as:
    ## ##   Blazanin, Michael. 2023. gcplyr: an R package for microbial growth
    ## ##   curve data analysis. bioRxiv doi: 10.1101/2023.04.30.538883
    ## ##

``` r
read.csv("./data/4_WTGrowthCurve.csv", header=TRUE) %>%
  group_by(well, replicate) %>%
  mutate(Time = cycle / 4,
         smoothed = smooth_data(x = Time, y = value,
                                sm_method = "moving-average", window_width_n = 9),
         smderivpercap_3 = calc_deriv(x = Time, y = smoothed,
                                       percapita = TRUE, blank = 0, window_width_n = 3)) %>%
  summarize(max_percap = max_gc(smderivpercap_3, na.rm = TRUE),
            max_percap_time = extr_val(Time, which_max_gc(smderivpercap_3)),
            doub_time = doubling_time(y = max_percap)) %>% 
  ungroup() %>%
  summarise(mean = mean(max_percap),
            sd = sd(max_percap)) %>% kable()
```

    ## `summarise()` has grouped output by 'well'. You can override using the
    ## `.groups` argument.

|      mean |        sd |
|----------:|----------:|
| 0.5412246 | 0.0553524 |

Results in an estimate:

- alpha = 0.54 ± 0.06.

Generate a table of values.

``` r
parameters <- data.frame(
  parameter = c("alpha"),
  value = c(0.5412246),
  sd = c(0.0553524)
)
```

This value is consistent with other unpublished measurements for SBW25
wild-type growth rate in the plate reader.

#### Plasmid fitness costs

The actual mean fitness values for the chromosomal and plasmid
mutations, from [Hall et al. 2021 PLoS
Biology](dx.doi.org/10.1371/journal.pbio.3001225) (see also
www.github.com/jpjh/COMPMUT), were:

``` r
source("https://raw.githubusercontent.com/jpjh/COMPMUT/main/functions/calculate_fitness.R")
id_cols <- c("experiment","replicate","host","plasmid")

chromosomal_compensation_marker <- calculate_fitness(
  read.csv("https://raw.githubusercontent.com/jpjh/COMPMUT/main/data/COMPMUT_exp_data_1.csv", header=TRUE), id_cols) %>%
  filter(host=="wild-type" & plasmid=="plasmid-free") %>% 
  summarise(mean = mean(W_gm)) %>% pull()    
```

    ## Warning: Using `all_of()` outside of a selecting function was deprecated in tidyselect
    ## 1.2.0.
    ## ℹ See details at
    ##   <https://tidyselect.r-lib.org/reference/faq-selection-context.html>
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
(chromosomal_compensation_fitness <- calculate_fitness(
  read.csv("https://raw.githubusercontent.com/jpjh/COMPMUT/main/data/COMPMUT_exp_data_1.csv", header=TRUE), id_cols) %>% 
    filter(host %in% c("wild-type", "delta_PFLU4242") &
                      plasmid %in% c("plasmid-free", "pQBR57")) %>%
    mutate(W = W_gm/chromosomal_compensation_marker) %>%
  group_by(host, plasmid) %>% summarise(fitness = mean(W), sd = sd(W)))
```

    ## `summarise()` has grouped output by 'host'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 4 × 4
    ## # Groups:   host [2]
    ##   host           plasmid      fitness     sd
    ##   <chr>          <chr>          <dbl>  <dbl>
    ## 1 delta_PFLU4242 pQBR57         0.969 0.0699
    ## 2 delta_PFLU4242 plasmid-free   0.989 0.0554
    ## 3 wild-type      pQBR57         0.824 0.0520
    ## 4 wild-type      plasmid-free   1     0.0332

``` r
id_colsp <- c("replicate","marker","plasmid")

plasmid_compensation_marker <- calculate_fitness(
  read.csv("https://raw.githubusercontent.com/jpjh/COMPMUT/main/data/COMPMUT_exp_data_2.csv", header=TRUE), id_colsp) %>%
  filter(plasmid=="plasmid-free") %>% 
  summarise(mean = mean(W_gm)) %>% pull()    

(plasmid_compensation_fitness <- calculate_fitness(
  read.csv("https://raw.githubusercontent.com/jpjh/COMPMUT/main/data/COMPMUT_exp_data_2.csv", header=TRUE), id_colsp) %>%
    mutate(W = ifelse(marker == "Sm", 1/(W_gm / plasmid_compensation_marker),
                               (W_gm / plasmid_compensation_marker)),
           host = "wild-type") %>%
  group_by(host, plasmid) %>% summarise(fitness = mean(W), sd = sd(W)))
```

    ## `summarise()` has grouped output by 'host'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 3 × 4
    ## # Groups:   host [1]
    ##   host      plasmid      fitness     sd
    ##   <chr>     <chr>          <dbl>  <dbl>
    ## 1 wild-type pQBR57_V100A   0.945 0.0413
    ## 2 wild-type pQBR57_anc     0.812 0.0338
    ## 3 wild-type plasmid-free   1     0.0574

``` r
(compensations <- bind_rows(chromosomal_compensation_fitness, plasmid_compensation_fitness) %>%
  filter(plasmid!="plasmid-free") %>%
  bind_cols(data.frame(strain = c("chromosome compensation","no compensation","plasmid compensation","no compensation")))) %>%
  kable()
```

| host           | plasmid      |   fitness |        sd | strain                  |
|:---------------|:-------------|----------:|----------:|:------------------------|
| delta_PFLU4242 | pQBR57       | 0.9690664 | 0.0698719 | chromosome compensation |
| wild-type      | pQBR57       | 0.8238168 | 0.0520151 | no compensation         |
| wild-type      | pQBR57_V100A | 0.9449607 | 0.0413181 | plasmid compensation    |
| wild-type      | pQBR57_anc   | 0.8122002 | 0.0337998 | no compensation         |

i.e.

- beta has two estimates: 0.82 ± 0.05; 0.81 ± 0.03
- beta_chrCM = 0.97 ± 0.07
- beta_plaCM = 0.94 ± 0.04

``` r
betas <- data.frame(
  parameter = c("beta","beta_chrCM","beta_plaCM"),
  value = c(mean(c(0.8238168,0.8122002)), 0.9690664, 0.9449607),
  sd = c(mean(c(0.0520151,0.0337998)), 0.0698719, 0.0413181)
)

parameters <- rbind(parameters, betas)
```

#### Carrying capacity

Taken from the endpoint of the plasmid-free competitions in the PLoS
Biology paper.

``` r
read.csv("https://raw.githubusercontent.com/jpjh/COMPMUT/main/data/COMPMUT_exp_data_1.csv", header=TRUE) %>% 
  filter(timepoint=="end" & plasmid=="plasmid-free" & host=="wild-type") %>%
  summarise(cfu_ml = mean((count_white + count_blue) * (1000/spread) * 10^dilution)/1e6,
            sd = sd((count_white + count_blue) * (1000/spread) * 10^dilution)/1e6) %>% kable()
```

|  cfu_ml |       sd |
|--------:|---------:|
| 5716.25 | 1590.994 |

- K = 5.7e9 ± 1.6e9

``` r
parameters[5,] <- c("K",5.71625e9,1.590994e9)
```

#### Gamma

Taken from the PLoS Biology paper, with estimate from the uncompensated
strains.

``` r
d7 <- read.csv("https://raw.githubusercontent.com/jpjh/COMPMUT/main/data/COMPMUT_exp_data_3-1.csv", header=TRUE)

d7_start_D <- d7 %>% 
  filter(timepoint == "start_donor") %>% select(!recipient & !timepoint) %>%
  mutate(D0 = (10^dilution) * (1000/spread) * count_donor * 0.5 * 0.01)

d7_start_R_pf <- d7 %>%
  filter(timepoint == "start_pf_recipient") %>% select(!donor & !plasmid & !timepoint) %>%
  mutate(R0 = (10^dilution) * (1000/spread) * count_recipient * 0.5 * 0.01)

d7_start_pf <- full_join(d7_start_D, d7_start_R_pf, by=c("replicate"),
                         suffix=c("donor","recipient")) %>%
  select(replicate, donor, plasmid, recipient, D0, R0)

d7_end_DR <- d7 %>% filter(timepoint == "end_DR") %>%
  mutate(D = (10^dilution) * (1000/spread) * count_donor,
         R = (10^dilution) * (1000/spread) * count_recipient) %>%
  select(replicate, donor, plasmid, recipient, D, R)

d7_end_T <- d7 %>% filter(timepoint == "end_T") %>%
  mutate(TC = (10^dilution) * (1000/spread) * count_recipient) %>%
  select(replicate, donor, plasmid, recipient, TC)

left_join(d7_start_pf, d7_end_DR) %>%
  left_join(d7_end_T) %>%
  mutate(N0 = D0 + R0,
         N = D + R,
         gamma = (log(N/N0)/24) * log(1+((TC/R)*(N/D))) * (1/(N-N0))) %>%
  filter(donor == "wild-type" & plasmid == "pQBR57_anc") %>%
  summarise(mean = mean(gamma),
            sd = sd(gamma),
            log10mean = mean(log10(gamma)),
            log10sd = sd(log10(gamma))) %>% kable(digits=20)
```

    ## Joining with `by = join_by(replicate, donor, plasmid, recipient)`
    ## Joining with `by = join_by(replicate, donor, plasmid, recipient)`

|         mean |           sd | log10mean |   log10sd |
|-------------:|-------------:|----------:|----------:|
| 4.573199e-12 | 1.013806e-12 | -11.34845 | 0.1022044 |

Estimates:

- gamma = 10^-11.35±0.1 i.e. 4.5e-12 ± 1e-12

``` r
parameters[6,] <- c("gamma",4.573199e-12,1.013806e-12)
```

#### Mu

Mu is the washout/death rate.

Populations were transferred every 24h losing 99% of their population
each transfer.

Therefore estimate mu at:

``` r
(mu <- 0.99/24)
```

    ## [1] 0.04125

``` r
parameters[7,] <- c("mu", mu, 0)
```

### Table of parameters

``` r
parameters %>% kable()
```

| parameter  | value        | sd           |
|:-----------|:-------------|:-------------|
| alpha      | 0.5412246    | 0.0553524    |
| beta       | 0.8180085    | 0.04290745   |
| beta_chrCM | 0.9690664    | 0.0698719    |
| beta_plaCM | 0.9449607    | 0.0413181    |
| K          | 5716250000   | 1590994000   |
| gamma      | 4.573199e-12 | 1.013806e-12 |
| mu         | 0.04125      | 0            |

``` r
write.table(parameters, file = c("./data/3_Parameters.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)
```

------------------------------------------------------------------------

**[Back to index.](../README.md)**
