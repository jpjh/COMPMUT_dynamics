# COMPMUT_dynamics
Ecological dynamics of plasmid compensatory mutations: data and analysis.

Associated with the preprint manuscript ***Superiority of chromosomal compared to plasmid-encoded compensatory mutations***, [doi: 10.1101/2024.01.15.575717](https://www.biorxiv.org/content/10.1101/2024.01.15.575717v1)

Rosanna C.T. Wright, A. Jamie Wood, Michael J. Bottery, Katie J. Muddiman, Steve Paterson, Ellie Harrison, Michael A. Brockhurst, James P.J. Hall.

---

### RMarkdown scripts describing figure generation, statistical analysis, numerical simulations, etc.

1. Figure generation: ([`1_ExperimentalPlots.Rmd`](./docs/1_ExperimentalPlots.md))
2. Statistical analysis: ([`2_ExperimentalAnalysis.Rmd`](./docs/2_ExperimentalAnalysis.md))
3. Parameter estimation: ([`3_ParameterEstimation.Rmd`](./docs/3_ParameterEstimation.md))
4. Numerical simulations of plaCM dynamics: ([`4_PlaCMSimulations.Rmd`](./docs/4_PlaCMSimulations.md))
5. Additional experiments on plaCM conjugation and acquisition ([`5_PlaCMExperiments.Rmd`](./docs/5_PlaCMExperiments.md))
6. Numerical simulations of CM dynamics: ([`6_CMDynamicsSimulations.Rmd`](./docs/6_CMDynamicsSimulations.md))
7. Plots of FP control experiments: ([`7_ControlPlots.Rmd`](./docs/7_ControlPlots.md))

---

#### Shiny App

1. Code for Shiny app: ([`app.R`](./shiny_app/app.R))
2. Link to Shiny app server: https://jpjh.shinyapps.io/COMPMOD_shiny/

---

#### Raw data

1. Relative fitness data: ([`1_RelativeFitnessData.csv`](./data/1_RelativeFitnessData.csv)). A comma-delimited text file in which columns indicate: an arbitrary treatment code (`ttmt`) describing the combination of test and reference strains, an arbitrary ID code (`ID`) which further describes whether the combination was subject to mercury selection (‘H’) or not, a replicate number (`rep`) coded 1-10, the concentration of mercury chloride used for the competition in µM (`mercury`), the counts of yellow-fluorescing competitors at the start (`Count.Y0`), the counts of yellow-fluorescing competitors at the end (`Count.Y24`), the count of red-fluorescing competitors at the start (`Count.R0`), the count of red-fluorescing competitors at the end (`Count.R24`), the relative fitness calculated as the ratio of Malthusian parameters (`RF_w`), the relative fitness calculated as the difference of Malthusian parameters (`RF_r`), the test strain identity (`strain`), the test strain chromosome genotype (`host`), the test strain fluorescent label (`host.label`), the test strain plasmid genotype (`plasmid`), the reference strain chromosome genotype (`reference_host`), the reference strain fluorescent label (`reference_host_label`), the reference plasmid genotype (`reference_plasmid`) and the reference strain identity (`reference`). 
2. Dynamics (fraction) data: ([`2_FractionData.csv`](./data/2_FractionData.csv)). A comma-delimited text file in which columns indicate: an arbitrary treatment code (`Treatment`) describing the experimental conditions, a replicate number (`Replicate`) coded 1-10, the transfer number (`transfer`) in days from 0 (start of the experiment) to 8, the total count (`Total_count`), the fractions of the Yellow, Red, Orange, and Blue populations (`Fraction_Yellow`, `Fraction_Red`, `Fraction_Orange`, `Fraction_Blue`, respectively), an arbitrary coding for Experiment (`Experiment`) in which ‘1’ refers to the experiments presented in Figures 3 and 4, and ‘2’ refers to the experiments presented in Figures 5-7, and the starting ratio of plasmid-free competitors (`Ratio`, provided only for experiment 2 in which Ratio was being varied).
   - For Experiment 1, treatment codes were as follows:
     - A: SBW25 ; SBW25(pQBR57∆PQBR57_0059::GFP) ; SBW25∆PFLU4242::dTomato(pQBR57)
     - AM: SBW25 ; SBW25(pQBR57∆PQBR57_0059::GFP) ; SBW25∆PFLU4242::dTomato(pQBR57) + Mercury
     - C: SBW25.GmR ; SBW25(pQBR57∆PQBR57_0059::tdTomato) ; SBW25∆PFLU4242::YFP(pQBR57)
     - CM: SBW25.GmR ; SBW25(pQBR57∆PQBR57_0059::tdTomato) ; SBW25∆PFLU4242::YFP(pQBR57) + Mercury
     - Therefore: For A and AM, chrCM = Fraction_Yellow and plaCM = Fraction_Red. For C and CM, chrCM = Fraction_Red and plaCM = Fraction_Yellow. In both cases, Fraction_Blue is SBW25 i.e. plasmid-free uncompensated competitors.
     - Treatments B, BM, D, and DM were as A, AM, C, and CM respectively, for Fractions Yellow and Red. However, here Fraction_Blue is SBW25(pQBR103).
   - For Experiment 2, treatment codes were as follows:
     - A: SBW25 ; SBW25(pQBR57∆pQBR57_0059::GFP) ; SBW25∆PFLU4242::dTomato(pQBR57)
     - B: SBW25 ; SBW25(pQBR57∆pQBR57_0059::GFP) ; SBW25∆PFLU4242(pQBR57::dTomato)
     - C: SBW25 ; SBW25(pQBR57∆pQBR57_0059::GFP) ; SBW25(pQBR57::dTomato)
3. Measured parameters: 
   - [`3_Parameters.txt`](./data/3_Parameters.txt), a tab-delimited text file indicating calculated parameters that is also represented in Table S1
   - [`4_WTGrowthCurve.csv`](./data/4_WTGrowthCurve.csv), associated growth curve data used to calculate alpha for `3_Parameters.txt` and Table S1. A comma-delimited text file in which columns indicate the well of the 96-well plate being measured (`well`), the replicate (`replicate`) a-k, the measurement cycle (`cycle`, measurements were taken every 15 minutes), and the optical density (`value`) as determined by a Tecan Infinite 200 Pro Nano.
4. Data used to calculate conjugation rates ([`5_ConjugationData.csv`](./data/5_ConjugationData.csv)). A comma-delimited text file in which columns indicate the strain being tested (`strain`, either pQBR57 [pQ57] or pQBR57∆PQBR57_0059 [delta_59]), the time of measurement in hours (`time_h`), the replicate a-f (`rep`), the media on which samples were plated (`media`, which can be used to determine whether counts are of donors, recipients, or transconjugants), the –log~10~ dilution factor of the plated sample (`dilution`) , the volume spread (`spread`), the counts of white or blue colonies (`count_white`, `count_blue`, which refer to recipients and donors respectively), notes on what the counts refer to (`notes`), and whether the data were discarded due to a more accurate measurement from the same culture (`crop`: measurements were cropped by inserting a ‘x’ in this column if more accurate CFU counts were available for this sample). 
5. Data from experiments investigating plaCM strains:

     - [`6_TransconjugantConjugation.csv`](./data/6_TransconjugantConjugation.csv), a comma-delimited text file in which columns indicate the strain being tested (`strain`, either pQBR57 or pQBR57∆PQBR570059 [pBR57_delta_59]), the sample taken and relevant timepoint and media used (`trt`, indicating whether the samples target donors, recipients, or transconjugants; `time` indicating whether the samples were taken at the start or end of the experiment; `media` indicating the media used for selective plating), the –log~10~ dilution factor of the plated sample (`dilution`), the volume spread (`spread`), the counts of white or blue colonies (`count_white`, `count_blue`). This table includes the counts for the 100-fold-excess-recipient experiment ('expt_2'), and the transconjugant-growth-rate experiment ('expt_1'), indicated in column `expt`.

     - [`7_AcquisitionCosts.csv`](./data/7_AcquisitionCosts.csv), a comma-delimited text file in which columns indicate the time of measurement in hours (`time`), the OD600 measurement (`Measurements`), the corresponding well and row and column (`Well`, `Row`, `Col`), the biological replicate (`Rep`), the strain used (`Strain`), and the treatment (`Trt`; 'T1_trt' are de novo transconjugants, 'T1old_trt' are established transconjugants, and 'D_trt' and 'R1_trt' are donors and recipient only, respectively). 

     - [`8_ControlCompetitionData.csv`](./data/8_ControlCompetitionData.csv), a comma-delimited text file with columns labelled as with `1_RelativeFitnessData.csv`.

6. Worksheet of model analysis:

   - [`RHWorkings.nb`](./data/RHWorkings.nb), a Wolfram Mathematica worksheet containing analyses from the model as described in the Supplementary Appendix.
---

Questions, comments, suggestions: please [get in touch](mailto:j.p.j.hall@liverpool.ac.uk)!
