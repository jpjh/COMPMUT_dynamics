######################################################################################
### Simulation and Shiny Application of COMPMOD (Compensatory Mutations Modelling) ###
######################################################################################


library(shiny)
library(deSolve)
library(tidyverse)
library(ggplot2)

ui <- fluidPage(

  headerPanel("COMPMOD: Plasmid Compensatory Mutations Modelling"),
  
  wellPanel(
    fluidRow(
      tags$body("Interactive models to investigate the fates of plasmid and chromosomally-encoded plasmid compensatory mutations.")
      ),
    fluidRow(
      tags$body("Notes:"),
      tags$div(
        tags$ul(
          tags$li("To set up a model that simulates discrete (batch) trasfers, set 'batch transfer dilution factor' to <1 (usually 0.01) and 'dilution rate' to zero."),
          tags$li("To set up a 'chrCM' type compensation, set the relative fitness of the corresponding transconjugant to that of the 'wild-type plasmid bearer'"),
          tags$li("To set up a 'plaCM' type compensation, set the relative fitness of the corresponding transconjugant to the same as the CM"),
          tags$li("To simulate a 'head-to-head' CM competition, set 'proportion with plasmid starting with CM' to 1 and 'starting proportion of CM 2' to 0.5")
        )
      )
      
      )),
  
  # Sidebar with sliders that enable the user to define parameter values
  fluidRow(
    column(4, wellPanel(
      h3("Simulation time"),
      sliderInput("transferTime", "time between transfers (hours):", 
                  min=1, max=192, value=24, step=1, ticks=FALSE),
      sliderInput("numTransfers", "number of transfers:", 
                  min=1, max=192, value=16, step=1, ticks=FALSE),
      numericInput("dilFac", "batch transfer dilution factor:", 
                  min=0, max=1, value=1, step=0.000001),
      numericInput("mu", "dilution rate:", 
                  min=0, max=1, value=0.04125, step=.01),
      h3("Densities and starting proportions"),
      numericInput("K", "maximal attainable cell density (K)",
                   min=1, max=1e11, value=5e9, step=1),
      p(textOutput("K_scinot")),
      numericInput("y0", "initial cell density",
                  min=0, max=1e11, value=5e7, step=1),
      p(textOutput("y0_scinot")),
      sliderInput("pcomp", "proportion with plasmid starting with CM", 
                  min=0, max=1, value=.5, step=0.05),
      sliderInput("pcm2", "starting proportion of CM 2", 
                  min=0, max=1, value=0, step=0.05)
    )),
    
    wellPanel(
      h3("Parameterization"),
      fluidRow(
        column(4,
               h4("Growth Rates"),
               sliderInput("alpha_Z_f", "maximum plasmid-free growth rate", 
                           min=0, max=3, value=0.5, step=.01), 
               
               sliderInput("beta_X_0", "relative fitness of wild-type plasmid-bearer", 
                           min=0, max=1.2, value=0.82, step=.01), 
               
               sliderInput("beta_X_1", "relative fitness of CM 1", 
                           min=0, max=1.2, value=0.94, step=.01),
               
               sliderInput("beta_X_1t", "relative fitness of CM 1 transconjugant", 
                           min=0, max=1.2, value=0.94, step=.01),
               
               sliderInput("beta_X_2", "relative fitness of CM 2", 
                           min=0, max=1.2, value=0.97, step=.01),

               sliderInput("beta_X_2t", "relative fitness of CM 2 transconjugant", 
                           min=0, max=1.2, value=0.82, step=.01),
        ),
        column(4,
               h4("Conjugation rates"),
               numericInput("gamma_X_0", "wild-type donor conjugation rate",
                            min = 0, value = 5e-12, step=0),
               p(textOutput("gamma_X_0_scinot")),
               
               numericInput("gamma_X_1", "CM 1 donor conjugation rate",
                            min = 0, value = 5e-12, step=.0),
               p(textOutput("gamma_X_1_scinot")),
               
               numericInput("gamma_X_1t", "CM 1 transconjugant conjugation rate",
                            min = 0, value = 5e-12, step=.0),
               p(textOutput("gamma_X_1t_scinot")),
               
               numericInput("gamma_X_2", "CM 2 donor conjugation rate",
                            min = 0, value = 5e-12, step=.0),
               p(textOutput("gamma_X_2_scinot")),

               numericInput("gamma_X_2t", "CM 2 transconjugant conjugation rate",
                            min = 0, value = 5e-12, step=.0),
               p(textOutput("gamma_X_2t_scinot"))
        ),
        column(4,
               h4("Selection"),
               sliderInput("eta_Z_f", "effect of selection on plasmid-free", 
                           min=0, max=1, value=0, step=.01)
        )
      )
    )
  ),
  
  # Show the simulation results
  mainPanel(
    plotOutput("graph", width = 1600, height=540)),
  
  
  # Show copyright
  wellPanel(
    fluidRow(
      column(12,
             tags$body("Created by "),
             tags$a("James P. J. Hall", 
                    href="http://www.jpjhall.net")
      ),
      column(12,
             tags$body("Inspired by "),
             tags$a("Zwanzig et al. 2019 doi:10.1128/mSystems.00186-18",
                    href="http://dx.doi.org/10.1128/mSystems.00186-18"))))
)

server <- function(input, output) {
  
  simulation <- reactive({
    
    # Collect parameters from sliders
    
    transferTime <- input$transferTime
    numTransfers <- input$numTransfers
    dilFac <- input$dilFac
    mu <- input$mu
    K <- input$K
    y0 <- input$y0
    pcomp <- input$pcomp
    pcm2 <- input$pcm2
    alpha_Z_f <- input$alpha_Z_f
    beta_X_0 <- input$beta_X_0
    beta_X_1 <- input$beta_X_1
    beta_X_2 <- input$beta_X_2
    beta_X_1t <- input$beta_X_1t
    beta_X_2t <- input$beta_X_2t
    gamma_X_0 <- input$gamma_X_0
    gamma_X_0t <- input$gamma_X_0
    gamma_X_1 <- input$gamma_X_1
    gamma_X_2 <- input$gamma_X_2
    gamma_X_1t <- input$gamma_X_1t
    gamma_X_2t <- input$gamma_X_2t
    eta_Z_f <- input$eta_Z_f
    
    # produce parameter set
    parms <- c(transferTime = transferTime,
               numTransfers = numTransfers,
               dilFac = dilFac,
               mu = mu,
               K = K,
               y0 = y0,
               pcomp = pcomp,
               pcm2 = pcm2,
               alpha_Z_f = alpha_Z_f,
               alpha_X_0 = alpha_Z_f * beta_X_0,
               alpha_X_1 = alpha_Z_f * beta_X_1,
               alpha_X_2 = alpha_Z_f * beta_X_2,
               alpha_X_0t = alpha_Z_f * beta_X_0,
               alpha_X_1t = alpha_Z_f * beta_X_1t,
               alpha_X_2t = alpha_Z_f * beta_X_2t,
               gamma_X_0 =  gamma_X_0,
               gamma_X_1 =  gamma_X_1,
               gamma_X_2 =  gamma_X_2,
               gamma_X_1t = gamma_X_1t,
               gamma_X_2t = gamma_X_2t,
               eta_Z_f = eta_Z_f,
               eta_X_0 = 0,
               eta_X_1 = 0,
               eta_X_2 = 0,
               eta_X_0t = 0,
               eta_X_1t = 0,
               eta_X_2t = 0)
    
    # defining initial state
    
    initial <- c(Z_f = 1, X_0 = (1-pcomp), X_1 = pcomp * (1-pcm2), X_2 = pcomp*pcm2,
                 X_0t = 0, X_1t = 0, X_2t = 0)
    
    initial_a <- c(0.9, 0.1, 0.1, 0.1, 0, 0, 0)
    initial_b <- c(0.5, 0.5, 0.5, 0.5, 0, 0, 0)
    initial_c <- c(0.1, 0.9, 0.9, 0.9, 0, 0, 0)
    initial_d <- c(0, 0.5, 0.5, 0.5, 0, 0, 0)
    
    y0_a <- as.numeric(y0) * initial * initial_a
    y0_b <- as.numeric(y0) * initial * initial_b
    y0_c <- as.numeric(y0) * initial * initial_c
    y0_d <- as.numeric(y0) * initial * initial_d
    
    # define simulation period and timings
    times <- seq(0, transferTime*numTransfers, by = 0.1) 
    
    ### Set up list of 'transfer events'
    
    event <- data.frame(var = rep(c("Z_f", "X_0", "X_1", "X_2", "X_0t", "X_1t", "X_2t"), numTransfers),
                        time = rep(seq(transferTime,transferTime*numTransfers,by=transferTime),each=7), 
                        value = c(dilFac),
                        method = c("mult"))

    ### Define ODEs
    
    COMPMOD_model<- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        # Subpopulations:
        # Z_f = no plasmid
        # X_0 = plasmid-carrier type 0 (e.g. no CM)
        # X_1 = plasmid-carrier type 1 (e.g. plaCM)
        # X_2 = plasmid-carrier type 2 (e.g. chrCM)
        # X_0t = transconjugant from X_0
        # X_1t = transconjugant from X_1
        # X_2t = transconjugant from X_2
        
        sigma = 1-((Z_f + X_0 + X_1 + X_2 + X_0t + X_1t + X_2t)/K)
        
        #         GROWTH AND DEATH      CONJUGATION               SELECTION     
        
        dZ_f <- (alpha_Z_f * Z_f) * (sigma) - 
          (mu * Z_f) - 
          (gamma_X_0 * Z_f * X_0) - 
          (gamma_X_1 * Z_f * X_1) - 
          (gamma_X_2 * Z_f * X_2) - 
          (gamma_X_0t * Z_f * X_0t) -
          (gamma_X_1t * Z_f * X_1t) -
          (gamma_X_2t * Z_f * X_2t) - (eta_Z_f * Z_f) 
        
        dX_0 <- (alpha_X_0 * X_0) * (sigma) - 
          (mu * X_0) -                                    (eta_X_0 * X_0)
        
        dX_1 <- (alpha_X_1 * X_1) * (sigma) - 
          (mu * X_1) -                                    (eta_X_1 * X_1) 
        
        dX_2 <- (alpha_X_2 * X_2) * (sigma) - 
          (mu * X_2) -                                    (eta_X_2 * X_2) 
        
        dX_0t <- (alpha_X_0t * X_0t) * (sigma) - 
          (mu * X_0t) + 
          (gamma_X_0  * Z_f * X_0) + 
          (gamma_X_0t * Z_f * X_0t) - 
          (eta_X_0t * X_0t)
        
        dX_1t <- (alpha_X_1t * X_1t) * (sigma) - 
          (mu * X_1t) + 
          (gamma_X_1 * Z_f * X_1) + 
          (gamma_X_1t * Z_f * X_1t) -
          (eta_X_1t * X_1t) 
        
        dX_2t <- (alpha_X_2t * X_2t) * (sigma) - 
          (mu * X_2t) + 
          (gamma_X_2 * Z_f * X_2) + 
          (gamma_X_2t * Z_f * X_2t) - 
          (eta_X_2t * X_2t) 
        
        list(c(dZ_f, dX_0, dX_1, dX_2, dX_0t, dX_1t, dX_2t))
      })
    }
    
    ### Run simulations
    
    mod_a <- data.frame(ode(y0_a, times, COMPMOD_model, parms, events=list(data=event))) %>%
        pivot_longer(cols=-time, names_to = "subpop", values_to = "density") %>%
        mutate(subpop = factor(subpop, levels=c('X_0','X_0t','X_1','X_1t','X_2','X_2t','Z_f'))) %>%
        group_by(time, subpop) %>% summarise(tot = sum(density), .groups='drop_last') %>%
        mutate(fraction = tot/sum(tot),
               subpanel = "10:1 plasmid-free")
    
    mod_b <- data.frame(ode(y0_b, times, COMPMOD_model, parms, events=list(data=event))) %>%
      pivot_longer(cols=-time, names_to = "subpop", values_to = "density") %>%
      mutate(subpop = factor(subpop, levels=c('X_0','X_0t','X_1','X_1t','X_2','X_2t','Z_f'))) %>%
      group_by(time, subpop) %>% summarise(tot = sum(density), .groups='drop_last') %>%
      mutate(fraction = tot/sum(tot),
             subpanel = "1:1 plasmid-free")
    
    mod_c <- data.frame(ode(y0_c, times, COMPMOD_model, parms, events=list(data=event))) %>%
      pivot_longer(cols=-time, names_to = "subpop", values_to = "density") %>%
      mutate(subpop = factor(subpop, levels=c('X_0','X_0t','X_1','X_1t','X_2','X_2t','Z_f'))) %>%
      group_by(time, subpop) %>% summarise(tot = sum(density), .groups='drop_last') %>%
      mutate(fraction = tot/sum(tot),
             subpanel = "1:10 plasmid-free")
    
    mod_d <- data.frame(ode(y0_d, times, COMPMOD_model, parms, events=list(data=event))) %>%
      pivot_longer(cols=-time, names_to = "subpop", values_to = "density") %>%
      mutate(subpop = factor(subpop, levels=c('X_0','X_0t','X_1','X_1t','X_2','X_2t','Z_f'))) %>%
      group_by(time, subpop) %>% summarise(tot = sum(density), .groups='drop_last') %>%
      mutate(fraction = tot/sum(tot),
             subpanel = "no plasmid-free")
    
    mod_output <- rbind(mod_a, mod_b, mod_c, mod_d) %>%
      mutate(subpanel = factor(subpanel, levels=c("10:1 plasmid-free",
                                            "1:1 plasmid-free",
                                            "1:10 plasmid-free",
                                            "no plasmid-free")),
             subpop = factor(subpop, levels=c("X_0","X_0t","X_1","X_1t",
                                              "X_2","X_2t","Z_f")))
    
    return(mod_output)
  })
  
  output$graph <- renderPlot({
    
    simres <- simulation()
    
    transferTime <- input$transferTime
    numTransfers <- input$numTransfers
    dilFac <- input$dilFac
    
    makePlot <- function(x, title.name){
      x %>% filter(time %in% seq(0,transferTime*numTransfers,transferTime)) %>%
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
        # scale_x_continuous(breaks = seq(0,24*transfers,24), labels = seq(0,transfers,1)) +
        facet_grid(.~subpanel) +
        theme_bw(18) + 
        theme(legend.position="bottom")
    }
    
    makePlot(simres, "TEST")
    
  }, height = 480, width=960)
   
  output$K_scinot <- renderText({
    formatC(input$K, format = "e", digits = 1)
  })

  output$y0_scinot <- renderText({
    formatC(input$y0, format = "e", digits = 1)
  })

  output$gamma_X_0_scinot <- renderText({
    formatC(input$gamma_X_0, format = "e", digits = 1)
  })
  
  output$gamma_X_1_scinot <- renderText({
    formatC(input$gamma_X_1, format = "e", digits = 1)
  })
  
  output$gamma_X_2_scinot <- renderText({
    formatC(input$gamma_X_2, format = "e", digits = 1)
  })
  
  output$gamma_X_1t_scinot <- renderText({
    formatC(input$gamma_X_1t, format = "e", digits = 1)
  })
  
  output$gamma_X_2t_scinot <- renderText({
    formatC(input$gamma_X_2t, format = "e", digits = 1)
  })
}

### Create the Shiny app

shinyApp(ui, server)