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
