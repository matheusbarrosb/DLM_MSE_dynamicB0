run_mse = function(n_sims, nyears, nages, init_nya,
                   waa, harvest_rate, rec_regime_length,
                   rec_type, sr_params,
                   maturity, threshold, vb_params, plot = FALSE) {
  
  # make containers
  abs_biomass_mat  = matrix(NA, nrow = nyears, ncol = n_sims)
  rel_biomass_mat  = matrix(NA, nrow = nyears, ncol = n_sims)
  rec_mat          = matrix(NA, nrow = nyears, ncol = n_sims)
  catch_mat        = matrix(NA, nrow = nyears, ncol = n_sims)
  mean_lengths_mat = matrix(NA, nrow = nyears, ncol = n_sims)

  # progress bar
  pb  = txtProgressBar(min = 0, max = n_sims, style = 3)
  sim = list()
  for (i in 1:n_sims) {
    
    setTxtProgressBar(pb, i)
    
    sim[[i]] =
      run_simulation(nyears            = nyears, 
                     nages             = nages,
                     init_nya          = init_nya,
                     waa               = waa, 
                     rec_regime_length = rec_regime_length,
                     sr_params         = sr_params,
                     rec_type          = rec_type,
                     survival_mean     = 0.8, 
                     survival_sd       = 0.1,
                     harvest_rate      = harvest_rate,
                     maturity          = maturity,
                     threshold         = threshold,
                     vb_params         = vb_params,
                     bin_size          = 1)
    
    abs_biomass_mat[,i]  = sim[[i]]$biomass
    rel_biomass_mat[,i]  = sim[[i]]$biomass / sim[[i]]$ghost_biomass
    catch_mat[,i]        = sim[[i]]$catch
    rec_mat[,i]          = sim[[i]]$nya_mat[,1]
    mean_lengths_mat[,i] = format_lya_mat(sim[[i]]$length_mat, vb_params)$mean_lengths
    
  }
  close(pb)
  
  # summarize outputs
  abs_biomass_mean = rowMeans(abs_biomass_mat)
  abs_biomass_sd   = apply(abs_biomass_mat, 1, sd)
  rel_biomass_mean = rowMeans(rel_biomass_mat)
  rel_biomass_sd   = apply(rel_biomass_mat, 1, sd)
  catch_mean       = rowMeans(catch_mat)
  catch_sd         = apply(catch_mat, 1, sd)
  recruitment_mean = rowMeans(rec_mat)
  recruitment_sd   = apply(rec_mat, 1, sd)
  mean_lenghts_mu  = rowMeans(mean_lengths_mat)
  mean_lengths_sd  = apply(mean_lengths_mat, 1, sd)

  # plotting
  if (plot == TRUE) {
    

    par(mfrow = c(2,3), # 2 rows, 3 columns 
        mar = c(3,3,2,1), # smaller margins for each plot 
        oma = c(1,1,1,1), # outer margins 
        mgp = c(2,0.7,0), # axis title, labels, line spacing 
        tcl = -0.3) 
    
    plot(abs_biomass_mean, type = "l", xlab = "Year", ylab = "Spawning biomass", ylim = c(0, max( (abs_biomass_mean + abs_biomass_sd)*1.2)))
    lines(abs_biomass_mean + abs_biomass_sd*1.92, lty = 2)
    lines(abs_biomass_mean - abs_biomass_sd*1.96, lty = 2)
    for(i in 1:10) {
      lines(abs_biomass_mat[,i], col = rgb(0,0,0,0.2))
    }
    
    plot(rel_biomass_mean, type = "l", xlab = "Year", ylab = expression(B/B[unfished]), ylim = c(0, max( (rel_biomass_mean + rel_biomass_sd)*1.2)))
    lines(rel_biomass_mean + rel_biomass_sd*1.92, lty = 2)
    lines(rel_biomass_mean - rel_biomass_sd*1.96, lty = 2)
    for(i in 1:10) {
      lines(rel_biomass_mat[,i], col = rgb(0,0,0,0.2))
    }
    abline(h = 0.5, lty = 2, col = "red")
    
    plot(recruitment_mean, type = "l", xlab = "Year", ylab = "Recruitment (N)", ylim = c(0, max( (recruitment_mean + recruitment_sd)*1.2)))
    lines(recruitment_mean + recruitment_sd*1.92, lty = 2)
    lines(recruitment_mean - recruitment_sd*1.96, lty = 2)
    for(i in 1:10) {
      lines(rec_mat[,i], col = rgb(0,0,0,0.2))
    }
    
    plot(catch_mean, type = "p", pch = 20, cex = 0.8,
         xlab = "Year", ylab = "Catch (t)", ylim = c(0, max( (catch_mean + catch_sd)*1.2)))
    arrows(1:nyears, catch_mean - catch_sd*1.96,
           1:nyears, catch_mean + catch_sd*1.96,
           angle = 90, code = 3, length = 0, col = rgb(0,0,0,0.5))

    plot(mean_lenghts_mu, type = "l", pch = 20, cex = 0.8,
         xlab = "Age", ylab = "Mean length (cm)", ylim = c(0, max( (mean_lenghts_mu + mean_lengths_sd)*1.2)))

  }
  
  output = 
    list(
      df = 
      data.frame(
        abs_biomass_mean = abs_biomass_mean,
        abs_biomass_sd   = abs_biomass_sd,
        rel_biomass_mean = rel_biomass_mean,
        rel_biomass_sd   = rel_biomass_sd,
        catch_mean       = catch_mean,
        catch_sd         = catch_sd,
        mean_lenghts_mu  = mean_lenghts_mu,
        mean_lengths_sd  = mean_lengths_sd
      ),
      sims = sim
    )

  
  return(output)
  
}





