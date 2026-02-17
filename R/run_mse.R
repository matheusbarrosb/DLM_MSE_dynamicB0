run_mse = function(n_sims, nyears, nages, init_nya,
                   waa, max_harvest_rate, rec_regime_length,
                   rec_type, sr_params, selectivity,
                   survival_mean, survival_sd, mcmc_setup,
                   burn_in_length, hist_harvest_rate, sca_model_path,
                   maturity, thresholds, vb_params, plot = FALSE,
                   estimation = FALSE) {

  sca_model_obj = NULL
  if (estimation == TRUE) {
    if (is.null(sca_model_path)) {
      stop("Estimation is TRUE but 'sca_model_path' is not provided.")
    }
    message("Compiling Stan model inside run_mse...")
    # Compile the model once before the simulation loop
    sca_model_obj = rstan::stan_model(file = sca_model_path)
  }
  
  total_time = nyears + burn_in_length

  # make containers
  abs_biomass_mat  = matrix(NA, nrow = total_time, ncol = n_sims)
  rel_biomass_mat  = matrix(NA, nrow = total_time, ncol = n_sims)
  est_biomass_mat  = matrix(NA, nrow = nyears, ncol = n_sims)
  rec_mat          = matrix(NA, nrow = total_time, ncol = n_sims)
  catch_mat        = matrix(NA, nrow = total_time, ncol = n_sims)
  mean_lengths_mat = matrix(NA, nrow = total_time, ncol = n_sims)

  # progress bar
  pb  = txtProgressBar(min = 0, max = n_sims, style = 3)
  sim = list()
  for (i in 1:n_sims) {

    setTxtProgressBar(pb, i)

    sim[[i]] =
      run_simulation(nyears            = nyears,
                     nages             = nages,
                     burn_in_length    = burn_in_length,
                     init_nya          = init_nya,
                     waa               = waa,
                     selectivity       = selectivity,
                     rec_regime_length = rec_regime_length,
                     sr_params         = sr_params,
                     rec_type          = rec_type,
                     survival_mean     = survival_mean,
                     survival_sd       = survival_sd,
                     max_harvest_rate  = max_harvest_rate,
                     hist_harvest_rate = hist_harvest_rate,
                     maturity          = maturity,
                     mcmc_setup        = mcmc_setup,
                     threshold         = thresholds,
                     vb_params         = vb_params,
                     bin_size          = 1,
                     estimation        = estimation,
                     sca_model         = sca_model_obj)

    abs_biomass_mat[,i]  = sim[[i]]$biomass
    rel_biomass_mat[,i]  = sim[[i]]$biomass / sim[[i]]$ghost_biomass
    est_biomass_mat[,i]  = sim[[i]]$est_biomass
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
  
  est_biomass_mean = rowMeans(est_biomass_mat)
  est_biomass_mean = c(rep(NA, burn_in_length), est_biomass_mean)
  est_biomass_sd   = apply(est_biomass_mat, 1, sd)
  est_biomass_sd = c(rep(NA, burn_in_length), est_biomass_sd)
  
  catch_mean       = rowMeans(catch_mat)
  catch_sd         = apply(catch_mat, 1, sd)
  recruitment_mean = rowMeans(rec_mat)
  recruitment_sd   = apply(rec_mat, 1, sd)
  mean_lenghts_mu  = rowMeans(mean_lengths_mat)
  mean_lengths_sd  = apply(mean_lengths_mat, 1, sd)

  # plotting
  if (plot == TRUE && estimation == FALSE) {

    par(mfrow = c(2,3), # 2 rows, 3 columns
        mar = c(3,3,2,1), # smaller margins for each plot
        oma = c(1,1,1,1), # outer margins
        mgp = c(2,0.7,0), # axis title, labels, line spacing
        tcl = -0.3
        )

    plot(abs_biomass_mean, type = "l", xlab = "Year", ylab = "Spawning biomass", ylim = c(0, max( (abs_biomass_mean + abs_biomass_sd)*1.2)))
    lines(abs_biomass_mean + abs_biomass_sd*1.92, lty = 2)
    lines(abs_biomass_mean - abs_biomass_sd*1.96, lty = 2)
    for(i in 1:10) {
      lines(abs_biomass_mat[,i], col = rgb(0,0,0,0.2))
    }
    abline(v = burn_in_length, lty = 2, col = "blue")
    text(x = burn_in_length/2, y = max( (abs_biomass_mean + abs_biomass_sd)*1.2)*0.8, labels = "Burn-in", col = "blue")
    text(x = burn_in_length + (total_time - burn_in_length)/2, y = max( (abs_biomass_mean + abs_biomass_sd)*1.2)*0.8, labels = "Projection", col = "blue")
    box(lwd = 2)

    plot(rel_biomass_mean, type = "l", xlab = "Year", ylab = expression(B/B[unfished]), ylim = c(0, max( (rel_biomass_mean + rel_biomass_sd)*1.2)))
    lines(rel_biomass_mean + rel_biomass_sd*1.92, lty = 2)
    lines(rel_biomass_mean - rel_biomass_sd*1.96, lty = 2)
    for(i in 1:10) {
      lines(rel_biomass_mat[,i], col = rgb(0,0,0,0.2))
    }
    abline(h = 0.5, lty = 2, col = "red")
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)

    plot(recruitment_mean, type = "l", xlab = "Year", ylab = "Recruitment (N)", ylim = c(0, max( (recruitment_mean + recruitment_sd)*1.2)))
    lines(recruitment_mean + recruitment_sd*1.92, lty = 2)
    lines(recruitment_mean - recruitment_sd*1.96, lty = 2)
    for(i in 1:5) {
      lines(rec_mat[,i], col = rgb(0,0,0,0.2))
    }
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)

    plot(catch_mean, type = "p", pch = 20, cex = 0.8,
         xlab = "Year", ylab = "Catch (t)", ylim = c(0, max( (catch_mean + catch_sd)*1.2)))
    lines(catch_mean + catch_sd*1.92, lty = 2)
    lines(catch_mean - catch_sd*1.96, lty = 2)
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)

    plot(mean_lenghts_mu, type = "l", pch = 20, cex = 1,
         xlab = "Age", ylab = "Mean length (cm)", ylim = c(0, max( (mean_lenghts_mu + mean_lengths_sd)*1.2)))
    lines(mean_lenghts_mu + mean_lengths_sd*1.92, lty = 2)
    lines(mean_lenghts_mu - mean_lengths_sd*1.96, lty = 2)
    abline(v = burn_in_length, lty = 2, col = "blue")
    box(lwd = 2)

  }

  output =
    list(
      df =
      data.frame(
        abs_biomass_mean = abs_biomass_mean,
        abs_biomass_sd   = abs_biomass_sd,
        rel_biomass_mean = rel_biomass_mean,
        rel_biomass_sd   = rel_biomass_sd,
        est_biomass_mean = est_biomass_mean,
        est_biomass_sd   = est_biomass_sd,
        catch_mean       = catch_mean,
        catch_sd         = catch_sd,
        mean_lenghts_mu  = mean_lenghts_mu,
        mean_lengths_sd  = mean_lengths_sd
      ),
      sims = sim
    )

  return(output)

}
