source("util.R")
library(furrr)
library(progressr)
library(patchwork)
library(modelr)
library(ggstance)

calc_r_adj <- function(estimate.rct, V.rct, estimate.pred) {
  ## Estimates correlation with adjustment for measurement uncertainty
  # Input:
  #   estimate.rct: vector of measured treatment effects
  #   V.rct: covariance matrix
  #   estimate.pred: vector of predicted treatment effects
  #
  # The correlation coefficient estimated using the model: (true effect, prediction) ~ N(mu, Sigma)
  
  data = bind_rows(
    tibble(name="rct", estimate=estimate.rct) %>% mutate(row_id = row_number()),
    tibble(name="pred", estimate=estimate.pred) %>% mutate(row_id = row_number()),
  ) 
  V = bdiag(V.rct, diag(0, nrow=length(estimate.rct))) * 1
  
  options(warn=-1) # Metafor warns by default that predictions have no uncertainty
  fit <- rma.mv(estimate~name, data=data, V=V, random=list(~estimator|row_id), sparse=T, struct="UN")
  options(warn=0)
  
  return(fit$rho)
}

run_analysis <- function() {
  
  # > 1. Load the human data ####
  
  # Load RCT results (choosing a random control condition per study)
  results_rct <- load_tess_rct_results()
  # Choose one outcome per study
  results_rct <- filter_rct_results(results_rct, random_outcome=T);
  df_rct <- results_rct$df
  V_rct <- results_rct$V
  
  
  # > 2. Join with predictions from forecasters and LLMs ####
  df <- df_rct
  V <- V_rct
  hypotheses <- results_rct$hypotheses
  
  models = c("gpt-4", "forecast", "combined")
  # Load predictions from each llm model
  for(model in c("gpt-4")) {
    df_model <- load_tess_llm_results(df_rct, model=model, hypotheses=hypotheses)  %>% select(-model)
    df <- df %>% left_join(df_model)
  }
  
  # Load predictions from human forecasters
  df_forecasts <- load_tess_forecast_results(df_rct, hypotheses=hypotheses)
  df <- df %>% left_join(df_forecasts)
  
  # Create a combined "forecasters + gpt-4" prediction
  df <- df %>% mutate(
    estimate.combined = (estimate.forecast + `estimate.gpt-4`)/2
  )
  
  # Remove any studies which were excluded from collecting GPT-4 predictions (due to e.g. outcome scale)
  for(model in models) {
    estimate_name = paste0('estimate.', model)
    keep_idx <- df %>% with(!is.na(get(estimate_name)))
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  
  
  # > 3. Calculate the correlations for each model
  out <- tibble()
  for(model in models) {
    estimate_name = paste0('estimate.', model)
    
    # Correlation, r
    r <- cor(df$estimate.rct, df[[estimate_name]])
    
    # Disattenuated correlation, r_adj
    # Note: calculating the CI is computationally intensive so not returning that here for ease of use. Full code is in util.R
    r_adj <- calc_r_adj(df$estimate.rct, V, df[[estimate_name]])
    
    # Proportion significant correct sign
    prop <- with(df %>% filter(p.value.rct<0.05), mean(sign(estimate.rct)==sign(get(estimate_name))))
    
    out <- bind_rows(out, tibble(model, r, r_adj, prop))
  }
  
  out
}




# Since our analysis is slightly stochastic (due to random selection of control group and outcome variable per study)
# we repeat 16 times. Approx runtime: 5 mins

n_runs = 16
plan(multisession, workers = 4)
handlers(handler_progress(format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta"))

df_runs <- with_progress({ p<-progressr::progressor(n_runs); future_map_dfr(1:n_runs, .options=furrr_options(seed=T), .f=function(sim_idx) {
  out <- run_analysis()
  p()
  out
})})


# Calculate the median 
df_runs %>%
  group_by(model) %>%
  summarize(n_runs=n(), across(everything(), median))
