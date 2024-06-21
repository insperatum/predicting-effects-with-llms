source("util.R")
library(furrr)
library(progressr)

run_analysis <- function(group_labels, filt1, filt2, filt_llm_1, filt_llm_2) {
  # group_labels: e.g. c('Black', 'White')
  # filt1, filt2: functions for filtering RCT data to each subgroup
  # filt_llm_1, filt_llm_2: functions for filtering LLM data to each subgroup
  
  tryCatch({
    # Load RCT results for group 1 (choosing a random reference-condition and outcome per study)
    results_1 <- load_tess_rct_results(filt = filt1, use_contrasts=T)
    results_1 <- filter_rct_results(results_1, random_outcome=T)
    # Load RCT results for group 2, for the same studies/reference0conditions/outcomes as group 1
    results_2 <- load_tess_rct_results(filt = filt2, use_contrasts=T, reference_conditions = results_1$df)
    
    # Ensure that we're looking at the exact same set of studies, conditions and outcomes for group 1 and group 2
    keep_1 <- results_1$df %>% left_join(results_2$df %>% select(study, outcome.name, condition.name) %>% mutate(keep_idx=T)) %>% with(!is.na(keep_idx))
    keep_2 <- results_2$df %>% left_join(results_1$df %>% select(study, outcome.name, condition.name) %>% mutate(keep_idx=T)) %>% with(!is.na(keep_idx))
    results_1$df <- results_1$df[keep_1,]
    results_1$V <- results_1$V[keep_1,keep_1]
    results_2$df <- results_2$df[keep_2,]
    results_2$V <- results_2$V[keep_2,keep_2]
    
    # Load LLM predictions
    df_llm_1 <- load_tess_llm_results(results_1$df, restricted=F, model='gpt-4', filt=filt_llm_1) %>% select(-model)
    df_llm_2 <- load_tess_llm_results(results_2$df, restricted=F, model='gpt-4', filt=filt_llm_2) %>% select(-model)
    
    # Calculate LLM-RCT correlation for group 1
    r_1 = cor(results_1$df$estimate.rct, df_llm_1$`estimate.gpt-4`, use="complete.obs")
    r_adj_1 <- calc_correlation(results_1$df$estimate.rct, results_1$V, df_llm_1$`estimate.gpt-4`, return_ci=F)[['estimate']]
    
    # Calculate LLM-RCT correlation for group 2
    r_2 = cor(results_2$df$estimate.rct, df_llm_2$`estimate.gpt-4`, use="complete.obs")
    r_adj_2 <- calc_correlation(results_2$df$estimate.rct, results_2$V, df_llm_2$`estimate.gpt-4`, return_ci=F)[['estimate']]
    
    # Calculate LLM-RCT correlation for difference between groups
    interaction = results_1$df$estimate.rct - results_2$df$estimate.rct
    V_interaction = results_1$V + results_2$V
    llm_interaction = df_llm_1$`estimate.gpt-4` - df_llm_2$`estimate.gpt-4`
    r_interaction = cor(interaction, llm_interaction, use="complete.obs")
    r_adj_interaction = calc_correlation(interaction, V_interaction, llm_interaction, return_ci=F)[['estimate']]
      
    tribble(
      ~subgroup,         ~r,            ~r_adj,
      group_labels[[1]], r_1,           r_adj_1,
      group_labels[[2]], r_2,           r_adj_2,
      "Interaction",     r_interaction, r_adj_interaction 
    )
  }, error=function(cond) {
    tibble(error = list(cond))
  })
}



#List of demographics to run the analysis with
df_groups <- tribble(
  ~demo, ~group_labels, ~filt1, ~filt_llm_1, ~filt2, ~filt_llm_2,
  "Ethnicity",
    c("White", "Black"),
    \(d) filter(d, race_4=="White"), \(x) filter(x, spec_template_group=="white"),
    \(d) filter(d, race_4=="Black"), \(x) filter(x, spec_template_group=="black"),
  "Gender",
    c("Male", "Female"),
    \(d) filter(d, GENDER=="Male"), \(x) filter(x, spec_template_group=="male"),
    \(d) filter(d, GENDER=="Female"), \(x) filter(x, spec_template_group=="female"),
  "Party",
    c("Democrat", "Republican"),
    \(d) filter(d, pid_3=="Democrat"), \(x) filter(x, spec_template_group=="democrat"),
    \(d) filter(d, pid_3=="Republican"), \(x) filter(x, spec_template_group=="republican"),
)


# Since our analysis is slightly stochastic (due to random selection of control group and outcome variable per study)
# we repeat 16 times. Approx runtime: 10 mins
n_runs = 16
handlers(handler_progress(format   = ":spin :current/:total [:bar] :percent in :elapsed ETA: :eta"))
plan(multisession, workers = 4)

df_runs <- scenarios %>% pmap_dfr(function(demo, group_labels, filt1, filt2, filt_llm_1, filt_llm_2) {
  cat("Running analyses by ", demo, "\n")
  with_progress({ p<-progressr::progressor(n_runs); future_map_dfr(1:n_runs, .options=furrr_options(seed=T), .f = function(sim_idx){
    run_analysis(group_labels, filt1, filt2, filt_llm_1, filt_llm_2) %>% mutate(demo)
  })})
})


# Calculate the median 
df_runs %>%
  group_by(demo, subgroup) %>%
  summarize(n_runs=n(), across(everything(), median))
