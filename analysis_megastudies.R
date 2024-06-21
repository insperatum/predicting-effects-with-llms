source("util.R")

df_studies_outcomes <- readRDS("data/megastudies.rds")

results <- df_studies_outcomes %>% pmap_dfr(function(dataset, outcome, df, V, N) {
  cat(str_glue("\n\n**** Running {dataset} {replace_na(outcome, '')}... ****\n\n"))
  
  # Calculate correlation between RCT data and GPT-4
  out <- tibble(dataset, outcome,
                r_gpt = cor(df$estimate.rct, df$`prediction.gpt-4`),
                r_adj_gpt = calc_correlation(df$estimate.rct, V, df$`prediction.gpt-4`, return_ci = F)[['estimate']])
  
  # Calculate correlation between RCT and expert forecasts
  if("prediction.expert" %in% colnames(df)) {
    out <- out %>% mutate(
      r_expert = cor(df$estimate.rct, df$`prediction.expert`),
      r_adj_expert = calc_correlation(df$estimate.rct, V, df$`prediction.expert`, return_ci = F)[['estimate']],
    )
  }
  
  out
})

results