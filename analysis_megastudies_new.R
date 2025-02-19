source("util.R")

df_studies_outcomes <- readRDS("data/megastudies_new.rds")
df_predictions_all = readRDS("data/individual_expert_predictions.rds")

d = df_studies_outcomes %>%
  left_join(df_predictions_all) %>%
  rowwise() %>% filter(!is.null(df_predictions)) %>% ungroup()


df_results = d %>% pmap_dfr(function(dataset, outcome, df, V, df_predictions, ...) {
  get_w = function(pred) {
    (pred==max(pred))/sum(pred==max(pred))
  }
  # Note: we have to use GPT to select the best treatment _for each forecaster_
  # separately, as forecasters (in some studies) only see a subset of treatments
  df_w = df_predictions %>%
    left_join(df) %>%
    group_by(forecaster_id) %>%
    filter(n()>1) %>% # If a forecaster only rated one condition, drop them
    mutate(
      w.random = 1/n(),
      w.forecaster = get_w(prediction),
      w.gpt = get_w(`prediction.gpt-4`),
    ) %>%
    group_by(condition.name) %>%
    summarize(across(c(w.forecaster, w.gpt, w.random), mean))
  
  df = df %>% left_join(df_w)

  f <- function(w) {tibble(
      estimate = sum(w * df$estimate.rct),
      se = as.numeric(sqrt(t(w) %*% V %*% w))
  )}
  
  bind_rows(
    f(df$w.gpt - df$w.random) %>% mutate(method="gpt"),
    f(df$w.forecaster - df$w.random) %>% mutate(method="expert"),
  ) %>%
    mutate(dataset, outcome)
})

df_results

df_results %>%
  pivot_wider(id_cols=c(dataset, outcome), names_from=method, values_from=c(estimate, se)) %>%
  mutate(ratio_gpt = estimate_gpt / estimate_expert)

         