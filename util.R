library(tidyverse)
library(metafor)
library(lme4)
library(rjson)
library(scales)
library(wkb)
library(systemfit)

cached_files <- list()
cached_read_csv <- function(x) {
  if(! (x %in% names(cached_files))) {
    cached_files[[x]] <- read_csv(x)
  }
  return(cached_files[[x]]) 
}
cached_read_rds <- function(x) {
  if(! (x %in% names(cached_files))) {
    cached_files[[x]] <- read_rds(x)
  }
  return(cached_files[[x]]) 
}

calc_correlation <- function(estimate.rct, V.rct, estimate.pred, V.pred=NULL, moderators=list(), return_ci=F, return_profile=F, return_proportion_vs_control=F, return_proportion_between_treatments=F) {
  data <- tibble(
    est1 = estimate.rct,
    est2 = estimate.pred,
    !!!moderators
  ) %>%
    mutate(row_id = row_number())
  
  data <- bind_rows(
    data %>% mutate(name="est1", value=est1),
    data %>% mutate(name="est2", value=est2),
  )
  
  form <- "value ~ 0 + name"
  if(length(moderators)>0) {
    form <- paste0(c(form, paste0("name*", names(moderators))), collapse=" + ")
  }
  
  if(is.null(V.pred)) {
    V.pred <- diag(0, nrow=length(estimate.rct))
  }
  
  V = bdiag(
    V.rct,
    V.pred
  ) * 1
  V = (V + t(V))/2
  
  options(warn=-1) # Mute metafor warning about non-positive sampling variances for predictions
  fit <- rma.mv(
    data = data,
    as.formula(form),
    random=list(~name|row_id),
    V = bdiag(
      V.rct,
      V.pred
    ) * 1,
    sparse=T,
    struct="UN"
  )
  options(warn=0)
  
  fit
  
  if(return_ci) {
    ci <- confint(fit, rho=1, verbose=T) # Takes 5 minutes...
    out <- ci$random["rho",]
    
    if(return_profile) {
      prof <- profile(fit, rho=1, parallel="multicore", ncpus=8, steps=201, xlim=c(-1, 1))
      out$profile <- list(tibble(rho=prof$rho, ll=prof$ll))
    }
  } else{
    out <- tibble(estimate=fit$rho, ci.lb=NA, ci.ub=NA)
  }
  
  if(return_proportion_vs_control) {
    keep_idx = estimate.rct/sqrt(diag(V.rct)) > 2
    prop = mean(sign(estimate.rct[keep_idx]) == sign(estimate.pred[keep_idx]), na.rm=T)
    out$prop = prop
  }
  
  if(return_proportion_between_treatments) {
    stopifnot(return_proportion_vs_control==F)
    diff.pred = outer(estimate.pred, estimate.pred, "-")
    diff.rct = outer(estimate.rct, estimate.rct, "-")
    V_diff = outer(diag(V.rct), diag(V.rct), "+") - 2*V.rct # Variance on difference
    
    sig = abs(diff.rct) / sqrt(V_diff) > 2 # Which pairwise differences are significant
    sig = (!is.na(sig) & sig)
    correct = sign(diff.pred) == sign(diff.rct)
    
    out$prop = sum(sig*correct) / sum(sig)
  }

  out$cor_raw = cor(estimate.rct, estimate.pred)
  
  return(out)
}

calc_rmse <- function(estimate.rct, V.rct, estimate.pred) {
  sqrt(mean((estimate.rct - estimate.pred)^2) - mean(diag(V.rct)))
}

calc_metaregression <- function(estimate.rct, V.rct, estimate.pred) {
  data <- tibble(
    est1 = estimate.rct,
    est2 = estimate.pred,
  ) %>%
    mutate(row_id = row_number())
  
  
  fit0 <- rma.mv(
    data = data,
    est1 ~ 1,
    random=list(~1|row_id),
    V = V.rct,
    sparse=T
  )
  
  fit <- rma.mv(
    data = data,
    est1 ~ est2,
    random=list(~1|row_id),
    V = V.rct,
    sparse=T
  )
  
  
  var.rct = fit0$sigma2
  resid_var.rct = fit$sigma2
  sampling_var.rct = mean(diag(V.rct))
  b = broom::tidy(fit) %>% filter(term=="est2") %>% pull(estimate)
  b_conf.low = broom::tidy(fit, conf.int=T) %>% filter(term=="est2") %>% pull(conf.low)
  b_conf.high = broom::tidy(fit, conf.int=T) %>% filter(term=="est2") %>% pull(conf.high)
  
  out = tibble(b, var.rct,
               b_conf.low,
               b_conf.high,
               resid_var.rct, sampling_var.rct)
  return(out)
}


calc_metaregression2 <- function(estimate.rct, V.rct, estimate.pred1, estimate.pred2) {
  data <- tibble(
    est1 = estimate.rct,
    est2 = estimate.pred1,
    est3 = estimate.pred2
  ) %>%
    mutate(row_id = row_number())
  
  fit <- rma.mv(
    data = data,
    est1 ~ est2 + est3,
    random=list(~1|row_id),
    V = V.rct,
    sparse=T
  )
  
  # var.rct = fit0$sigma2
  resid_var.rct = fit$sigma2
  sampling_var.rct = mean(diag(V.rct))
  b1 = broom::tidy(fit) %>% filter(term=="est2") %>% pull(estimate)
  b2 = broom::tidy(fit) %>% filter(term=="est3") %>% pull(estimate)
  b1_conf.low = broom::tidy(fit, conf.int=T) %>% filter(term=="est2") %>% pull(conf.low)
  b1_conf.high = broom::tidy(fit, conf.int=T) %>% filter(term=="est2") %>% pull(conf.high)
  b2_conf.low = broom::tidy(fit, conf.int=T) %>% filter(term=="est3") %>% pull(conf.low)
  b2_conf.high = broom::tidy(fit, conf.int=T) %>% filter(term=="est3") %>% pull(conf.high)
  
  out = tibble(b1, b2, b1_conf.low, b1_conf.high, b2_conf.low, b2_conf.high, #var.rct,
               resid_var.rct, sampling_var.rct)
  return(out)
}


compare_with_reference <- function(df_pred, df_rct) {
  df_pred <- df_rct %>%
    select(study, outcome.name, condition.name, reference_condition) %>%
    left_join(df_pred %>% rename(reference_condition = condition.name, estimate.reference = estimate)) %>%
    left_join(df_pred) %>%
    mutate(estimate = estimate - estimate.reference) %>%
    select(-estimate.reference)
  
  return(df_pred)
}


load_tess_forecast_results <- function(df_rct, individual=F, hypotheses=NULL) {
  df_forecasting <- cached_read_rds("data/forecasting_responses.RDS")
  
  if(individual) {
    return(
      df_forecasting %>%
        rename(estimate=value) %>%
        compare_with_reference(df_rct) %>%
        rename(estimate.forecast_individual = estimate)
    )  
  }
  
  if(!is.null(hypotheses)) {
    df_forecasting <- df_forecasting %>% left_join(hypotheses) %>%
      filter(!is.na(hypothesis)) %>%
      mutate(condition.name = paste0(hypothesis, "_", t_hypothesis))
  }
  
  ## combine
  df_forecasting_estimates <- df_forecasting %>%
    group_by(PROLIFIC_PID, outcome.name) %>%
    mutate(value = (value - mean(value, na.rm=T))/100) %>%
    mutate(x = str_glue("{study}__{outcome.name}__{condition.name}")) %>%
    lm(value ~ 0 + x, .) %>%
    broom::tidy(conf.int=T) %>%
    transmute(
      study = str_match(term, "x(.*)__(.*)__(.*)")[,2],
      outcome.name = str_match(term, "x(.*)__(.*)__(.*)")[,3],
      condition.name = str_match(term, "x(.*)__(.*)__(.*)")[,4],
      estimate,
      #conf.low,
      #conf.high
    )
  
  # Compare with reference condition
  df_forecasting_estimates <- df_forecasting_estimates %>%
    compare_with_reference(df_rct) %>%
    rename(estimate.forecast = estimate)
  return(df_forecasting_estimates)
}


load_tess_llm_results <- function(df_rct, restricted=F, model=NA, n_prompts=NA, keep_spec_group=F, filt=NULL, hypotheses=NULL, use_contrasts=TRUE) {
  df_llm <- cached_read_rds("data/llm_responses.RDS")
  if(!is.na(model)) {
    df_llm <- df_llm %>% filter(model == .env$model)
  }
  
  if(!is.null(hypotheses)) {
    df_llm <- df_llm %>% left_join(hypotheses) %>%
      filter(!is.na(hypothesis)) %>%
      mutate(condition.name = paste0(hypothesis, "_", t_hypothesis))
  }
  
  if(!is.null(filt)) df_llm <- filt(df_llm)
  
  df_llm_predictions <- df_llm %>%
    mutate(expectation = ifelse(scale_flip, outcome_scale_max - (expectation-outcome_scale_min), expectation)) %>%
    mutate(expectation_rescaled = (expectation - outcome_scale_min)/(outcome_scale_max- outcome_scale_min)) %>%
    group_by(model, study, outcome.name, spec_group, spec_template_group, weight, scale.direction.intuitive, scale_flip) %>%
    mutate(expectation_rel = expectation_rescaled - mean(expectation_rescaled)) %>%
    ungroup() 
  
  df_llm_predictions <- df_llm_predictions %>%
    filter((scale_flip==0 & scale.direction.intuitive %in% c(NA, "-", 1)) | (scale_flip==1 & scale.direction.intuitive== 0))
  
  if(!is.na(n_prompts)) {
    keep_spec_groups <- df_llm_predictions %>%
      select(spec_group) %>%
      distinct() %>%
      slice_sample(n=n_prompts) 
    df_llm_predictions <- df_llm_predictions %>% right_join(keep_spec_groups)
  }
  
  if(keep_spec_group) {
    df_ates_llm <- df_llm_predictions %>%
      group_by(model, study, outcome.name, condition.name, spec_group) %>%
      summarize(estimate = mean(expectation_rel)) %>%
      ungroup()
  } else {
    df_ates_llm <- df_llm_predictions %>%
      group_by(model, study, outcome.name) %>%
      group_modify(function(d,k) {
        d %>%
          lm(expectation_rel ~ 0 + condition.name, weights = weight, .) %>%
          broom::tidy(conf.int=T) %>%
          mutate(condition.name = str_match(term, "condition.name(.*)")[,2], .before=1) %>% select(-term)
      }) %>%
      select(model, study, outcome.name, condition.name, estimate)
  }
  
  # Compare with reference condition
  if(use_contrasts) df_ates_llm <- df_ates_llm %>%
      compare_with_reference(df_rct)
  
  if(is.na(model)) {
    df_ates_llm <- df_ates_llm %>%
      rename(estimate.llm = estimate)
  } else {
    name = paste0("estimate.", model)
    df_ates_llm <- df_ates_llm %>%
      rename({{name}} := estimate)
  }
  
  return(df_ates_llm)
}

get_crossval_estimates <- function(df_rct, df_llm_by_spec) {
  df_augmented <- df_rct %>% left_join(df_llm_by_spec)
  study_folds <- crossv_kfold(df %>% select(study) %>% distinct, k=5)
  df_crossval <- study_folds %>% pmap_dfr(function(train, test, .id) {
    keep_spec_groups <- df_augmented %>%
      right_join(train$data[train$idx,]) %>%
      group_by(spec_group, study) %>%
      summarize(score=cor(estimate.rct, estimate.gpt4prompt)) %>%
      summarize(score=mean(score, na.rm=T)) %>%
      slice_max(score, prop=0.5) %>%
      pull(spec_group)
    df_llm_by_spec %>%
      right_join(test$data[test$idx,]) %>%
      filter(spec_group %in% keep_spec_groups) %>%
      group_by(study, outcome.name, condition.name, reference_condition) %>%
      summarize(estimate.crossval = mean(estimate.gpt4prompt)) %>%
      ungroup()
  })
    
}

load_tess_hypotheses <- function() {
  df_hypotheses <- cached_read_csv("data/hypotheses.csv")%>%
    pivot_longer(c(hyp1, hyp2, hyp3, hyp4), names_to="hypothesis", values_to="t_hypothesis") %>%
    filter(!is.na(t_hypothesis)) %>%
    group_by(study, outcome.name, hypothesis) %>%
    group_modify(function(d,k){tibble(hypothesis_data=list(d))})
  
  foo <- df_hypotheses %>%
    group_by(study) %>%
    slice_sample(n=1) %>%
    group_by(study, outcome.name, hypothesis) %>%
    reframe(bind_rows(hypothesis_data)) %>%
    ungroup()
  
  return(foo)
}

load_tess_rct_results <- function(use_contrasts = T, filt=NULL, reference_conditions=NULL, use_hypotheses=F) {
  df_responses <- cached_read_rds("data/rct_responses.RDS")
  hypotheses <- if(use_hypotheses) load_tess_hypotheses() else NULL
  
  df_regression_by_study <- df_responses %>%
    rowwise() %>%
    mutate(data = list(data %>% ungroup() %>% mutate(
      participant_id = as.character(row_number()),
      outcome.name,
      y = rescale(as.numeric(y), to=c(0,1), from=c(outcome.min, outcome.max)),
      y = y - mean(y, na.rm=T)))) %>%
    group_by(study) %>%
    summarise(data=list(bind_rows(data))) %>%
    ungroup() %>%
    pmap_dfr(function(study, data, ...) {
      #cat(study)
      # making the data long and then wide again such that 1 column = 1 outcome
      d_long <- data %>%
        # inner_join(df_llm %>% filter(study ==.env$study) %>% select(outcome.name, condition.name) %>% distinct, by =c("outcome.name", "condition.name")) %>%
        group_by(participant_id) %>%
        filter(!any(is.na(y))) %>%
        ungroup()
      if(nrow(d_long)==0) {
        cat(study)
        print(table(data$outcome.name, data$condition.name))
        
        ## TODO: This happens for HamiltonS31 because outcomes are measued on distinct participants.
        ## In this case, we should just revert to doing independent regressions
        ## But dropping for now to simplify
        return(tibble())
      }
      
      ###
      if(!is.null(hypotheses)) {
        foo <- hypotheses %>% filter(study==.env$study)
        if(nrow(foo) == 0) {
          return(tibble())
        }
        d_long <- d_long %>% left_join(foo) %>% filter(!is.na(t_hypothesis)) %>% mutate(condition.name = paste0(hypothesis, "_", t_hypothesis)) %>% select(-t_hypothesis)
      }
      
      d_wide <- d_long %>%
        pivot_wider(names_from="outcome.name", values_from="y")
      
      outcomes <- unique(d_long$outcome.name)
      
      if(use_contrasts) {
        if(is.null(reference_conditions)) {
          reference_condition = sample(d_wide$condition.name, 1)  
        } else {
          reference_condition = reference_conditions %>% filter(study==.env$study) %>% pull(reference_condition) %>% unique
        }
        
        d_wide <- d_wide %>% mutate(condition.name = fct_relevel(condition.name, reference_condition))
        equations <- outcomes %>% map(~as.formula(str_glue("`{.}` ~ condition.name")))  
      } else{
        reference_condition = NA
        equations <- outcomes %>% map(~as.formula(str_glue("`{.}` ~ 0 + condition.name")))  
      }
      
      
      if(!is.null(filt)) {
        tryCatch({
            d_wide <- filt(d_wide)
            if(nrow(d_wide)==0) return(tibble())    
          },
          error = function(e) {
            return(tibble())
          }
        )
        
      }
      
      fit <- systemfit(equations, data=d_wide, method = "SUR")
      
      term_to_outcome <- function(term) {outcomes[as.integer(str_match(term, "eq([0-9])*_")[,2])]}
      term_to_condition <- function(term) {str_match(term, "condition.name(.*)")[,2]}
      df_estimates <- broom::tidy(fit) %>%
        transmute(
          study,
          outcome.name = term_to_outcome(term),
          condition.name = term_to_condition(term),
          estimate.rct = estimate,
          std.error.rct = std.error,
          p.value.rct = p.value,
          N_per_condition.rct = nrow(d_wide)/n_distinct(d_wide$condition.name),
          N_total.rct = nrow(d_wide),
          N_conditions.rct = n_distinct(d_wide$condition.name),
          N_outcomes.rct = length(outcomes),
          N_effects.rct = (N_conditions.rct - 1) * N_outcomes.rct
        )
      
      if(any(is.na(df_estimates$std.error.rct))) {
        cat(study, " failing!!\n") 
        return(tibble())
      } else{
        # df_estimates <- df_estimates %>% select(-std.error)
      }
      
      
      V <- vcov(fit)
      colnames(V) <- colnames(V) %>%
        {paste0(term_to_outcome(.), ": ", term_to_condition(.))}
      rownames(V) <- rownames(V) %>%
        {paste0(term_to_outcome(.), ": ", term_to_condition(.))}
      
      if(use_contrasts) {
        df_estimates <- df_estimates %>%
          mutate(reference_condition,
                 keep_idx = !is.na(condition.name)) 
        V <- V[df_estimates$keep_idx, df_estimates$keep_idx]
        df_estimates <- df_estimates %>% filter(keep_idx) %>% select(-keep_idx)
      }
      tibble(study, reference_condition, estimates=list(df_estimates), V=list(V))
    })
  
  df_regression <- bind_rows(df_regression_by_study$estimates)
  V_regression <- df_regression_by_study %>%
    pull(V) %>%
    bdiag %>%
    as.matrix
  
  list("df"=df_regression, "V"=V_regression, "hypotheses"=hypotheses)
}

filter_rct_results <- function(
    results_rct,
    predicate=NA,
    published=NA, criterion3=NA, small_effects=NA, short_text=NA, significant=NA, first_outcome=NA,
    random_outcome=F, random_conditions=NA, low_std_error=NA, large_N=NA, max_effect_size=NA, min_effect_size=NA,
    gender_related=NA, race_related=NA, partisanship_related=NA,
    gpt_recognizes=NA, field=NA) {
  
  df <- results_rct$df
  V <- results_rct$V
  hypotheses <- results_rct$hypotheses
  
  if(!is.na(predicate)) {
    keep_idx <- predicate(df)
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(gender_related)) {
    df_bentable <- cached_read_csv("data/study_coding/TESSArchive_StudyCoding - Ben_table.csv")
    is_gender_related <- df_bentable %>% mutate(keep=(gender_relatedness + gender_relatedness_tommy)==2) %>% select(study, keep)
    keep_idx <- df %>% left_join(is_gender_related) %>% with(keep == gender_related)

    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(race_related)) {
    df_bentable <- cached_read_csv("data/study_coding/TESSArchive_StudyCoding - Ben_table.csv")
    is_race_related <- df_bentable %>% mutate(keep=(race_relatedness + race_relatedness_tommy)==2) %>% select(study, keep)
    keep_idx <- df %>% left_join(is_race_related) %>% with(keep == race_related)
    
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(partisanship_related)) {
    df_bentable <- cached_read_csv("data/study_coding/TESSArchive_StudyCoding - Ben_table.csv")
    is_partisanship_related <- df_bentable %>% mutate(keep=(partisanship_relatedness + partisanship_relatedness_tommy)==2) %>% select(study, keep)
    keep_idx <- df %>% left_join(is_partisanship_related) %>% with(keep == partisanship_related)
    
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(published)) {
    df_bentable <- cached_read_csv("data/study_coding/TESSArchive_StudyCoding - Ben_table.csv")
    is_published <- df_bentable %>% mutate(keep= (published=="yes" | published=="yes, dissertation") & year <=2021) %>% select(study, keep)
    keep_idx <- df %>% left_join(is_published) %>% with(keep == published)
    
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(criterion3)) {
    df_outcomes <- cached_read_csv("processed_data/outcomes.csv")
    is_criterion3 <- df_outcomes %>% mutate(keep=fails_criterion3==0) %>% select(study, outcome.name, keep)
    keep_idx <- df %>% left_join(is_criterion3) %>% with(keep == criterion3)
    
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(small_effects)) {
    keep_idx <- (abs(df$estimate.rct) < median(abs(df$estimate.rct))) == small_effects
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(max_effect_size)) {
    keep_idx <- abs(df$estimate.rct) < max_effect_size
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(min_effect_size)) {
    keep_idx <- abs(df$estimate.rct) > min_effect_size
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(short_text)) {
    df_conditions <- cached_read_csv("processed_data/conditions.csv")
    df_study_len <- df_conditions %>% mutate(len = str_length(paste0(condition.text, blurb))) %>% group_by(study) %>% summarize(len=mean(len))
    is_short <- df_study_len %>% mutate(keep=len<median(len)) %>% select(study, keep)
    
    keep_idx <- df %>% left_join(is_short) %>% with(keep == short_text)
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(significant)) {
    keep_idx <- (df$p.value.rct < 0.05) == significant
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(first_outcome)) {
    df_outcomes <- cached_read_csv("processed_data/outcomes.csv")
    is_first_outcome <- df_outcomes %>% mutate(keep = outcome.number==1) %>% select(study, outcome.name, keep) %>% distinct
    
    keep_idx <- df %>% left_join(is_first_outcome) %>% with(replace_na(keep, T) == first_outcome)
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(random_outcome) {
    keep_outcome <- df %>% select(study, outcome.name) %>% distinct() %>% group_by(study) %>% slice_sample(n=1) %>% ungroup() %>% mutate(keep=T)
    keep_idx <- df %>% left_join(keep_outcome) %>% with(replace_na(keep, F))
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(random_conditions)) {
    keep_condition <- df %>% select(study, condition.name) %>% distinct() %>% group_by(study) %>% slice_sample(n=random_conditions) %>% ungroup() %>% mutate(keep=T)
    keep_idx <- df %>% left_join(keep_condition) %>% with(replace_na(keep, F))
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(low_std_error)) {
    keep_idx <- (df$std.error.rct < median(df$std.error.rct)) == low_std_error
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(large_N)) {
    keep_idx <- (df$N_per_condition.rct > median(df$N_per_condition.rct)) == large_N
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(gpt_recognizes)) {
    df_gpt_author_recognition <- cached_read_csv("processed_data/gpt_author_recognition.csv") %>%
      transmute(study, gpt_recognizes_authors = proportion_yes>0.5)
    keep_idx <- df %>% left_join(df_gpt_author_recognition) %>% with(gpt_recognizes_authors == gpt_recognizes)
    
    keep_idx = replace_na(keep_idx, F)
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  
  if(!is.na(field)) {
    df_fields <- cached_read_csv("data/study_coding/TESSArchive_StudyCoding - Ben_table.csv") %>%
      select(study, sociology, political_science, communication, psychology, social_policy)
    
    keep_idx <- df %>% left_join(df_fields) %>% with(get(field)==1)
    df <- df[keep_idx,]
    V <- V[keep_idx, keep_idx]
  }
  return(list("df"=df, "V"=V, "hypotheses"=hypotheses))
}
