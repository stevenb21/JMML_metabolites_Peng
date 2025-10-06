library(tidyverse)
library(caret)
library(randomForest)
library(xgboost)
library(e1071)
library(pROC)
library(PRROC)
library(yardstick)

# 1) Read data once

full_data <- readRDS("../res/datamine/full_combined_matrix.rds") |> 
  select(-meta_bd_blood_collection_time_x)

# 2) Define feature‐selection scenarios
scenarios <- list(
  bc_bin      = function(df) df %>% select(matches("^(bc_|bin_)")),
  bc_only     = function(df) df %>% select(starts_with("bc_")),
  bc_bin_meta = function(df) df %>% select(matches("^(bc_|bin_|meta_)")),
  bc_meta     = function(df) df %>% select(matches("^(bc_|meta_)"))
)

# 3) Data‐prep helper
prepare_data <- function(df, select_fn) {
  # 1) pull out the outcome
  y <- factor(df$meta_sample_group, levels = c("Control", "JMML"))
  
  # 2) select your raw predictors
  raw <- select_fn(df)
  
  # 3) drop the label and sample ID if they exist
  raw <- raw %>% select(-any_of(c("meta_sample_group", "meta_SampleID")))
  
  # 4) one-hot encode factors + leave numeric alone
  X <- model.matrix(~ . -1, data = raw) %>% as.data.frame()
  
  list(X = X, y = y)
}



# helper to flip anything under 0.5
flip_half <- function(x) if_else(x < 0.5, 1 - x, x)

# 4) CV + model‐training helper
run_cv <- function(X, y, k = 5, base_seed = 123) {
  set.seed(base_seed)
  folds <- createFolds(y, k = k, returnTrain = FALSE)
  
  cv_list <- map2(folds, seq_along(folds), function(test_idx, fold_id) {
    X_train <- X[-test_idx,]; X_test  <- X[test_idx,]
    y_train <- factor(y[-test_idx], levels = levels(y))
    y_test  <- factor(y[ test_idx], levels = levels(y))
    
    # — models —
    set.seed(base_seed + fold_id)
    rf_mod  <- randomForest(x = X_train, y = y_train)
    rf_p    <- predict(rf_mod, X_test, type = "prob")[, "JMML"]
    
    dtrain  <- xgb.DMatrix(data = as.matrix(X_train), label = as.numeric(y_train == "JMML"))
    dtest   <- xgb.DMatrix(data = as.matrix(X_test),  label = as.numeric(y_test  == "JMML"))
    xgb_mod <- xgb.train(
      params   = list(objective = "binary:logistic", eval_metric = "auc"),
      data     = dtrain, nrounds = 100, verbose = 0
    )
    xgb_p    <- predict(xgb_mod, dtest)
    
    svm_mod  <- svm(x = X_train, y = y_train, probability = TRUE)
    svm_p    <- attr(predict(svm_mod, X_test, probability = TRUE), "probabilities")[, "JMML"]
    
    glm_mod  <- glm(y_train ~ ., data = as.data.frame(X_train), family = binomial)
    glm_p    <- predict(glm_mod, newdata = as.data.frame(X_test), type = "response")
    
    # — metrics & preds —
    rocs <- list(
      RandomForest = roc(y_test, rf_p,  levels = levels(y), direction = "auto"),
      XGBoost      = roc(y_test, xgb_p, levels = levels(y), direction = "auto"),
      SVM          = roc(y_test, svm_p, levels = levels(y), direction = "auto"),
      Logistic     = roc(y_test, glm_p, levels = levels(y), direction = "auto")
    )
    prs  <- list(
      RandomForest = pr.curve(scores.class0 = rf_p[y_test=="Control"],  scores.class1 = rf_p[y_test=="JMML"], curve=FALSE),
      XGBoost      = pr.curve(scores.class0 = xgb_p[y_test=="Control"], scores.class1 = xgb_p[y_test=="JMML"], curve=FALSE),
      SVM          = pr.curve(scores.class0 = svm_p[y_test=="Control"], scores.class1 = svm_p[y_test=="JMML"], curve=FALSE),
      Logistic     = pr.curve(scores.class0 = glm_p[y_test=="Control"], scores.class1 = glm_p[y_test=="JMML"], curve=FALSE)
    )
    
    models <- names(rocs)
    metrics <- map_dfr(models, function(m) {
      tibble(
        Fold     = fold_id,
        Model    = m,
        Accuracy = mean(
          if (m == "RandomForest")    predict(rf_mod, X_test) == y_test else
            if (m == "XGBoost")         ifelse(xgb_p > 0.5, "JMML", "Control") == y_test else
              if (m == "SVM")             predict(svm_mod, X_test) == y_test else
                ifelse(glm_p    > 0.5, "JMML", "Control") == y_test
        ),
        AUROC = as.numeric(rocs[[m]]$auc),
        PRAUC = prs[[m]]$auc.integral
      )
    })
    
    preds <- map_dfr(models, function(m) {
      tibble(
        Fold  = fold_id,
        Model = m,
        truth = y_test,
        prob  = switch(
          m,
          RandomForest = rf_p,
          XGBoost      = xgb_p,
          SVM          = svm_p,
          Logistic     = glm_p
        )
      )
    })
    
    list(metrics = metrics, preds = preds)
  })
  
  # combine, summarize, pooled
  cv_results <- map_dfr(cv_list, "metrics")
  cv_preds   <- map_dfr(cv_list, "preds")
  
  cv_summary <- cv_results %>%
    group_by(Model) %>%
    summarise(
      Accuracy = mean(Accuracy),
      AUROC    = mean(AUROC),
      PRAUC    = mean(PRAUC),
      .groups  = "drop"
    ) %>%
    mutate(Fold = "Mean") %>%
    select(Fold, Model, Accuracy, AUROC, PRAUC)
  
  pooled_summary <- cv_preds %>%
    group_by(Model) %>%
    summarise(
      Accuracy = mean((prob > .5 & truth == "JMML") | (prob <= .5 & truth == "Control")),
      AUROC    = as.numeric(roc(truth, prob, levels = levels(y), direction = "auto")$auc),
      PRAUC    = pr.curve(scores.class0 = prob[truth=="Control"], scores.class1 = prob[truth=="JMML"], curve = FALSE)$auc.integral,
      .groups  = "drop"
    ) %>%
    mutate(Fold = "Pooled") %>%
    select(Fold, Model, Accuracy, AUROC, PRAUC)
  
  
  # — after your existing cv_summary —  
  cv_summary <- cv_summary %>%
    mutate(
      AUROC = flip_half(AUROC)
    )
  
  # — after your existing pooled_summary —  
  pooled_summary <- pooled_summary %>%
    mutate(
      AUROC = flip_half(AUROC)
    )
  
  final_table <- bind_rows(
    cv_results %>% mutate(Fold = as.character(Fold)),
    pooled_summary
  ) %>%
    select(-PRAUC) %>% 
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    mutate(
      Fold  = factor(Fold, levels = c(as.character(1:5), "Pooled")),
      Model = factor(Model, levels = c("Logistic","SVM","RandomForest","XGBoost"))
    ) %>%
    arrange(Model, Fold)
  
  list(
    per_fold = cv_results,
    pooled   = pooled_summary,
    final    = final_table,
    preds    = cv_preds
  )
}

# 5) Run all scenarios
results <- map(scenarios, ~ {
  dat <- prepare_data(full_data, .x)
  run_cv(dat$X, dat$y)
})



fold_res <- imap_dfr(results, ~ .x$per_fold %>% mutate(Scenario = .y)) %>%
  select(Scenario, everything()) %>%
  select(-PRAUC) %>% 
  mutate(AUROC = flip_half(AUROC)) %>% 
  mutate(Model = factor(Model, levels = c("Logistic", "SVM", "RandomForest", "XGBoost")))


write_csv(fold_res, "../res/datamine/fold-summaries.csv")

# 6) Aggregate pooled summaries across scenarios
pooled_all <- imap_dfr(results, ~ .x$pooled %>% mutate(Scenario = .y)) %>%
  select(Scenario, everything()) %>% select(-PRAUC) %>% 
  mutate(Model = factor(Model, levels = c("Logistic", "SVM", "RandomForest", "XGBoost")))

# Write out a CSV of pooled results
write_csv(pooled_all, "../res/datamine/pooled_summaries_all_scenarios.csv")

# 7) Plot & save ROC for each scenario

# ──  helper  ──
fix_roc_for_plot <- function(truth, prob, levels = c("Control","JMML")) {
  r <- roc(truth, prob, levels = levels, direction = "auto")
  if (as.numeric(r$auc) < 0.5) {
    r <- roc(truth, 1 - prob, levels = levels, direction = "auto")
  }
  r
}

for (sc in names(results)) {
  preds <- results[[sc]]$preds
  
  roc_list <- preds %>%
    group_by(Model) %>%
    summarize(
      roc_obj = list(
        fix_roc_for_plot(truth, prob)
      ),
      .groups = "drop"
    ) %>%
    { set_names(.$roc_obj, .$Model) }
  
  aucs <- map_dbl(roc_list, ~ as.numeric(.$auc))
  
  # build legend labels
  legend_labels <- paste0(names(aucs), " (AUC=", sprintf("%.3f", aucs), ")")
  
  p <- ggroc(
    roc_list,
    aes     = "colour",     # ← here!
    legacy.axes = TRUE
  ) +
    geom_abline(
      slope     = 1, 
      intercept = 0,
      linetype  = "dashed",
      color     = "gray70"
    ) +
    geom_line(size = 1.2) +
    scale_colour_brewer(       # ← note scale_**colour**_*
      palette = "Set2",
      breaks  = names(aucs),
      labels  = paste0(names(aucs), " (AUC=", sprintf("%.3f", aucs), ")"),
      name    = "",
      guide   = guide_legend(
        ncol   = 2,       # 2 columns
        nrow   = 2,       # 2 rows
        byrow  = TRUE     # fill across rows first
      )
    ) +
    coord_equal() +
    theme_bw(base_size = 14) +
    labs(
      x = "False Positive Rate (1 – Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme(
      panel.grid.minor     = element_blank(),
      panel.grid.major     = element_line(color = "gray90"),
      axis.title           = element_text(face = "bold"),
      axis.text            = element_text(color = "black"),
      legend.position      = "bottom",
      legend.direction     = "horizontal",
      legend.key.width     = unit(1.5, "cm"),
      legend.background    = element_blank()
    )
  
  
  ggsave(
    filename = file.path("../res/datamine", paste0(sc, "_pooledROC.png")),
    plot     = p,
    width    = 6, height = 6,
    units    = "in",
    dpi      = 300
  )
  
}

