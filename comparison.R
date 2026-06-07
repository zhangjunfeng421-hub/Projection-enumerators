library(lhs)
library(SLHD)
library(MaxPro)
library(DiceDesign)
library(foreach)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(patchwork)

source("GSOA_opt.R")

# ============================================================
# 1. criteria
# ============================================================

euclidean_distance <- function(D) {
  n <- nrow(D)
  min_dist <- Inf
  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      dist <- sum((D[i,] - D[j,])^2)
      if(dist < min_dist) min_dist <- dist
    }
  }
  return(min_dist)
}

WD <- function(P) {
  n <- nrow(P)
  m <- ncol(P)
  
  constant <- -(4/3)^m
  
  pairs <- expand.grid(a = 1:n, b = 1:n)
  
  pair_products <- apply(pairs, 1, function(idx) {
    a <- idx[1]; b <- idx[2]
    diffs <- abs(P[a, ] - P[b, ])
    prod(3/2 - diffs + diffs^2)
  })
  
  pair_term <- (1/n^2) * sum(pair_products)
  
  return((constant + pair_term)^0.5)
}

CD <- function(D) {
  n <- nrow(D)
  m <- ncol(D)
  
  term1 <- (13/12)^m
  
  term2 <- 0
  for(a in 1:n) {
    prod_term <- prod(1 + 0.5*abs(D[a,] - 0.5) - 0.5*(D[a,] - 0.5)^2)
    term2 <- term2 + prod_term
  }
  term2 <- (2/n) * term2
  
  term3 <- 0
  for(a in 1:n) {
    for(b in 1:n) {
      prod_term <- prod(1 + 0.5*abs(D[a,] - 0.5) + 0.5*abs(D[b,] - 0.5) - 0.5*abs(D[a,] - D[b,]))
      term3 <- term3 + prod_term
    }
  }
  term3 <- (1/n^2) * term3
  
  CD <- term1 - term2 + term3
  return(CD)
}

corr_aveabs <- function(D) {
  cor_matrix <- abs(cor(D))
  rho_ave <- mean(cor_matrix[upper.tri(cor_matrix)])
  return(rho_ave)
}

# ============================================================
# 2. 投影评估函数
# ============================================================

evaluate_projections <- function(D, k, criteria_funs) {
  m <- ncol(D)
  
  combos <- combn(m, k, simplify = FALSE)
  n_combos <- length(combos)
  
  results <- list()
  for(crit_name in names(criteria_funs)) {
    results[[crit_name]] <- numeric(n_combos)
  }
  
  for(idx in 1:n_combos) {
    D_proj <- D[, combos[[idx]], drop = FALSE]
    for(crit_name in names(criteria_funs)) {
      results[[crit_name]][idx] <- criteria_funs[[crit_name]](D_proj)
    }
  }
  
  summary_stats <- list()
  for(crit_name in names(criteria_funs)) {
    if(crit_name == "Euclidean") {
      summary_stats[[crit_name]] <- min(results[[crit_name]])
    } else {
      summary_stats[[crit_name]] <- max(results[[crit_name]])
    }
  }
  
  return(summary_stats)
}

# ============================================================
# 3. main fuction
# ============================================================

compare_designs <- function(n = 27, m = 13, n_cores = NULL) {
  
  cat("生成实验设计...\n")
  
  # 生成设计
  D_gsoa <- GSOA_opt(c(1,0,2,1), 3, K = 1, add = c(1,1,1))
  D_gsoa <- (2 * D_gsoa + 1)/(2 * 27)
  
  X <- lhsDesign(n, m)$design
  D_unif <- discrepESE_LHS(X)$design
  
  D_maximin <- maximinSLHD(1, n, m)$Design
  D_maximin <- (2 * D_maximin - 1)/(2 * max(D_maximin))
  
  D_maxpro <- MaxProLHD(n, m)$Design
  
  designs <- list(
    GSOA = D_gsoa,
    Uniform = D_unif,
    Maximin = D_maximin,
    MaxPro = D_maxpro
  )
  
  criteria <- list(
    Euclidean = euclidean_distance,
    WD = WD,
    CD = CD,
    rho_ave = corr_aveabs
  )
  
  k_values <- 1:(m-1)
  
  # 设置并行核心数
  if(is.null(n_cores)) {
    n_cores <- detectCores() - 1
  }
  
  # 创建所有任务组合
  all_tasks <- expand.grid(
    design_idx = 1:length(designs),
    k_idx = 1:length(k_values)
  )
  
  n_tasks <- nrow(all_tasks)
  cat(sprintf("总任务数: %d 个设计 × %d 个k值 = %d\n", 
              length(designs), length(k_values), n_tasks))
  cat(sprintf("使用 %d 个核心进行并行计算...\n", n_cores))
  
  # 创建并行集群
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 导出必要的函数和变量到集群
  clusterExport(cl, c("designs", "criteria", "k_values", "evaluate_projections",
                      "euclidean_distance", "WD", "CD", "corr_aveabs"),
                envir = environment())
  
  # 在每个worker中加载必要的包
  clusterEvalQ(cl, {
    library(foreach)
    library(lhs)
  })
  
  # 并行执行
  cat("开始并行计算...\n")
  results_list <- foreach(task_idx = 1:n_tasks, .combine = 'rbind',
                          .errorhandling = 'stop') %dopar% {
                            
                            design_idx <- all_tasks[task_idx, "design_idx"]
                            k_idx <- all_tasks[task_idx, "k_idx"]
                            k <- k_values[k_idx]
                            
                            current_design <- designs[[design_idx]]
                            proj_stats <- evaluate_projections(current_design, k, criteria)
                            
                            # 构建结果数据框
                            result_row <- data.frame(
                              design_idx = design_idx,
                              design_name = names(designs)[design_idx],
                              k_idx = k_idx,
                              k = k,
                              stringsAsFactors = FALSE
                            )
                            
                            # 动态添加各个准则的结果
                            for(crit_name in names(proj_stats)) {
                              result_row[[crit_name]] <- proj_stats[[crit_name]]
                            }
                            
                            result_row
                          }
  
  # 关闭集群
  stopCluster(cl)
  cat("并行计算完成！\n")
  
  # 整理投影结果为矩阵格式
  criteria_names <- names(criteria)
  projection_results <- list()
  for(crit_name in criteria_names) {
    proj_mat <- matrix(NA, nrow = length(designs), ncol = length(k_values))
    rownames(proj_mat) <- names(designs)
    colnames(proj_mat) <- paste0("k=", k_values)
    
    for(idx in 1:nrow(results_list)) {
      proj_mat[results_list$design_idx[idx], results_list$k_idx[idx]] <- 
        results_list[idx, crit_name]
    }
    projection_results[[crit_name]] <- proj_mat
  }
  
  return(list(
    projection_results = projection_results,
    k_values = k_values,
    designs = designs,
    raw_results = results_list,
    criteria_names = criteria_names
  ))
}

# ============================================================
# 4. 运行并行计算
# ============================================================

result <- compare_designs(n_cores = 6)

save(result,file = "comparison_result.Rdata")
# ============================================================
# 5. plot
# ============================================================
plot_comparison <- function(results) {

  if (!is.null(results$full_results)) {
    design_names <- rownames(results$full_results)
  } else if (!is.null(names(results$designs))) {
    design_names <- names(results$designs)
  } else {
    design_names <- c("GSOA", "Uniform", "Maximin", "MaxPro")
  }
  design_names <- as.character(design_names)  
  
  k_values <- results$k_values
  
  prepare_relative_data <- function(data_mat, measure_name) {
    rel_mat <- sweep(data_mat, 2, data_mat[1, ], "-")
    rel_mat[1, ] <- 0
    colnames(rel_mat) <- k_values
    rownames(rel_mat) <- design_names
    as.data.frame(rel_mat) %>%
      mutate(design = rownames(.)) %>%
      pivot_longer(cols = -design, names_to = "k", values_to = "value") %>%
      mutate(k = as.numeric(k),
             measure = measure_name)
  }
  
  euclidean_rel <- prepare_relative_data(results$projection_results$Euclidean, "Euclidean")
  wd_rel       <- prepare_relative_data(results$projection_results$WD, "WD")
  cd_rel       <- prepare_relative_data(results$projection_results$CD, "CD")
  
  rho_mat <- results$projection_results$rho_ave
  colnames(rho_mat) <- k_values
  rownames(rho_mat) <- design_names
  rho_long <- as.data.frame(rho_mat) %>%
    mutate(design = rownames(.)) %>%
    pivot_longer(cols = -design, names_to = "k", values_to = "value") %>%
    mutate(k = as.numeric(k),
           measure = "rho_ave")
  
  base_theme <- theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      legend.position = "none",       
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 10)
    )
  
  p1 <- euclidean_rel %>%
    ggplot(aes(x = k, y = value, color = design, linetype = design, shape = design)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "k", y = "Relative Euclidean distance",
         title = "Euclidean distance (Relative to GSOA)") +
    base_theme +
    scale_color_manual(values = rep("black", length(design_names))) +   
    scale_linetype_manual(values = 1:length(design_names)) +
    scale_shape_manual(values = 15:(15+length(design_names)-1))+
    guides(linetype = guide_legend(title = "Design", nrow = 1),
           shape = guide_legend(title = "Design", nrow = 1),
           color = "none")  
  
  p2 <- wd_rel %>%
    ggplot(aes(x = k, y = value, color = design, linetype = design, shape = design)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "k", y = "Relative WD", title = "WD (Relative to GSOA)") +
    base_theme +
    scale_color_manual(values = rep("black", length(design_names))) +
    scale_linetype_manual(values = 1:length(design_names)) +
    scale_shape_manual(values = 15:(15+length(design_names)-1))+
    guides(linetype = guide_legend(title = "Design", nrow = 1),
           shape = guide_legend(title = "Design", nrow = 1),
           color = "none")  
  
  p3 <- cd_rel %>%
    ggplot(aes(x = k, y = value, color = design, linetype = design, shape = design)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = "k", y = "Relative CD", title = "CD (Relative to GSOA)") +
    base_theme +
    scale_color_manual(values = rep("black", length(design_names))) +
    scale_linetype_manual(values = 1:length(design_names)) +
    scale_shape_manual(values = 15:(15+length(design_names)-1))+
    guides(linetype = guide_legend(title = "Design", nrow = 1),
           shape = guide_legend(title = "Design", nrow = 1),
           color = "none")  
  

  p4 <- rho_long %>%
    ggplot(aes(x = k, y = value, color = design, linetype = design, shape = design)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    labs(x = "k", y = "Correlation", title = "Average absolute correlation") +
    base_theme +
    theme(legend.position = "bottom") +  
    scale_color_manual(values = rep("black", length(design_names))) +
    scale_linetype_manual(values = 1:length(design_names)) +
    scale_shape_manual(values = 15:(15+length(design_names)-1)) +
    guides(linetype = guide_legend(title = "Design", nrow = 1),
           shape = guide_legend(title = "Design", nrow = 1),
           color = "none")   
  
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(theme = theme(plot.title = element_blank())) &
    theme(legend.position = "bottom",legend.title = element_text(size = 12),legend.text = element_text(size = 12))     
  
  return(combined)
}

P <- plot_comparison(result)
ggsave("comparison_result.png", P,
       width = 16, height = 16, dpi=300, bg="white")
