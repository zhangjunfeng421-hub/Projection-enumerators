library(lhs)           # random Latin hypercube
library(SLHD)          # maximin Latin hypercube
library(MaxPro)        # maximum projection Latin hypercube
library(DiceDesign)   
library(DiceKriging)   # GP model
library(DoE.base)
library(ggplot2)
library(gridExtra)
source("GSOA_opt.R")
source("functions_test.R")

n_runs <- 27
n_factors <- 13
s <- 3

# ============================================================
# plot function
# ============================================================ 
design_labels <- c(
  "ud" = "Uniform ",
  "gsoa" = "GSOA", 
  "rlhd" = "Random LHD",
  "maximin" = "Maximin",
  "maxpro" = "MaxPro"
)

Plot_function <- function(results_list, name_title, Ylim = c(0.3,1.5)){
  plot_data <- data.frame()
  for (name in names(results_list)) {
    plot_data <- rbind(plot_data, 
                       data.frame(
                         Design = name,
                         MSE = results_list[[name]]
                       ))
  }
  
  plot_data$Design <- factor(plot_data$Design, 
                             levels = names(design_labels),
                             labels = design_labels)
  
  p <- ggplot(plot_data, aes(x = Design, y = MSE)) +
    geom_boxplot(fill = "gray70", color = "black", alpha = 0.7, width = 0.7, fatten = 1.2) +
    stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.5, color = "black") +
    coord_cartesian(ylim = Ylim) +
    labs(
      title = name_title,
      y = "Normalized RMSE",
      x = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "gray85", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
  return(p)
}

# ============================================================
# Define MSE
# ============================================================
MSE <- function(y_pred, y_true,y_bar) {
  mse <- mean((y_pred - y_true)^2)
  dominator <- mean((y_pred-y_bar)^2)
  return((mse/dominator)^0.5)
}

# ============================================================
# transform the design
# ============================================================
transform_to_original <- function(X,ranges) {
  n_rf=nrow(ranges)
  n <- nrow(X)
  d <- ncol(X)
  X_orig <- matrix(0, n, d)
  
  for (j in 1:min(n_rf,d)) {
    X_orig[, j] <- ranges[j, 1] + X[, j] * (ranges[j, 2] - ranges[j, 1])
  }
  
  if(n_rf < d){
    X_orig[, (n_rf+1):d] <- X[,(n_rf+1):d]
  }
  
  return(X_orig)
}

# ============================================================
# main function for the simulations
# ============================================================
run_simulation <- function(design_name, design_matrix,  n_reps = 100, test_data = X_test, 
                           gp_type ="gauss", model = function(x) {
                             
                             x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
                             x5 <- x[5]; x6 <- x[6]; x7 <- x[7]; x8 <- x[8]
                             
                             numerator <- 2 * pi * x3 * (x4 - x6)
                             denominator <- log(x2 / x1) * (1 + (2 * pi * x3) / (log(x2 / x1) * x1^2 * x8) + x3 / x5)
                             
                             y <- numerator / denominator
                             return(y)
                           },Ranges= matrix(c(
                             0.05, 0.15,      # x1
                             100, 50000,      # x2
                             63070, 115600,   # x3
                             990, 1110,       # x4
                             63.1, 116,       # x5
                             700, 820,        # x6
                             1120, 1680,      # x7
                             9855, 12045      # x8
                           ), nrow = 8, ncol = 2, byrow = T)) {
  
  n <- nrow(design_matrix)
  m <- ncol(design_matrix)
  results <- numeric(n_reps)
  test_data <- transform_to_original(test_data, ranges = Ranges)
  
  # compute the real responses on test set
  y_test <- apply(test_data, 1, model)
  #y_test_mean <- mean(y_test)
  
  for (rep in 1:n_reps) {
    D_1 <- design_matrix
    if (rep > 1) {
      col_order <- sample(1:m,m)
      D_perm <- D_1[, col_order]
      dimnames(D_perm) <- dimnames(D_1)
    }else{
      D_perm <- D_1
    }
    
    D_train_orig <- transform_to_original(D_perm, ranges = Ranges)
    
    y_train <- apply(D_train_orig, 1, model)
    y_train_log <- log(y_train)
    y_train_mean <- mean(y_train)
    
    # fit model
    gp_model <- km(
      formula = ~1,
      design = D_train_orig,
      response = y_train_log,
      nugget = 1e-6,
      covtype = gp_type,
      control = list(pop.size = 100, trace = FALSE)
    )
    
    if(any(coef(gp_model)$range < 1e-5)){
      gp_model <- km(
        design = D_train_orig,
        response = y_train_log,
        nugget = 1e-8,
        covtype = gp_type,
        control = list(pop.size = 100, trace = FALSE)
      )
    }
    
    
    # predict response
    y_pred_log <- predict(gp_model, newdata = test_data, type = "UK")$mean
    y_pred <- exp(y_pred_log)
    
    # calculate mse
    mse <- MSE(y_pred, y_test,y_train_mean)
    
    if(mse>4){
      gp_model <- km(
        design = D_train_orig,
        response = y_train_log,
        nugget = 1e-8,
        covtype = gp_type,
        control = list(pop.size = 100, trace = FALSE)
      )
      
      # predict response
      y_pred_log <- predict(gp_model, newdata = test_data, type = "UK")$mean
      y_pred <- exp(y_pred_log)
      
      # calculate mse
      mse <- MSE(y_pred, y_test,y_train_mean)
    }
    
    results[rep] <- mse
  }
  
  return(results)
}

# ============================================================
# Generate test data (random Latin hypercube, N = 10000)
# ============================================================

randomLHD <- function(n,k)
{ # generate a random LHD
  x <- matrix(0,n,k)
  for(j in 1:k) x[,j] = sample(0:(n-1), n)
  return(x)
}
N_test <- 10000
X_test <- randomLHD(N_test, n_factors)
X_test <- (2 * X_test + 1)/(2 * max(X_test)+2)
# ============================================================
# Generate designs
# ============================================================
cat("Generate the designs...\n")
designs <- list()

# 1.Generate UD
cat("  1. UD...\n")
X <- lhsDesign(n_runs,n_factors)$design
designs$ud <- discrepESE_LHS(X)$design

# 2.Generate optimal GSOA
cat("  2. GSOA...\n")
X <- GSOA_opt(c(1,0,2,1), s, K = 1, add = c(1,1,1))
X <- (2 * X + 1)/(2 * max(X)+2)
designs$gsoa <- X

# 3.random LHD
cat("  3. Random LHD...\n")
X <- randomLHD(n_runs, n_factors)
X <- (2 * X + 1)/(2 * max(X)+2)
designs$rlhd <- X

# 4. Maximin Latin hypercube
cat("  4. Maximin LHD...\n")
X <- maximinSLHD(1,n_runs, n_factors)$Design
X <- (2 * X - 1)/(2 * max(X))
designs$maximin <- X

# 5. Maximum projection Latin hypercube
cat("  5. MaxPro LHD...\n")
designs$maxpro <- MaxProLHD(n_runs, n_factors)$Design


# ============================================================
# matern5 (borehole)
# ============================================================


cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")

results_list_matern5_borehole <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_matern5_borehole[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "matern5_2"
  )
}

save(results_list_matern5_borehole,file = "matern5_borehole.Rdata")
p_matern5_borehole <- Plot_function(results_list_matern5_borehole,"Matern 5/2 kernel and Borehole", Ylim = c(0.3,1.4))
# ============================================================
# gauss (borehole)
# ============================================================

cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")
results_list_gauss_borehole <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_gauss_borehole[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "gauss"
  )
}

save(results_list_gauss_borehole,file = "gauss_borehole.Rdata")
p_gauss_borehole <- Plot_function(results_list_gauss_borehole,"Gaussian kernel and Borehole",Ylim = c(0.1,0.9))
# ============================================================
# marten3 (borehole)
# ============================================================

cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")
results_list_matern3_borehole <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_matern3_borehole[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "matern3_2"
  )
}

save(results_list_matern3_borehole,file = "matern3_borehole.Rdata")
p_matern3_borehole <- Plot_function(results_list_matern3_borehole,"Matern 3/2 kernel and Borehole")
# ============================================================
# matern5 (Wing)
# ============================================================


cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")

results_list_matern5_wing <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_matern5_wing[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "matern5_2",
    model = Wing,
    Ranges = Ranges_Wing
  )
}

save(results_list_matern5_wing,file = "matern5_Wing.Rdata")
p_matern5_wing <- Plot_function(results_list_matern5_wing,"Matern 5/2 kernel and Wing",Ylim = c(0.1,1.4))
# ============================================================
# gauss (Wing)
# ============================================================

cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")
results_list_gauss_wing <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_gauss_wing[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "gauss",
    model = Wing,
    Ranges = Ranges_Wing
  )
}

save(results_list_gauss_wing,file = "gauss_Wing.Rdata")
p_gauss_wing <- Plot_function(results_list_gauss_wing,"Gaussian kernel and Wing",Ylim = c(0.1,1.1))
# ============================================================
# marten3 (Wing)
# ============================================================

cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")
results_list_matern3_wing <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_matern3_wing[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "matern3_2",
    model = Wing,
    Ranges = Ranges_Wing
  )
}

save(results_list_matern3_wing,file = "matern3_Wing.Rdata")
p_matern3_wing <- Plot_function(results_list_matern3_wing,"Matern 3/2 kernel and Wing",Ylim = c(0.2,1.7))
# ============================================================
# matern5 (OTL)
# ============================================================


cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")

results_list_matern5_OTL <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_matern5_OTL[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "matern5_2",
    model = OTL,
    Ranges = Ranges_OTL
  )
}

save(results_list_matern5_OTL,file = "matern5_OLT.Rdata")
p_matern5_OTL <- Plot_function(results_list_matern5_OTL,"Matern 5/2 kernel and OTL",Ylim = c(0.1,0.6))
# ============================================================
# gauss (OLT)
# ============================================================

cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")
results_list_gauss_OTL <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_gauss_OTL[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "gauss",
    model = OTL,
    Ranges = Ranges_OTL
  )
}

save(results_list_gauss_OTL,file = "gauss_OLT.Rdata")
p_gauss_OTL <- Plot_function(results_list_gauss_OTL,"Gaussian kernel and OTL",Ylim = c(0.1,0.4))
# ============================================================
# marten3 (OTL)
# ============================================================

cat("Start the simulation...\n")

cat("\n excute simulation (100 repetitions)...\n")
results_list_matern3_OTL <- list()

for (name in names(designs)) {
  cat(sprintf(" working in progress: %s...\n", name))
  results_list_matern3_OTL[[name]] <- run_simulation(
    design_name = name,
    design_matrix = designs[[name]],
    gp_type = "matern3_2",
    model = OTL,
    Ranges = Ranges_OTL
  )
}

save(results_list_matern3_OTL,file = "matern3_OLT.Rdata")
p_matern3_OTL <- Plot_function(results_list_matern3_OTL,"Matern 3/2 kernel and OTL",Ylim = c(0.1,0.7))
cat("\n The simulation is completed\n")
P <- grid.arrange( 
  p_gauss_OTL, p_matern5_OTL, p_matern3_OTL,
  p_gauss_borehole, p_matern5_borehole, p_matern3_borehole,
  p_gauss_wing, p_matern5_wing, p_matern3_wing,
  ncol = 3
)

ggsave("simulation_results.png", P,
       width = 12, height = 12, dpi=300, bg="white")
