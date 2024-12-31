# 1) source the external R file
source("basic_functions.R")

# 2) Now all the functions are available
#    We can call them directly:

###### Usage 1: Initial Parameters Setup and Simulation 1 Setup
epss=1e-8
p=3
v_vec= c(0.4,0.6)
eta=0.5
k_B=length(v_vec)
par_test <- init_par_test(p=p, k_B=k_B, seed = 2025) # c(theta_1, ..., theta_p, phi_1, ..., phi_{k_B})
par_test[(p+1):(p+k_B)]=c(-1.5,-0.1)
par_test[1:p]=par_test[1:p]-mean(par_test[1:p])

t(construct_U(3))%*%par_test[1:p]

connected_pair_matrix <-init_edge_generate_RA(p=p, method="full", reg=2) #"full","regular","min" 
k=nrow(connected_pair_matrix)
H_test=build_soft_ranking_H(par_test, p=p, eta=eta, delta=1)
#build_H(p, k_B, eta_theta = 1, eta = 0.5)
pi_init <- init_pi_mat(k=k, v_vec=v_vec)
# solve
res <- solve_pi_nloptr(
  par_test, connected_pair_matrix, k_B=k_B, H=H_test,
  v_vec= v_vec,
  pi_init= pi_init,
  maxeval=2000
);
pi_soft_limit=res$pi




###### Usage 2: Initial Parameters Setup and Simulation 1 perform
set.seed(2025)
# initial data round-robin
out_gen <- data_random_generator_roundrobin_BTD(
  p, connected_pair_matrix, k_B,
  theta_min=-2,theta_max=2,
  phi_min=-3, phi_max=3,par_true=par_test
)
df_data_init <- out_gen$data_random
vec_b_init <- build_vec_b_from_df(df_data_init, k_B=k_B)
T_final <- 50
par_true     <- par_test  # from your usage

# we define a small function that, given par_mle, returns H:
H_dynamic_fun <- function(par_mle) {
  # e.g. build_soft_ranking_H( par_true= par_mle, p= p, eta=0.5, delta=1 )
  build_soft_ranking_H(par_true= par_mle, p= p, eta=0.5, delta=1)
}

res_final_dyn <- multiple_rounds_adaptive_dynamicH(
  T              = T_final,
  df_data_init   = df_data_init,
  pi_n_init      = NULL,
  connected_pair_matrix= connected_pair_matrix,
  p= p, 
  k_B= k_B,
  R= 3,
  phi_lower= -3, phi_upper= 3,
  init_par= NULL,  # random or pass something
  method= "NLOPT_LD_SLSQP",
  maxeval=2000,
  verbose=FALSE,
  v= v_vec,
  vec_b_init= vec_b_init,
  par_true= par_true,
  H_dynamic_fun= H_dynamic_fun  # supply the function that builds H
)



par_test-res_final_dyn$final_par_mle

###### plot

# 1) Create a numeric vector to store the differences
total_point=T_final-k+1
diff_vec <- numeric(total_point)

# 2) Loop (or use sapply) over i from 1..998
for (i in 1:total_point) {
  val_i <- F_pi_par(
    res_final_dyn$pi_tensor[,, i], 
    par_test, 
    connected_pair_matrix,
    k_B, 
    H_test
  )
  val_soft <- F_pi_par(
    pi_soft_limit, 
    par_test, 
    connected_pair_matrix,
    k_B, 
    H_test
  )
  diff_vec[i] <- val_i - val_soft
}

# 3) Plot the vector
# a) Open a graphics device with high resolution
png(
  filename = "my_figure.png",
  width = 1200,    # in pixels
  height = 900,    # in pixels
  res = 300        # DPI (dots per inch)
)

# b) Suppose you have diff_vec, total_point already defined
#    Create a line plot with a red horizontal line at 0
plot(
  x = 1:total_point +2,    # or 1:total_point + 2
  y = diff_vec,
  type = "l",
  xlab = "Sample Size n",
  # ylab with subscript => expression(F[n]) for F_n
  ylab = expression(F[n]),
  # optionally remove main title or set your own
  main = NULL
)

abline(h = 0, col = "red", lwd = 2)

# c) Close the device
dev.off()








