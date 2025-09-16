#RABVITG for estimating parameters knowing the state variable
# Load necessary libraries
library(stats)
library(ggplot2)
library(cowplot)
library(pracma)
library(nloptr)
library(plotly)
install.packages("tidyr")
library(tidyr)
library(randtoolbox)
install.packages("gridExtra")
library(gridExtra)
 
rm(list = ls())
# Define Constants


xmin <- 0
xmax <- 5
Num <-7   #número de puntos discretización espacial #7
dt <- 0.01  #discretización temporal
rho1 <- 0.1
rho2 <- 0.1
r <- 2
a1 <- 1
a2 <- 1

A<-8
B<-4
K<-4
eps<-0.25
sigma<-0  #si es estocástico, distinto de cero
s<-2
u_star<-0.5
ustar<-0.5
Tfin<-30000




# Define the RBF interpolation function
rbf_interpolate <- function(X, y, s, xnew) {
  # Define the Multiquadric RBF function
  gaussian_rbf <- function(r, s) {
    exp(-(r^2) / (2 * s^2))
  }
  # Calculate the distance matrix for the input points
  dist_matrix <- as.matrix(dist(X))
  
  # Compute the RBF matrix for the input points
  K <- gaussian_rbf(dist_matrix, s)
  
  # Solve for the weights
  weights <- solve(K, y)
  
  # Function to calculate the RBF for new points
  interpolate_point <- function(xnew_row) {
    dist_new <- sqrt(rowSums((t(t(X) - xnew_row))^2))
    rbf_new <- gaussian_rbf(dist_new, s)
    sum(rbf_new * weights)
  }
  
  
  # Ensure xnew is a matrix with at least one row
  xnew <- matrix(xnew, nrow = 1)
  
  # Apply the interpolation to each row in xnew
  interpolated_values <- apply(xnew, 1, interpolate_point)
  
  
  return(interpolated_values)
}



# Funciones Valor por jugador
valfun1 <- function(k, x01, u01, u02, x1, v0, ustar,r,a1,a2,sigma) {
  xnew1 <- x01 + dt * (-r * x01 + a1 * k + a2 * u02) + sqrt(dt)*sigma
  xnew1 <- pmax(xmin*0.90, pmin(xmax*1.10, xnew1))
  xnew2 <- x01 + dt * (-r * x01 + a1 * k + a2 * u02) - sqrt(dt)*sigma
  xnew2 <- pmax(xmin*0.90, pmin(xmax*1.10, xnew2))
  
  X <- matrix(x1)
  # Define output values
  y <- matrix(v0[,1])
  # Define shape parameter
  
  
  
  # Perform interpolation
  val1_interp1  <- rbf_interpolate(X, y, s, xnew1)
  val1_interp2  <- rbf_interpolate(X, y, s, xnew2)
  

  
  val1_interp<- mean(c( val1_interp1,val1_interp2))
  
  
  
    
    val1 <- dt * (A*(x01+1)^alpha1 -  (beta1* (k - ustar)^2+ ((1-beta1)*(k-ustar)^2)/(1+(k-ustar)^2))) + (1 - dt * rho1) * val1_interp
    
    
    
  
  
  val <- -val1
  return(val)
}

valfun2 <- function(k, x01, u01, u02, x1, v0, ustar,r,a1,a2,sigma) {
  xnew1 <- x01 + dt * (-r * x01 + a1 * u01 + a2 * k) + sqrt(dt)*sigma
  xnew1 <- pmax(xmin*0.90, pmin(xmax*1.10, xnew1))
  xnew2 <- x01 + dt * (-r * x01 + a1 * u01 + a2 * k) - sqrt(dt)*sigma
  xnew2 <- pmax(xmin*0.90, pmin(xmax*1.10, xnew2))
  
  
  X <- matrix(x1)
  # Define output values
  y <- matrix(v0[,2])
  # Define shape parameter
  
  # Perform interpolation
  val2_interp1  <- rbf_interpolate(X, y, s, xnew1)
  val2_interp2  <- rbf_interpolate(X, y, s, xnew2)
  
  
  
  
  
  
  val2_interp<- mean(c( val2_interp1,val2_interp2))
  
  
    
  
  val2 <- dt * (A*(x01+1)^alpha2 -  (beta2* (k - ustar)^2+ ((1-beta2)*(k-ustar)^2)/(1+(k-ustar)^2))) + (1 - dt * rho2) * val2_interp
  
  
  
  
  
  
  
  val <- -val2
  
  return(val)
}

#Función principal RaBVItG

LOVE2D <- function(Num, ustar, itmax,r,a1,a2,alpha1,alpha2,beta1,beta2,sigma) {
  tol <- 0.00001
  
  # Initialize variables
  N <- Num
  difference <- numeric(itmax)
  
  x1 <- seq(xmin, xmax, length.out = N)
  v0 <- cbind(dt*A*(x1+1)^alpha1, dt*A*(x1+1)^alpha2)
  dif <- c(tol + 1000, 2)
  its <- 1
  
  # Initialize arrays to store results
  k11 <- numeric(N)
  k22 <- numeric(N)
  v1 <- matrix(0, nrow = N, ncol = 2)
  
  while (max(dif) > tol & its < itmax) {
    uold1 <- k11
    uold2 <- k22
    
    for (i in 1:N) {
      x01 <- x1[i]
      
      u01 <- uold1[i]
      u02 <- uold2[i]
      
      tolg <- 0.00001
      errorg <- tolg + 1
      itsg <- 1
      
      
      #Game Iteration
      
      while (errorg > tolg & itsg < 25) {
        itsg <- itsg + 1
        
        k1 <- optimise(valfun1, interval = c(ustar, 100), x01 = x01, u01 = u01, u02 = u02, x1 = x1, v0 = v0,ustar=ustar,r=r,a1=a1,a2=a2, sigma=sigma)$minimum
        k2 <- optimise(valfun2, interval = c(ustar, 100), x01 = x01, u01 = u01, u02 = u02, x1 = x1, v0 = v0,ustar=ustar,r=r,a1=a1,a2=a2, sigma=sigma)$minimum
        
        policystar <- c(k1, k2)
        policyold <- c(u01, u02)
        errorg <- max(abs(policystar - policyold))
        u01 <- 0.75 * k1 + 0.25 * u01
        u02 <- 0.75 * k2 + 0.25 * u02
      }
      
      k11[i] <- u01
      k22[i] <- u02
      
      v1[i, 1] <- (-valfun1(u01, x01, u01, u02, x1, v0,ustar,r,a1,a2,sigma))
      v1[i, 2] <- (-valfun2(u02, x01, u01, u02, x1, v0,ustar,r,a1,a2,sigma))
    }
    
    unew <- cbind(k11, k22)
    uold <- cbind(uold1,uold2)
    dif <- max(abs(unew - uold)) #RECUERDA CAMBIARLO
    #dif <- max(abs(v0 - v1))
    difference[its]<- dif
    dif <- tail(dif, 1)
    
    #Value Iteration
    v0[,1] <- v1[,1]
    v0[,2] <- v1[,2]
    #par(mfrow = c(1, 2))
    #plot(k11)
    #plot(v1[,1])
    #par(mfrow = c(1, 2))
    #print(its)
    
    
    its <- its + 1
  }
  
  return(list(v0_0 = v0, k11 = k11, k22 = k22))
}



############################################################################



alpha1_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
alpha2_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# Crear data frame vacío para almacenar resultados
C1<-matrix(NA, nrow = length(alpha1_values), ncol = length(alpha2_values))
C2<-matrix(NA, nrow = length(alpha1_values), ncol =  length(alpha2_values))
V1<-matrix(NA, nrow = length(alpha1_values), ncol = length(alpha2_values))
V2<-matrix(NA, nrow = length(alpha1_values), ncol =  length(alpha2_values))
D12<-matrix(NA, nrow = length(alpha1_values), ncol =  length(alpha2_values))
Xstar<-matrix(NA, nrow = length(alpha1_values), ncol =  length(alpha2_values))
val1<-0
val2<-0
value1<-0
value2<-0

# Bucle para calcular la dinámica para cada valor de a2
for (k in 1:length(alpha1_values)) {
  
  for (l in 1:length(alpha2_values)) {
    
    alpha1 <- alpha1_values[k]
    alpha2 <- alpha2_values[l]
    beta1<-1
    beta2<-1
    
    # Calcular la política de control óptimo usando LOVE2D
    LOVE2D_result <- LOVE2D(Num, ustar, 400, r, a1, a2,alpha1,alpha2,beta1,beta2, sigma)
    U_fin1 <- LOVE2D_result$k11
    U_fin2 <- LOVE2D_result$k22
    V_fin1 <- LOVE2D_result$v0_0[,1]
    V_fin2 <- LOVE2D_result$v0_0[,2]
    
    
    # Inicialización de las trayectorias
    x1 <- seq(xmin, xmax, length.out = Num)
    
    x1 <- seq(xmin, xmax, length.out = Num)
    
    x<-0
    x[1]<-4.5
    u1<-0
    u1[1]<-cubicspline(x1, U_fin1, xi = x[1])
    u2<-0
    u2[1]<-cubicspline(x1, U_fin2, xi = x[1])
    
    
    for (i in 1:Tfin){
      
      
      x[i+1]<-x[i] + dt * (-r * x[i] + a1 *u1[i] + a2 *u2[i])
      
      
      
      X <- matrix(x1)
      # Define output values
      y1 <- matrix(U_fin1)
      # Define shape parameter
      y2 <- matrix(U_fin2)
      # Perform interpolation
      
      
      
      
      u1[i+1] <-rbf_interpolate(X, y1, s, x[i+1])
      u1[i+1] <-pmax(ustar,pmin(50,u1[i+1]))
      u2[i+1] <-rbf_interpolate(X, y2, s, x[i+1])
      u2[i+1] <-pmax(ustar,pmin(50,u2[i+1]))
      
     
    
     
        val1 <- dt * (A*(x[i]+1)^alpha1 -  (beta1* (u1[i] - ustar)^2+ ((1-beta1)*(u1[i]-ustar)^2)/(1+(u1[i]-ustar)^2)))
        
        
  
      
      val1[i]<-val1
      
      
      
        
      val2 <- dt * (A*(x[i]+1)^alpha2 -  (beta2* (u2[i] - ustar)^2+ ((1-beta2)*(u2[i]-ustar)^2)/(1+(u2[i]-ustar)^2)))
      
  
      val2[i]<-val2
      
      value1[i]<-exp(-rho1*dt*(i-1))*val1[i]
      value2[i]<-exp(-rho2*dt*(i-1))*val2[i]
      
      
    }
    
    
    
    # Guardar resultados en formato largo
    
    C1[k,l]<-(u1[Tfin])
    C2[k,l]<- (u2[Tfin])
    V1[k,l]<- sum(value1)
    V2[k,l]<- sum(value2)
    D12[k,l]<-(u1[Tfin]-u2[Tfin])
    Xstar[k,l]<-(x[Tfin])
    
    print(c(k,l))
    
  }
  
}

save(C1, C2, V1, V2, D12,Xstar, file = "MODELOB_FIJO_BETA.RData")







#############################CAMBIAR NOMBRE ARCHIVO





beta1_values <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
beta2_values <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
# Crear data frame vacío para almacenar resultados
C1<-matrix(NA, nrow = length(beta1_values), ncol = length(beta2_values))
C2<-matrix(NA, nrow = length(beta1_values), ncol =  length(beta2_values))
V1<-matrix(NA, nrow = length(beta1_values), ncol = length(beta2_values))
V2<-matrix(NA, nrow = length(beta1_values), ncol =  length(beta2_values))
D12<-matrix(NA, nrow = length(beta1_values), ncol =  length(beta2_values))
Xstar<-matrix(NA, nrow = length(beta1_values), ncol =  length(beta2_values))
val1<-0
val2<-0
value1<-0
value2<-0

# Bucle para calcular la dinámica para cada valor de a2
for (k in 1:length(beta1_values)) {
  
  for (l in 1:length(beta2_values)) {
    
    alpha1 <- 0.5
    alpha2 <- 0.5
    beta1<-beta1_values[k]
    beta2<-beta2_values[l]
    
    # Calcular la política de control óptimo usando LOVE2D
    LOVE2D_result <- LOVE2D(Num, ustar, 500, r, a1, a2,alpha1,alpha2,beta1,beta2, sigma)
    U_fin1 <- LOVE2D_result$k11
    U_fin2 <- LOVE2D_result$k22
    V_fin1 <- LOVE2D_result$v0_0[,1]
    V_fin2 <- LOVE2D_result$v0_0[,2]
    
    
    # Inicialización de las trayectorias
    x1 <- seq(xmin, xmax, length.out = Num)
    
    x1 <- seq(xmin, xmax, length.out = Num)
    
    x<-0
    x[1]<-4.5
    u1<-0
    u1[1]<-cubicspline(x1, U_fin1, xi = x[1])
    u2<-0
    u2[1]<-cubicspline(x1, U_fin2, xi = x[1])
    
    
    for (i in 1:Tfin){
      
      
      x[i+1]<-x[i] + dt * (-r * x[i] + a1 *u1[i] + a2 *u2[i])
      
      
      
      X <- matrix(x1)
      # Define output values
      y1 <- matrix(U_fin1)
      # Define shape parameter
      y2 <- matrix(U_fin2)
      # Perform interpolation
      
      
      
      
      u1[i+1] <-rbf_interpolate(X, y1, s, x[i+1])
      u1[i+1] <-pmax(ustar,pmin(50,u1[i+1]))
      u2[i+1] <-rbf_interpolate(X, y2, s, x[i+1])
      u2[i+1] <-pmax(ustar,pmin(50,u2[i+1]))
      
      
      val1 <- dt * (A*(x[i]+1)^alpha1 -  (beta1* (u1[i] - ustar)^2+ ((1-beta1)*(u1[i]-ustar)^2)/(1+(u1[i]-ustar)^2)))
      
      
      
      
      val1[i]<-val1
      
      
      
      
      val2 <- dt * (A*(x[i]+1)^alpha2 -  (beta2* (u2[i] - ustar)^2+ ((1-beta2)*(u2[i]-ustar)^2)/(1+(u2[i]-ustar)^2)))
      
      
      val2[i]<-val2
      
      value1[i]<-exp(-rho1*dt*(i-1))*val1[i]
      value2[i]<-exp(-rho2*dt*(i-1))*val2[i]
      
      
    }
    
    
    
    # Guardar resultados en formato largo
    
    C1[k,l]<-(u1[Tfin])
    C2[k,l]<- (u2[Tfin])
    V1[k,l]<- sum(value1)
    V2[k,l]<- sum(value2)
    D12[k,l]<-(u1[Tfin]-u2[Tfin])
    Xstar[k,l]<-(x[Tfin])
    
    print(c(k,l))

  }
  
}

C1_a<-C1
C2_a<-C2
V1_a<-V1
V2_a<-V2
D12_a<-D12
Xstar_a<-Xstar

save(C1_a, C2_a, V1_a, V2_a, D12_a,Xstar_a, file = "MODELOB_FIJO_ALPHA.RData")




#################### Joint variation





alpha1_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
beta1_values <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
# Crear data frame vacío para almacenar resultados
C1<-matrix(NA, nrow = length(alpha1_values), ncol = length(beta1_values))
C2<-matrix(NA, nrow = length(alpha1_values), ncol =  length(beta1_values))
V1<-matrix(NA, nrow = length(alpha1_values), ncol = length(beta1_values))
V2<-matrix(NA, nrow = length(alpha1_values), ncol =  length(beta1_values))
D12<-matrix(NA, nrow = length(alpha1_values), ncol =  length(beta1_values))
Xstar<-matrix(NA, nrow = length(alpha1_values), ncol =  length(beta1_values))
val1<-0
val2<-0
value1<-0
value2<-0

# Bucle para calcular la dinámica para cada valor de a2
for (k in 1:length(alpha1_values)) {
  
  for (l in 1:length(beta1_values)) {
    
    alpha1 <- alpha1_values[k]
    alpha2 <- 0.5
    beta1<-beta1_values[l]
    beta2<-2
    
    # Calcular la política de control óptimo usando LOVE2D
    LOVE2D_result <- LOVE2D(Num, ustar, 500, r, a1, a2,alpha1,alpha2,beta1,beta2, sigma)
    
   
     
     
    U_fin1 <- LOVE2D_result$k11
    U_fin2 <- LOVE2D_result$k22
    V_fin1 <- LOVE2D_result$v0_0[,1]
    V_fin2 <- LOVE2D_result$v0_0[,2]
    
    
    # Inicialización de las trayectorias
    x1 <- seq(xmin, xmax, length.out = Num)
    
    x1 <- seq(xmin, xmax, length.out = Num)
    
    x<-0
    x[1]<-4.5
    u1<-0
    u1[1]<-cubicspline(x1, U_fin1, xi = x[1])
    u2<-0
    u2[1]<-cubicspline(x1, U_fin2, xi = x[1])
    
    
    for (i in 1:Tfin){
      
      
      x[i+1]<-x[i] + dt * (-r * x[i] + a1 *u1[i] + a2 *u2[i])
      
      
      
      X <- matrix(x1)
      # Define output values
      y1 <- matrix(U_fin1)
      # Define shape parameter
      y2 <- matrix(U_fin2)
      # Perform interpolation
      
      
      
      
      u1[i+1] <-rbf_interpolate(X, y1, s, x[i+1])
      u1[i+1] <-pmax(ustar,pmin(50,u1[i+1]))
      u2[i+1] <-rbf_interpolate(X, y2, s, x[i+1])
      u2[i+1] <-pmax(ustar,pmin(50,u2[i+1]))
      
      
      
      val1 <- dt * (A*(x[i]+1)^alpha1 -  (beta1* (u1[i] - ustar)^2+ ((1-beta1)*(u1[i]-ustar)^2)/(1+(u1[i]-ustar)^2)))
      
      
      
      
      val1[i]<-val1
      
      
      
      
      val2 <- dt * (A*(x[i]+1)^alpha2 -  (beta2* (u2[i] - ustar)^2+ ((1-beta2)*(u2[i]-ustar)^2)/(1+(u2[i]-ustar)^2)))
      
      
      val2[i]<-val2
      
      value1[i]<-exp(-rho1*dt*(i-1))*val1[i]
      value2[i]<-exp(-rho2*dt*(i-1))*val2[i]
      
      
    }
    
    
    
    # Guardar resultados en formato largo
    
    C1[k,l]<-(u1[Tfin])
    C2[k,l]<- (u2[Tfin])
    V1[k,l]<- sum(value1)
    V2[k,l]<- sum(value2)
    D12[k,l]<-(u1[Tfin]-u2[Tfin])
    Xstar[k,l]<-(x[Tfin])
    
    print(c(k,l))
   
    
  }
  
}

C1_b<-C1
C2_b<-C2
V1_b<-V1
V2_b<-V2
D12_b<-D12
Xstar_b<-Xstar

save(C1_b, C2_b, V1_b, V2_b, D12_b,Xstar_b, file = "ModeloB_JOINT.RData")







# =========================================================

############################################################
# LIBRERÍAS
############################################################
library(ggplot2)
library(dplyr)
library(patchwork)
library(metR)      # geom_text_contour
library(viridis)
library(MBA)
library(akima)     # interp() con extrapolación
library(cowplot)   # ggdraw, draw_plot, wrap_elements
library(grid)      # unit(), arrow()
library(scales)    # number()
library(showtext)  # fuentes TTF/OTF embebidas
library(sysfonts)  # font_add_google

# ---- Fuente con letras + ♥ (evita cuadrados) ----
sysfonts::font_add_google("Noto Serif", "paperfont")
showtext::showtext_auto()

# ---- Tema auxiliar: rótulos en cursiva (sin tocar títulos) ----
theme_italics <- theme(
  axis.title   = element_text(face = "italic"),
  axis.text    = element_text(face = "italic"),
  legend.title = element_text(face = "italic"),
  legend.text  = element_text(face = "italic")
)

# Helper para guardar EPS SIN cairo:
#   1) genera PDF vectorial (useDingbats = FALSE para contornos fiables)
#   2) convierte PDF -> EPS con Ghostscript (eps2write)
save_eps <- function(plot, eps_file, width, height) {
  pdf_file <- sub("\\.eps$", ".pdf", eps_file)
  grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  print(plot)
  dev.off()
  gs <- Sys.which("gs"); if (gs == "") gs <- Sys.which("gswin64c.exe"); if (gs == "") gs <- Sys.which("gswin32c.exe")
  if (gs != "") {
    cmd <- sprintf('"%s" -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=eps2write -sOutputFile="%s" "%s"', gs, eps_file, pdf_file)
    system(cmd)
  } else {
    warning("Ghostscript no encontrado. Se dejó el PDF en: ", pdf_file, ". Convierte a EPS manualmente.")
  }
}

# Transparencias
alpha_fill <- 0.85
tile_alpha <- 0.90

############################################################
# FIGURE 3  ->  New_F3.eps
# Utility & Disutility
############################################################
U <- function(x, alpha) 8 * (x + 1)^alpha
D_B <- function(c, beta) beta * (c-0.5)^2 + ((1 - beta) * (c-0.5)^2) / (1 + (c-0.5)^2)

alphas <- c(0.1, 0.3, 0.5, 0.7, 0.9)
betas  <- c(1,1.50,2,2.50,3)
x_vals <- seq(0.001, 4.999, length.out = 600)
c_vals <- seq(0.001, 2.999, length.out = 600)

df_U <- expand.grid(x = x_vals, alpha = alphas) |> mutate(y = U(x, alpha), alpha = factor(alpha))
df_DB <- expand.grid(c = c_vals, beta = betas)   |> mutate(y = D_B(c, beta), beta = factor(beta))

pU <- ggplot(df_U, aes(x, y, color = alpha)) +
  geom_line(size = 1) +
  labs(title = "Utility", x = "x", y = "U(x)", color = expression(alpha)) +
  scale_color_brewer(palette = "Blues") +
  theme_minimal(base_family = "paperfont") + theme(legend.position = "right") + theme_italics

pDB <- ggplot(df_DB, aes(c, y, color = beta)) +
  geom_line(size = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  labs(title = "Disutility", x = "c ", y = "D(c)", color = expression(beta)) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1, 2, 3),
    labels = c(expression(0), expression(c^"*"==0.5), expression(1), expression(2), expression(3))
  ) +
  scale_color_brewer(palette = "Reds") +
  theme_minimal(base_family = "paperfont") + theme(legend.position = "right") + theme_italics

fig3 <- pU + pDB
save_eps(fig3, "New_F3.eps", width = 7.0, height = 3.6)

############################################################
# FIGURE 4  ->  New_F4.eps
# MODELOB_FIJO_BETA
############################################################
# --- Librerías ---
library(ggplot2)
library(akima)
library(metR)
library(patchwork)
library(scico)
library(scales)

# --- Datos / modelo ---
load("MODELOB_FIJO_BETA.RData")

alpha1_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
alpha2_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# --- Helper: contornos α–α con 1:1 (heat + contorno + labels) ---
make_contour_plot_alpha <- function(z_vals, breaks_contour, fill_lab, main_title) {
  df <- expand.grid(alpha1 = alpha1_values, alpha2 = alpha2_values)
  df$z <- as.vector(z_vals); df <- na.omit(df)
  
  ip <- akima::interp(
    x = df$alpha1, y = df$alpha2, z = df$z,
    xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
    yo = seq(min(alpha2_values), max(alpha2_values), length.out = 220),
    linear = TRUE, extrap = TRUE, duplicate = "mean"
  )
  
  interp_df <- expand.grid(alpha1 = ip$x, alpha2 = ip$y)
  interp_df$z <- as.vector(ip$z)
  
  ggplot(interp_df, aes(x = alpha1, y = alpha2)) +
    geom_raster(aes(fill = z), interpolate = TRUE, na.rm = TRUE, alpha = tile_alpha) +
    geom_contour(aes(z = z), color = "grey30", linewidth = 0.3,
                 breaks = breaks_contour, na.rm = TRUE) +
    metR::geom_text_contour(aes(z = z), breaks = breaks_contour, size = 3.0, stroke = 0,
                            label.placer = metR::label_placer_n(1)) +
    # --- Escala azul → rojo elegante (scico "vik") continua ---
    scale_fill_scico(palette = "vik", direction = 1, oob = scales::squish) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    coord_equal(expand = FALSE) +
    labs(title = main_title, fill = fill_lab,
         x = expression(alpha[1]), y = expression(alpha[2])) +
    theme_minimal(base_family = "paperfont") +
    theme(
      legend.position = "none",
      plot.title   = element_text(size = 10, face = "bold", margin = margin(b = 2)),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.x  = element_text(size = 10),
      axis.text.y  = element_text(size = 10),
      plot.margin  = margin(2, 2, 2, 2)
    ) + theme_italics
}

# --- Renombrados/alias ---
C1_B1 <- C1; C2_B1 <- C2; V1_B1 <- V1; V2_B1 <- V2

c1.B1 <- make_contour_plot_alpha(
  C1_B1, c(0.6,0.7,0.8,1,1.2,1.3,1.5,1.7),
  expression(c[1]^scriptstyle("\u2665")), "P.1"
)

c2.B1 <- make_contour_plot_alpha(
  C2_B1, c(0.6,0.7,0.8,1,1.2,1.3,1.5,1.7),
  expression(c[2]^scriptstyle("\u2665")), "P.2"
)

# --- Feeling (Xstar) α–α con 1:1 ---
df_x <- expand.grid(alpha1 = alpha1_values, alpha2 = alpha2_values)
df_x$z <- as.vector(Xstar); df_x <- na.omit(df_x)

ipx <- akima::interp(
  x = df_x$alpha1, y = df_x$alpha2, z = df_x$z,
  xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
  yo = seq(min(alpha2_values), max(alpha2_values), length.out = 220),
  linear = TRUE, extrap = TRUE, duplicate = "mean"
)

interp_x <- expand.grid(alpha1 = ipx$x, alpha2 = ipx$y)
interp_x$z <- as.vector(ipx$z)

p3 <- ggplot(interp_x, aes(x = alpha1, y = alpha2, z = z)) +
  geom_contour_filled(aes(fill = after_stat(level)),
                      bins = 15, show.legend = FALSE,
                      alpha = alpha_fill, na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = c(0.7, 0.85, 1.5, 1.75, 1.9),
                          size = 3.6, stroke = 0, skip = 0,
                          check_overlap = FALSE,
                          label.placer = metR::label_placer_n(1)) +
  labs(title = NULL, x = expression(alpha[1]), y = expression(alpha[2])) +
  # --- Escala DISCRETA azul → rojo (scico "vik") para los niveles binned ---
  scale_fill_manual(values = scico::scico(15, palette = "vik"),
                    drop = FALSE, guide = "none") +
  coord_equal(expand = FALSE) +
  theme_minimal(base_family = "paperfont") +
  theme(
    legend.position = "none",
    axis.title.x  = element_text(size = 10),
    axis.title.y  = element_text(size = 10),
    axis.text.x   = element_text(size = 10),
    axis.text.y   = element_text(size = 10),
    plot.margin   = margin(2, 2, 2, 2)
  ) + theme_italics +
  annotate("segment", x = 0.45, y = 0.45, xend = 0.65, yend = 0.60,
           arrow = arrow(type = "closed", angle = 20, length = unit(0.15, "inches")),
           linewidth = 0.6, color = "black") +
  annotate("text",
           x = 0.40, y = 0.40,
           label = 'italic(x)^scriptstyle("\u2665")', parse = TRUE,
           size = 5.2, color = "black", hjust = 0.5, vjust = 0.5)

# --- Títulos de columnas ---
col_title <- function(txt) {
  ggplot() + theme_void() + ggtitle(txt) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.margin = margin(0,0,0,0))
}
t1 <- col_title("Equilibrium Effort")
t2 <- col_title("Equilibrium Feeling")
t3 <- col_title("Well-being")

# --- Layout 3x3 ---
areas_sq <- c(
  patchwork::area(t = 1, l = 1, b = 1, r = 1),
  patchwork::area(t = 1, l = 2, b = 1, r = 2),
  patchwork::area(t = 1, l = 3, b = 1, r = 3),
  patchwork::area(t = 2, l = 1, b = 2, r = 1),
  patchwork::area(t = 2, l = 2, b = 2, r = 2),
  patchwork::area(t = 2, l = 3, b = 2, r = 3),
  patchwork::area(t = 3, l = 1, b = 3, r = 1),
  patchwork::area(t = 3, l = 2, b = 3, r = 2),
  patchwork::area(t = 3, l = 3, b = 3, r = 3)
)

# --- Composición final ---
final_plot_4 <-
  t1 + t2 + t3 +
  c1.B1 + p3 + make_contour_plot_alpha(V1_B1, c(85,95,100,120,140,150,170), expression(v[1]), "P.1") +
  c2.B1 + patchwork::plot_spacer() + make_contour_plot_alpha(V2_B1, c(85,95,100,120,140,150,170), expression(v[2]), "P.2") +
  plot_layout(design = areas_sq, heights = c(0.06, 1, 1), widths = c(1, 1, 1)) &
  theme(
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 10),
    plot.margin = margin(2, 2, 2, 2)
  )

# --- Export ---
save_eps(final_plot_4, "New_F4.eps", width = 7.2, height = 6.2)




############################################################
# FIGURE 5  ->  New_F5.eps
# Dinámica con alpha1=0.5 y alpha2 variable + scatter
############################################################
# Requiere: LOVE2D(), Num, ustar, r, a1, a2, sigma, xmin, xmax
alpha1_values <- 0.5
alpha2_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

theme_italics <- theme(
  axis.title = element_text(face = "italic"),
  axis.text  = element_text(face = "italic"),
  plot.title = element_text(face = "italic")
)

# Paletas: azul (bajo) -> rojo (alto)
cols_main <- c(  # para los plots principales
  "#08306B", "#1858A8", "#2C7FB8", "#4BA3D1", "#74C0E3",
  "#F6A385", "#F16D5B", "#D94841", "#A50F15"
)
cols_zoom_dark <- c(  # versión más fuerte para el zoom
  "#05244F", "#0E3E82", "#1E63B5", "#2F7CC1", "#5AAFD6",
  "#F07F61", "#E24D3A", "#C3322B", "#800A0A"
)

results_a2 <- data.frame()
for (k in seq_along(alpha2_values)) {
 

  dt=0.01;   alpha1 <- 0.5; alpha2 <- alpha2_values[k]; beta1 <- 1; beta2 <- 1
  LOVE2D_result <- LOVE2D(Num, ustar, 4000, r, a1, a2, alpha1, alpha2, beta1, beta2, sigma)
  U_fin1 <- LOVE2D_result$k11; U_fin2 <- LOVE2D_result$k22
  V_fin1 <- LOVE2D_result$v0_0[,1]; V_fin2 <- LOVE2D_result$v0_0[,2]
  x1 <- seq(xmin, xmax, length.out = Num)
  results_a2 <- rbind(results_a2, data.frame(
    Feeling = x1,
    Effort_Partner1 = U_fin1,
    Effort_Partner2 = U_fin2,
    Wellbeing_partner1 = V_fin1,
    Wellbeing_partner2 = V_fin2,
    a2 = as.factor(alpha2)
  ))
}

# Zoom (bajado para no tapar rótulos)
x_zoom <- c(1, 2)
tmp_zoom <- subset(results_a2, Feeling >= x_zoom[1] & Feeling <= x_zoom[2])
y_zoom <- range(tmp_zoom$Effort_Partner1, na.rm = TRUE)
pad <- diff(y_zoom) * 0.08; if (!is.finite(pad) || pad == 0) pad <- 0.03
y_zoom <- c(y_zoom[1] - pad, y_zoom[2] + pad)

plot1_E_a_base <- ggplot(results_a2, aes(x = Feeling, y = Effort_Partner1, color = a2, group = a2,
                                         linetype = ifelse(a2 == "0.5", "dashed", "solid"))) +
  geom_line(linewidth = 0.25) +
  geom_hline(yintercept = ustar, linetype = "dashed", color = "black") +
  labs(x = "feeling", y = expression(c[1]), title = "Effort Partner 1") +
  scale_color_manual(values = cols_main, guide = "none") +              # <- paleta azul->rojo
  scale_linetype_identity(guide = "none") +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.25),
                     labels = function(x) ifelse(x == ustar, expression(c[1]^"*"), x),
                     limits = c(0.45, 2)) +
  xlim(0, 5) +
  theme_classic(base_family = "paperfont") + theme(legend.position = "none") + theme_italics

plot_zoom <- ggplot(results_a2, aes(x = Feeling, y = Effort_Partner1, color = a2, group = a2,
                                    linetype = ifelse(a2 == "0.5", "dashed", "solid"))) +
  geom_line(linewidth = 0.25) +
  geom_hline(yintercept = ustar, linetype = "dashed", color = "black") +
  scale_color_manual(values = cols_zoom_dark, guide = "none") +         # <- paleta más fuerte
  scale_linetype_identity(guide = "none") +
  coord_cartesian(xlim = x_zoom, ylim = y_zoom, expand = FALSE) +
  theme_classic(base_family = "paperfont") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.background  = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border     = element_blank())

zoom_x <- 0.63; zoom_y <- 0.52; zoom_w <- 0.40; zoom_h <- 0.40  # ↓ bajado
xlim_main <- c(0, 5); ylim_main <- c(0.45, 2)
x0_data <- 2; y0_data <- 1.25; x1_data <- 2.9; y1_data <- 1.5
x0 <- (x0_data - xlim_main[1]) / diff(xlim_main)
y0 <- (y0_data - ylim_main[1]) / diff(ylim_main)
x1 <- (x1_data - xlim_main[1]) / diff(xlim_main)
y1 <- (y1_data - ylim_main[1]) / diff(ylim_main)

final_plot_a2 <- ggdraw() +
  draw_plot(plot1_E_a_base) +
  draw_plot(plot_zoom, x = zoom_x, y = zoom_y, width = zoom_w, height = zoom_h) +
  draw_line(x = c(x0, x1), y = c(y0, y1),
            arrow = arrow(type = "closed", length = unit(0.16, "cm")),
            linewidth = 0.4, color = "black")

plot1_E_a <- wrap_elements(full = final_plot_a2)

plot2_E_a <- ggplot(results_a2, aes(x = Feeling, y = Effort_Partner2, color = a2, group = a2,
                                    linetype = ifelse(a2 == "0.5", "dashed", "solid"))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = ustar, linetype = "dashed", color = "black") +
  labs(x = "feeling", y = expression(c[2]), title = "Effort Partner 2") +
  scale_color_manual(values = cols_main, name = expression(alpha[2])) + # <- paleta azul->rojo
  scale_linetype_identity(guide = "none") +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.25),
                     labels = function(x) ifelse(x == ustar, expression(c[2]^"*"), x),
                     limits = c(0.45, 2)) +
  xlim(0, 5) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(linetype = "solid", linewidth = 2))) +
  theme_classic(base_family = "paperfont") + theme_italics

plot1_V_a <- ggplot(results_a2, aes(x = Feeling, y = Wellbeing_partner1, color = a2, group = a2,
                                    linetype = ifelse(a2 == "0.5", "dashed", "solid"))) +
  geom_line(linewidth = 0.5) +
  labs(x = "feeling", y = expression(v[1]), title = "Well-being Partner 1") +
  scale_color_manual(values = cols_main, guide = "none") +              # <- paleta azul->rojo
  scale_linetype_identity(guide = "none") +
  xlim(0, 5) + ylim(80, 170) +
  theme_classic(base_family = "paperfont") + theme(legend.position = "none") + theme_italics

plot2_V_a <- ggplot(results_a2, aes(x = Feeling, y = Wellbeing_partner2, color = a2, group = a2,
                                    linetype = ifelse(a2 == "0.5", "dashed", "solid"))) +
  geom_line(linewidth = 1) +
  labs(x = "feeling", y = expression(v[2]), title = "Well-being Partner 2") +
  scale_color_manual(values = cols_main, guide = "none") +              # <- paleta azul->rojo
  scale_linetype_identity(guide = "none") +
  xlim(0, 5) + ylim(80, 170) +
  theme_classic(base_family = "paperfont") + theme(legend.position = "none") + theme_italics

left_grid_a2 <- ((plot1_E_a + plot2_E_a) / (plot1_V_a + plot2_V_a)) + plot_layout(guides = "collect")
left_grid_a2 <- left_grid_a2 &
  theme(legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.title      = element_text(size = 9, face = "italic"),
        legend.text       = element_text(size = 8, face = "italic"),
        legend.key.width  = unit(0.9, "lines"),
        legend.key.height = unit(0.7, "lines"),
        legend.margin     = margin(t = 2, b = -2))

# Scatter lateral (usa Xstar de MODELOB_FIJO_BETA)
load("MODELOB_FIJO_BETA.RData")
dats_xstar_a2 <- data.frame(Xs = Xstar[,5], alpha2 = alpha2_values)

scatter_right_a2 <- ggplot(dats_xstar_a2, aes(x = alpha2, y = Xs)) +
  geom_point(size = 1.3, alpha = 0.7) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = expression(alpha[2]), y = bquote(italic(x)^scriptstyle("\u2665"))) +
  theme_classic(base_family = "paperfont") +
  theme(axis.title = element_text(size = 10, face = "italic"),
        axis.text  = element_text(size = 9, face = "italic"),
        plot.margin = margin(6, 6, 6, 6))

final_with_scatter_a2 <- (left_grid_a2 | scatter_right_a2) + plot_layout(widths = c(9, 3))
save_eps(final_with_scatter_a2, "New_F5.eps", width = 7.2, height = 4.8)




############################################################
# FIGURE 6  ->  New_F6.eps
# MODELOB_FIJO_ALPHA
############################################################


# --- Librerías ---
library(ggplot2)
library(akima)
library(metR)
library(patchwork)
library(scico)
library(scales)

theme_italics <- theme(
  axis.title = element_text(face = "italic"),
  axis.text  = element_text(face = "italic"),
  plot.title = element_text(face = "italic")
)

# --- Datos / modelo ---
load("MODELOB_FIJO_ALPHA.RData")

beta1_values <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
beta2_values <- beta1_values

# --- Helper contornos β–β con 1:1 (heat + contorno + labels) ---
make_contour_plot_beta <- function(z_vals, breaks_contour, fill_lab, main_title) {
  df <- expand.grid(beta1 = beta1_values, beta2 = beta2_values)
  df$z <- as.vector(z_vals); df <- na.omit(df)
  
  ip <- akima::interp(
    x = df$beta1, y = df$beta2, z = df$z,
    xo = seq(min(beta1_values), max(beta1_values), length.out = 220),
    yo = seq(min(beta2_values), max(beta2_values), length.out = 220),
    linear = TRUE, extrap = TRUE, duplicate = "mean"
  )
  
  interp_df <- expand.grid(beta1 = ip$x, beta2 = ip$y)
  interp_df$z <- as.vector(ip$z)
  
  ggplot(interp_df, aes(x = beta1, y = beta2)) +
    geom_raster(aes(fill = z), interpolate = TRUE, na.rm = TRUE, alpha = tile_alpha) +
    geom_contour(aes(z = z), color = "grey30", linewidth = 0.3,
                 breaks = breaks_contour, na.rm = TRUE) +
    metR::geom_text_contour(aes(z = z), breaks = breaks_contour,
                            size = 3.0, stroke = 0,
                            label.placer = metR::label_placer_n(1)) +
    # --- Escala continua azul → rojo elegante (scico "vik") ---
    scale_fill_scico(palette = "vik", direction = 1, oob = scales::squish) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = 0)) +
    coord_equal(expand = FALSE) +
    labs(title = main_title, fill = fill_lab,
         x = expression(beta[1]), y = expression(beta[2])) +
    theme_minimal(base_family = "paperfont") +
    theme(
      legend.position = "none",
      plot.title   = element_text(size = 10, face = "bold", margin = margin(b = 2)),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.x  = element_text(size = 10),
      axis.text.y  = element_text(size = 10),
      plot.margin  = margin(2, 2, 2, 2)
    ) + theme_italics
}

# --- Títulos de columnas (por si no estaba definido antes) ---
col_title <- function(txt) {
  ggplot() + theme_void() + ggtitle(txt) +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.margin = margin(0,0,0,0))
}

# --- Renombrados/alias ---
C1_B1a <- C1_a;  C2_B1a <- C2_a;  V1_B1a <- V1_a;  V2_B1a <- V2_a

c1.B1b <- make_contour_plot_beta(
  C1_B1a, c(0.95,1,1.05),
  expression(c[1]^scriptstyle("\u2665")), "P.1"
)

c2.B1b <- make_contour_plot_beta(
  C2_B1a, c(0.95,1,1.05),
  expression(c[2]^scriptstyle("\u2665")), "P.2"
)

# --- Feeling (Xstar_a) β–β con 1:1 ---
df_xb <- expand.grid(beta1 = beta1_values, beta2 = beta2_values)
df_xb$z <- as.vector(Xstar_a); df_xb <- na.omit(df_xb)

ipxb <- akima::interp(
  x = df_xb$beta1, y = df_xb$beta2, z = df_xb$z,
  xo = seq(min(beta1_values), max(beta1_values), length.out = 220),
  yo = seq(min(beta2_values), max(beta2_values), length.out = 220),
  linear = TRUE, extrap = TRUE, duplicate = "mean"
)

interp_xb <- expand.grid(beta1 = ipxb$x, beta2 = ipxb$y)
interp_xb$z <- as.vector(ipxb$z)

p3b <- ggplot(interp_xb, aes(x = beta1, y = beta2, z = z)) +
  geom_contour_filled(aes(fill = after_stat(level)),
                      bins = 15, show.legend = FALSE,
                      alpha = alpha_fill, na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = c(0.99,1.02,1.05,1.10),
                          size = 3.6, stroke = 0, skip = 0,
                          check_overlap = FALSE,
                          label.placer = metR::label_placer_n(1)) +
  labs(title = NULL, x = expression(beta[1]), y = expression(beta[2])) +
  # --- Escala DISCRETA azul → rojo (scico "vik") para niveles binned ---
  scale_fill_manual(values = scico::scico(15, palette = "vik"),
                    drop = FALSE, guide = "none") +
  coord_equal(expand = FALSE) +
  theme_minimal(base_family = "paperfont") +
  theme(
    legend.position = "none",
    axis.title.x  = element_text(size = 10),
    axis.title.y  = element_text(size = 10),
    axis.text.x   = element_text(size = 10),
    axis.text.y   = element_text(size = 10),
    plot.margin   = margin(2, 2, 2, 2)
  ) + theme_italics +
  annotate("segment", x = 2.46, y = 2.46, xend = 2.10, yend = 2.15,
           arrow = arrow(type = "closed", angle = 20, length = unit(0.15, "inches")),
           linewidth = 0.4, color = "black") +
  annotate("text",
           x = 2.52, y = 2.52,
           label = 'italic(x)^scriptstyle("\u2665")', parse = TRUE,
           size = 5.2, color = "black", hjust = 0.5, vjust = 0.5)

# --- Títulos y layout ---
t1 <- col_title("Equilibrium Effort")
t2 <- col_title("Equilibrium Feeling")
t3 <- col_title("Well-being")

areas_b_sq <- c(
  patchwork::area(t = 1, l = 1, b = 1, r = 1),
  patchwork::area(t = 1, l = 2, b = 1, r = 2),
  patchwork::area(t = 1, l = 3, b = 1, r = 3),
  patchwork::area(t = 2, l = 1, b = 2, r = 1),
  patchwork::area(t = 2, l = 2, b = 2, r = 2),
  patchwork::area(t = 2, l = 3, b = 2, r = 3),
  patchwork::area(t = 3, l = 1, b = 3, r = 1),
  patchwork::area(t = 3, l = 2, b = 3, r = 2),
  patchwork::area(t = 3, l = 3, b = 3, r = 3)
)

# --- Composición final ---
final_plot_6 <-
  t1 + t2 + t3 +
  c1.B1b + p3b + make_contour_plot_beta(V1_B1a, c(113,113.5,114,114.5), expression(v[1]), "P.1") +
  c2.B1b + patchwork::plot_spacer() + make_contour_plot_beta(V2_B1a, c(113,113.5,114,114.5), expression(v[2]), "P.2") +
  plot_layout(design = areas_b_sq, heights = c(0.06, 1, 1), widths  = c(1, 1, 1)) &
  theme(
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 10),
    plot.margin = margin(2, 2, 2, 2)
  )

# --- Export ---
save_eps(final_plot_6, "New_F6.eps", width = 7.2, height = 6.2)



############################################################
# EXTRA  ->  New_F7.eps
# Dinámica con beta2 variable + scatter


############################################################

############################################################
## DATOS
############################################################
beta2_values <- c(1, 1.5, 1.75, 2, 2.5, 2.75, 3)
results_b2 <- data.frame()

for (k in seq_along(beta2_values)) {
  

  alpha1 <- 0.5; alpha2 <- 0.5; beta1 <- 1; beta2 <- beta2_values[k]
  LOVE2D_result <- LOVE2D(Num, ustar, 4000, r, a1, a2, alpha1, alpha2, beta1, beta2, sigma)
  U_fin1 <- LOVE2D_result$k11; U_fin2 <- LOVE2D_result$k22
  V_fin1 <- LOVE2D_result$v0_0[,1]; V_fin2 <- LOVE2D_result$v0_0[,2]
  x1 <- seq(xmin, xmax, length.out = Num)
  results_b2 <- rbind(results_b2, data.frame(
    Feeling = x1,
    Effort_Partner1    = U_fin1,
    Effort_Partner2    = U_fin2,
    Wellbeing_partner1 = V_fin1,
    Wellbeing_partner2 = V_fin2,
    b2 = as.character(beta2)
  ))
}

## --- LIBRERÍAS ---
library(ggplot2)
library(patchwork)  # inset_element, plot_layout
library(cowplot)    # ggdraw, draw_plot

## --- PALETA ÚNICA (β baja = azul, β alta = rojo) + THEME ---
lvl_b2 <- as.character(beta2_values)
results_b2$b2 <- factor(results_b2$b2, levels = lvl_b2, ordered = TRUE)

# 9 colores (azul -> rojo) para panel principal
cols_beta_main <- c(
  "#08306B", "#1858A8", "#2C7FB8", "#4BA3D1", "#74C0E3",
  "#F6A385", "#F16D5B", "#D94841", "#A50F15"
)
# 9 colores más “fuertes” para el zoom
cols_beta_zoom <- c(
  "#05244F", "#0E3E82", "#1E63B5", "#2F7CC1", "#5AAFD6",
  "#F07F61", "#E24D3A", "#C3322B", "#800A0A"
)
cols_beta_main_named <- stats::setNames(cols_beta_main, lvl_b2)
cols_beta_zoom_named <- stats::setNames(cols_beta_zoom, lvl_b2)

theme_italics <- theme(
  axis.title = element_text(face = "italic"),
  axis.text  = element_text(face = "italic"),
  plot.title = element_text(face = "italic")
)

## --- RANGOS ---
EFF_Y_LIMITS <- c(0, 1.28)
WB_MIN <- 105

## -------- Effort Partner 1 (con zoom) --------
x_zoom_E  <- c(1.45, 1.75)
tmp_E     <- subset(results_b2, Feeling >= x_zoom_E[1] & Feeling <= x_zoom_E[2])
y_zoom_Er <- range(tmp_E$Effort_Partner1, na.rm = TRUE)
pad_E     <- max(0.005, 0.02 * diff(y_zoom_Er[is.finite(y_zoom_Er)]))
y_zoom_E  <- c(max(EFF_Y_LIMITS[1], y_zoom_Er[1] - pad_E),
               min(EFF_Y_LIMITS[2], y_zoom_Er[2] + pad_E))
if (!all(is.finite(y_zoom_E)) || diff(y_zoom_E) <= 1e-6) y_zoom_E <- c(0.60, 0.70)

plot1_E_b_base <- ggplot(results_b2, aes(x = Feeling, y = Effort_Partner1,
                                         color = b2, group = b2,
                                         linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = ustar, linetype = "dashed", color = "black") +
  labs(x = "feeling", y = expression(c[1]), title = "Effort Partner 1") +
  scale_color_manual(values = cols_beta_main_named, breaks = lvl_b2,
                     name = expression(beta[2]), guide = "none") +   # sin leyenda aquí
  scale_linetype_identity(guide = "none") +
  scale_x_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = EFF_Y_LIMITS,
                     breaks = seq(0, 1.28, by = 0.25), expand = c(0, 0),
                     labels = function(x) ifelse(x == ustar, expression(c[1]^"*"), x)) +
  theme_classic(base_family = "paperfont") +
  theme(plot.margin = margin(6, 6, 10, 6)) + theme_italics

plot1_E_b_zoom <- ggplot(results_b2, aes(x = Feeling, y = Effort_Partner1,
                                         color = b2, group = b2,
                                         linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = ustar, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_color_manual(values = cols_beta_zoom_named, breaks = lvl_b2, guide = "none") +
  scale_linetype_identity(guide = "none") +
  coord_cartesian(xlim = x_zoom_E, ylim = y_zoom_E, expand = FALSE) +
  theme_classic(base_family = "paperfont") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0))

plot1_E_b <- plot1_E_b_base +
  inset_element(plot1_E_b_zoom, left = 0.51, bottom = 0.00, right = 1.00, top = 0.75, align_to = "panel")

## -------- Effort Partner 2 (deja la ÚNICA leyenda) --------
plot2_E_b <- ggplot(results_b2, aes(x = Feeling, y = Effort_Partner2,
                                    color = b2, group = b2,
                                    linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = ustar, linetype = "dashed", color = "black") +
  labs(x = "feeling", y = expression(c[2]), title = "Effort Partner 2") +
  scale_color_manual(values = cols_beta_main_named, breaks = lvl_b2,
                     name = expression(beta[2])) +  # esta mantiene la leyenda
  scale_linetype_identity(guide = "none") +
  scale_x_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = EFF_Y_LIMITS,
                     breaks = seq(0, 1.28, by = 0.25), expand = c(0, 0),
                     labels = function(x) ifelse(x == ustar, expression(c[2]^"*"), x)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE,
                              override.aes = list(linetype = "solid", linewidth = 3))) +
  theme_classic(base_family = "paperfont") +
  theme(plot.margin = margin(6, 6, 10, 6)) + theme_italics

## -------- Well-being Partner 1 (con zoom) --------
wb1_max <- max(results_b2$Wellbeing_partner1, na.rm = TRUE)
if (!is.finite(wb1_max) || wb1_max <= WB_MIN) wb1_max <- WB_MIN + 0.1

x_zoom_W  <- c(1.45, 1.75)
tmp_W     <- subset(results_b2, Feeling >= x_zoom_W[1] & Feeling <= x_zoom_W[2])
y_zoom_Wr <- range(tmp_W$Wellbeing_partner1, na.rm = TRUE)
pad_W     <- max(0.05, 0.02 * diff(y_zoom_Wr[is.finite(y_zoom_Wr)]))
y_zoom_W  <- c(max(WB_MIN, y_zoom_Wr[1] - pad_W), min(wb1_max, y_zoom_Wr[2] + pad_W))
if (!all(is.finite(y_zoom_W)) || diff(y_zoom_W) <= 1e-6) y_zoom_W <- c(WB_MIN + 0.2, WB_MIN + 0.4)

plot1_V_b_base <- ggplot(results_b2, aes(x = Feeling, y = Wellbeing_partner1,
                                         color = b2, group = b2,
                                         linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.5) +
  labs(x = "feeling", y = expression(v[1]), title = "Well-being Partner 1") +
  scale_color_manual(values = cols_beta_main_named, breaks = lvl_b2,
                     name = expression(beta[2]), guide = "none") +
  scale_linetype_identity(guide = "none") +
  scale_x_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(WB_MIN, wb1_max), expand = c(0, 0)) +
  theme_classic(base_family = "paperfont") +
  theme(plot.margin = margin(6, 6, 12, 6)) + theme_italics

plot1_V_b_zoom <- ggplot(results_b2, aes(x = Feeling, y = Wellbeing_partner1,
                                         color = b2, group = b2,
                                         linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = cols_beta_zoom_named, breaks = lvl_b2, guide = "none") +
  scale_linetype_identity(guide = "none") +
  coord_cartesian(xlim = x_zoom_W, ylim = y_zoom_W, expand = FALSE) +
  theme_classic(base_family = "paperfont") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0))

plot1_V_b <- plot1_V_b_base 
  

## -------- Well-being Partner 2 (con zoom) --------
wb2_max <- max(results_b2$Wellbeing_partner2, na.rm = TRUE)
if (!is.finite(wb2_max) || wb2_max <= WB_MIN) wb2_max <- WB_MIN + 0.1

x_zoom_W2  <- c(1.45, 1.75)
tmp_W2     <- subset(results_b2, Feeling >= x_zoom_W2[1] & Feeling <= x_zoom_W2[2])
y_zoom_W2r <- range(tmp_W2$Wellbeing_partner2, na.rm = TRUE)
pad_W2     <- max(0.05, 0.02 * diff(y_zoom_W2r[is.finite(y_zoom_W2r)]))
y_zoom_W2  <- c(max(WB_MIN, y_zoom_W2r[1] - pad_W2), min(wb2_max, y_zoom_W2r[2] + pad_W2))
if (!all(is.finite(y_zoom_W2)) || diff(y_zoom_W2) <= 1e-6) y_zoom_W2 <- c(WB_MIN + 0.2, WB_MIN + 0.4)

plot2_V_b_base <- ggplot(results_b2, aes(x = Feeling, y = Wellbeing_partner2,
                                         color = b2, group = b2,
                                         linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.5) +
  labs(x = "feeling", y = expression(v[2]), title = "Well-being Partner 2") +
  scale_color_manual(values = cols_beta_main_named, breaks = lvl_b2,
                     name = expression(beta[2]), guide = "none") +
  scale_linetype_identity(guide = "none") +
  scale_x_continuous(limits = c(0, 5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(WB_MIN, wb2_max), expand = c(0, 0)) +
  theme_classic(base_family = "paperfont") +
  theme(plot.margin = margin(6, 6, 12, 6)) + theme_italics

plot2_V_b_zoom <- ggplot(results_b2, aes(x = Feeling, y = Wellbeing_partner2,
                                         color = b2, group = b2,
                                         linetype = ifelse(b2 == "1","dashed","solid"))) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = cols_beta_zoom_named, breaks = lvl_b2, guide = "none") +
  scale_linetype_identity(guide = "none") +
  coord_cartesian(xlim = x_zoom_W2, ylim = y_zoom_W2, expand = FALSE) +
  theme_classic(base_family = "paperfont") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0))

plot2_V_b <- plot2_V_b_base +
  inset_element(plot2_V_b_zoom, left = 0.59, bottom = 0.05, right = 1.05, top = 0.59, align_to = "panel")

## -------- GRID 2x2 con UNA sola leyenda (collect) --------
left_grid_b2 <- ((plot1_E_b + plot2_E_b) / (plot1_V_b + plot2_V_b)) + plot_layout(guides = "collect")

left_grid_b2 <- left_grid_b2 &
  theme(legend.position   = "bottom",
        legend.direction  = "horizontal",
        legend.title      = element_text(size = 9, face = "italic"),
        legend.text       = element_text(size = 8, face = "italic"),
        legend.key.width  = unit(0.9, "lines"),
        legend.key.height = unit(0.7, "lines"),
        legend.margin     = margin(t = 2, b = -2))

## -------- SCATTER y EXPORT (ajusta si cambiaste Xstar) --------
# (mantengo tu bloque de 5 puntos; si tienes 9 Xstar, dime cómo indexarlos)
load("MODELOB_FIJO_ALPHA.RData")
Xstar_new <- c(Xstar_a[1,1], Xstar_a[1,6], Xstar_a[1,11], Xstar_a[1,16], Xstar_a[1,21])
beta2_scatter <- c(1, 1.5, 2.0, 2.5, 3.0)
dats_xstar_b2 <- data.frame(Xs = Xstar_new, beta2 = beta2_scatter)

scatter_right_b2 <- ggplot(dats_xstar_b2, aes(x = beta2, y = Xs)) +
  geom_point(size = 1.6, alpha = 0.8) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(x = expression(beta[2]), y = bquote(italic(x)^scriptstyle("\u2665"))) +
  theme_classic(base_family = "paperfont") +
  theme(axis.title = element_text(size = 10, face = "italic"),
        axis.text  = element_text(size = 9, face = "italic"),
        plot.margin = margin(6, 6, 6, 6))

final_with_scatter_b2 <- (left_grid_b2 | scatter_right_b2) + plot_layout(widths = c(9, 3))
save_eps(final_with_scatter_b2, "New_F7.eps", width = 7.2, height = 4.8)


############################################################
# PLOT 8  ->  New_F8.eps
# Tres columnas (alpha): Effort diff, Well-being diff, Feeling (α–α con 1:1)
############################################################
# --- Datos del modelo ---
load("MODELOB_FIJO_BETA.RData")

alpha1_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
alpha2_values <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# --- Librerías usadas explícitamente ---
library(ggplot2)
library(grid)    # para grid::unit (flecha en p3)

# Opacidad del relleno en mapas
alpha_fill <- 0.92

# Rangos (por si los necesitas en otros ajustes)
x_rng <- range(alpha1_values); y_rng <- range(alpha2_values)

# --- Interpolaciones ---
d1 <- expand.grid(alpha1 = alpha1_values, alpha2 = alpha2_values)
d1$z <- as.vector(C1 - C2); d1 <- na.omit(d1)
ip1 <- akima::interp(d1$alpha1, d1$alpha2, d1$z,
                     xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
                     yo = seq(min(alpha2_values), max(alpha2_values), length.out = 220),
                     extrap = TRUE, duplicate = "mean")
df1 <- expand.grid(alpha1 = ip1$x, alpha2 = ip1$y); df1$z <- as.vector(ip1$z)

d2 <- expand.grid(alpha1 = alpha1_values, alpha2 = alpha2_values)
d2$z <- as.vector(V1 - V2); d2 <- na.omit(d2)
ip2 <- akima::interp(d2$alpha1, d2$alpha2, d2$z,
                     xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
                     yo = seq(min(alpha2_values), max(alpha2_values), length.out = 220),
                     extrap = TRUE, duplicate = "mean")
df2 <- expand.grid(alpha1 = ip2$x, alpha2 = ip2$y); df2$z <- as.vector(ip2$z)

d3 <- expand.grid(alpha1 = alpha1_values, alpha2 = alpha2_values)
d3$z <- as.vector(Xstar); d3 <- na.omit(d3)
ip3 <- akima::interp(d3$alpha1, d3$alpha2, d3$z,
                     xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
                     yo = seq(min(alpha2_values), max(alpha2_values), length.out = 220),
                     extrap = TRUE, duplicate = "mean")
df3 <- expand.grid(alpha1 = ip3$x, alpha2 = ip3$y); df3$z <- as.vector(ip3$z)

# --- Paletas y tema ---
pos_col <- "#D73027"   # positivo (1>2)
neg_col <- "#1858A8"   # negativo (2>1)
neu_col <- "#F0F0F0"   # neutro

theme_compacto <- theme_minimal(base_family = "paperfont") +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12, face = "italic"),
        axis.text.x  = element_text(size = 10, face = "italic"),
        axis.text.y  = element_text(size = 10, face = "italic"),
        plot.margin  = margin(4, 4, 4, 4))

# IMPORTANTE: no tocar el título aquí; solo cursiva en ejes
theme_italics <- theme(
  axis.title = element_text(face = "italic"),
  axis.text  = element_text(face = "italic")
)

# --- Helper para situar rótulos dentro de la región correcta ---
# sign_target = +1 para z>0 (">"),  -1 para z<0 ("<")
# xf,yf son fracciones del rectángulo [0..1] (0.25 ~ 25% desde la esquina inferior-izquierda)
pick_in_region <- function(df, sign_target = 1, xf = 0.75, yf = 0.10) {
  xr <- range(df$alpha1); yr <- range(df$alpha2)
  x0 <- xr[1] + xf * diff(xr);  y0 <- yr[1] + yf * diff(yr)
  d2 <- (df$alpha1 - x0)^2 + (df$alpha2 - y0)^2
  s  <- sign(df$z)
  idx <- which(s == sign_target)
  if (length(idx) == 0) {  # fallback si no hay estrictamente >0 o <0
    idx <- which(abs(df$z) > 0)
  }
  j <- idx[which.min(d2[idx])]
  c(x = df$alpha1[j], y = df$alpha2[j])
}

# ========= p1: Δ Effort (C1 - C2) =========
z1_rng <- range(df1$z, na.rm = TRUE)
max_abs1 <- max(abs(z1_rng))
breaks1  <- pretty(c(-max_abs1, max_abs1), n = 14)
if (!any(abs(breaks1) < .Machine$double.eps)) breaks1 <- sort(c(breaks1, 0))
pal1 <- grDevices::colorRampPalette(c(neg_col, neu_col, pos_col))(length(breaks1) - 1)

# posiciones ~25% de cada esquina, ajustadas al signo de la región
p1_pos_lt <- pick_in_region(df1, sign_target = -1, xf = 0.10, yf = 0.75)  # c1 < c2
p1_pos_rb <- pick_in_region(df1, sign_target = +1, xf = 0.75, yf = 0.10)  # c1 > c2

p1 <- ggplot(df1, aes(alpha1, alpha2, z = z)) +
  geom_contour_filled(breaks = breaks1, show.legend = FALSE, alpha = alpha_fill) +
  # contorno 0 bien marcado
  geom_contour(aes(z = z), breaks = 0, color = "white", linewidth = 1.6) +
  geom_contour(aes(z = z), breaks = 0, color = "black", linewidth = 0.9) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = breaks1, min.size = 1, skip = 0, rotate = TRUE,
                          size = 4.2, color = "black",
                          label.placer = metR::label_placer_n(1)) +
  scale_fill_manual(values = pal1, drop = FALSE) +
  labs(title = "Difference in Efforts", x = expression(alpha[1]), y = expression(alpha[2])) +
  coord_equal(expand = FALSE) +
  # Rótulos dentro de cada región (sin fondo, ligera inclinación)
  annotate("text",
           x = p1_pos_lt["x"], y = p1_pos_lt["y"],
           label = 'c[1]^scriptstyle("\\u2665") < c[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5) +
  annotate("text",
           x = p1_pos_rb["x"], y = p1_pos_rb["y"],
           label = 'c[1]^scriptstyle("\\u2665") > c[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5) +
  theme_compacto

# ========= p2: Δ Well-being (V1 - V2) =========
z2_rng <- range(df2$z, na.rm = TRUE)
max_abs2 <- max(abs(z2_rng))
breaks2  <- pretty(c(-max_abs2, max_abs2), n = 14)
if (!any(abs(breaks2) < .Machine$double.eps)) breaks2 <- sort(c(breaks2, 0))
pal2 <- grDevices::colorRampPalette(c(neg_col, neu_col, pos_col))(length(breaks2) - 1)

p2_pos_lt <- pick_in_region(df2, sign_target = -1, xf = 0.10, yf = 0.75)  # v1 < v2
p2_pos_rb <- pick_in_region(df2, sign_target = +1, xf = 0.75, yf = 0.10)  # v1 > v2

p2 <- ggplot(df2, aes(alpha1, alpha2, z = z)) +
  geom_contour_filled(breaks = breaks2, show.legend = FALSE, alpha = alpha_fill) +
  geom_contour(aes(z = z), breaks = 0, color = "white", linewidth = 1.6) +
  geom_contour(aes(z = z), breaks = 0, color = "black", linewidth = 0.9) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = breaks2, min.size = 1, skip = 0, rotate = TRUE,
                          size = 3.8, color = "black",
                          label.placer = metR::label_placer_n(1)) +
  scale_fill_manual(values = pal2, drop = FALSE) +
  labs(title = "Difference in Well-being", x = expression(alpha[1]), y = expression(alpha[2])) +
  coord_equal(expand = FALSE) +
  annotate("text",
           x = p2_pos_lt["x"], y = p2_pos_lt["y"],
           label = 'v[1]^scriptstyle("\\u2665") < v[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5) +
  annotate("text",
           x = p2_pos_rb["x"], y = p2_pos_rb["y"],
           label = 'v[1]^scriptstyle("\\u2665") > v[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5) +
  theme_compacto

# ========= p3: Feeling (azul -> rojo como antes) =========
p3 <- ggplot(df3, aes(x = alpha1, y = alpha2, z = z)) +
  geom_contour_filled(aes(fill = after_stat(level)),
                      bins = 15, show.legend = FALSE,
                      alpha = alpha_fill, na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = c(0.7, 0.85, 1.5, 1.75, 1.9),
                          size = 3.6, stroke = 0, skip = 0,
                          check_overlap = FALSE,
                          label.placer = metR::label_placer_n(1)) +
  labs(title = "Equilibrium Feeling", x = expression(alpha[1]), y = expression(alpha[2])) +
  scale_fill_manual(values = scico::scico(15, palette = "vik"),
                    drop = FALSE, guide = "none") +
  coord_equal(expand = FALSE) +
  theme_compacto +   # mantiene el título en negrita como en p1 y p2
  theme_italics +    # cursiva en ejes, NO cambia el título
  annotate("segment", x = 0.45, y = 0.45, xend = 0.65, yend = 0.60,
           arrow = arrow(type = "closed", angle = 20, length = grid::unit(0.15, "inches")),
           linewidth = 0.6, color = "black") +
  annotate("text",
           x = 0.40, y = 0.40,
           label = 'italic(x)^scriptstyle("\\u2665")', parse = TRUE,
           size = 5.2, color = "black", hjust = 0.5, vjust = 0.5)

# --- Export ---
final_row_alpha <- cowplot::plot_grid(p1, p2, p3, ncol = 3, align = "hv", rel_widths = c(1,1,1))
save_eps(final_row_alpha, "New_F8.eps", width = 14.0, height = 4.8)





############################################################
# FIGURA 9  ->  New_F9.eps
# Tres columnas (beta): Effort diff, Well-being diff, Feeling (β–β con 1:1)
############################################################
## --- LIBRERÍAS ---
library(ggplot2)
library(akima)
library(metR)
library(scico)
library(cowplot)
library(scales)

## --- DATOS β–β ---
load("ModeloB_FIJO_ALPHA.RData")   # cambia el nombre si tu .RData es otro


## --- DATOS β–β ---

## --- DATOS β–β ---
# --- Datos ---
load("ModeloB_FIJO_ALPHA.RData")   # según tu snippet inicial; si fuera otro .RData, cámbialo aquí

beta1_values <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
                  2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
beta2_values <- beta1_values
df_beta <- expand.grid(beta1 = beta1_values, beta2 = beta2_values)

## --- Interpolaciones (C1_a - C2_a ; V1_a - V2_a ; Xstar_a) ---
db1 <- df_beta; db1$z <- as.vector(C1_a - C2_a); db1 <- subset(db1, is.finite(z))
ipb1 <- akima::interp(db1$beta1, db1$beta2, db1$z,
                      xo = seq(min(beta1_values), max(beta1_values), length.out = 240),
                      yo = seq(min(beta2_values), max(beta2_values), length.out = 240),
                      extrap = TRUE, duplicate = "mean")
dfb1 <- expand.grid(beta1 = ipb1$x, beta2 = ipb1$y); dfb1$z <- as.vector(ipb1$z)

db2 <- df_beta; db2$z <- as.vector(V1_a - V2_a); db2 <- subset(db2, is.finite(z))
ipb2 <- akima::interp(db2$beta1, db2$beta2, db2$z,
                      xo = seq(min(beta1_values), max(beta1_values), length.out = 240),
                      yo = seq(min(beta2_values), max(beta2_values), length.out = 240),
                      extrap = TRUE, duplicate = "mean")
dfb2 <- expand.grid(beta1 = ipb2$x, beta2 = ipb2$y); dfb2$z <- as.vector(ipb2$z)

db3 <- df_beta; db3$z <- as.vector(Xstar_a); db3 <- subset(db3, is.finite(z))
ipb3 <- akima::interp(db3$beta1, db3$beta2, db3$z,
                      xo = seq(min(beta1_values), max(beta1_values), length.out = 240),
                      yo = seq(min(beta2_values), max(beta2_values), length.out = 240),
                      extrap = TRUE, duplicate = "mean")
dfb3 <- expand.grid(beta1 = ipb3$x, beta2 = ipb3$y); dfb3$z <- as.vector(ipb3$z)

## --- Paletas y tema ---
pos_col <- "#D73027"   # positivo (1>2)
neg_col <- "#1858A8"   # negativo (2>1)
neu_col <- "#F0F0F0"   # neutro

if (!exists("alpha_fill")) alpha_fill <- 0.95  # transparencia por defecto si no definida

theme_compacto <- theme_minimal(base_family = "paperfont") +
  theme(legend.position = "none",
        plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12, face = "italic"),
        axis.text.x  = element_text(size = 10, face = "italic"),
        axis.text.y  = element_text(size = 10, face = "italic"),
        plot.margin  = margin(4, 4, 4, 4))

# IMPORTANTE: no tocar plot.title aquí para no sobreescribir la negrita del tema principal
theme_italics <- theme(
  axis.title = element_text(face = "italic"),
  axis.text  = element_text(face = "italic")
)

## ========= p1b: Δ Effort (C1_a - C2_a) =========
z1_rng   <- range(dfb1$z, na.rm = TRUE)
max_abs1 <- max(abs(z1_rng))
breaks1  <- pretty(c(-max_abs1, max_abs1), n = 14)
if (!any(abs(breaks1) < .Machine$double.eps)) breaks1 <- sort(c(breaks1, 0))
pal1     <- grDevices::colorRampPalette(c(neg_col, neu_col, pos_col))(length(breaks1) - 1)

p1b <- ggplot(dfb1, aes(x = beta1, y = beta2, z = z)) +
  geom_contour_filled(breaks = breaks1, show.legend = FALSE, alpha = alpha_fill) +
  # Contorno 0 muy marcado
  geom_contour(aes(z = z), breaks = 0, color = "white", linewidth = 1.6, na.rm = TRUE) +
  geom_contour(aes(z = z), breaks = 0, color = "black", linewidth = 0.9,  na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = breaks1, min.size = 1, skip = 0, rotate = TRUE,
                          size = 4.2, color = "black",
                          label.placer = metR::label_placer_n(1)) +
  scale_fill_manual(values = pal1, drop = FALSE) +
  labs(title = "Difference in Efforts", x = expression(beta[1]), y = expression(beta[2])) +
  coord_equal(expand = FALSE) +
  annotate("text",
           x = max(beta1_values) - 0.10 * diff(range(beta1_values)),
           y = min(beta2_values) + 0.15 * diff(range(beta2_values)),
           label = expression(c[1]^scriptstyle("\u2665") < c[2]^scriptstyle("\u2665")),
           size = 6, hjust = 1, vjust = 0, fontface = "bold", color = "black") +
  annotate("text",
           x = max(beta1_values) - 0.90 * diff(range(beta1_values)),
           y = min(beta2_values) + 0.90 * diff(range(beta2_values)),
           label = expression(c[1]^scriptstyle("\u2665") > c[2]^scriptstyle("\u2665")),
           size = 6, hjust = 0, vjust = 1, fontface = "bold", color = "black") +
  theme_compacto

## ========= p2b: Δ Well-being (V1_a - V2_a) =========
z2_rng   <- range(dfb2$z, na.rm = TRUE)
max_abs2 <- max(abs(z2_rng))
breaks2  <- pretty(c(-max_abs2, max_abs2), n = 14)
if (!any(abs(breaks2) < .Machine$double.eps)) breaks2 <- sort(c(breaks2, 0))
pal2     <- grDevices::colorRampPalette(c(neg_col, neu_col, pos_col))(length(breaks2) - 1)

p2b <- ggplot(dfb2, aes(x = beta1, y = beta2, z = z)) +
  geom_contour_filled(breaks = breaks2, show.legend = FALSE, alpha = alpha_fill) +
  geom_contour(aes(z = z), breaks = 0, color = "white", linewidth = 1.6, na.rm = TRUE) +
  geom_contour(aes(z = z), breaks = 0, color = "black", linewidth = 0.9,  na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = breaks2, min.size = 1, skip = 0, rotate = TRUE,
                          size = 3.8, color = "black",
                          label.placer = metR::label_placer_n(1)) +
  scale_fill_manual(values = pal2, drop = FALSE) +
  labs(title = "Difference in Well-being", x = expression(beta[1]), y = expression(beta[2])) +
  coord_equal(expand = FALSE) +
  annotate("text",
           x = max(beta1_values) - 0.10 * diff(range(beta1_values)),
           y = min(beta2_values) + 0.10 * diff(range(beta2_values)),
           label = expression(v[1]^scriptstyle("\u2665") > v[2]^scriptstyle("\u2665")),
           size = 6, hjust = 1, vjust = 0, fontface = "bold", color = "black") +
  annotate("text",
           x = max(beta1_values) - 0.85 * diff(range(beta1_values)),
           y = min(beta2_values) + 0.85 * diff(range(beta2_values)),
           label = expression(v[1]^scriptstyle("\u2665") < v[2]^scriptstyle("\u2665")),
           size = 6, hjust = 0, vjust = 1, fontface = "bold", color = "black") +
  theme_compacto

## ========= p3b2: Feeling (azul → rojo con scico "vik") =========
p3b2 <- ggplot(dfb3, aes(x = beta1, y = beta2, z = z)) +
  geom_contour_filled(aes(fill = after_stat(level)),
                      bins = 15, show.legend = FALSE,
                      alpha = alpha_fill, na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = c(0.99, 1.02, 1.05, 1.10),
                          size = 3.6, stroke = 0, skip = 0,
                          check_overlap = FALSE,
                          label.placer = metR::label_placer_n(1)) +
  labs(title = "Equilibrium Feeling", x = expression(beta[1]), y = expression(beta[2])) +
  scale_fill_manual(values = scico::scico(15, palette = "vik"),
                    drop = FALSE, guide = "none") +
  coord_equal(expand = FALSE) +
  theme_compacto +   # mantiene el título en negrita como en p1b y p2b
  theme_italics +    # cursiva en ejes, NO cambia el título
  annotate("segment", x = 2.46, y = 2.46, xend = 2.10, yend = 2.15,
           arrow = arrow(type = "closed", angle = 20, length = grid::unit(0.15, "inches")),
           linewidth = 0.4, color = "black") +
  annotate("text",
           x = 2.52, y = 2.52,
           label = 'italic(x)^scriptstyle("\\u2665")', parse = TRUE,
           size = 5.2, color = "black", hjust = 0.5, vjust = 0.5)

## --- Export ---
final_row_beta <- cowplot::plot_grid(p1b, p2b, p3b2, ncol = 3, align = "hv", rel_widths = c(1,1,1))
save_eps(final_row_beta, "New_F9.eps", width = 14.0, height = 4.8)



############################################################
# FIGURA 10  ->  New_F10.eps
# Joint (alpha1, beta1) (ejes distintos, NO 1:1)
############################################################

## --- LIBRERÍAS ---
library(ggplot2)
library(akima)
library(metR)
library(scico)
library(cowplot)
library(scales)

## --- DATOS JOINT α1–β1 ---
## --- DATOS JOINT α1–β1 ---

## --- DATOS JOINT α1–β1 ---
load("MODELOB_JOINT.RData")

alpha1_values <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta1_values  <- c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3)
df_joint <- expand.grid(alpha1 = alpha1_values, beta1 = beta1_values)

## --- LIBRERÍAS ---
library(ggplot2)
library(grid)      # para grid::unit en la flecha
# metR / scico / cowplot se usan con ::
# akima se llama como akima::interp

## --- PARÁMETROS DE DIBUJO ---
alpha_fill <- 0.92  # opacidad de los rellenos

## --- Interpolaciones (C1_b - C2_b ; V1_b - V2_b ; Xstar_b) ---
jj1 <- df_joint; jj1$z <- as.vector(C1_b - C2_b); jj1 <- subset(jj1, is.finite(z))
ipj1 <- akima::interp(jj1$alpha1, jj1$beta1, jj1$z,
                      xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
                      yo = seq(min(beta1_values),  max(beta1_values),  length.out = 260),
                      extrap = TRUE, duplicate = "mean")
dfj1 <- expand.grid(alpha1 = ipj1$x, beta1 = ipj1$y); dfj1$z <- as.vector(ipj1$z)

jj2 <- df_joint; jj2$z <- as.vector(V1_b - V2_b); jj2 <- subset(jj2, is.finite(z))
ipj2 <- akima::interp(jj2$alpha1, jj2$beta1, jj2$z,
                      xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
                      yo = seq(min(beta1_values),  max(beta1_values),  length.out = 260),
                      extrap = TRUE, duplicate = "mean")
dfj2 <- expand.grid(alpha1 = ipj2$x, beta1 = ipj2$y); dfj2$z <- as.vector(ipj2$z)

jj3 <- df_joint; jj3$z <- as.vector(Xstar_b); jj3 <- subset(jj3, is.finite(z))
ipj3 <- akima::interp(jj3$alpha1, jj3$beta1, jj3$z,
                      xo = seq(min(alpha1_values), max(alpha1_values), length.out = 220),
                      yo = seq(min(beta1_values),  max(beta1_values),  length.out = 260),
                      extrap = TRUE, duplicate = "mean")
dfj3 <- expand.grid(alpha1 = ipj3$x, beta1 = ipj3$y); dfj3$z <- as.vector(ipj3$z)

## --- Paletas y tema (mismo estilo que el script anterior) ---
pos_col <- "#D73027"   # positivo (1>2)
neg_col <- "#1858A8"   # negativo (2>1)
neu_col <- "#F0F0F0"   # neutro

theme_compacto <- theme_minimal(base_family = "paperfont") +
  theme(legend.position = "none",
        plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 12, face = "italic"),
        axis.title.y = element_text(size = 12, face = "italic"),
        axis.text.x  = element_text(size = 10, face = "italic"),
        axis.text.y  = element_text(size = 10, face = "italic"),
        plot.margin  = margin(4, 4, 4, 4))

# Panel cuadrado (misma relación de aspecto en los tres gráficos)
theme_panel <- theme_compacto + theme(aspect.ratio = 1)

# Sólo cursiva en ejes; no toca títulos
theme_italics <- theme(
  axis.title = element_text(face = "italic"),
  axis.text  = element_text(face = "italic")
)

## --- Helper: elegir un punto dentro de la región con signo deseado ---
## sign_target = +1 (z>0) para ">",  -1 (z<0) para "<"
## xf,yf son fracciones 0..1 del rectángulo (0.15 ~ 15% desde el borde inferior/izquierdo)
pick_in_region_ab <- function(df, sign_target = 1, xf = 0.75, yf = 0.15) {
  xr <- range(df$alpha1); yr <- range(df$beta1)
  x0 <- xr[1] + xf * diff(xr);  y0 <- yr[1] + yf * diff(yr)
  d2 <- (df$alpha1 - x0)^2 + (df$beta1 - y0)^2
  s  <- sign(df$z)
  idx <- which(s == sign_target)
  if (length(idx) == 0) idx <- which(abs(df$z) > 0)  # fallback
  j <- idx[which.min(d2[idx])]
  c(x = df$alpha1[j], y = df$beta1[j])
}


## ========= p1j: Δ Effort (C1_b - C2_b) =========
z1_rng   <- range(dfj1$z, na.rm = TRUE)
max_abs1 <- max(abs(z1_rng))
breaks1  <- pretty(c(-max_abs1, max_abs1), n = 14)
if (!any(abs(breaks1) < .Machine$double.eps)) breaks1 <- sort(c(breaks1, 0))
pal1     <- grDevices::colorRampPalette(c(neg_col, neu_col, pos_col))(length(breaks1) - 1)

## ← añade este límite para las ETIQUETAS (no para el relleno):
LABEL_MAX <- 0.6
label_breaks1 <- breaks1[breaks1 >= -LABEL_MAX & breaks1 <= LABEL_MAX]
if (length(label_breaks1) == 0) label_breaks1 <- c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)

# Posiciones interiores (~15–25% desde esquinas) según signo
p1j_pos_neg <- pick_in_region_ab(dfj1, sign_target = -1, xf = 0.15, yf = 0.75)  # c1 < c2
p1j_pos_pos <- pick_in_region_ab(dfj1, sign_target = +1, xf = 0.75, yf = 0.15)  # c1 > c2

p1j <- ggplot(dfj1, aes(x = alpha1, y = beta1, z = z)) +
  geom_contour_filled(breaks = breaks1, show.legend = FALSE, alpha = alpha_fill) +
  geom_contour(aes(z = z), breaks = 0, color = "white", linewidth = 1.6, na.rm = TRUE) +
  geom_contour(aes(z = z), breaks = 0, color = "black", linewidth = 0.9,  na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = label_breaks1, min.size = 1, skip = 0, rotate = TRUE,
                          size = 4.2, color = "black",
                          label.placer = metR::label_placer_n(1)) +
  scale_fill_manual(values = pal1, drop = FALSE) +
  labs(title = "Difference in Efforts", x = expression(alpha[1]), y = expression(beta[1])) +
  # Panel cuadrado (sin coord_equal, que te comprimía)
  theme_panel +
  # Rótulos interiores (sin fondo, ligera inclinación)
  annotate("text",
           x = p1j_pos_neg["x"], y = p1j_pos_neg["y"],
           label = 'c[1]^scriptstyle("\\u2665") < c[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5) +
  annotate("text",
           x = p1j_pos_pos["x"], y = p1j_pos_pos["y"],
           label = 'c[1]^scriptstyle("\\u2665") > c[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5)

## ========= p2j: Δ Well-being (V1_b - V2_b) =========
z2_rng   <- range(dfj2$z, na.rm = TRUE)
max_abs2 <- max(abs(z2_rng))
breaks2  <- pretty(c(-max_abs2, max_abs2), n = 14)
if (!any(abs(breaks2) < .Machine$double.eps)) breaks2 <- sort(c(breaks2, 0))
pal2     <- grDevices::colorRampPalette(c(neg_col, neu_col, pos_col))(length(breaks2) - 1)

p2j_pos_neg <- pick_in_region_ab(dfj2, sign_target = -1, xf = 0.15, yf = 0.75)  # v1 < v2
p2j_pos_pos <- pick_in_region_ab(dfj2, sign_target = +1, xf = 0.75, yf = 0.15)  # v1 > v2

p2j <- ggplot(dfj2, aes(x = alpha1, y = beta1, z = z)) +
  geom_contour_filled(breaks = breaks2, show.legend = FALSE, alpha = alpha_fill) +
  geom_contour(aes(z = z), breaks = 0, color = "white", linewidth = 1.6, na.rm = TRUE) +
  geom_contour(aes(z = z), breaks = 0, color = "black", linewidth = 0.9,  na.rm = TRUE) +
  metR::geom_text_contour(aes(label = after_stat(level), z = z),
                          breaks = breaks2, min.size = 1, skip = 0, rotate = TRUE,
                          size = 3.8, color = "black",
                          label.placer = metR::label_placer_n(1)) +
  scale_fill_manual(values = pal2, drop = FALSE) +
  labs(title = "Difference in Well-being", x = expression(alpha[1]), y = expression(beta[1])) +
  theme_panel +
  annotate("text",
           x = p2j_pos_neg["x"], y = p2j_pos_neg["y"],
           label = 'v[1]^scriptstyle("\\u2665") < v[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5) +
  annotate("text",
           x = p2j_pos_pos["x"], y = p2j_pos_pos["y"],
           label = 'v[1]^scriptstyle("\\u2665") > v[2]^scriptstyle("\\u2665")',
           parse = TRUE, size = 5.2, color = "black", angle = 6,
           hjust = 0.5, vjust = 0.5)

## ========= p3j: Feeling (azul → rojo con scico "vik") =========
p3j <- ggplot(dfj3, aes(x = alpha1, y = beta1, z = z)) +
  geom_contour_filled(aes(fill = after_stat(level)), bins = 15,
                      show.legend = FALSE, alpha = alpha_fill, na.rm = TRUE) +
  metR::geom_text_contour(
    aes(label = after_stat(level), z = z),
    breaks = c(0.99, 1.10, 1.20, 1.30),
    min.size = 1, skip = 0, rotate = TRUE, size = 3.6, color = "black",
    label.placer = metR::label_placer_n(1)
  ) +
  labs(title = "Equilibrium Feeling", x = expression(alpha[1]), y = expression(beta[1])) +
  scale_fill_manual(values = scico::scico(15, palette = "vik"),
                    drop = FALSE, guide = "none") +
  theme_panel +     # panel cuadrado
  theme_italics +   # cursiva sólo en ejes
  annotate("segment", x = 0.31, y = 2.20, xend = 0.60, yend = 1.75,
           arrow = arrow(type = "closed", angle = 20, length = grid::unit(0.15, "inches")),
           linewidth = 0.6, color = "black") +
  annotate("text", x = 0.30, y = 2.25,
           label = 'italic(x)^scriptstyle("\\u2665")', parse = TRUE,
           size = 5.2, color = "black", hjust = 0.5, vjust = 0.5)

## --- Export ---
final_row_joint <- cowplot::plot_grid(p1j, p2j, p3j, ncol = 3, align = "hv", rel_widths = c(1,1,1))
save_eps(final_row_joint, "New_F10.eps", width = 14.0, height = 4.8)
