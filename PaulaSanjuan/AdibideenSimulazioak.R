#BERNOULLI
library(scales)

set.seed(123)
n <- 20
p <- 0.4
num_sim <- 10
base_colors <- rainbow(num_sim)

ruta <- "C:/Users/paula/OneDrive/Documentos/bernoulli04.png"

png(filename = ruta, width = 1200, height = 900, res = 150)

ylim_max <- 2

par(pty="s") #Grafiko karratua

plot(0:n, rep(NA, n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Bernoulli"~p==.(p)),
     xlim=c(0,8), ylim=c(0, ylim_max),
     axes = FALSE)

axis(1)
axis(2, at = seq(0, ylim_max, 1), labels = seq(0, ylim_max, 1), las = 1)
box()

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rbinom(Z[i-1], 1, p)))
  }
  
  jitter_offset <- runif(1, -0.02, 0.02)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()


#GEOMETRIKOA1

library(scales)  

set.seed(123)
m <- 2
n <- 20
p <- 0.4
num_sim <- 10
base_colors <- rainbow(num_sim)

ruta_archivo <- "C:/Users/paula/OneDrive/Documentos/ofigeom04.png"

png(filename = ruta_archivo, width = 1200, height = 900, res = 150)

par(pty="s")  #grafiko karratua

plot(0:n, rep(NA,n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Geometrikoa"~p==.(p)),
     xlim=c(0,15), ylim=c(0,10))

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rgeom(Z[i-1], p)))
  }
  
  jitter_offset <- runif(1, -0.05, 0.05)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()




#GEOMETRIKOA2

library(scales)  

set.seed(123)
m <- 2
n <- 20
p <- 0.7
num_sim <- 10
base_colors <- rainbow(num_sim)

ruta_archivo <- "C:/Users/paula/OneDrive/Documentos/ofigeom07.png"

png(filename = ruta_archivo, width = 1200, height = 900, res = 150)  

par(pty="s")  #Grafiko karratua

plot(0:n, rep(NA,n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Geometrikoa"~p==.(p)),
     xlim=c(0,12), ylim=c(0,4))

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rgeom(Z[i-1], p)))
  }
  
  jitter_offset <- runif(1, -0.05, 0.05)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()




#BINOMIAL1

library(scales)  

set.seed(123)
m <- 2
n <- 20
p <- 0.4
num_sim <- 10
base_colors <- rainbow(num_sim)


ruta_archivo <- "C:/Users/paula/OneDrive/Documentos/binomial04.png"

png(filename = ruta_archivo, width = 1200, height = 900, res = 150)

par(pty="s")  #Grafiko karratua

plot(0:n, rep(NA,n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Binomiala"~p==.(p)),
     xlim=c(0,18), ylim=c(0,4))

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rbinom(Z[i-1], m, p)))
  }
  
  jitter_offset <- runif(1, -0.05, 0.05)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()






#BINOMIAL2

library(scales)

set.seed(123)
m <- 2
n <- 20
p <- 0.7
num_sim <- 10
base_colors <- rainbow(num_sim)

ruta <- "C:/Users/paula/OneDrive/Documentos/binomial07.png"

png(filename = ruta, width = 1200, height = 900, res = 150)

par(pty="s") #Grafiko karratua

plot(0:n, rep(NA,n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Binomiala"~p==.(p)),
     xlim=c(0,18), ylim=c(0,10))

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rbinom(Z[i-1], m, p)))
  }
  
  jitter_offset <- runif(1, -0.05, 0.05)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()


#POISSON1

library(scales)

set.seed(123)
n <- 20
lambda <- 0.7
num_sim <- 10
base_colors <- rainbow(num_sim)

ruta <- "C:/Users/paula/OneDrive/Documentos/poisson07.png"

png(filename = ruta, width = 1200, height = 900, res = 150)
 
par(pty="s") #Grafiko karratua

plot(0:n, rep(NA,n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Poisson"~lambda==.(lambda)),
     xlim=c(0,10), ylim=c(0,5))

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rpois(Z[i-1], lambda)))
  }
  
  jitter_offset <- runif(1, -0.05, 0.05)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()


#POISSON2

library(scales)

set.seed(123)
n <- 20
lambda <- 1.5
num_sim <- 10
base_colors <- rainbow(num_sim)

ruta <- "C:/Users/paula/OneDrive/Documentos/poisson15.png"

png(filename = ruta, width = 1200, height = 900, res = 150)

par(pty="s") #Grafiko karratua

plot(0:n, rep(NA,n+1), type="n",
     xlab="n", ylab=expression(Z[n]),
     main=bquote("Poisson"~lambda==.(lambda)),
     xlim=c(0,15), ylim=c(0,8))

for (j in 1:num_sim) {
  
  Z <- numeric(n+1)
  Z[1] <- 1
  
  for (i in 2:(n+1)) {
    Z[i] <- ifelse(Z[i-1]==0, 0, sum(rpois(Z[i-1], lambda)))
  }
  
  jitter_offset <- runif(1, -0.05, 0.05)
  col_trans <- alpha(base_colors[j], 0.55)
  
  lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
}

dev.off()
