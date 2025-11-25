  #BERNOULLI
  library(scales)
  
  set.seed(123)
  n <- 20
  p <- 0.4
  num_sim <- 10
  Z0_values <- c(1, 10, 100)
  base_colors <- rainbow(num_sim)
  
  ruta_carpeta <- "C:/Users/paula/OneDrive/Documentos/"
  
  for (Z0 in Z0_values) {
    
    archivo <- paste0(ruta_carpeta, "bernoulli_Z0_", Z0, ".png")
    
    png(filename = archivo, width = 1500, height = 1200, res = 150)
    
    par(pty="s")   #Grafiko karratua
    par(mar=c(4,4,4,1))
    
    # Ardatzen tamainak Z0-ren arabera
    if (Z0 == 1) {
      xlim_vals <- c(0, 10)
      ylim_vals <- c(0, 2)
      y_ticks <- seq(0, 2, by=1)  # Zenbaki osoak soilik
    } else if (Z0 == 10) {
      xlim_vals <- c(0, 12)
      ylim_vals <- c(0, 12)
      y_ticks <- seq(0, 12, by=2)
    } else {
      xlim_vals <- c(0, 12)
      ylim_vals <- c(0, 110)
      y_ticks <- seq(0, 110, by=10)
    }
    
    
    plot(0:n, rep(NA, n+1), type="n",
         xlab="n", ylab=expression(Z[n]),
         main=bquote("Bernoulli"~p==.(p)~","~Z[0]==.(Z0)),
         xlim=xlim_vals, ylim=ylim_vals,
         axes=FALSE,
         cex.main=2, font.main=2)  
    
    axis(1)
    axis(2, at=y_ticks, labels=y_ticks, las=1)  
    box()
    
    
    for (j in 1:num_sim) {
      Z <- numeric(n+1)
      Z[1] <- Z0
      for (i in 2:(n+1)) {
        Z[i] <- ifelse(Z[i-1]==0, 0, sum(rbinom(Z[i-1], 1, p)))
      }
      
      jitter_offset <- runif(1, -0.02, 0.02)
      col_trans <- alpha(base_colors[j], 0.55)
      
      lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
    }
    
    dev.off()
  }



#GEOMETRIKOA1


library(scales)

set.seed(123)
n <- 20
p <- 0.4
num_sim <- 10
Z0_values <- c(1, 10, 100)
base_colors <- rainbow(num_sim)
ruta_carpeta <- "C:/Users/paula/OneDrive/Documentos/"

for (Z0 in Z0_values) {
  archivo <- paste0(ruta_carpeta, "Geometriko_p", p, "_Z0_", Z0, ".png")
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s") #Grafiko karratua
  par(mar=c(4,4,4,1))
  
  # Ardatzen tamainak Z0-ren arabera
  if (Z0 == 1) { ylim_vals <- c(0,10) }
  else if (Z0 == 10) { ylim_vals <- c(0,50) }
  else { ylim_vals <- c(0,200) }
  
  if(Z0 == 1) {
    xlim_vals <- c(0,15)
  } else if(Z0 == 10) {
    xlim_vals <- c(0,20)
  } else {
    xlim_vals <- c(0,10)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Geometrikoa"~p==.(p)~","~Z[0]==.(Z0)),
       xlim=xlim_vals, ylim=ylim_vals,
       axes=FALSE, cex.main=2, font.main=2)
  
  axis(1)
  axis(2, at=seq(0, ylim_vals[2], by=ifelse(Z0==1,1,ifelse(Z0==10,10,20))), las=1)
  box()
  
  for (j in 1:num_sim) {
    Z <- numeric(n+1)
    Z[1] <- Z0
    for (i in 2:(n+1)) { Z[i] <- if(Z[i-1]==0) 0 else sum(rgeom(Z[i-1], p)) }
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}



#GEOMETRIKOA2

library(scales)

set.seed(123)
n <- 20
p <- 0.4
num_sim <- 10
Z0_values <- c(1, 10, 100)
base_colors <- rainbow(num_sim)
p <- 0.7

for (Z0 in Z0_values) {
  archivo <- paste0(ruta_carpeta, "Geometriko_p", p, "_Z0_", Z0, ".png")
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s") #Grafiko karratua
  par(mar=c(4,4,4,1))
  
  if (Z0 == 1) { ylim_vals <- c(0,4) }
  else if (Z0 == 10) { ylim_vals <- c(0,10) }
  else { ylim_vals <- c(0,100) }
  
  if(Z0 == 1) {
    xlim_vals <- c(0,10)
  } else if(Z0 == 10) {
    xlim_vals <- c(0,6)
  } else {
    xlim_vals <- c(0,6)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Geometrikoa"~p==.(p)~","~Z[0]==.(Z0)),
       xlim=xlim_vals, ylim=ylim_vals,
       axes=FALSE, cex.main=2, font.main=2)
  
  axis(1)
  axis(2, at=seq(0, ylim_vals[2], by=ifelse(Z0==1,1,ifelse(Z0==10,2,10))), las=1)
  box()
  
  for (j in 1:num_sim) {
    Z <- numeric(n+1)
    Z[1] <- Z0
    for (i in 2:(n+1)) { Z[i] <- if(Z[i-1]==0) 0 else sum(rgeom(Z[i-1], p)) }
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}




#BINOMIAL1
library(scales)

set.seed(123)
n <- 20
num_sim <- 10
Z0_values <- c(1, 10, 100)
base_colors <- rainbow(num_sim)
p <- 0.4
n_trial <- 2

for (Z0 in Z0_values) {
  archivo <- paste0(ruta_carpeta, "Binomial_p", p, "_n", n_trial, "_Z0_", Z0, ".png")
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s") #Grafiko karratua
  par(mar=c(4,4,4,1))
  
  if (Z0 == 1) { ylim_vals <- c(0,4) }
  else if (Z0 == 10) { ylim_vals <- c(0,12) }
  else { ylim_vals <- c(0,100) }
  
  if(Z0 == 1) {
    xlim_vals <- c(0,12)
  } else if(Z0 == 10) {
    xlim_vals <- c(0,20)
  } else {
    xlim_vals <- c(0,20)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Binomiala"~p==.(p)~","~Z[0]==.(Z0)),
       xlim=xlim_vals, ylim=ylim_vals,
       axes=FALSE, cex.main=2, font.main=2)
  
  axis(1)
  axis(2, at=seq(0, ylim_vals[2], by=ifelse(Z0==1,1,ifelse(Z0==10,2,10))), las=1)
  box()
  
  for (j in 1:num_sim) {
    Z <- numeric(n+1)
    Z[1] <- Z0
    for (i in 2:(n+1)) { 
      Z[i] <- if(Z[i-1]==0) 0 else sum(rbinom(Z[i-1], n_trial, p)) 
    }
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}





#BINOMIAL2

p <- 0.7
n_trial <- 2

for (Z0 in Z0_values) {
  archivo <- paste0(ruta_carpeta, "Binomial_p", p, "_n", n_trial, "_Z0_", Z0, ".png")
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s") #Grafiko karratua
  par(mar=c(4,4,4,1))
  
  if (Z0 == 1) { ylim_vals <- c(0,8) }
  else if (Z0 == 10) { ylim_vals <- c(0,50) }
  else { ylim_vals <- c(0,200) }
  
  if(Z0 == 1) {
    xlim_vals <- c(0,12)
  } else if(Z0 == 10) {
    xlim_vals <- c(0,10)
  } else {
    xlim_vals <- c(0,10)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Binomiala"~p==.(p)~","~Z[0]==.(Z0)),
       xlim=xlim_vals, ylim=ylim_vals,
       axes=FALSE, cex.main=2, font.main=2)
  
  axis(1)
  axis(2, at=seq(0, ylim_vals[2], by=ifelse(Z0==1,1,ifelse(Z0==10,10,20))), las=1)
  box()
  
  for (j in 1:num_sim) {
    Z <- numeric(n+1)
    Z[1] <- Z0
    for (i in 2:(n+1)) { 
      Z[i] <- if(Z[i-1]==0) 0 else sum(rbinom(Z[i-1], n_trial, p)) 
    }
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}







#POISSON1

lambda <- 0.7

for (Z0 in Z0_values) {
  archivo <- paste0(ruta_carpeta, "Poisson_lambda", lambda, "_Z0_", Z0, ".png")
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s") #Grafiko karratua
  par(mar=c(4,4,4,1))
  
  if (Z0 == 1) { ylim_vals <- c(0,6) }
  else if (Z0 == 10) { ylim_vals <- c(0,12) }
  else { ylim_vals <- c(0,100) }
  
  if(Z0 == 1) {
    xlim_vals <- c(0,10)
  } else if(Z0 == 10) {
    xlim_vals <- c(0,12)
  } else {
    xlim_vals <- c(0,18)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Poisson"~lambda==.(lambda)~","~Z[0]==.(Z0)),
       xlim=xlim_vals, ylim=ylim_vals,
       axes=FALSE, cex.main=2, font.main=2)
  
  axis(1)
  axis(2, at=seq(0, ylim_vals[2], by=ifelse(Z0==1,1,ifelse(Z0==10,2,10))), las=1)
  box()
  
  for (j in 1:num_sim) {
    Z <- numeric(n+1)
    Z[1] <- Z0
    for (i in 2:(n+1)) { Z[i] <- if(Z[i-1]==0) 0 else sum(rpois(Z[i-1], lambda)) }
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}






#POISSON2

lambda <- 1.5

for (Z0 in Z0_values) {
  archivo <- paste0(ruta_carpeta, "Poisson_lambda", lambda, "_Z0_", Z0, ".png")
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s") #Grafiko karratua
  par(mar=c(4,4,4,1))
  
  if (Z0 == 1) { ylim_vals <- c(0,10) }
  else if (Z0 == 10) { ylim_vals <- c(0,52) }
  else { ylim_vals <- c(0,200) }
  
  if(Z0 == 1) {
    xlim_vals <- c(0,15)
  } else if(Z0 == 10) {
    xlim_vals <- c(0,10)
  } else {
    xlim_vals <- c(0,10)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Poisson"~lambda==.(lambda)~","~Z[0]==.(Z0)),
       xlim=xlim_vals, ylim=ylim_vals,
       axes=FALSE, cex.main=2, font.main=2)
  
  axis(1)
  axis(2, at=seq(0, ylim_vals[2], by=ifelse(Z0==1,1,ifelse(Z0==10,10,20))), las=1)
  box()
  
  for (j in 1:num_sim) {
    Z <- numeric(n+1)
    Z[1] <- Z0
    for (i in 2:(n+1)) { Z[i] <- if(Z[i-1]==0) 0 else sum(rpois(Z[i-1], lambda)) }
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}


