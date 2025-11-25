# GEOMETRICO + BINOMIAL TXANDAKATUZ (MISTOA)
library(scales)

set.seed(123)
n <- 20
p_geo <- 0.4
p_bin <- 0.7
num_sim <- 10
Z0_values <- c(1, 10, 100)
base_colors <- rainbow(num_sim)

ruta_carpeta <- "C:/Users/paula/OneDrive/Documentos/"

for (Z0 in Z0_values) {
  
  archivo <- paste0(ruta_carpeta, "mixto_Z0_", Z0, ".png")
  
  png(filename = archivo, width = 1500, height = 1200, res = 150)
  
  par(pty="s")   
  par(mar=c(4,4,4,1))
  
  if (Z0 == 1) {
    xlim_vals <- c(0, 15)
    ylim_vals <- c(0, 15)
    y_ticks <- seq(0, 15, by=3)
  } else if (Z0 == 10) {
    xlim_vals <- c(0, 10)
    ylim_vals <- c(0, 25)
    y_ticks <- seq(0, 25, by=5)
  } else {
    xlim_vals <- c(0, 10)
    ylim_vals <- c(0, 200)
    y_ticks <- seq(0, 200, by=20)
  }
  
  plot(0:n, rep(NA, n+1), type="n",
       xlab="n", ylab=expression(Z[n]),
       main=bquote("Geom(0.4)~Bin(2,0.7),"~Z[0]==.(Z0)),
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
      if (Z[i-1] == 0) {
        Z[i] <- 0
      } else {
        if ((i-1) %% 2 == 0) {
          ondorengoak <- rgeom(Z[i-1], p_geo)
        } else {
          ondorengoak <- rbinom(Z[i-1], 2, p_bin)
        }
        Z[i] <- sum(ondorengoak)
      }
    }
    
    jitter_offset <- runif(1, -0.05, 0.05)
    col_trans <- alpha(base_colors[j], 0.55)
    
    lines(0:n, Z + jitter_offset, col=col_trans, lwd=1.7)
  }
  
  dev.off()
}
