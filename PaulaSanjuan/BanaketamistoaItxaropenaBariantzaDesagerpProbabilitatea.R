library(scales)

set.seed(123)

ruta_carpeta <- "C:/Users/paula/OneDrive/Documentos/"

simulatu_prozesua <- function(Z0, n_belaunaldi, banaketa) {
  Z <- numeric(n_belaunaldi + 1)
  Z[1] <- Z0
  for (n in 1:n_belaunaldi) {
    if (Z[n] == 0) {
      Z[n+1] <- 0
    } else {
      if (banaketa == "mistoa") {
        if (n %% 2 == 0) {
          ondorengoak <- rgeom(Z[n], prob = 0.4)
        } else {
          ondorengoak <- rbinom(Z[n], size = 2, prob = 0.7)
        }
      } else if (banaketa == "binomiala") {
        ondorengoak <- rbinom(Z[n], size = 2, prob = 0.7)
      } else if (banaketa == "geometrikoa") {
        ondorengoak <- rgeom(Z[n], prob = 0.4)
      }
      Z[n+1] <- sum(ondorengoak)
    }
  }
  return(Z)
}

# Parametroak
N <- 1000
num_sim_vis <- 10
n_belaunaldi <- 20
Z0 <- 1
banaketak <- c("binomiala","geometrikoa","mistoa")
nombres <- c("Binomial","Geometrikoa","Mistoa")
colors <- c("blue","red","darkgreen")
ltys <- c(1,2,3)
x_vals <- 0:n_belaunaldi

estimaciones_lista <- list()
for (i in 1:3) {
  sim <- replicate(N, simulatu_prozesua(Z0, n_belaunaldi, banaketak[i]))
  estimaciones_lista[[i]] <- list(
    media = apply(sim, 1, mean),
    varianza = apply(sim, 1, var),
    extincion = apply(sim == 0, 1, mean)
  )
}

# Banaketa mistoaren simulazio independente batzuk
sim_mixta_vis <- replicate(num_sim_vis, simulatu_prozesua(Z0, n_belaunaldi, "mistoa"))

escala_var <- 1000

# --- 1. Itxaropena ---
archivo1 <- paste0(ruta_carpeta, "Comparacion_Esperanza.png")
png(filename = archivo1, width = 1500, height = 1500, res = 150)
par(pty="s", mar=c(4,4,4,1))  #Grafiko karratua
plot(x_vals, rep(NA, n_belaunaldi+1), type="n",
     xlab="n", ylab=expression(E[Z[n]]),
     main="E[Z_n]", axes=FALSE,
     cex.main=2, font.main=2,
     xlim=c(0,n_belaunaldi),
     ylim=c(0,max(sapply(estimaciones_lista,function(e) max(e$media)))))
axis(1); axis(2, las=1); box()
for(i in 1:3) lines(x_vals, estimaciones_lista[[i]]$media, col=colors[i], lty=ltys[i], lwd=2)
for(j in 1:num_sim_vis) lines(x_vals, sim_mixta_vis[,j], col=rgb(0,0.5,0,0.3), lwd=1)
dev.off()

# --- 2. Bariantza ---
archivo2 <- paste0(ruta_carpeta, "Comparacion_Varianza.png")
png(filename = archivo2, width = 1500, height = 1500, res = 150)
par(pty="s", mar=c(4,4,4,1))
plot(x_vals, rep(NA, n_belaunaldi+1), type="n",
     xlab="n", ylab=bquote(Var[Z[n]]~"(x10^3)"),
     main="Var[Z_n]", axes=FALSE,
     cex.main=2, font.main=2,
     xlim=c(0,n_belaunaldi),
     ylim=c(0, max(sapply(estimaciones_lista,function(e) quantile(e$varianza/escala_var,0.95)))))
axis(1); axis(2, las=1); box()
for(i in 1:3) lines(x_vals, estimaciones_lista[[i]]$varianza/escala_var, col=colors[i], lty=ltys[i], lwd=2)
for(j in 1:num_sim_vis) lines(x_vals, sim_mixta_vis[,j]^2/escala_var, col=rgb(0,0.5,0,0.3), lwd=1)
dev.off()

# --- 3. Desagerpenerako probabilitatea ---
archivo3 <- paste0(ruta_carpeta, "Comparacion_Extincion.png")
png(filename = archivo3, width = 1500, height = 1500, res = 150)
par(pty="s", mar=c(4,4,4,1))
plot(x_vals, rep(NA, n_belaunaldi+1), type="n",
     xlab="n", ylab=expression(e[n]),
     main="e_n", axes=FALSE,
     cex.main=2, font.main=2,
     xlim=c(0,n_belaunaldi), ylim=c(0,1))
axis(1); axis(2, las=1); box()
for(i in 1:3) lines(x_vals, estimaciones_lista[[i]]$extincion, col=colors[i], lty=ltys[i], lwd=2)
for(j in 1:num_sim_vis) lines(x_vals, sim_mixta_vis[,j]==0, col=rgb(0,0.5,0,0.3), lwd=1)
dev.off()
