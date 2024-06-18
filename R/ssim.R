# Funzione per calcolare l'indice di similarità strutturale (SSIM) tra due vettori
ssim <- function(x, y, K1 = 0.01, K2 = 0.03, L = 1) {

  
  are_matrix <- is.matrix(x) & is.matrix(y)
  if(are_matrix)
    # Assicurati che i vettori siano della stessa lunghezza
    if (nrow(x) != nrow(y) | ncol(x)!=ncol(y)) {
      stop("I vettori devono avere la stessa lunghezza.")
    }
  else
    if (length(x) != length(y)) {
      stop("I vettori devono avere la stessa lunghezza.")
    }
  # Calcola le medie
  mu_x <- mean(x)
  mu_y <- mean(y)


  if(are_matrix)
  {
    n_item <- nrow(x)*ncol(x)
    # Calculate the variance using the formula
    sigma_x2 <- sum((x - mu_x)^2) / ( n_item - 1)
    sigma_y2 <- sum((y - mu_y)^2) / (n_item - 1)
    # Calculate the covariance using the formula
    sigma_xy <- sum((x - mu_x) * (y - mu_y)) / (n_item - 1)
  }
  else
  {
    # Calcola le varianze
    sigma_x2 <- var(x)
    sigma_y2 <- var(y)
    # Calcola la covarianza
    sigma_xy <- cov(x, y)
  }



  # Costanti per stabilizzare la divisione con denominatori piccoli
  C1 <- (K1 * L)^2
  C2 <- (K2 * L)^2

  # Calcolo dei tre termini del SSIM
  luminance <- (2 * mu_x * mu_y + C1) / (mu_x^2 + mu_y^2 + C1)
  contrast <- (2 * sqrt(sigma_x2) * sqrt(sigma_y2) + C2) / (sigma_x2 + sigma_y2 + C2)
  structure <- (sigma_xy + C2 / 2) / (sqrt(sigma_x2) * sqrt(sigma_y2) + C2 / 2)

  # Calcola SSIM
  ssim_value <- luminance * contrast * structure
  return(ssim_value)
}
