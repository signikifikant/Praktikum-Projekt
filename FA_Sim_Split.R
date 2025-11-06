install.packages("psych")
library(psych)

#Liste erstellen
Ergebnis <- list(
  Simulierte_Werte = NULL,
  KorrMatrix = NULL,
  Anzahl_Faktoren = NULL,
  Faktorladungen = NULL,
  Kommunalität = NULL,
  Erklärte_Varianz = NULL,
  Anzahl_Items = NULL, 
  Fehlermatrix = NULL,
  Chi_Quadrat_Test = NULL,
  RMSEA = NULL,
  TLI = NULL,
  SRMR = NULL
)


#Anzahl der VPn pro Durchgang
n <- 200
###Korrelationsmatrix (Faktor x Item)
#Maximale Anzahl Faktoren = max_factors
max_factors <- 10
#Anzahl der Faktoren ziehen
n_factors <- sample(1:max_factors, size = 1, replace = TRUE)
#Maximale Anzahl Items = max_Items
max_items <- 10

##Alle Fakotren unterschiedliche Item-Anzahl
#Anzahl der Items pro Faktor (Anzahl nach Faktor variabel)
items_per_factor <- sample(1:max_items, size = n_factors, replace = TRUE)
#Anzahl der Items gesamt = nrow
n_items <- sum(items_per_factor)

#Hauptladungen Wertebereich
main_loading_range <- c(0.5, 0.9)
#Nebenladungen Wertebereich
cross_loading_range <- c(0, 0.3)

Korrelationsmatrix <- function(n_items, n_factors, items_per_factor, main_loading_range,
                               cross_loading_range) {
  #leere Korrleationsmatrix
  lambda <- matrix(0, nrow = n_items, ncol = n_factors)
  start <- 1
  for(f in 1:n_factors){
    num_items <- items_per_factor[f]
      if(num_items > 0){
        high_items <- start:(start + num_items - 1)
        lambda[high_items, f] <- runif(num_items, min=main_loading_range[1], 
                                       max=main_loading_range[2])
    start <- start + num_items
    }
  }
  # Nebenladungen zufällig für alle Items/Faktoren, wo noch 0 ist
  lambda[lambda == 0] <- runif(sum(lambda==0), min=cross_loading_range[1], 
                             max=cross_loading_range[2])

  lambda <- round(lambda, 4)
  #Ergebnisse in Liste speichern
  Ergebnis$Faktorladungen <<- lambda
}

Korrelationsmatrix(n_items, n_factors, items_per_factor, main_loading_range,
                   cross_loading_range)


##Alle Faktoren gleich viele Item-Anzahl
#Anzahl der Items pro Faktor (Alle Faktoren gleich viele Items)
items_per_factor_con <- sample(1:max_items, size = 1, replace = TRUE)
#Anzahl der Items gesamt = nrow
n_items_con <- items_per_factor_con * n_factors


Korrelationsmatrix_con <- function(n_items_con, n_factors, items_per_factor_con, 
                                   main_loading_range, cross_loading_range) {
  #leere Korrleationsmatrix
  lambda <- matrix(0, nrow = n_items_con, ncol = n_factors)
  start <- 1
  for(f in 1:n_factors){
    num_items <- items_per_factor_con
    if(num_items > 0){
      high_items <- start:(start + num_items - 1)
      lambda[high_items, f] <- runif(num_items, min=main_loading_range[1], 
                                     max=main_loading_range[2])
      start <- start + num_items
    }
  }
  # Nebenladungen zufällig für alle Items/Faktoren, wo noch 0 ist
  lambda[lambda == 0] <- runif(sum(lambda==0), min=cross_loading_range[1], 
                               max=cross_loading_range[2])
  
  lambda <- round(lambda, 4)
  #Ergebnisse in Liste speichern
  Ergebnis$Faktorladungen <<- lambda
}

Korrelationsmatrix_con(n_items_con, n_factors, items_per_factor_con, 
                       main_loading_range, cross_loading_range) 


##Faktorausprägungen
X <- matrix(rnorm(n*n_factors), ncol = n_factors)

##Fehlermatrix
E <- matrix(rnorm(n*n_items), ncol = n_items)

Y <- X %*% t(Ergebnis$Faktorladungen) + E


Daten_simulieren <- function(n, items, sigma_e, lambda, fac) {
  #Argumente Check
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("n muss eine positive Zahl sein")
  }
  if (!is.numeric(items) || length(items) != 1 || items <= 0) {
    stop("items muss eine positive Zahl sein")
  }
  if (!is.numeric(sigma_e) || length(sigma_e) != 1 || sigma_e < 0) {
    stop("sigma_e muss eine nicht-negative Zahl sein")
  }
  if (!(is.matrix(lambdas) || is.numeric(lambdas)) || length(lambdas) != items) {
    stop("lambda muss ein numerischer Vektor oder eine Matrix mit Länge = items sein")
  }
  if (!is.numeric(fac) || length(fac) != 1 || fac <= 0) {
    stop("fac muss eine positive Zahl sein")
  }
  
  #Lambda-Matrix prüfen und in Matrix umwandeln
  lambda <- as.matrix(lambdas)
  if (fac == 1 && is.vector(lambda)) {
    lambda <- matrix(lambda, nrow=1)
  }
  if (!all(dim(lambda) == c(fac, items))) {
    stop("lambda muss Dimension fac x items haben")
  }
  
  #Ausprägung des Faktors (1 Faktor)
  X <- matrix(rnorm(n*fac), ncol = fac)
  #Fehlermatrix
  E <- matrix(rnorm(n*items, 0, sigma_e), nrow = n, ncol = items)
  
  # Items simulieren
  Y <- X %*% lambda + E
  #In Data Frame für Output
  Y_dataf <- as.data.frame(Y)
  #Spalten umbenennen 
  colnames(Y_dataf) <- paste0("Item", 1:items)
  
  #Liste erstellen
  List <- list(
    Anzahl_Rep = n,
    Items = items,
    SD_Fehler = sigma_e,
    Faktorladungen = lambdas,
    Faktorwerte = X,
    Fehlermatrix = E,
    Itemmatrix = Y_dataf
  )
  return(List)
}



