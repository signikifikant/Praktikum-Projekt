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
#Erstellen einer Liste ("Listenvektor") der Länge n_sim, um später dort die 
#Ergebnisse der einzelnen Durchläufe zu speichern 
sim_results <- vector("list", n_sim)


Faktorenanalyse <- function(n, items, betas, sigma_e, nfactors = NULL) {
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
  if (!(is.matrix(betas) || is.numeric(betas)) || length(betas) != items) {
    stop("betas muss ein numerischer Vektor oder eine Matrix mit Länge = items sein")
  }
  
 
  #Prädiktor simulieren
  X <- matrix(rnorm(n), ncol=1)
  #Fehlermatrix
  E <- matrix(rnorm(n*items, 0, sigma_e), nrow = n, ncol = items)
  #Beta Formatierung 
  beta <- as.matrix(betas)         
  if (ncol(beta) == 1) {
    Y <- X %*% t(beta) + E         
  } else if (nrow(beta) == 1) {
    Y <- X %*% beta + E             
  } else {
    stop("betas muss entweder ein Zeilen- oder Spaltenvektor sein")
  }
  #In Data Frame für Output
  Y_dataf <- as.data.frame(Y)
  #Spalten umbenennen 
  colnames(Y_dataf) <- paste0("Item", 1:items)
  #Korrelationen berechnen
  CorY <- cor(Y)
  #Eigenwert
  eigenvalues <- eigen(CorY)$values
  #Anzahl der Faktoren
  if (is.null(nfactors)) {
    #Parallel-Analysis
    pa <- psych::fa.parallel(Y_dataf, fa="fa", n.iter=100, plot=FALSE)
    nfactors <- pa$nfact
  }
  #Faktoranzahl Kontrolle
  if (nfactors > 0) {
    #EFA rechnen
    efa_result <- fa(Y_dataf, nfactors = nfactors, rotate = "varimax")
    #Liste für Ergebnisse
    return(list(
      Data = Y_dataf,
      Correlation = CorY,
      Eigenvalues = eigenvalues,
      nFactors = nfactors,
      Loadings = efa_result$loadings,
      VarExplained = efa_result$Vaccounted,
      Communalities = efa_result$communality
    ))
  } else {
    warning("Keine Faktoren extrahiert")
    return(list(
      Data = Y_dataf,
      Correlation = CorY,
      Eigenvalues = eigenvalues,
      nFactors = 0,
      Loadings = NULL,
      VarExplained = NULL,
      Communalities = NULL
    ))
  }
}

#n = Anzahl der Personen/Beobachtungen
#itmes = Anzahl der Items
#sigma_e = Standardabweichung des Fehlers 
#lambas = Faktorladungen 
#fac = Anzahl der Faktoren


#Lambda simulieren

Daten_simulieren <- function(n, items, sigma_e, lambdas, fac) {
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


###Korrelationsmatrix (Faktor x Item)
#Maximale Anzahl Faktoren = mF
max_factors <- 10
#Anzahl der Faktoren ziehen
n_factors <- sample(1:max_factors, size = 1, replace = TRUE)
#Maximale Anzahl Items = mI
max_items <- 10

##Alle Fakotren unterschiedliche Item-Anzahl
#Anzahl der Items pro Faktor (Anzahl nach Faktor variabel)
items_per_factor <- sample(1:max_items, size = n_factors, replace = TRUE)
#Anzahl der Items gesamt = nrow
n_items <- sum(items_per_factor)

##Alle Faktoren gleich viele Item-Anzahl
#Anzahl der Items pro Faktor (Alle Faktoren gleich viele Items)
items_per_factor_con <- sample(1:max_items, size = 1, replace = TRUE )
#Anzahl der Items gesamt = nrow
n_items_con <- items_per_factor_con * n_factors

#leere Korrleationsmatrix
lambda <- matrix(0, nrow = n_items, ncol = n_factors)

#Hauptladungen Wertebereich
main_loading_range <- c(0.5, 0.9)
#Nebenladungen Wertebereich
cross_loading_range <- c(0, 0.3)

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

##Erklärung Liste
#Anzahl_Rep => Enthält n, die Anzahl der simulierten Beobachtungen (VPn-Anzahl)
#Items => Anzahl der simulierten Items
#SD_Fehler => Standardabweichung der Fehlermatrix 
#Faktorladungen => Enthält die Faktorladungen (betas)
#Predictors = latente Variablen (Faktoren)
X <- matrix(rnorm(20), ncol=2)
E <- matrix(rnorm(20*5, 0, 0.3, nrow = 20, ncol = 5)
Daten_simulieren(50, 6, 0.3, c(0.3, 0.4, 0.2, 0.5, 0.7, 0.1))
