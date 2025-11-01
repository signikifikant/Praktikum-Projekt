#Packages
install.packages("psych")
library(psych)


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


Faktorenanalyse(500, 5, c(0.3, 0.4, 0.5, 0.1, 0.4), 0.8)


###Funktion für verschiedene Potenzen 
Faktorenanalyse_Pot <- function(n, items, betas, sigma_e, nfactors = NULL, nPot, Pot, betaPot) {
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
  if (!is.numeric(nPot) || length(nPot) != 1 || nPot < 0 || nPot > items) {
    stop("nPot muss zwischen 0 und items liegen")
  }
  if (nPot > 0 && length(betaPot) != nPot) {
    stop("Länge von betaPot muss nPot entsprechen")
  }
  
  #Prädiktor simulieren
  X <- matrix(rnorm(n), ncol=1)
  #Fehlermatrix
  E <- matrix(rnorm(n*items, 0, sigma_e), nrow = n, ncol = items)
  #Beta Formatierung 
  beta <- as.matrix(betas)         
  if (ncol(beta) == 1) {
    Y <- X %*% t(beta) + E         
  } else {
    Y <- X %*% beta + E             
  } 
  #Potenz für i-Items
  if (nPot > 0) {
    for(i in 1:nPot) {
    Y[,i] <- Y[,i] + betaPot[i] * (X^Pot)
    }
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
    
    #Liste erstellen
    res <- list(
      Stichprobengröße = n,
      Itemanzahl = items,
      Regressionskoeffizienten = betas,
      SD_Fehler = sigma_e,
      Data = Y_dataf,
      Correlation = CorY,
      Eigenvalues = eigenvalues,
      Anzahl_an_Faktoren = nfactors,
      Loadings = efa_result$loadings,
      VarExplained = efa_result$Vaccounted,
      Communalities = efa_result$communality,
      Anzahl_der_Potenzen = nPot,
      Potenz = Pot,
      Koeffizient_Potenz = betaPot
    )
  
    class(res) <- "FaktorSim"
    return(res)
    
  } else {
    warning("Keine Faktoren extrahiert")
    
    res <- list(
      Data = Y_dataf,
      Correlation = CorY,
      Eigenvalues = eigenvalues,
      nFactors = 0,
      Loadings = NULL,
      VarExplained = NULL,
      Communalities = NULL
    )
    
    class(res) <- "FaktorSim"
    return(res)
  }
  print.FaktorSim <- function(x, ...) {
    cat("FaktorSim Objekt\n")
    cat("----------------\n")
    cat("Stichprobengröße:", x$Stichprobengröße, "\n")
    cat("Anzahl Items:", x$Itemanzahl, "\n")
    cat("Anzahl extrahierter Faktoren:", x$Anzahl_an_Faktoren, "\n")
    cat("Anzahl der Potenzen:", x$Anzahl_der_Potenzen, "\n")
    if (!is.null(x$Data)) {
      cat("Data: (", nrow(x$Data), "Zeilen ×", ncol(x$Data), "Spalten, zum Anzeigen x$Data aufrufen)\n")
    }
    if (!is.null(x$Eigenvalues)) {
      cat("Eigenwerte:", round(x$Eigenvalues, 3), "\n")
    }
    if (!is.null(x$Loadings)) {
      cat("Loadings: (zum Anzeigen x$Loadings aufrufen)\n")
    }
  }
  }
R.version

Faktorenanalyse_Pot(200, 8, c(0.4, 0.3, 0.5, 0.4, 0.6, 0.7, 0.1, 0.6), 1, nfactors = NULL, 3, 3, c(0.2, 0.5, 0.3) )


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

