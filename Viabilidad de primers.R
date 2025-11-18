#####################################################################
##   Código para A.M.P.L.I.F.Y. Simulador in silico de una PCR     ##
##            Amplification of proteins for MicroBaby´s            ##
##  Por: Pacheco Estrella María Fernanda y Puente Rivera Diego     ##
#####################################################################

#Cargar paqueterías necesarias#
library(Biostrings)


cat("Bienvenidx a A.M.P.L.I.F.Y. - Simulador PCR in silico\n") 
#Mensaje inicial#
fw <- toupper(readline("Ingresa tu primer forward (5' -> 3'):")) 
#primer foward#
rv <- toupper(readline("Ingresa tu primer reverse (3' -> 5'):")) 
#primer reverse#
ciclos <- as.numeric(readline("¿Cuántos ciclos quieres simular?"))
#cuantos ciclos necesita el estiduante simular#
sec <- readDNAStringSet(readline("Coloca tu archivo en formato FNA:"))
#secuencia blanco#


# función que normaliza la secuencia (porque a veces es DNAStringSet)
normalizar <- function(sec){
  if (class(sec) == "DNAStringSet") sec <- sec [[1]]
  if (class(sec) == "character") sec <- DNAString(sec)
  return(sec)
}
#Se crea función que idique que hacer dependiendo de los datos que ingrese el usuario#

## Análisis de los primers ##

# reverso complementario
revcomp <- function(x) {
  as.character(reverseComplement(DNAString(x)))
}

# Porcentaje de GC#
gc_porcent <- function(primer){
  bases <- strsplit(primer, "")[[1]]
  gc <- sum(bases %in% c("G","C"))
  return(round(100 * gc / length(bases), 1))
}

# Homodimero #
homodimero <- function(primer){
  seq1 <- strsplit(primer, "")[[1]]
  seq2 <- strsplit(revcomp(primer), "")[[1]]
  n1 <- length(seq1)
  n2 <- length(seq2)
  
  max_match <- 0
  
  for(i in 1:n1){
    for(j in 1:n2){
      k <- 0
      while(i+k <= n1 && j+k <= n2 && seq1[i+k] == seq2[j+k]){
        k <- k + 1
      }
      if(k > max_match) max_match <- k
    }
  }
  
  return(max_match)
}

# Heterodimero #
heterodimero <- function(fw, rv){
  seq1 <- strsplit(fw, "")[[1]]
  seq2 <- strsplit(revcomp(rv), "")[[1]]
  n1 <- length(seq1)
  n2 <- length(seq2)
  
  max_match <- 0
  
  for(i in 1:n1){
    for(j in 1:n2){
      k <- 0
      while(i+k <= n1 && j+k <= n2 && seq1[i+k] == seq2[j+k]){
        k <- k + 1
      }
      if(k > max_match) max_match <- k
    }
  }
  
  return(max_match)
}

# Temperatura (Tm)#
tm <- function(primer) {
  bases <- strsplit(primer, "")[[1]]
  A <- sum(bases=="A")
  T <- sum(bases=="T")
  G <- sum(bases=="G")
  C <- sum(bases=="C")
  Tm <- 2*(A+T)+ 4*(G+C)
  return(round(Tm,1))
}

#GC al final del primer#
gc_final <- function(primer){
  bases <- strsplit(primer, "")[[1]]
  n <- length(bases)
  end <- paste0(bases[(n-1):n], collapse = "")
  sum(end == "GC")
}

#Homopolímeros (AAA, CCC, GGG, TTT)#
hopol<- function(primer, k=3){
  bases <- strsplit(primer, "")[[1]]
  count <- 1
  for (i in 2:length(bases)) {
    if(bases[i]==bases[i-1]){
      count <- count + 1
      if (count>=k) return(T)
    } else {
      count <- 1
    }
  }
  return(F)
}

# Detectar horquillas #
horquillas <- function(primer, min_match=4){
  bases <- strsplit(primer, "")[[1]]
  nucl <- length(bases)
  if (nucl < min_match) return(F)
  
  for (i in 1: (nucl-min_match+1)) {
    fragmento <- paste0(bases[i:(i+min_match-1)], collapse = "")
    
    frag_rvcm <- as.character(revcomp(fragmento))
    
    for (j in 1: (nucl-min_match+1)) {
      horq <- paste0(bases[j:(j+min_match-1)], collapse = "")
      
      if(horq==frag_rvcm) return(T)
    }
  }
  return(F)
}

# Evitar que el primer se una consigo mismo o con otro #

dim_prim <- function(seq1, seq2, min_match=4){
  bases1 <- strsplit(seq1, "")[[1]]
  bases2 <- strsplit(revcomp(seq2), "")[[1]]
  n1 <- length(bases1)
  n2 <- length(bases2)
  for (i in 1: (n1-min_match+1)) {
    for (j in 1: (n2-min_match+1)) {
      if (all(bases1[i:(i+min_match-1)] == bases2[j:(j+min_match-1)])) return(T)
    }
  }
  return(F)
}

## Mismo primer##
mism_prim <- function(primer, min_match=4){
  dim_prim(primer, primer, min_match=4)
}

## foward y reverse##
cross_dim <- function(fw, rv){
  dim_prim(fw, rv)
}

#Evaluación del primer# 
eval_primer <- function(primer, tipo) {
  long <- nchar(primer)
  gc_percent <- gc_porcent(primer)
  homod <- homodimero(primer)
  tmel <- tm(primer)
  clamp <- gc_final(primer)
  polim <- hopol(primer)
  hqlla<- horquillas(primer)
  self_prim <- mism_prim(primer)
  
  cat(tipo, "-", primer, "\n")
  cat("Longitud:", long, "pb. ", 
              ifelse(long %in% 18:24, "OK", "No óptimo"), "\n")
  cat("%GC:", gc_percent, "%. ", 
              ifelse(gc_percent >= 40 & gc_percent<=60, "OK", "No óptimo"), "\n")
  cat("Tm:", tmel, "°C. ", 
              ifelse(tmel>=55 & tmel<=65, "OK", "No óptimo"), "\n")
  cat("Homodímeros máximos:", homod, "pb. ", 
              ifelse(homod >= 3, "No óptimo", "OK"), "\n")
  cat("Horquillas:", 
              ifelse(hqlla, "Sí, no óptimo", "No"), "\n\n")
}

# Viabilidad de ambos
eval_pair <- function(fw, rv) {
  dif <- abs(tm(fw) - tm(rv))
  cross <- heterodimero(fw, rv)
  
  cat("Compatibilidad del par:\n")
  cat("ΔTm: ", dif, " °C", 
              ifelse(dif<=2, "OK (ideal <=2°C",
                     ifelse(dif<=5, "OK (<=5°C)", "No óptimo")), "\n")
  
  cat("Heterodimero: ", cross, "pb. ", 
              ifelse(cross >= 3, "Es aceptable", "Demasiados heterodímeros, no funcional"), "\n\n")
  
  return(list(difTm=dif, cross=cross))
}

## Búsqueda de primers en la secuencia ##
# Se unen o no a la secuencia
union_sec <- function(fw, rv, sec) {
  rc <- revcomp(rv)
  
  uni_fw <- start(matchPattern(fw, sec))
  uni_rv <- start(matchPattern(rc, sec))
  
  if (length(uni_fw) && length(uni_rv)) {
    ini <- min(uni_fw)
    fin <- max(uni_rv) + nchar(rc) - 1
    size <- fin - ini + 1
    
    cat("Ambos primers se unen. \n")
    cat("Forward en:", ini, "\n")
    cat("Reverse en:", fin - nchar(rc) + 1, "\n")
    cat("Amplicón:", size, "pb \n\n")
    
    return(list(ok=TRUE, size=size))
  }
  
  cat("Problemas de unión: \n")
  if (!length(uni_fw)) cat("Forward no se une.\n")
  if (!length(uni_rv)) cat("Reverse no se une.\n")
  cat("\n")
  return(list(ok=FALSE, size=0))
}

# VIABILIDAD COMPLETA (RESUMEN)
evaluar_primers <- function(fw, rv, sec) {
  sec <- normalizar(sec)
  
  cat("Evaluación de primers\n")
  cat("------------------------------\n")
  cat("1) Individual\n")
  eval_primer(fw, "Forward")
  eval_primer(rv, "Reverse")
  cat("------------------------------\n")
  cat("2) Compatibilidad de ambos\n")
  eval_pair(fw, rv)
  cat("------------------------------\n")
  cat("3) Unión a la secuencia\n")
  res <- union_sec(fw, rv, sec)
  cat("------------------------------\n")
  cat("4) Resultados\n")
  if (res$ok) {
    if (res$size %in% 100:500) cat("VIABLE\n")
    else if (res$size < 100)  cat("Amplicón muy pequeño\n")
    else                      cat("Amplicón muy grande\n")
    
    cat("Tamaño:", res$size, "pb \n")
  } else {
    cat("NO VIABLE\n")
  }
  
  return(res)
}

blaOxy <-readDNAStringSet("C:/Users/diego/OneDrive/Documentos/GitHub/Proyecto-final/blaOXY.fna")
resultado <- evaluar_primers(
  fw = "AATTGATGATGGAATTCCAT",
  rv = "GGTCCGCAGACGGCATGAA",
  sec = blaOxy
)




#Eficiencia
calcular_eficiencia <- function(fw, rv){
  bases_fw <- strsplit(fw, "")[[1]]
  gc_fw <- gc_dimero(fw)
  bases_rv <- strsplit(rv, "")[[1]]
  gc_rv <- gc_dimero(rv)
  tm_fw <- tm(fw)
  tm_rv <- tm(rv)
  dif_tm <- abs(tm_fw - tm_rv)
  pen_gc <- ifelse(gc_fw <0.40 | gc_fw > 0.60, 0.85, 1)
  pen_gc <- ifelse(gc_rv <0.40 | gc_rv > 0.60, 0.85, 1)
  pen_tm <- ifelse(dif_tm > 2, 0.80, 1)
  eficiencia <- 1 * pen_gc * pen_tm
  return(eficiencia)
}

calcular_eficiencia (fw = "AATTGATGATGGAATTCCAT",
rv = "GGTCCGCAGACGGCATGAA")

#############################



#Simulación PCR 
pcr <- function(ciclos, eficiencia){
  copias <- numeric(ciclos)
  copias [1] <- 1
  for(i in 2:ciclos){
    copias [i] <- copias[i-1] * (1 + eficiencia)
  }
  return(copias)
}

pcr(20, )


#Sugerencia de mejores primers
sug_primer <- function(sec){
  sec <- as.character(sec)
  largo <- nchar(sec)
  tamaño <- 18:24
  print("Buscando mejores primers dentro de tu secuencia...")
  print("Esto puede tardar unos segundos...")
  for (len_)
} 
