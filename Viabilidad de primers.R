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
#Se crea función que idnique que hacer dependiendo de los datos que ingrese el usuario#

## Análisis de los primers ##

# reverso complementario
revcomp <- function(x) {
  reverseComplement(DNAString(x))
}

#GC individual#
gc_dimero <- function(primer){
  bases <- strsplit(primer, "")[[1]]
  count <- 0
  for (i in 1:(length(bases)-1)){
    if(bases[i]=="G"&& bases[i+1]=="C"){
      count <- count +1
    }
  }
  
  return(count)
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
  end <- paste0(bases[(length(bases)-2):(length(bases)-1)],
                bases[(length(bases)-1):length(bases)])
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
  for (i in 1: (nucl-min_match)) {
    fragmento <- bases[i:(i+min_match-1)]
    frag_rvcm <- as.character(revcomp(fragmento))
    for (i in 1: (nucl-min_match+1)) {
      horq <- paste(bases[j:(j+min_match-1)], collapse = "")
      if(all(horq==frag_rvcm)) return(T)
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
 ##Mismo primer##
mism_prim <- function(primer){
  dim_prim(primer, primer)
}

 ##foward y reverse##
cross_dim <- function(fw, rv){
  dim_prim(fw, rv)
}

#Evaluación del primer# 
eval_primer <- function(primer, tipo) {
  long <- nchar(primer)
  gc   <- gc_dimero(primer)
  tmel <- tm(primer)
  clamp <- gc_final(primer)
  polim <- hopol(primer)
  hqlla<- horquillas(primer)
  self_prim <- mism_prim(primer)
  
  print(paste(tipo, "-", primer))
  print(paste("Longitud:", long, "pb. ", ifelse(long %in% 18:24, "OK", "No óptimo")))
  print(paste("Dímeros GC:", gc, "%. ", ifelse(gc >=40 & gc<=60, "OK", "No óptimo")))
  print(paste("Tm:", tmel, "°C. ", ifelse(tmel>=55 & tmel<=65, "OK", "No óptimo")))
}

# Viabilidad de ambos
eval_pair <- function(fw, rv) {
  dif <- abs(tm(fw) - tm(rv))
  cross <- cross_dim(fw, rv)
  print(paste("Compatibilidad del par:"))
  print(paste("ΔTm: ", dif, " °C", ifelse(dif<=2, "OK (ideal <=2°C",
                                          ifelse(dif<=5, "OK (<=5°C)", "No óptimo"))))
  print(paste("Cross dímero >=4pb: ", ifelse(cross, "Sí", "No")))
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
    
   print("Ambos primers se unen.")
   print(paste("Forward en:", ini))
   print(paste("Reverse en:", fin - nchar(rc) + 1))
   print(paste("Amplicón:", size, "pb"))
    
    return(list(ok=TRUE, size=size))
  }
  
  print("Problemas de unión:")
  if (!length(uni_fw)) print("Forward no se une.")
  if (!length(uni_rv)) print("Reverse no se une.")
  
  return(list(ok=FALSE, size=0))
}

# VIABILIDAD COMPLETA (RESUMEN)
evaluar_primers <- function(fw, rv, sec) {
  sec <- as.character(normalizar(sec))
  
  print("Evaluación de primers")
  print("--------------------------")
  print(" 1) Individual")
  eval_primer(fw, "Forward")
  eval_primer(rv, "Reverse")
  print("--------------------------")
  print(" 2) Compatibilidad de ambos")
  eval_pair(fw, rv)
  
  print(" 3) Unión a la secuencia")
  res <- union_sec(fw, rv, sec)
  
  print(" 4) Resultados")
  if (res$ok) {
    if (res$size %in% 100:500) print("VIABLE")
    else if (res$size < 100)  print("Amplicón muy pequeño")
    else                      print("Amplicón muy grande")
    
    print(paste("Tamaño:", res$size, "pb"))
  } else {
    print("NO VIABLE")
  }
  
  return(res)
}

blaOxy <-readDNAStringSet("C:/Users/HP/Documents/Fer/Bioinformática/Proyecto_Final/blaOXY.fna")
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

#Simulación PCR 
pcr <- function(ciclos, eficiencia){
  copias <- numeric(ciclos)
  copias [1] <- 1
  for(i in 2:ciclos){
    copias [i] <- copias[i-1] * (1 + eficiencia)
  }
  return(copias)
}

#Sugerencia de mejores primers
sug_primer <- function(sec){
  sec <- as.character(sec)
  largo <- nchar(sec)
  tamaño <- 18:24
  print("Buscando mejores primers dentro de tu secuencia...")
  print("Esto puede tardar unos segundos...")
  for (len_)
 
} 
