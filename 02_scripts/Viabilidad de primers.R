#####################################################################
##   Código para A.M.P.L.I.F.Y. Simulador in silico de una PCR     ##
##            Amplification of proteins for MicroBaby´s            ##
##  Por: Pacheco Estrella María Fernanda y Puente Rivera Diego     ##
#####################################################################

#Cargar paqueterías necesarias#
library(Biostrings)

######## VIABILIDAD DE PCR ###################

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
              ifelse(cross >= 3, "Demasiados heterodímeros, no funcional", "Es aceptable"), "\n\n")
  
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


########### PCR
sim_pcr <- function(E, cycles = cycles, N0 = 1) {
  copias <- numeric(cycles + 1)
  copias[1] <- N0
  
  for (i in 2:(cycles + 1)) {
    copias[i] <- copias[i - 1] * (1 + E)
  }
  
  return(copias)
}

# VIABILIDAD COMPLETA (RESUMEN)
evaluar_primers <- function(fw, rv, sec, cycles=30, N0=1) {
  
  sec <- normalizar(sec)
  
  cat("[Evaluación de primers]\n")
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
  
  cat("4) Viabilidad final y simulación PCR")
  
#Si los primers no se unen bien#
  if (res$ok) {
    cat("Los primers no funcionan porque no se pudieron alinear bien con la 
        secuencia.\n")
    cat("Por esta razón no se hace la simulación PCR.\n")
    return(list(union_sec=res, copias = NULL))
  }
  
#Si el tamaño del amplicón no es útil#
  if (res$size < 100 || res$size > 2000) {
    cat("Los primers si se unen, sin embargo el amplicón no queda en un tamaño
        que sea adecuado.\n")
    cat("Tamaño del amplicón: ", res$size, " pb\n\n")
    cat("No se hace simulación PCR .\n")
    return(list(union_sec = res, copias = NULL))
  }
  
#Si llegó hasta acá, todo esta bien#
  cat("Los primers son viables .\n")
  cat("Tamaño del amplicón: ", res$size, " pb\n\n")
  
#Función para estimar eficiencia segpun la diferencia de la Tm#
  eficiencia <- function(fw, rv){
    dif_tm <- abs(tm(fw) - tm(rv))
    if (dif_tm <=2) {
      return(0.95)
    } else if (dif_tm <=5) {
      return(0.85)
    } else {
      return(0.60)
    }
  }

  
#Calcular eficiencia y simular PCR#
efi <- eficiencia(fw, rv)
copias <- sim_pcr(efi, cycles = cycles, N0 = N0)
cat("Eficiencia astimada: ", round(efi, 3), "\n")
cat("Copias después de ", cycles, " ciclos:", format(copias[length(copias)], 
                                                     scientific = T), "\n")
return(list(union_sec = res, eficiencia = efi, copias = copias))
} 
blaOxy <-Biostrings::readDNAStringSet("blaOXY.fna")
#### No funcionan estos
resultado <- evaluar_primers(
  fw = "AATTGATGATGGAATTCCAT",
  rv = "GGTCCGCAGACGGCATGAA",
  sec = blaOxy,
  cycles = 45
)

#### Estos primers sí funcionan
resultado <- evaluar_primers(
  fw = "AGAGCGGAATGACGCTGGCT",
  rv = "CGATCGAGACGAAAGGTGGC",
  sec = blaOxy,
  cycles = 45
)


#Si los primer ingresados no funcionan, se sugieren mejores primers#
sug_primer <- function(sec, fw, rv){
  #Aquí vamos a revisar si los primers que ingresó lx usuarix si sirven#
  rev_seq <- reverseComplement(DNAString(sec))
  pos_fw <- matchPattern(fw, sec)
  pos_rv <- matchPattern(rv, sec)
  if (length(pos_fw) > 0 && length(pos_rv) > 0) {
    #significa que si se pudo alinear no recomendemos nada#
    return(NULL)
  }


# Si llega a esta parte, significa que no sirven sus primers y buscamos nuevos :(#

  cat("\n --- Buscando primers alternativos (los ingresados por lx alumnx no
      funcionaron) ---\n")
  
  sec_len <- nchar(sec)
  sec_vec <- as.character(sec)
  sug_fw <- c()
  sug_rv <- c()
  
  #Buscar primer foward:#
  for (i in 1: (sec_len - 24)) {
    for (l in 18:24) {
      subn <- substr(sec_vec, i, i + 1 - 1)
      gc <- letterFrequency(DNAString(subn), "GC") / nchar(subn)
      tm_v <- tm(subn)
      if (tm_v >=55 && tm_v <=65 && gc >=0.40 && gc <=0.60) {
        sug_fw <- rbind(sug_fw, data.frame(sec = subn, pos = i, tm = tm_v))
      }
    }
  }
#Buscar primers reverse:#
  rev_vec <- as.character(reverseComplement(DNAString(sec)))
  for (i in 1:(sec_len-24)){
  for (l in 18:24){
    subn <- substr (rev_vec, i, i + l - 1)
    gc <- letterFrequency(DNAString(subn), "GC") / nchar(subn)
    tm_v <- tm(subn)
    if (tm_v >=55 && tm_v<=65 && gc >= 0.40 && gc <=0.60){
      sug_rv <- rbind(sug_rv, data.frame(sec = subn, pos = i, tm=tm_v))
    }
  }
  }
  
  if(nrow(sug_fw) == 0 || nrow(sug_rv) == 0){
    cat("No pude encontrar primers adecuados en la secuencia.\n")
    return(NULL)
  }
  
  #Probar parejas de primers hasta encontrar una que genere un amplicón válido#
  for (i in 1:nrow(sug_fw)){
    for (j in 1:nrow(sug_rv)) {
      inicio <- sug_fw$pos[i]
      fin <-seq_len - sug_rv$pos[j]
      tamaño <- fin - inicio + 1
      if (tamaño >=100 && tamaño <= 2000){
        cat("Se encontró un par recomendable:\n")
        cat("Primer foward: ", sug_fw$sec[i], "\n")
        cat("Primer reverse: ", sug_rv$sec[j], "\n")
        cat("Tamaño del amplicón: ", tamaño, " pb\n")
        return(list(
          foward = sug_fw[i,],
          reverse = sug_rv[j,],
          amplicon = tamaño
        ))
      }
    }
  }
  cat ("Se encontraron primers, pero ninguno forma un amplicón que cumpla con 
       el rango establecido.\n")
  return(NULL)

}
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
