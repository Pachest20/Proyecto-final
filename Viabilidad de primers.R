library(Biostrings)

# reverso complementario
revcomp <- function(x) reverseComplement(DNAString(x))

# temperatura
tm <- function(x) {
  gc <- sum(strsplit(x, "") %in% c("G", "C")) 
  round(64.9 + 41 * (gc - 16.4) / nchar(x), 1)
}

# Viabilidad de cada uno
eval_primer <- function(primer, tipo) {
  long <- nchar(primer)
  gc   <- sum(strsplit(primer, "") %in% c("G", "C")) / long * 100
  tmel <- tm(primer)
  
  print(paste(tipo, "-", primer))
  print(paste("Longitud:", long, "pb. ", ifelse(long %in% 18:24, "OK", "No óptimo")))
  print(paste("GC:", round(gc,1), "%. ", ifelse(gc >=40 & gc<=60, "OK", "No óptimo")))
  print(paste("Tm:", tmel, "°C. ", ifelse(tmel>=55 & tmel<=65, "OK", "No óptimo")))
}

# Viabilidad en conjunto
eval_pair <- function(fw, rv) {
  diff <- abs(tm(fw) - tm(rv))
  print(paste("Diferencia de Tm:", diff,  "°C. ", ifelse(diff<=2, "OK", "Muy alta")))
}

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
  sec <- as.character(sec)
  
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

blaOxy <-readDNAStringSet("C:/Users/diego/OneDrive/Desktop/Bioinformatica 2025/blaOXY.fna")
resultado <- evaluar_primers(
  fw = "AATTGATGATGGAATTCCAT",
  rv = "GGTCCGCAGACGGCATGAA",
  sec = blaOxy
)

