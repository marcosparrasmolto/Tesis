# Versión 1.01
# Poner todas los archivos multifasta de los cuales se quieran analizar sus frecuencias nucleotídicas en una misma carpeta. Modificar los valores de correlación y frecuencia nucleotídica al gusto. 

source("https://bioconductor.org/biocLite.R") #Instalamos el paquete necesario para leer las secuencias de nuestros contigs
biocLite("Biostrings") 

library(Biostrings) #Cargamos la librería

nucltot <- 4 #Indicamos la frecuencia nucleotídica que queramos porcesar. En neustro caso es 4 porque evaluamos la frecuencia tetranucleotídica
corr <- 0.96 #Aquí podemos modificar la correlación que establecemos como corte. Por defecto 

dir_files<-choose.dir() #Seleccionamos la carpeta que contiene nuestros archivos multifasta de nuestros contigs
files <- list.files(path=dir_files) #Cargamos los nombres de los archivos que están en esta carpeta

tabla_freq_muestras_2<-data.frame(matrix(0, (4^nucltot)+1)) # Tabla para guardar los resultados de nuestra frecuencia

for (arch_tot in 1:length(files))
{
  fastaFile <- readDNAStringSet(paste(dir_files,files[arch_tot],sep="\\")) #Leemos la secuencia de ADN de un contig concreto 
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)

  tabla_freq_muestras<-data.frame(matrix(0, ncol = (length(fastaFile)), (4^nucltot)+1)) # Creamos una tabla de resultados para las frecuencias nucleotídicas dependiendo de la frecuencia nucleotídica en cuestión
  
  for(i in 1:length(fastaFile))
  {
    colnames(tabla_freq_muestras)[i]<-as.character(df[i,1])
    
    contig <- DNAString(df[i,2]) #Convertimos el objeto que contiene la información de la secuencia a objeto tipo cadena de ADN  
    
    contig<-DNAString(paste(contig,reverseComplement(contig),sep="")) #Calculamos la reversa complementaria
    
    obs_df <- (oligonucleotideFrequency(contig,nucltot)) #Calculamos la frecuencia observada de cada nucleotido
    obs_df <- obs_df/(length(contig)-(nucltot-1)) #Calculamos el % de cada nucleotido
    exp_df <- (expectedDinucleotideFrequency(contig,nucltot)) #Calculamos la frencuencia observada de cada combinación según la frecuencia nucleotídica que queramos evaluar en cada caso
    exp_df <- exp_df/(length(contig)) #Calculamos el % de cada nucleotido de la esperada
    
    #####Hacemos la media de las frecuencias. Santos F. 2010. The metavirome of a hypersaline environment
    
    res_freq<-obs_df/exp_df #La frecuencia viene dada como la frecuencia observada entre la esperada
    
    tabla_freq_muestras[1:4^nucltot,i]<-as.vector(as.numeric(as.matrix(res_freq)))
  }
  
  tabla_freq_muestras[(4^nucltot)+1,]<-files[arch_tot]
  
  tabla_freq_muestras_2<-cbind(tabla_freq_muestras_2,tabla_freq_muestras)
  
}  

rownames(tabla_freq_muestras_2)[1:4^nucltot]<-rownames(as.matrix(oligonucleotideFrequency(contig,nucltot)))
rownames(tabla_freq_muestras_2)[(4^nucltot)+1]<-"Archivo de procedencia"
tabla_freq_muestras_2<-tabla_freq_muestras_2[, -1]
tabla_freq_muestras<-tabla_freq_muestras_2
trans_muestras<-t(tabla_freq_muestras)
familias<-as.factor(t(tabla_freq_muestras[(4^nucltot)+1,]))

####Transformación en data frame de numeros

tabla_matrix<-as.data.frame(matrix(as.numeric(trans_muestras[,-((4^nucltot)+1)]),ncol=(4^nucltot)))
rownames(tabla_matrix)<-rownames(trans_muestras)

####Transformación logaritmo10

tabla_dats_log<-log10(tabla_matrix+1) #Transformamos los números resultantes en logaritmo en base 10

############## Sacara los que tengan correlacion por encima de 0.96, para Cytoscape

tabla_cor<-cor(t(tabla_dats_log)) #Calculamos la correlación de Pearson entre las frecuencias nucleotídicas calculadas de todos los contigs

resultados_cor<-data.frame(matrix(0, ncol=3))

for(i in 1:length(tabla_cor[1,]))
{
  for(o in 1:length(tabla_cor[1,]))
  {
    if(tabla_cor[i,o]>corr&tabla_cor[i,o]<1)
    {
      found_cor<-data.frame(rownames(tabla_cor)[i],colnames(tabla_cor)[o],tabla_cor[i,o])
      colnames(found_cor)<-(c("X1","X2","X3"))
      resultados_cor<-rbind(resultados_cor,found_cor)
    }
  }
}

resultados_cor<-resultados_cor[-1,]

expectedDinucleotideFrequency <- function(x,num) #Con esta función calculamos la frecuencia esperada según el modelo matemático Hidden Markov Chain model. Hay que cargarla antes de cargar el resto del script
{
  vec_resultados_mono<-oligonucleotideFrequency(x,1)/length(x) #Frecuencia de mononucleótido
  vec_resultados_anterior<-vec_resultados_mono
  
  for(freq_nuc in 2:num)
  {
    vec_resultados<-rep(0,(4^freq_nuc))
    pos<-1
    
    for(i in 1:length(vec_resultados_mono))
    {
      for(o in 1:length(vec_resultados_anterior))
      {
        vec_resultados[pos]<-vec_resultados_mono[i]*vec_resultados_anterior[o]
        pos<-pos+1
      }
    }
    
    vec_resultados_anterior<-vec_resultados
    
  }
  
  data.frame(matrix(vec_resultados*length(x),ncol=4^(num-1),nrow=4))#Esperada tri
}
