#!/usr/bin/perl

use strict;

my $PROGNAME = "Trifonov_Complex.pl";
my $VERSION = "1.01";
my $USAGE=<<"USAGE" ;

#############################################################

 $PROGNAME $VERSION

 Este script evalua la complejidad de las secuencias a través del modelo de complejidad lingüistica de Trifonov
 
 Usage

 1- Modificar la ruta del archivo (path) y el nombre de archivo de entra y salida (file_contigs y file_Salida respectivamente)
 
 2- Indicar por argumentos el tamaño de ventana de nucleótidos a analizar (primer parámetro) y el tamaño del paso a dar (segundo parámetro).
 
	Ejemplo: "perl $PROGNAME 50 20"

 Autor:    

 Marcos Parras

#############################################################

USAGE

#############################################################

unless ($ARGV[0] =~ /[0-9]/ && $ARGV[1] =~ /[0-9]/) { print $USAGE;  exit; }

my $cont_Trif=$ARGV[0]-1;
my $cont_found_Trif;
my $contig;
my $file_Salida="Resultados_complex_Trifonov2.txt";
my $file_contigs="contigs.fasta_OUT.txt_12_OUT.txt";
my $path="/ngs/alopezbueno/Marcos/Shared/Procesado_Metagenomas/MapVir1_ALL/";
my $pos_Trif=0;
my $pos_final=0;
my $Result_Trif;
my $step_size=$ARGV[1]; # Posicion sumada a partir de la cual analizar una ventana de nucléotidos
my $sub_Trif;
my $sub_cad;
my $Tot_sub;
my $win_size=$ARGV[0]; # Tamaño de la ventana a evalua

open(contigs,$path.$file_contigs) || die; # Ruta del archivo de entrada
open(res,">".$path.$file_Salida) || die; # Ruta del archivo de salida

print res "Contig\tPosicion\tComplejidad medida en ventana de $win_size pb\n";

while(<contigs>)
{
	chomp $_;
	
	if($_ =~ /^>/)
	{
		$contig=(split(/ /,(split /^>/,$_)[1]))[0]; # Divide el nombre de la secuencia para obtener un nombre único. Puede ser el nombre o número del contig. Modificar la función de split para coger uno u otro campo
	}
	else 
	{		
		while ($pos_final<=(length($_)) && (length($_)-$pos_final)>($win_size)) # La secuencia se evalua a partir de los parámetros de entrada 
		{
			$sub_cad=substr $_,$pos_final,$win_size; # Se extrae una subcadena desde la secuencia original dependiendo de los parámetros de entrada
						
			### Trifonov. Este modelo evalua la complejidad lingüistica de una cadena dada
			
			$pos_Trif=0;
			
			while ($pos_Trif<=$win_size) # Recorremos la secuencia completa hasta que el paso supere la longitud total de la secuencia
			{
				
				$sub_Trif=substr $sub_cad, $pos_Trif, 2; # Generamos una subcadena de la longitud del tamaño de ventana
				
				if($Tot_sub !~ $sub_Trif)
				{
					$cont_found_Trif++;
					$Tot_sub.=$sub_Trif." ";
				}
				
				$pos_Trif++;
			}
			
			$Result_Trif=($cont_found_Trif/$cont_Trif);
			
			print res $contig."\t".$pos_final."\t".$Result_Trif."\n"; # Los resultados se imprimer a un archivo delimitado por tabuladores
			
			$pos_final=$pos_final+$step_size; # La nueva posición de inicio se determina según el tamaño del paso indicado en los parámetros de entrada
			
			$Tot_sub=undef;
			$cont_found_Trif=0;
			
		}
		
		$pos_final=0;
		$pos_Trif=0;
	}
	
}

close contigs;
close res;
