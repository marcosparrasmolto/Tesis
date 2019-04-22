#!/usr/bin/perl

use strict;

my $PROGNAME = "Busqueda_Primers_Mapeo_SISPA.pl";
my $VERSION = "1.01";
my $USAGE=<<"USAGE" ;

#############################################################

 $PROGNAME $VERSION

 Este script busca regiones en cada contig que contengan una secuencia similar a los primers SISPA utilizados (En nuestro caso FRV20, K y 454) 
 
 Usage

 1- Modificar la ruta del archivo (path) y el nombre de archivo de entrada y salida (file_contigs y file_Salida respectivamente)
 
 2- En caso de querer probar otros primers en lugar de los que vienen por defecto, modificar la secuencia dentro del script

 3- Para ejecutar el script, utilizar el argumento "V"
 
	Ejemplo: "perl $PROGNAME V"
	
 Autor:    

 Marcos Parras

#############################################################

USAGE

#############################################################

unless ($ARGV[0] eq 'V') { print $USAGE;  exit; }

my $A454=('ATCGTCGTCGTAGGCTGCTC','GAGCAGCCTACGACGACGAT');
my $FRV20=('GCCGGAGCTCTGCAGATATC','GATATCTGCAGAGCTCCGGC');
my $K=('GACCATCTAGCGACCTCCAC','GTGGAGGTCGCTAGATGGTC');
my $cont=0;
my $file_Salida="Posibles_localizaciones_primers.txt";
my $file_contigs="contigs.fasta_OUT.txt_12_OUT.txt";
my $item;
my $nombre_contig;
my $path="/ngs/alopezbueno/Marcos/Shared/Procesado_Metagenomas/EstudioContigsDesplazados_Vir1_ALL/MapeoTodosContigs/";
my %Primers_cont=(1 => "FRV20", 2 => "K", 3 => "454", 4 => "FRV20", 5 => "K", 6 => "454");
my @Primers=('AGCTCTGCAGATATC','TCTAGCGACCTCCAC','CGTCGTAGGCTGCTC','CTGCAGAGCTCCGGC','GGTCGCTAGATGGTC','GCCTACGACGACGAT'); # Últimos 15 pares de bases de cada primer, tanto de la secuencia Forward como Reverse
 
open(contigs,$path.$file_contigs) || die; # Ruta del archivo de entrada
open(salida,">".$path.$file_Salida) || die; # Ruta del archivo de salida

print salida "Secuencia encontrada"."\t"."Contig"."\t"."Posición del contig"."\t"."Dirección"."\t"."Longitud de la secuencia"."\t"."Primer"."\n";

while(<contigs>)
{
	if($_ =~ "^>")
	{
		$nombre_contig=(split(/>/,(split(/ /,$_))[0]))[1];
	}else
	{
		foreach $item (@Primers)
		{
			$cont++;

			for(my $i=0;$i<8;$i++)
			{
				for(my $e=8;$e<15;$e++)
				{
			  		my $string = $_;
			  		my $char = substr $item, $i, $e;
			  		my $offset=0;

			  		my $result = index($string, $char, $offset);

					while ($result != -1) 
					{
						if($cont<3)
						{
							print salida $char."\t".$nombre_contig."\t".$result."\t"."F"."\t".length($char)."\t".$item."\t".$Primers_cont{$cont}."\n"; # Forward
						}else
						{
							print salida $char."\t".$nombre_contig."\t".$result."\t"."R"."\t".length($char)."\t".$item."\t".$Primers_cont{$cont}."\n"; # Reverse
						}

						$offset = $result + 1;
						$result = index($string, $char, $offset);
					}
				}
			}
		}

		$cont=0;	
	}
}

close contigs;
close salida;
