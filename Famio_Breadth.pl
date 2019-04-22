#!/usr/bin/perl

use strict;

my $PROGNAME = "Famio_Breadth.pl";
my $VERSION = "1.00";
my $USAGE=<<"USAGE" ;

##############################
##############

# Este script calcula el "Dice" entre dos genomas, teniendo en cuenta el breadth en el numerador del más corto, y el breadth contra si mismo del más corto de cada choque. Ident 35 y longitud 45 por defecto.

 Usage

 1- Modificar la ruta del archivo (path) para indicar la carpeta donde se alojan las secuencias fasta que se van a analizar
 
 3- Para ejecutar el script, utilizar el argumento "V"
 
	Ejemplo: "perl $PROGNAME V"
	
##############
##############################

USAGE

unless ($ARGV[0] eq 'V') { print $USAGE;  exit; }

my $path="./";
my $file;
my $file2;
my %self_hit;
my %lista_scores=undef;
my $i=0;
my $item;
my $key;
my $cont_files_tot=0;
my $izq;
my $der;
my $band=0;
my $ident=0;
my $aa_l=0;
my $pos1=0;
my $pos2=0;
my $pos3=0;
my %posOcupadas=undef;
my $contPos=0;
my $conteoTot=0;
my $seqCont=0;
my $l_contig=0;
my %lista_l;
my $Log_num1;
my $Log_num2;
my $Log_denizq;
my $Log_dender;

opendir(DIR, $path) or die $!; #se abre el directorio
my @files = grep(!/^\.fasta/,readdir(DIR));
closedir(DIR);

open(salida,">Resultados_score.txt"); #Nombre del archivo de saliva con los resultados

print salida "\t";

foreach $file (sort @files)
{
	if($file =~ /\.fasta/)
	{
		print salida $file."\t";
	}
}

print salida $file."\n";


foreach $file (@files)
{
	if($file =~ /\.fasta/)
	{
		system("/ngs/software/bin/makeblastdb -dbtype nucl -in $path$file -out $path$file.DB"); #Modificar la ruta absoluta de la aplicación makeblastdb

		open(filecontigs,$path.$file);

		$seqCont=undef;

		while (<filecontigs>)
		{
			chomp $_;
	
			if($_ =~ "^>")
			{
			}else
			{
				$seqCont.=$_;
			}
		}

		$l_contig=length($seqCont);

		$lista_l{$file}=$l_contig;

		foreach $file2 (@files)
		{
			if($file2 =~ /\.fasta/)
			{			
				system("/ngs/software/bin/tblastx -evalue 0.00001 -query $path$file2 -db $path$file.DB -matrix BLOSUM45 -outfmt 6 -num_threads 12 -out $path$file$file2.blast"); #Modificar la ruta absoluta de la aplicación tblastx

				open(fasta_leido, $path.$file.$file2.'.blast') or die;

				while(<fasta_leido>)
				{
					chomp $_;

					$ident=(split(/\t/,$_))[2];
					$aa_l=(split(/\t/,$_))[3];

					$pos1=(split(/\t/,$_))[6];
					$pos2=(split(/\t/,$_))[7];

					if($pos2<$pos1)
					{
						$pos3=$pos1;
						$pos1=$pos2;
						$pos2=$pos3;
					}	


					if($ident>=35 && $aa_l>=45) # Modificar los parámetros $ident para el mínimo de identidas y $aa_l para el mínimo de longitud de aminoácidos
					{
						for($contPos=$pos1;$contPos<=$pos2;$contPos++)
						{
							$posOcupadas{$contPos}=1;
						}

					}

					$band=0;
				}

				foreach $item (sort keys %posOcupadas)
				{
					if($item =~ /[0-9]/)
					{
						$conteoTot++;
					}
				}

				%posOcupadas=undef;
				
				print salida2 "\n";
				print $file." ".$file2."\n";

				close fasta_leido;

				if($file eq $file2)
				{
					$self_hit{$file}=$conteoTot;
				}	

				$lista_scores{$file."?".$file2}=$conteoTot;

				$conteoTot=0;
				system("rm -rf *blast"); #Se eliminan progresivamente los blast que se van generando

			}
		}

		system("rm -rf $path$file.DB*"); #Se eliminan las bases DB creadas mediante makeblastdb una vez que todas las secuencias se han alineado contra ella

		$cont_files_tot++;
		
	}
}
close salida2;
print $cont_files_tot."\n";

foreach $key (sort keys %lista_scores)
{
	if($key =~ /./)#Solo imprime las keys que tengas información
	{

		$izq=(split(/\?/,$key))[0];
		$der=(split(/\?/,$key))[1];		

		if($band==0)
		{
			print salida $izq;
			$band++;
		}


		#####Vamos a usar para cada interacción el menor de las dos secuencias, para intentas evitar así numeros negativos

		if($lista_scores{$key}>0)
		{
			$Log_num1=log10($lista_scores{$key});		
		}else
		{
			$Log_num1=$lista_scores{$key};	
		}
	
		if($lista_scores{$der."?".$izq}>0)
		{
			$Log_num2=log10($lista_scores{$der."?".$izq})
		}else
		{
			$Log_num2=$lista_scores{$der."?".$izq}
		}

		if($self_hit{$izq}>0)
		{
			$Log_denizq=log10($self_hit{$izq})
		}else
		{
			$Log_denizq=$self_hit{$izq}
		}

		if($self_hit{$der}>0)
		{
			$Log_dender=log10($self_hit{$der})
		}else
		{
			$Log_dender=$self_hit{$der}
		}

		if($Log_denizq<=$Log_dender)
		{
			print salida "\t".(1-(($Log_num2)/($Log_denizq)));
		}else
		{
			print salida "\t".(1-(($Log_num1)/($Log_dender)));
		}

		$i++;

		if($cont_files_tot==$i)
		{
			print salida "\n"; #Se imprime progresivamente la información de distancia para cada par de archivos fasta
			$i=0;
			$band=0;

		}
	}	
}

close salida;

sub log10 { #Función logaritmo 10
    my $n = shift;
    return log($n)/log(10);
}
