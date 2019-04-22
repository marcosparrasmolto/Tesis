#!/usr/bin/perl

use strict;

my $PROGNAME = "ParseadorBlast.pl";
my $VERSION = "1.21";
my $USAGE=<<"USAGE" ;

#############################################################

 $PROGNAME $VERSION

 Este script analiza los resultados de Blast de los Orfs de todos los contigs lanzados por Blastx contras las bases de datos nr y Phast. Categoriza a los contigs en función de los cinco primero hits de cada ORF y lo asigna a Virus, Bacteria, No hit y De Sastre.
 
 Usage

 * Modificar en el /path/ y poner la ruta de los archivos blast en formato de salida tipo 0.
 * Son necesarios 3 archivos para la taxonomia de virus y su host. Estos archivos deben estar en la misma ruta que el script:
	
	-Resultados_GeneBank20.txt
	-Taxa_ProtVir_2016.txt
	-ICTV-Master-Species-List-2014.csv

 * Los archivos con los ORF sacados por Prodigal se colocarán en una carpeta a un nivel superior de los archivos blast (../), conservando el mismo codigo de nombre para el metagenoma que en el archivo blast.
 * Los archivos con los contigs de CLC se colocarán en una carpeta a dos niveles superiores de los archivos blast (../../), conservando el mismo codigo de nombre para el metagenoma que en el archivo blast.

	Ejemplo:

	Blast: ALH6_CLC_assembly-nr.blast
	FileORF: ../ALH6_CLC_assembly.fa
	FileContigs: ../../ALH6_CLC_assembly.fa 

 1- Si no se quiere realizar una criba por longitud y cobertura de los contigs, indicarlo con el parámetro "N".

 	Ejemplo: "$PROGNAME N"

 2- Si se quiere realizar una criba por longitud y cobertura de los contigs con un doble criterio, indicar los cuatro parámetros del modo: longitud 1, cobertura 1, longitud 2, cobertura 2.

 	Ejemplo: "$PROGNAME 3000 15 10000 4"

 Autor:    

 Marcos Parras

#############################################################

USAGE

#############################################################

unless ($ARGV[0] eq 'N' || ($ARGV[0] =~ /[0-9]/ && $ARGV[1] =~ /[0-9]/ && $ARGV[2] =~ /[0-9]/ && $ARGV[3] =~ /[0-9]/)) { print $USAGE;  exit; }

my $band_info_organ_Primera=0;
my $cinco_primeros_hits=0;
my $concat_prueba_ocho;
my $cont=0;
my $cont2=0;
my $cont3=0;
my $cont4=0;
my $cont5=0;
my $cont_hits_bacterias=0;
my $cont_hits_virus=0;
my $contigAnt;
my $contigNumero1;
my $contigNumero2;
my $contigeste;
my $coverCribaReads;
my $coverageCont;
my $criba; #Variable queguarda la información sobre si queremos o no realizar una criba por Longitud y Cobertura
my $criba_coverage;
my $criba_coverage2;
my $criba_len;
my $criba_len2;
my $der;
my $evalue;
my $evalue_cinco_hits=0;
my $f1=0;
my $f2=0;
my $file;
my $file2;
my $file_reads;
my $hit;
my $hit_best;
my $hit_ident;
my $hits_totales=0;
my $hits_totales_ORF=0;
my $ident;
my $ident_cinco_hits=0;
my $ident_saved;
my $izq;
my $mayor=0;
my $nom_orga1;
my $nom_orga2;
my $nom_orga;
my $nom_orga_best;
my $nombreCont;
my $nombreContigActual;
my $orga;
my $orga_cinco_hits=0;
my $orgas_cinco_best_hit;
my $path="/ngs/alopezbueno/Marcos/Shared/CLC/Contigs/Prot/Blast_Orfs/";
my $ReadsBact=0;
my $ReadsBactCriba=0;
my $ReadsDeSastre=0;
my $ReadsDeSastreCriba=0;
my $ReadsVir=0;
my $ReadsVirCriba=0;
my $relacionTamOrgaCont;
my $seqCont;
my $seqOrf;
my $spe;
my $spe_tun;
my $temp;
my $temp2;
my $tipo;
my $tipo_cinco_hits=0;
my $tipos_cinco_best_hit;
my $todoHits;
my $todoHitsBact;
my $todoHitsDeSastre;
my $todoHitsViruses;

my @evalue=undef;
my @evalue_cinco=undef;
my @ident_cinco=undef;
my @nom_orga=undef;
my @nom_orga1=undef;
my @nom_orga2=undef;
my @nom_orga_best=undef;
my @nombreCont=undef;
my @orgas_cinco=undef;
my @spe=undef;
my @tax=undef;
my @taxa_prot=undef;
my @tipos_cinco=undef;
my @todoHitsDeSastre=undef;
my @todo_el_hit=undef;
my @todosHitsBact=undef;
my @todosHitsViruses=undef;

my %averageCont=undef;
my %cad_organ=undef;
my %cont_orgas_totales=undef;
my %fam_organ=undef;
my %genero_organ=undef;
my %host_organ=undef;
my %ident=undef;
my %ident_virus=undef;
my %ord_organ=undef;
my %orga_hit_ORF=undef;
my %orga_hit_ORF_virus=undef;
my %readsContigs=undef;
my %seqCont=undef;
my %seqOrf=undef;
my %sub_organ=undef;
my %tam_organ=undef;
my %taxa_prot_virales=undef;
my %taxhost_organ=undef;
my %taxonomia_ORF=undef;
my %tipcad_organ=undef;

open (indice, "Resultados_GeneBank20.txt") || die "No se puede abrir el archivo de virus\n"; 
open (taxa, "Taxa_ProtVir_2016.txt") || die "No se puede abrir el archivo de virus\n"; 
open (ictv,"ICTV-Master-Species-List-2014.csv") || die "No se puede abrir el archivo de virus\n";

$criba=$ARGV[0];

if($criba ne 'N') #Si se requiere, se piden dos criterios de criba por longitud y cobertura
{
	$criba_len=$ARGV[0];
	$criba_coverage=$ARGV[1];
	$criba_len2=$ARGV[2];
	$criba_coverage2=$ARGV[3];
}

while(<ictv>)
{
	chomp $_;
	
	if((split(/\;/,$_))[9] =~ /\(/)
	{
		$tipcad_organ{(split(/\(/,(split(/\;/,$_))[9]))[0]}=1;
	}elsif((split(/\;/,$_))[9] =~ /\-/)
	{
		$tipcad_organ{(split(/\-/,(split(/\;/,$_))[9]))[0]}=1;
	}else
	{
		$tipcad_organ{(split(/\;/,$_))[9]}=1;
	}	

	$ord_organ{(split(/\;/,$_))[0]}=1;
	$fam_organ{(split(/\;/,$_))[1]}=1;
	$sub_organ{(split(/\;/,$_))[2]}=1;
	$genero_organ{(split(/\;/,$_))[3]}=1;
}
	
close ictv;
	
while (<indice>)
{
	chomp $_;
		
	$nom_orga1=(split(/\t/,$_))[1];
	$nom_orga2=(split(/\t/,$_))[0];
	
	@nom_orga1=split(/[-'"._ ]/,$nom_orga1);
	@nom_orga2=split(/[-'"._ ]/,$nom_orga2);


	foreach my $item (@nom_orga1)
	{
		$cont++;
	}
	
	$nom_orga1=undef;

	for(my $i=0;$i<$cont;$i++)
	{
		$nom_orga1.=$nom_orga1[$i];
	}

	foreach my $item (@nom_orga2)
	{
		$cont++;
	}
	
	$nom_orga2=undef;

	for(my $i=0;$i<$cont;$i++)
	{
		$nom_orga2.=$nom_orga2[$i];
	}

	$tipcad_organ{$nom_orga1}=(split(/\t/,$_))[2];
	$ord_organ{$nom_orga1}=(split(/\t/,$_))[3];
	$fam_organ{$nom_orga1}=(split(/\t/,$_))[4];
	$genero_organ{$nom_orga1}=(split(/\t/,$_))[6];
	$tam_organ{$nom_orga1}=(split(/\t/,$_))[7];
	$host_organ{$nom_orga1}=(split(/\t/,$_))[8];
	$taxhost_organ{$nom_orga1}=(split(/\t/,$_))[9];
	$taxhost_organ{$nom_orga1}.="\t".(split(/\t/,$_))[10];
	$taxhost_organ{$nom_orga1}.="\t".(split(/\t/,$_))[11];
	$taxhost_organ{$nom_orga1}.="\t".(split(/\t/,$_))[12];
	$taxhost_organ{$nom_orga1}.="\t".(split(/\t/,$_))[13];
	$taxhost_organ{$nom_orga1}.="\t".(split(/\t/,$_))[14];

	$tipcad_organ{$nom_orga2}=(split(/\t/,$_))[2];
	$ord_organ{$nom_orga2}=(split(/\t/,$_))[3];
	$fam_organ{$nom_orga2}=(split(/\t/,$_))[4];
	$genero_organ{$nom_orga2}=(split(/\t/,$_))[6];
	$tam_organ{$nom_orga2}=(split(/\t/,$_))[7];
	$host_organ{$nom_orga2}=(split(/\t/,$_))[8];
	$taxhost_organ{$nom_orga2}=(split(/\t/,$_))[9];
	$taxhost_organ{$nom_orga2}.="\t".(split(/\t/,$_))[10];
	$taxhost_organ{$nom_orga2}.="\t".(split(/\t/,$_))[11];
	$taxhost_organ{$nom_orga2}.="\t".(split(/\t/,$_))[12];
	$taxhost_organ{$nom_orga2}.="\t".(split(/\t/,$_))[13];
	$taxhost_organ{$nom_orga2}.="\t".(split(/\t/,$_))[14];

	$cont=0;
}

while(<taxa>)
{
	chomp $_;

	my $cont_taxa=0;
	my @taxa_prot_new;

	if($_ =~ /\;/)
	{
		@taxa_prot=split(/\;/,$_);
		
		$cont_taxa=scalar(@taxa_prot);

		for(my $i=0;$i<$cont_taxa;$i++)
		{
			if($taxa_prot[($cont_taxa-1)-$i] =~ /ssRNA viruses/ || $taxa_prot[($cont_taxa-1)-$i] =~ /dsRNA viruses/ || $taxa_prot[($cont_taxa-1)-$i] =~ /dsDNA viruses/ || $taxa_prot[($cont_taxa-1)-$i] =~ /ssDNA viruses/)
			{
				$taxa_prot[($cont_taxa-1)-$i]=(split(/ /,$taxa_prot[($cont_taxa-1)-$i]))[0];
			}		
				$taxa_prot_new[$i]=$taxa_prot[($cont_taxa-1)-$i];	
		}
	
		@taxa_prot=@taxa_prot_new;

		$taxa_prot_virales{$taxa_prot[0]}=join("\t", @taxa_prot);
		
		if($taxa_prot[0] =~ /./)
		{
			$cad_organ{$taxa_prot[0]}="Virus";
		}
	}else
	{

		@taxa_prot=split(/\t/,$_);

		$taxa_prot_virales{$taxa_prot[0]}=$_;
		
		if($taxa_prot[0] =~ /./)
		{
			$cad_organ{$taxa_prot[0]}="Virus";
		}
	}

}

close taxa;

opendir(DIR, $path) or die $!; #Se abre el directorio con los ficheros
my @files = grep(!/^\.blast/,readdir(DIR));
closedir(DIR);
mkdir($path.'/'."Clasificacion_contigs");

open (readsALL, ">".$path.'/'."Clasificacion_contigs/ReadsContigs.txt");
print readsALL "Metagenoma" . "\t". "Virus" . "\t" . "Bacterias" . "\t" . "DeSastre" . "\n";
open (especies, ">".$path.'/'."Clasificacion_contigs/PotencialesNuevasEspecies.txt");

if($criba ne 'N') #Si se ha requerido, se crea el fichero correspondiente al número de Reads contenidos en contigs que pasan la criba
{
	open (readsALLCriba, ">".$path.'/'."Clasificacion_contigs/ReadsContigsCriba.txt");
	print readsALLCriba "Metagenoma" . "\t". "Virus" . "\t" . "Bacterias" . "\t" . "DeSastre" . "\n";
}

foreach $file (@files)
{
	if($file =~ /\.blast/)
	{
		open (bacteria, ">".$path.'/'."Clasificacion_contigs".'/'."Bacterias_".(split(/\-nr.blast/,$file))[0].".fasta");
		open (contOrgas, ">".$path.'/'."Clasificacion_contigs".'/'."Conteo_Orga_".(split(/\-nr.blast/,$file))[0].".txt");
		open (contVirus, ">".$path.'/'."Clasificacion_contigs".'/'."Conteo_VIRUSES_Relacion".(split(/\-nr.blast/,$file))[0].".txt");
		open (DeSastre,">".$path.'/'."Clasificacion_contigs".'/'."DeSastre_".(split(/\-nr.blast/,$file))[0].".fasta");
		open (virusReads, ">".$path.'/'."Clasificacion_contigs".'/'."ContigReadsViruses_".(split(/\-nr.blast/,$file))[0].".txt");
		open (virusReadsCriba, ">".$path.'/'."Clasificacion_contigs".'/'."ContigReadsVirusesCriba_".(split(/\-nr.blast/,$file))[0].".txt");
		open (viruses, ">".$path.'/'."Clasificacion_contigs".'/'."Virus_".(split(/\-nr.blast/,$file))[0].".fasta");

		chdir $path;

		open (filecontigs, "../../".(split(/-nr.blast/,$file))[0].".fa") || die "El archivo no se ha podido abrir (contigs), se cerrara la aplicacion"; #El archivo que contiene los contigs debe tener una extensión .fa
		open (fileorfs, "../".(split(/-nr.blast/,$file))[0].".fa") || die "El archivo no se ha podido abrir (ORFs), se cerrara la aplicacion\n";

		$file_reads=(split(/CLC_assembly-nr.blast/,$file))[0];
		
		my $pathreads="/ngs/alopezbueno/Marcos/Shared/CLC/AssemblyData/";

		opendir(DIR, $pathreads) or die $!; #se abre el directorio
		my @files3 = grep(!/^\.blast/,readdir(DIR));
		closedir(DIR);		

		foreach my $item (@files3)
		{
			if($item =~ $file_reads)
			{
				$file_reads=$item;
				
			}
		}

		open (filereads, $pathreads.$file_reads) || die "El archivo no se ha podido abrir (reads), se cerrara la aplicacion\n"; #Con esta función se extrae la información de cuantos reads hay por compartimentos en nuestras muestras

		while (<filereads>)
		{
			chomp $_;
			
			($contigNumero1,$contigNumero2)=(split(/ /,(split(/\;/,$_))[0]))[1,2];
			$contigNumero1=$contigNumero1."_".$contigNumero2;

			$readsContigs{$contigNumero1}="Reads: ".(split(/\"/,(split(/\;/,$_))[2]))[1];
		}

		close filereads;

		while (<filecontigs>)
		{
			chomp $_;
			#chop $_; #Solo activar si el multifasta de salida de contigs se ha generado en CLC en Windows, ya que genera un caracter más al final de cada línea
	
			if($_ =~ "^>")
			{
				if($cont2==1)
				{
					$seqCont{$nombreCont}=$seqCont;
				
					$cont2=0;
					$cont=0;
					$seqCont=undef;
				}
				@nombreCont=split(/_/,$_);
		
				foreach my $keys (@nombreCont)
				{
					$cont++;
				}

				$nombreCont="contig_".(split(/ /,$nombreCont[$cont-1]))[0];
				$coverageCont="Coverage: ".(split(/ /,$nombreCont[$cont-1]))[3];			
				$averageCont{$nombreCont}=$coverageCont;
			
			}else
			{
				$seqCont.=$_;
				$cont2=1;
			}
		}

		$seqCont{$nombreCont}=$seqCont;

		$cont2=0;
		$cont=0;
		$seqCont=undef;

		close filecontigs;

		while (<fileorfs>)
		{
			chomp $_;

			if($_ =~ "^>")
			{

				if($cont2==1)
				{
					$seqOrf{$nombreCont}=$seqOrf;
					$cont2=0;
					$cont=0;
					$seqOrf=undef;
				}
				@nombreCont=split(/_/,(split,(/ #/,$_))[0]);
		
				foreach my $keys (@nombreCont)
				{
					$cont++;
				}

				$nombreCont="contig_".$nombreCont[$cont-2]."_".$nombreCont[$cont-1]; 
			}else
			{
				$seqOrf.=$_;
				$cont2=1;
			}
		}

		$seqOrf{$nombreCont}=$seqOrf;
				
		$cont=0;
		$cont2=0;
		$seqOrf=undef;

		close fileorfs;

		open (entrada,$path.$file) || die "El archivo no se ha podido abrir, se cerrara la aplicacion, 2\n";
		open (salida,">".$path.'/'."Clasificacion_contigs".'/'."OUTPUT_".(split(/\-nr.blast/,$file))[0].".txt"); #Salida para cara metagenoma donde se muestra la info de cada ORF según cada contig
		print salida "Contig" . "\t" . "Proteina" . "\t" . "TamProteina(aa)". "\t" ."Identidad" . "\t" . "Organismo" . "\t" . "TipoCad". "\t" . "\t" .  "Cadena" . "\t" ."Orden" . "\t" ."Familia" . "\t". "Genero"."\t"."Host"."\t"."Taxahost"."\n"; #Se imprimer en el archivo de salida una fila con los campos que se van a evaluar sobre cada ORF y contig

		while(<entrada>) #Aquí se empieza a recorrer cada archivo blast
		{
				chomp $_;

				if($_ =~ /^>/ && $cont==0 && $cinco_primeros_hits<=4) 
				{
					$concat_prueba_ocho=undef;
					$hit=(split(/>/,$_))[1];
					$hit_ident=$hit;
					$cont++;
					$f1=1;
					$cont2++;
					$cinco_primeros_hits++;
					$concat_prueba_ocho=$_;

				}elsif($cont>=1 && $cont<=5) #Desde aquí se empiza a recabar la información sobre los 5 primeros hits de cada ORF
				{	
					if($hit_ident !~ "Identities")#Recorre linea a linea hasta llegar al campo Identidad, para tener así todos los datos necesarios recogidos
					{
						$hit_ident.=$_;
						if($_ =~ /         /)
						{
							$hit.=(split(/         /,$_))[1];
						}else
						{
							$hit.=$_;
						}

					}else
					{					
						$cont++;
					}

					$concat_prueba_ocho.=$_;
				}elsif($cont==6)
				{
					$concat_prueba_ocho.=$_;
					
					if ($hit_ident !~ "Identities")
					{
						print "En el hit_ident: ".$hit_ident."\n"; #Linea control para controlar que obtenemos la identidad en toda la info que hemos cogido
					}

					$ident=(split(/\(/,(split(/\),/,(split(/ /,(split(/Identities = /,$hit_ident))[1]))[1]))[0]))[1];
					$evalue=(split(/\, /,(split(/Expect = /,$hit_ident))[1]))[0];
				
					if($evalue =~ /\./)
					{
						@evalue=split(/\./,$evalue);
						
						if($evalue[0]==0 && $evalue != "0.0") #Modificamos el formato del e-value para que sea compresible en cálculo
						{
							$evalue=substr($evalue[1],-1,1)."e-".(length($evalue[1]));
						}

					}elsif($evalue =~ /^e/)
					{
						$evalue="1".$evalue;
					}

					$cont++;
					
					$ident_cinco_hits.=$ident."?";

					$evalue_cinco_hits.=$evalue."?";

					@todo_el_hit=split(/ /,$hit);

					foreach my $item (@todo_el_hit)
					{
						$cont4++;
					}				
				
					$hit=undef;
				
					for (my $i=0;$i<$cont4;$i++)
					{
						if($todo_el_hit[$i] !~ / / && $todo_el_hit[$i] !~ /Score/ && $todo_el_hit[$i] !~ /Identities/ && $todo_el_hit[$i] !~ /Frame/)
						{
							$hit.=$todo_el_hit[$i]." ";
						}
					}			

					$cont4=0;	

					$orga=(split(/\[/,(split(/\]/,$hit))[0]))[1];

					if($orga !~ "[a-zA-Z]")
					{
						$orga=(split(/\[/,(split(/\]/,$hit))[0]))[2];
						$orga.=(split(/\]/,$hit))[1];
					}
				
					$tipo=(split(/\|/,(split(/\[/,$hit))[0]))[2];
			
					if($tipo =~ /gb/ || $tipo =~ /ref/ || $tipo =~ /emb/ || $tipo =~ /sp/ || $tipo =~ /dbj/)
					{
						$tipo=(split(/\|/,(split(/\[/,$hit))[0]))[4];		
					}					

					$tipo =~ s/^\s+//; #Elimina los espacios al principio
	 				$tipo =~ s/\s+$//; #Elimina los espacios al final
					$tipo_cinco_hits.=$tipo."?"; #Se concatena toda la información referente a los tipos

					$nom_orga=$orga;          

					@nom_orga=split(/[-'"._ ]/,$nom_orga);

					foreach my $item (@nom_orga)           
					{
						$cont3++;
					}
	
					$nom_orga=undef;

					for(my $i=0;$i<$cont3;$i++)
					{
						$nom_orga.=$nom_orga[$i];
					}

					if($cinco_primeros_hits==1)
					{
						$nom_orga_best=$orga;          

						@nom_orga_best=split(/[-'"._ ]/,$nom_orga_best);

						foreach my $item (@nom_orga_best)           
						{
							$cont3++;
						}
	
						$nom_orga_best=undef;

						for(my $i=0;$i<$cont3;$i++)
						{
							$nom_orga_best.=$nom_orga_best[$i]; #Se concatena toda la información referente al organismo
						}
					}
				
					$cont3=0;
					$orga =~ s/   / /;
					$orga =~ s/  / /;
					$orga_cinco_hits.=$orga."?"; #Se concatenan los primeros cinco organismos

				}elsif($cont==7)
				{
					$cont=0;
					$cont2=0;
					$hit=undef; #Se reinician variables cuando ya se ha obtenido toda la información significativa
					$ident=undef;
					$relacionTamOrgaCont=undef;

				}elsif($_ =~ /^Query=/ && $f1==1)
				{
					&ImprimirInfoContig;
							
					$cinco_primeros_hits=0;
					$contigeste=$_;
					$evalue_cinco_hits=undef;
					$f1=0;
					$f2=1;
					$ident_cinco_hits=undef;
					$orga_cinco_hits=undef;
					$tipo_cinco_hits=undef;

				}elsif($_ =~ /^Query=/ && $f1==0) #Aquí se entra en el caso que se analice un ORF pero no haya resultados (No hits)
				{
					$cinco_primeros_hits=0;
					$cont2=0;
					$evalue_cinco_hits=undef;
					$ident_cinco_hits=undef;
					$orga_cinco_hits=undef;
					$tipo_cinco_hits=undef;

					if($contigeste =~ "[a-zA-Z]")
					{
						$nombreContigActual=(split(/_/,$contigeste))[0]."_".(split(/_/,$contigeste))[1]; #Comprobamos si estamos analizando un contig diferente

						if($nombreContigActual ne $contigAnt)
						{
							foreach my $keys (sort keys %orga_hit_ORF)
							{								
								if($keys =~ "[a-zA-Z]")
								{												
									if ($orga_hit_ORF{$keys}>$mayor)
									{
										$mayor=$orga_hit_ORF{$keys};
										$spe=$keys;
										$ident_saved=$ident{$keys};
									}elsif($orga_hit_ORF{$keys}==$mayor && $ident{$keys}<$ident_saved)
									{
										$mayor=$orga_hit_ORF{$keys};
										$spe=$keys;
										$ident_saved=$ident{$keys};
									}
								}
							}

							if((split(/Coverage: /,$averageCont{$contigAnt}))[1] =~ /\./)
							{
								($izq,$der)=split(/\./,(split(/Coverage: /,$averageCont{$contigAnt}))[1]);
								$coverCribaReads=$izq.$der;
							}else
							{
								$coverCribaReads=(split(/Coverage: /,$averageCont{$contigAnt}))[1];
							}

							if($hits_totales_ORF>0 && ($cont_hits_bacterias/$hits_totales_ORF)>=0.7 && $todoHitsBact !~ $contigAnt." " && $coverCribaReads<100) #Se comprueba que un 70% de los ORFs dan contra bacteria
							{	
								$todoHitsBact.=$contigAnt." ";
							}elsif($hits_totales>0 && $hits_totales_ORF>0 && $hits_totales_ORF>0 && ($cont_hits_virus/$hits_totales)>=0.05 && $todoHitsViruses !~ $contigAnt." ") #Si el caso anterior no se cumple se evalua que un 5% de los hits sean virales
							{	
								$todoHitsViruses.=$contigAnt." ";
							}elsif($todoHitsBact !~ $contigAnt." " && $todoHitsViruses !~ $contigAnt." " && $todoHitsDeSastre !~ $contigAnt." ") #Si el contig no se categoriza dentro de ninguno de los anteriores grupos, entra en cajón de sastre
							{
								$todoHitsDeSastre.=$contigAnt." ";
							}
								
							if($band_info_organ_Primera==1)
							{
								@spe=split(/[-'"._ ]/,$spe);

								foreach my $item (@spe)
								{
									$cont++;
								}
	
								$spe_tun=undef;

								for(my $i=0;$i<$cont;$i++)
								{
									$spe_tun.=$spe[$i];
								}
			
								$cont=0;

								if($tam_organ{$spe_tun}>0)
								{
									$relacionTamOrgaCont=length($seqCont{$contigAnt})/$tam_organ{$spe_tun};

									print salida "Organismo:"." ".$spe." Tipo de cadena: ".$tipcad_organ{$spe_tun}." Longitud: ".$tam_organ{$spe_tun}." Relacion longitud: ".$relacionTamOrgaCont." Nº hits: ".$mayor."/".$hits_totales." Evalue best hit: ".$ident_saved."\n";
									
									if(length($seqCont{$contigAnt})/$tam_organ{$spe_tun}>=0.7) #En este punto se comprueba si estamos en el caso de un contig casi completo
									{
										print especies $file . "\t" . $contigAnt . "\t" . $spe. "\t" . length($seqCont{$contigAnt}) . "\t" . $tam_organ{$spe_tun}  . "\t" . $relacionTamOrgaCont . "\n";
									}
								}else
								{			
									print salida "Organismo:"." ".$spe." Nº hits: ".$mayor."/".$hits_totales." Evalue best hit: ".$ident_saved."\n";
								}					
					
								$taxonomia_ORF{$contigAnt}="Organismo:"." ".$spe." Tipo de cadena: ".$tipcad_organ{$spe_tun}." ".$mayor."/".$hits_totales." Evalue: ".$ident_saved; #Guarda la info para imprimila en los archivos de salida multifasta

							}
							
							$band_info_organ_Primera=1;

							print salida "\n".$nombreContigActual." "."Longitud ".length($seqCont{$nombreContigActual})." ".$averageCont{$nombreContigActual}."\n";

							$contigAnt=$nombreContigActual;
							$todoHits.=$nombreContigActual." ";
							
							print contVirus $file."\t".$contigAnt."\t".$hits_totales_ORF."\t".$hits_totales."\t".$cont_hits_virus."\n";

							$cont_hits_bacterias=0;
							$cont_hits_virus=0;
							$hits_totales=0;
							$hits_totales_ORF=0;
							$ident_saved=undef;
							$mayor=0;
							$spe=undef;
							%ident=undef;
							%ident_virus=undef;
							%orga_hit_ORF=undef;
							%orga_hit_ORF_virus=undef;
						}
					
					print salida $contigeste . "\t" . "\t" . "\t" .$ident . "\t" . "No hit" . "\n";
						
					$hits_totales_ORF++; #Contamos TODOS los ORFs incluyendo los NoHits. Si queremos solo contar los que son Hits conocidos debemos silenciar esta linea.

					}

					$contigeste=$_;
					$f2=1;

				}elsif($f2==1)
				{
					$cont5=0;
					$f2=0;
					$contigeste.=$_;
					$contigeste=(split(/ #/,(split(/Query= /, $contigeste))[1]))[0];
					@nombreCont=split(/_/,$contigeste);

					foreach my $keys (@nombreCont)
					{
						$cont5++;
					}
				
					$contigeste="contig_".$nombreCont[$cont5-2]."_".$nombreCont[$cont5-1];
					
					if($contigeste =~ /\#/)
					{
						$contigeste=(split(/\#/,$contigeste))[0];
					}
				}
		}

		&ImprimirInfoContig; #Al finalizar un archivo, el último ORF se analiza de la misma forma que el resto, pero sin requerir entrar a un nuevo Query	

		if((split(/Coverage: /,$averageCont{$contigAnt}))[1] =~ /\./)
		{
			($izq,$der)=split(/\./,(split(/Coverage: /,$averageCont{$contigAnt}))[1]);
			$coverCribaReads=$izq.$der;
		}else
		{
			$coverCribaReads=(split(/Coverage: /,$averageCont{$contigAnt}))[1];
		}

		#Se imprimen las listas de hits de virus, bacterias y no hits

		@todosHitsViruses=split(/ /,$todoHitsViruses); #Se imprime el multifasta de contigs virales
			
		foreach my $keys (@todosHitsViruses)
		{
			$temp=$keys." ";
			if($todoHitsBact !~ $temp)
			{
				print viruses ">".$keys." "."Longitud ".length($seqCont{$keys})." ".$averageCont{$keys}." ".$readsContigs{$keys}." ".$taxonomia_ORF{$keys}."\n".$seqCont{$keys}."\n";
				$ReadsVir+=(split(/Reads: /,$readsContigs{$keys}))[1];
		
				print virusReads ">".$keys." ".$readsContigs{$keys}."\n";
				if((split(/Coverage: /,$averageCont{$keys}))[1] =~ /\./)
				{
					($izq,$der)=split(/\./,(split(/Coverage: /,$averageCont{$keys}))[1]);
					$coverCribaReads=$izq.$der;
				}else
				{
					$coverCribaReads=(split(/Coverage: /,$averageCont{$keys}))[1];
				}

				if($criba ne 'N' && ((length($seqCont{$keys})>=$criba_len && $coverCribaReads>=$criba_coverage) || (length($seqCont{$keys})>=$criba_len2 && $coverCribaReads>=$criba_coverage2)))
				{	
					$ReadsVirCriba+=(split(/Reads: /,$readsContigs{$keys}))[1];
					print virusReadsCriba ">".$keys." ".$readsContigs{$keys}."\n";
				}
			}
		}
		
		@todosHitsBact=split(/ /,$todoHitsBact); #Se imprime el multifasta de contigs bacterianos
	
		foreach my $keys (@todosHitsBact)
		{
			$temp=$keys." ";
		
			print bacteria ">".$keys." "."Longitud ".length($seqCont{$keys})." ".$averageCont{$keys}." ".$readsContigs{$keys}." ".$taxonomia_ORF{$keys}."\n".$seqCont{$keys}."\n";
			$ReadsBact+=(split(/Reads: /,$readsContigs{$keys}))[1];
			$temp2=(split(/Reads: /,$readsContigs{$keys}))[1];
			
			if((split(/Coverage: /,$averageCont{$keys}))[1] =~ /\./)
			{
				($izq,$der)=split(/\./,(split(/Coverage: /,$averageCont{$keys}))[1]);
				$coverCribaReads=$izq.$der;
			}else
			{
				$coverCribaReads=(split(/Coverage: /,$averageCont{$keys}))[1];
			}			
					
			if($criba ne 'N' && ((length($seqCont{$keys})>=$criba_len && $coverCribaReads>=$criba_coverage) || (length($seqCont{$keys})>=$criba_len2 && $coverCribaReads>=$criba_coverage2)))
			{	
				$ReadsBactCriba+=(split(/Reads: /,$readsContigs{$keys}))[1];
			}
		}
		
		@todoHitsDeSastre=split(/ /,$todoHitsDeSastre); #Se imprime el multifasta de contigs sin categoría
			
		foreach my $keys (@todoHitsDeSastre)
		{
			$temp=$keys." ";
			if($todoHitsBact !~ $temp && $todoHitsViruses !~ $temp)
			{
				print DeSastre ">".$keys." "."Longitud ".length($seqCont{$keys})." ".$averageCont{$keys}." ".$readsContigs{$keys}." ".$taxonomia_ORF{$keys}."\n".$seqCont{$keys}."\n";
				$ReadsDeSastre+=(split(/Reads: /,$readsContigs{$keys}))[1];
			
				if((split(/Coverage: /,$averageCont{$keys}))[1] =~ /\./)
				{
					($izq,$der)=split(/\./,(split(/Coverage: /,$averageCont{$keys}))[1]);
					$coverCribaReads=$izq.$der;
				}else
				{
					$coverCribaReads=(split(/Coverage: /,$averageCont{$keys}))[1];
				}

				if($criba ne 'N' && ((length($seqCont{$keys})>=$criba_len && $coverCribaReads>=$criba_coverage) || (length($seqCont{$keys})>=$criba_len2 && $coverCribaReads>=$criba_coverage2)))
				{	
					$ReadsDeSastreCriba+=(split(/Reads: /,$readsContigs{$keys}))[1];
				}
			}
		}
	
		print readsALL $file . "\t". $ReadsVir . "\t" . $ReadsBact . "\t" . $ReadsDeSastre . "\n";

		if($criba ne 'N')
		{	
			print readsALLCriba $file . "\t". $ReadsVirCriba . "\t" . $ReadsBactCriba . "\t" . $ReadsDeSastreCriba . "\n";
		}
	
		foreach my $keys (sort keys %cont_orgas_totales)
		{
			print contOrgas $keys."\t".$cont_orgas_totales{$keys}."\n";
		} 
	
		$mayor=0;

		foreach my $keys (sort keys %orga_hit_ORF)
		{
			if($keys =~ "[a-zA-Z]")
			{										
				if ($orga_hit_ORF{$keys}>$mayor)
				{
					$mayor=$orga_hit_ORF{$keys};
					$spe=$keys;
					$ident_saved=$ident{$keys};
						
				}elsif($orga_hit_ORF{$keys}==$mayor && $ident{$keys}<$ident_saved)
				{
					$mayor=$orga_hit_ORF{$keys};
					$spe=$keys;
					$ident_saved=$ident{$keys};
				}
			}
		}

		@spe=split(/[-'"._ ]/,$spe);

		foreach my $item (@spe)
		{
			$cont++;
		}
		
		$spe_tun=undef;
		
		for(my $i=0;$i<$cont;$i++)
		{
			$spe_tun.=$spe[$i];
		}
	
		$cont=0;

		if($tam_organ{$spe_tun}>0)
		{
			$relacionTamOrgaCont=length($seqCont{$contigAnt})/$tam_organ{$spe_tun};
			print salida "Organismo:"." ".$spe." Tipo de cadena: ".$tipcad_organ{$spe_tun}." Longitud: ".$tam_organ{$spe_tun}." Relacion longitud: ".$relacionTamOrgaCont." Nº hits: ".$mayor."/".$hits_totales." Evalue best hit: ".$ident_saved."\n";
							
			if(length($seqCont{$contigAnt})/$tam_organ{$spe_tun}>=0.7)
			{
				print especies $file . "\t" . $contigAnt . "\t" . $spe. "\t" . length($seqCont{$contigAnt}) . "\t" . $tam_organ{$spe_tun}  . "\t" . $relacionTamOrgaCont . "\n";
			}
		}else
		{
			print salida "Organismo:"." ".$spe." Nº hits: ".$mayor."/".$hits_totales." Evalue best hit: ".$ident_saved."\n";
		}					
						
		$taxonomia_ORF{$contigAnt}="Organismo:"." ".$spe." Tipo de cadena: ".$tipcad_organ{$spe_tun}." ".$mayor."/".$hits_totales." ".$ident_saved;

		$band_info_organ_Primera=0;
		$cont_hits_bacterias=0;
		$cont_hits_virus=0;
		$contigAnt=undef;
		$contigeste=undef;
		$hits_totales=0;
		$hits_totales_ORF=0;
		$ident_saved=undef;
		$mayor=0;
		$nombreContigActual=undef;
		$ReadsBact=0;
		$ReadsBactCriba=0;
		$ReadsDeSastre=0;
		$ReadsDeSastreCriba=0;
		$ReadsVir=0;
		$ReadsVirCriba=0;
		$spe=undef;
		$spe=undef;
		$todoHits=undef;
		$todoHitsBact=undef;
		$todoHitsDeSastre=undef;
		$todoHitsViruses=undef;
		%averageCont=undef;
		%cont_orgas_totales=undef;
		%ident=undef;
		%ident_virus=undef;
		%orga_hit_ORF=undef;
		%orga_hit_ORF_virus=undef;
		%readsContigs=undef;
		%seqCont=undef;
		close bacteria;
		close contOrgas;
		close contVirus;
		close DeSastre;
		close entrada;
		close filecontigs;
		close salida;
		close viruses;

	}			
}

close readsALL;
close readsALLCriba;
close especies;

if($criba ne 'N')
{
	my $cont;
	my $contig;
	my $contigant;
	my $coverage_cal;
	my $file;
	my $path = $path.'/'."Clasificacion_contigs/";
	my $path3;
	my $seq;
	my $total_tam;

	opendir(DIR, $path) or die $!; #se abre el directorio
	my @files = grep(!/^\.fa/,readdir(DIR));
	closedir(DIR);
	mkdir $path3=$path.'/'."L".$criba_len."C".$criba_coverage."L".$criba_len2."C".$criba_coverage2.'/'; #Se crea una carpeta extra y se criban todos los resultados generados por longitud y cobertura

	foreach $file (@files)
	{
		if($file =~ /\.fasta/)
		{
			$file2 = $path.'/'.$file; #path absoluto del fichero o directorio
			open (entrada,$file2) || die "El archivo no se ha podido abrir, se cerrara la aplicacion\n";
			open (salida,">".$path3.$file."_".$criba_len."_".$criba_coverage."_".$criba_len2."_".$criba_coverage2.".fasta");#Silenciar si no es IDBA	

			while(<entrada>)
			{
				chomp $_;
				#chop $_;#Silenciar si no es CLC (Archivos de Windows)
			
					if($_ =~ /^>/ && $cont==0) 
					{
						$contig=$_;
						$cont=1;

					}elsif($_ !~ /^>/ && $cont==1)
					{
						$seq.=$_;
						$cont2++;
					}elsif($_ =~ /^>/ && $cont==1)
					{	
						$contigant=$contig;
						$coverage_cal=(split(/Coverage: /,$contig))[1];

						if($coverage_cal =~ /\./)
						{
							($izq,$der)=split(/\./,$coverage_cal);
							$coverage_cal=$izq.$der;
						}

						$contig=$_;
						$total_tam=length($seq);

						if(($total_tam>=$criba_len && $coverage_cal>=$criba_coverage) || ($total_tam>=$criba_len2 && $coverage_cal>=$criba_coverage2))
						{
							print salida $contigant . "\n" . $seq . "\n";#Silenciar si no es CLC. #Recorta el nombre por guiones bajos, si proviene directamente de una limpieza por prinseq, para poderlo subir a metavir y no se corte el final del nombre.
							$total_tam=undef;
							$seq=undef;
							$cont2=0;

						}else{$total_tam=undef;$seq=undef;$cont2=0;}
					}		
			}

			$contigant=$contig;
			$coverage_cal=(split(/Coverage: /,$contig))[1];

			if($coverage_cal =~ /\./)
			{
				($izq,$der)=split(/\./,$coverage_cal);
				$coverage_cal=$izq.$der;
			}

			$contig=$_;
			$total_tam=length($seq);

			if(($total_tam>=$criba_len && $coverage_cal>=$criba_coverage) || ($total_tam>=$criba_len2 && $coverage_cal>=$criba_coverage2))
			{
				print salida $contigant . "\n" . $seq . "\n";#Silenciar si no es CLC. #Recorta el nombre por guiones bajos, si proviene directamente de una limpieza por prinseq, para poderlo subir a metavir y no se corte el final del nombre.
				$total_tam=undef;
				$seq=undef;
				$cont2=0;
			}else{$total_tam=undef;$seq=undef;$cont2=0;}

			$cont2=0;
			$contigant=undef;
			$coverage_cal=undef;
			$seq=undef;
			$total_tam=undef;
			close entrada;
			close salida;

		}
	}
}

sub ImprimirInfoContig {

my $bandera_nueva_lista_hits=0;
my $eval_mayor_cinco_exp;
my $eval_probando_expo;
my $evalue_cinco_best_hit;
my $evalue_mayor_de_cinco=0;
my $fam_organ=undef;
my $genero_organ=undef;
my $ident_cinco_best_hit;
my $ident_mayor_de_cinco=0;
my $ident_shift=0;
my $ident_split=0;
my $ord_organ=undef;
my $probando_umbral=0;
my $tipcad_organ=undef;
my $u=0;
my $ya_contado_bact=0;
my $ya_contado_vir;

my @evalue_cinco_filt=undef;
my @ident_cinco_filt=undef;
my @orgas_cinco_filt=undef;
my @prop_orga;
my @taxa_este_orga;
my @tipos_cinco_filt=undef;

@evalue_cinco=split(/\?/,$evalue_cinco_hits);
@ident_cinco=split(/\?/,$ident_cinco_hits);
@orgas_cinco=split(/\?/,$orga_cinco_hits);
@tipos_cinco=split(/\?/,$tipo_cinco_hits);

$evalue_cinco_best_hit=$evalue_cinco[0];
$ident_cinco_best_hit=$ident_cinco[0];
$orgas_cinco_best_hit=$orgas_cinco[0];
$tipos_cinco_best_hit=$tipos_cinco[0];

#Evaluacion de cuales superan el 10% del best, dentro de los 5 primeros

$ident_mayor_de_cinco=(split(/%/,$ident_cinco[0]))[0];
						
$evalue_mayor_de_cinco=$evalue_cinco[0];
				
foreach my $evalue_mayor (@evalue_cinco) #Se comprueba que resultados pasan el corte del 90% con respecto al e-value del best hit
{
	if($evalue_mayor<$evalue_mayor_de_cinco)
	{
		$evalue_mayor_de_cinco=$evalue_mayor;	
	}
}

for(my $p=0;$p<=4;$p++)
{
	$probando_umbral=$evalue_cinco[$p];

	$eval_probando_expo=abs((split(/e/,$probando_umbral))[1]);
							
	$eval_mayor_cinco_exp=abs((split(/e/,$evalue_mayor_de_cinco))[1]);

	if($eval_probando_expo>=($eval_mayor_cinco_exp*0.9))
	{
		$evalue_cinco_filt[$bandera_nueva_lista_hits]=$evalue_cinco[$p];
		$ident_cinco_filt[$bandera_nueva_lista_hits]=$ident_cinco[$p];
		$orgas_cinco_filt[$bandera_nueva_lista_hits]=$orgas_cinco[$p];
		$tipos_cinco_filt[$bandera_nueva_lista_hits]=$tipos_cinco[$p];								
		$bandera_nueva_lista_hits++;								
	}
}

@evalue_cinco=@evalue_cinco_filt; #Se renombran las variables temporales
@ident_cinco=@ident_cinco_filt;
@orgas_cinco=@orgas_cinco_filt; 
@tipos_cinco=@tipos_cinco_filt;
		
$evalue_cinco_hits=undef;
$ident_cinco_hits=undef;		
$orga_cinco_hits=undef;
$tipo_cinco_hits=undef;

foreach my $refill (@orgas_cinco) #Guardamos los resultados de solo los organismos que han pasado el filtro
{
	$evalue_cinco_hits.=$evalue_cinco[$u]."?";
	$orga_cinco_hits.=$refill."?";
	$ident_cinco_hits.=$ident_cinco[$u]."?";	
	$tipo_cinco_hits.=$tipos_cinco[$u]."?";
	$u++;
}

if($contigeste =~ "[a-zA-Z]")
{
	$nombreContigActual=(split(/_/,$contigeste))[0]."_".(split(/_/,$contigeste))[1]; #Comprobamos si estamos analizando un contig diferente

	if($nombreContigActual ne $contigAnt)
	{
		foreach my $keys (sort keys %orga_hit_ORF)
		{								
			if($keys =~ "[a-zA-Z]")
			{												
				if ($orga_hit_ORF{$keys}>$mayor)
				{
					$mayor=$orga_hit_ORF{$keys};
					$spe=$keys;
					$ident_saved=$ident{$keys};
				}elsif($orga_hit_ORF{$keys}==$mayor && $ident{$keys}<$ident_saved)
				{
					$mayor=$orga_hit_ORF{$keys};
					$spe=$keys;
					$ident_saved=$ident{$keys};
				}
			}
		}

		if((split(/Coverage: /,$averageCont{$contigAnt}))[1] =~ /\./)
		{
			($izq,$der)=split(/\./,(split(/Coverage: /,$averageCont{$contigAnt}))[1]);
			$coverCribaReads=$izq.$der;
		}else
		{
			$coverCribaReads=(split(/Coverage: /,$averageCont{$contigAnt}))[1];
		}

		if($hits_totales_ORF>0 && ($cont_hits_bacterias/$hits_totales_ORF)>=0.7 && $todoHitsBact !~ $contigAnt." " && $coverCribaReads<100) #Se comprueba que un 70% de los ORFs dan contra bacteria
		{	
			$todoHitsBact.=$contigAnt." ";
		}elsif($hits_totales>0 && $hits_totales_ORF>0 && $hits_totales_ORF>0 && ($cont_hits_virus/$hits_totales)>=0.05 && $todoHitsViruses !~ $contigAnt." ") #Si el caso anterior no se cumple se evalua que un 5% de los hits sean virales
		{	
			$todoHitsViruses.=$contigAnt." ";
		}elsif($todoHitsBact !~ $contigAnt." " && $todoHitsViruses !~ $contigAnt." " && $todoHitsDeSastre !~ $contigAnt." ") #Si el contig no se categoriza dentro de ninguno de los anteriores grupos, entra en cajón de sastre
		{
			$todoHitsDeSastre.=$contigAnt." ";
		}
								
		if($band_info_organ_Primera==1)
		{
			@spe=split(/[-'"._ ]/,$spe);

			foreach my $item (@spe)
			{
				$cont++;
			}
	
			$spe_tun=undef;

			for(my $i=0;$i<$cont;$i++)
			{
				$spe_tun.=$spe[$i];
			}
			
			$cont=0;

			if($tam_organ{$spe_tun}>0)
			{
				$relacionTamOrgaCont=length($seqCont{$contigAnt})/$tam_organ{$spe_tun};

				print salida "Organismo:"." ".$spe." Tipo de cadena: ".$tipcad_organ{$spe_tun}." Longitud: ".$tam_organ{$spe_tun}." Relacion longitud: ".$relacionTamOrgaCont." Nº hits: ".$mayor."/".$hits_totales." Evalue best hit: ".$ident_saved."\n";
									
				if(length($seqCont{$contigAnt})/$tam_organ{$spe_tun}>=0.7) #En este punto se comprueba si estamos en el caso de un contig casi completo
				{
					print especies $file . "\t" . $contigAnt . "\t" . $spe. "\t" . length($seqCont{$contigAnt}) . "\t" . $tam_organ{$spe_tun}  . "\t" . $relacionTamOrgaCont . "\n";
				}
			}else
			{			
				print salida "Organismo:"." ".$spe." Nº hits: ".$mayor."/".$hits_totales." Evalue best hit: ".$ident_saved."\n";
			}					
					
			$taxonomia_ORF{$contigAnt}="Organismo:"." ".$spe." Tipo de cadena: ".$tipcad_organ{$spe_tun}." ".$mayor."/".$hits_totales." Evalue: ".$ident_saved; #Guarda la info para imprimila en los archivos de salida multifasta

		}
							
		$band_info_organ_Primera=1;

		print salida "\n".$nombreContigActual." "."Longitud ".length($seqCont{$nombreContigActual})." ".$averageCont{$nombreContigActual}."\n";

		$contigAnt=$nombreContigActual;
		$todoHits.=$nombreContigActual." ";
							
		print contVirus $file."\t".$contigAnt."\t".$hits_totales_ORF."\t".$hits_totales."\t".$cont_hits_virus."\n";

		$cont_hits_bacterias=0;
		$cont_hits_virus=0;
		$hits_totales=0;
		$hits_totales_ORF=0;
		$ident_saved=undef;
		$mayor=0;
		$spe=undef;
		%ident=undef;
		%ident_virus=undef;
		%orga_hit_ORF=undef;
		%orga_hit_ORF_virus=undef
	}

	#De aquí se saca la taxonomía completa a partir de todas las proteínas virales existentes en la DB

	@taxa_este_orga=split(/\t/,$taxa_prot_virales{$orgas_cinco_best_hit});
						
	foreach my $item (@taxa_este_orga)
	{
		if($tipcad_organ{$item}==1)
		{
			$tipcad_organ=$item;
		}elsif($ord_organ{$item}==1)
		{
			$ord_organ=$item;
		}elsif($fam_organ{$item}==1)
		{
			$fam_organ=$item;
		}elsif($genero_organ{$item}==1)
		{
			$genero_organ=$item;
		}
	}
					
	print salida $contigeste ."\t" . $tipos_cinco_best_hit . "\t" . (split(/ /,(split(/Length = /,$hit_best))[1]))[0] . "\t" .$evalue_cinco_best_hit . "\t" . $orgas_cinco_best_hit . "\t" . $tipcad_organ . "\t" . $cad_organ{$orgas_cinco_best_hit} . "\t" . $ord_organ . "\t" . $fam_organ . "\t" . $genero_organ . "\t" . $host_organ{$nom_orga_best} . "\t" . $taxhost_organ{$nom_orga_best} . "\n"; #Se imprime la información de cada ORF

	if($cont_orgas_totales{$orgas_cinco_best_hit}>=1)
	{
		$cont_orgas_totales{$orgas_cinco_best_hit}++;
	}else
	{
		$cont_orgas_totales{$orgas_cinco_best_hit}=1;
	}
		
	$hits_totales_ORF++;

}
		
	if($tipos_cinco_best_hit =~ "[a-zA-Z]")
	{
		$hits_totales++;

		foreach my $orgas_cinco (@orgas_cinco)
		{
			$ident_split=$evalue_cinco[$ident_shift];									
			$ident_shift++;
			@prop_orga=split(/\:/,$orga);

			if($cad_organ{$orgas_cinco_best_hit} =~ "virus" || $cad_organ{$orgas_cinco_best_hit} =~ "Virus" || $tipos_cinco_best_hit =~ "phage " || $tipos_cinco_best_hit =~ "phage-" || $tipos_cinco_best_hit =~ "phage_" || $tipos_cinco_best_hit =~ "Phage " || $tipos_cinco_best_hit =~ "Phage-" || $tipos_cinco_best_hit =~ "Phage_" || $orgas_cinco_best_hit =~ "phage " || $orgas_cinco_best_hit =~ "phage-" || $orgas_cinco_best_hit =~ "phage_" || $orgas_cinco_best_hit =~ "Phage " || $orgas_cinco_best_hit =~ "Phage-" || $orgas_cinco_best_hit =~ "Phage_" || $orgas_cinco_best_hit =~ "virus" || $tipos_cinco_best_hit =~ "capsid" || $tipos_cinco_best_hit =~ "Capsid") #Aquí se comprueba que el ORF sea o no viral
			{ 
				if($ya_contado_vir==0)
				{
					$cont_hits_virus++;
					$ya_contado_vir++;
				}

				$orga_hit_ORF_virus{$orgas_cinco_best_hit}++; #Suma uno al total de ORFs virales

				if($ident_virus{$orgas_cinco_best_hit}==undef || $ident_split<$ident_virus{$orgas_cinco_best_hit}) ### Cambiar signo si es identidad
				{
					$ident_virus{$orgas_cinco_best_hit}=$ident_split;
				}
			}else
			{
				if($ya_contado_bact==0)
				{
					$cont_hits_bacterias++;
					$ya_contado_bact++;									
				}
			}
			
			$orga_hit_ORF{$orgas_cinco}++;

			if($ident{$orgas_cinco}==undef || $ident_split<$ident{$orgas_cinco}) #Cambiar signo si es identidad
			{
				$ident{$orgas_cinco}=$ident_split;
			}

		}
			$ident_shift=0;
			$ya_contado_bact=0;
			$ya_contado_vir=0;
	}
}
