# Tarea1-Bioinformatica

Asunto: Tarea 1 Bioinformática

Programa de alineamiento de secuencias de pares utilizando BLOSUM62 y PAM250.

- Busca las secuencias en GenBank.
- Realiza alineamiento global y local de cada par.
- Muestra score y porcentaje de identidad. 

## Tabla de contenidos
- [Tarea1-Bioinformatica](#tarea1-bioinformatica)
  - [Tabla de contenidos](#tabla-de-contenidos)
  - [Requerimientos](#requerimientos)
  - [Uso](#uso)
  - [Información API GenBank](#información-api-genbank)
  - [Bases de datos disponibles](#bases-de-datos-disponibles)


## Requerimientos

Para instalar las bibliotecas requeridas

~~~
pip3 install -r requirements.txt
~~~

## Uso
Colocar el nombre de las secuencias en pares que se desean buscar en la base de datos de GenBank en un archivo separados por espacios con formato

~~~
acc acc
acc acc
acc acc
~~~

En donde:
    acc: Corresponde al accession number

Ejecutar programa con 

~~~python
python3 main.py -file nucleotide.txt -database nucleotide -format gb --export


python3 main.py -file protein.txt -database protein -format gb --export
~~~

options:
  -h, --help                show this help message and exit
  -database DATABASE        Database of alignment
  -file FILE                File containing accession numbers
  -format FORMAT            Fileformat fasta, fa, gb, etc
  -e, --export, --no-export Export sequences
  
## Información API GenBank

Link de la API de GenBank https://www.ncbi.nlm.nih.gov/home/develop/api/

## Bases de datos disponibles

Para consultar las bases de datos disponibles https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi
