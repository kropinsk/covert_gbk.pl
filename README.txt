Directions for using Perl Scripts on Windows 10 computer equipped with Strawberry Perl.  Please note on my computer the perl scripts are in a directory called c:\Max

Starting at "Type here to search" type:

cmd [return]
cd c:\   [return]
cd c:\Max   [return]
c:\Max   [return]
dir   [return]
perl (perl script name.pl) -i (input file name) -o .\ (output file to same directory)--xyl or ffn or ffa or ptt or -1 (upstream sequence length)

perl convert_gbk.pl -i T7.gbk -o .\--xls

fastq-splitter.pl

perl fastq-splitter.pl -i “Phage XYZ Illumina sequence_R1.fastq” -o .\ --n-parts 7

This would divide the fastq file into 7 parts of ~ equal size in bytes. The “ ” are used because the file name has spaces in it, any file name with spaces in it should be quoted on the windows command line. 
