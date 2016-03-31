miRNA Discovery Project : simple, easy-to-use, rule-based miRNA candidate discovery with multithread acceleration

Version : beta (core development complete, in optimization stage)
Author : Lee Sang-gil, dept. of Applied Biology & Chemistry, Seoul National University, South Korea
E-mail : tkdrlf9202@snu.ac.kr


[Quick Installation & Usage Guide]
Need Python 2.7.6 or newer, bowtie, RNAfold
Simply unzip and run main.py like the following:
python main.py -r "reference fasta format genome path" -i "input fasta format RNA-seq path" -o "output path"
For help screen:
python main.py -h


[Introduction]
This program discovers miRNA candidates, showing mature miRNA-miRNA* duplex and pre-miRNA sequence
It provides alignment information, simplified tabular information, and length distribution of genome-mapped RNA-seq reads
Some filtering procedure, I/O formats, and parameter schemes are inspired by miREAP (by Li Qibin <liqb@genomics.org.cn> Bioinformatics department, Beijing Genomics Institute)
It aims to minimize input file lists and software dependencies, it is written in python script with no installation hassle, and it doesn't require sudo (administrator) privilege to run, which is sometimes a problem when you deploy it at a server
It provides multithread implementation to maximize performance, which is scaled to the cpu threads you use
Albeit there are a number of similar softwares with machine learning techniques achieving high accuracy, this program can be used with in-house data with no extra data, and with simple and fast use
Scoring scheme is defined as simple own criteria, but users can modify the scheme by editing score_seq function in the source code


[Software Dependencies]
The program needs following software to be installed (in any path):
bowtie (http://bowtie-bio.sourceforge.net/index.shtml)
RNAfold (in ViennaRNA package) (https://www.tbi.univie.ac.at/RNA/)

If both software can be directly executed via terminal (bowtie-build, bowtie, RNAfold), you don't need to specify any parameters in the script
If you have installed these software but cannot directly be executed via terminal, you need to specify the path containing the executables in the script


[Parameters]
  -h, --help
  	show the help screen and exit


[Core]  	
  -r --reference
  	reference genome file location(.fa, .fasta) (default : 'ref.fa')
  	ex) /home/myname/(...)/my_ref.fa
  	if no path is specified, it tries to load "ref.fa" in the current folder

  -i --input
  	input small RNA file location(.fa, .fasta) (default : 'smrna.fa')
  	ex) /home/myname/(...)/my_smrna.fa
  	if no path is specified, it tries to load "smrna.fa" in the current folder

  	the file format should look like the following:
    	>t0000035 3234
    	GAATGGATAAGGATTAGCGATGATACA
    	>t0000072 1909
    	TTGCAGTATGTAGGAAATCAAAACGTTC
    	(i.e. >name (tab) reads (newline) seq)
  	most raw RNA-seq data formats are similar to this, but you may need to pre-process the file to match the format
  	i'm planning to implement automatic formatting in the future

  -o --output (default : 'result' subfolder)
  	output file location
  	ex) /home/myname/(...)/myoutput
  	the program generates output files in the specified output location
  	if the path do not exist, the program tries to make the directory
  	if no path is specified, results are generated in "result" subfolder in the current folder

  -w --bowtiepath
  	bowtie location (default: command "bowtie" and "bowtie-build")
  	ex) /home/myname/(...)/bowtie-1.1.2
  	you need to specify the path CONTAINING the executables "bowtie" and "bowtie-build"
  	if no path is specified, the program tries to execute with terminal command "bowtie" and "bowtie-build"

  -p --RNAfoldpath
  	RNAfold location (default : command "RNAfold")
  	ex) /home/myname/(...)/ViennaRNA/bin
  	you need to specify the path CONTAINING the executable "RNAfold"
  	if no path is specified, the program tries to execute with terminal command "RNAfold"


[pre-miRNA Discovery Criteria]
  -l --minlength (default : 18)
  	min. length of mature miRNA

  -L --maxlength (default : 26)
  	max. length of mature miRNA 

  -u --multloci (default : 20)
  	max. multiple loci of miRNA matches to reference genome
  	if the raw RNA-seq read matches to genome with too many locations, it should be considered as repeat, not miRNA candidate
  	this value is passed to bowtie, and only reads corresponding to this criterion will be mapped

  -d --distance (default : 35)
  	max. distance between miRNA-miRNA* of precursor hairpin loop
  	if you want to discover "bigger" pre-miRNAs (like plant miRNA or non-canonical miRNA), you can increase this value, but the program will tend to prefer "bigger" pre-miRNAs and "smaller" pre-miRNAs might be discarded

  -a --arm (default : 10)
    length of arms(both ends) of precursor from mature miRNA
    also called as "flank", this value specifies the length of potential drosha-cut sites
    the program assumes this value to calculate MFE of putative pre-miRNAs

  -f --mfe (default : 18)
  	min. abs. MFE for valid miRNA precursor
  	pre-miRNAs whose absolute MFE value are larger than this value (-18 kcal/mol) are considered as novel pre-miRNAs
  	the program "picks" pre-miRNA with the best normalized MFE near the searching location
  	pre-miRNAs not passing this criterion are considered "unstable" and discarded  


[Mature miRNA Scoring Criteria]
the site of pre-miRNA passing all the criteria with the best score is selected as putative mature miRNA-miRNA* duplex
if you want to allow more non-canonical miRNA-miRNA* duplex, you may want to increase these values

  -m --serialmismatch (default : 2)
  	max. serial mismatch of miRNA-miRNA* duplex
  	if any single mismatch site is "too large", it is discarded

  -M --multmismatch (default : 2)
  	max. multiple mismatches of miRNA-miRNA* duplex
  	if the number of mismatch sites is "too much", it is discarded

  -b --serialbulge (default : 2)
  	max. serial bulge of miRNA-miRNA* duplex
  	if any single bulge is "too large", it is discarded

  -B --multbulge (default : 2)
  	max. multiple bulges of miRNA-miRNA* duplex
  	if the number of bulges is "too much", it is discarded
    

[Performance & Memory Efficiency]
  -t --thread (default : automatic, all threads of the system)
  	number of CPU threads
  	if you want to run the program with fewer cpu threads (and less memory), specify the value

  -s --step (default : 2)
  	step size of RNAfold precursor extend loop
  	the program progressively calculates MFE of putative pre-miRNAs extending the length of the sequence with this specified value
  	for example, with default settings, the program starts with [ (minlength + arm) * 2 + distance ] and extends the length "left" or "right" automatically
  	if you increase this value, the program calculates with RNAfold more sparsely and can achieve faster performance, but it could miss the "optimal" choice of pre-miRNA
  	if you have set larger "distance" value, you may want to increase this value too  

  -c --mincount (default : automatic, number of unique RNA-seq reads / 500000)
  	min. read count of smRNA seq displayed at result_mature.txt
  	the program allocates copy of loaded map data to all threads to maximize performance, so if map data is too large, memory capacity could be an issue, especially if the map file exceeds several hundred MBs
  	also, if the RNA-seq file has many reads with small count (1~2), result_mature.txt might be hard to read, and result would show less probable miRNA candidates
  	by default, this value is defined in a heuristic way to achieve stable memory usage and reliable result
  	if you specify this value, genome-mapped reads with count lower than this value will be considered "insignificant", and will not be used for mature miRNA finding


[Miscellaneous & Additional Features] : recommended to keep default values
  --batch_size (default : number of CPU threads * 3)
  	size of input data "chunk" for internal processing, 
  	it can be arbitrary value, but (number of cpu threads) * n is recommended
  	if you increase the size, you COULD earn faster performance (but not always), but the memory usage will be increased
  	optimal value is depending on the system you are running, so you may want to experiment

  --plot (default:false)
    draw simple bar plot of length distribution of genome-aligned RNA sequence
    if you set this value to true, the program generates "length-distribution.png" file to output path
    need matplotlib python package to be installed
                        



