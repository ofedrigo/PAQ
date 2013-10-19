PAQ
===

Documentation for PAQ Version 1.0    August 24, 2007 

* Copyright 1999 by Prasith Baccam, Robert J. Thompson, Olivier 
* Fedrigo, James L. Cornette, and Susan Carpenter, Iowa State 
* University, Ames, IA 50011 
* This computer program and documentation may be freely copied 
* and used by anyone, provided no fee is charged for it. 



Contents: 

Overview 
Files in this Package 
Input File Format 
Running the Program 
Program Menu 
Literature Cited 


Overview

There are many tools available for studying molecular evolution. Phylogenetic programs like PHYLIP (Felsenstein, 1993) and PAUP (Swofford, 1999) have been increasingly utilized to analyze the evolutionary relationship among genetic sequences. Genetic conversions and recombinations, however, may not be modeled well under the bifurcating assumptions of phylogenetic tree construction. Progress is being made on inferring genetic relationships when recombination has occurred. Some examples of these methods include spectral analysis (Charleston, 1998), split decomposition (Dopazo et al., 1993), and median networks (Guenoche, 1986). Another method of finding similarity among sequences is to describe the multidimensional information in low dimensions. Examples of this method include principal components analysis (PCA) (Mardia et al., 1979), principal coordinates analysis (PCOORD) (Gower, 1966, Higgins, 1992), and multivariate statistical sequence analysis (van Heel, 1991). Some data sets, however, may not be analyzed well by this method as the best 3-dimensional representation of the data may not capture even the majority of information in the data set. In addition, the algorithms used to reduce the dimensionality of the data set may be non- intuitive to some users, and the results of the analysis may likewise be non- intuitive.  The clusters that are identified can be used to study the molecular evolution of viral quasispecies (Eigen, 1993) or other sequences that exhibit  a quasispecies behavior. 

PAQ is written in ANSI C++.  It should be possible to compile and run the program on a wide variety of platforms and operating systems. 

If you find PAQ useful, and if it is used for analyses that end up as part of a publication, we would appreciate it if you would cite the following reference. The program is mentioned in this paper. 

PAQ: partitioning analysis of quasispecies, Bioinformatics, 2001, 17(1):16-22. 


Files in this Package

UNIX
  PAQ.cpp	This code was written in C++ under the following operating systems:  Digital Unix 3.2x, DEC Ultrix 4.3x, SGI Irix 5.3 and compiled with GNU project C++ compiler, version 2.7.2 
Windows Operating Systems (98, NT, and XP tested)
  PAQ.exe 	This version of PAQ is a windows application.  The program has a graphical user interface (GUI) rather than being menu-driven. 

Macintosh 	We are having some difficulties with memory allocation with this version.  Please be patient, and check out the PAQ homepage for the release of the Macitosh version of the program.  Macintosh users running Mac OSX may be able to compile and run the unix version (I’ve been told, but have not tested).  Alternatively, Macintosh users may run the Windows version of PAQ by using the PC emulation program, Virtual PC.

Sample data: 
  hiv.paq	This data set was derived from the paper by Shapshak et al.  Briefly, HIV-1 gp120 envelope (V1-V5) sequences were isolated from different regions of the brain.  The data includes an outgroup sequence (HIVHXB3) and 63 sequences from Shapshak et al., with accession numbers M14100, and AF125810 - AF125874, respectively. The 63 sequences were coded as follows: P###A#G, where the first 3 numbers indicate the patient number, A represents the region of the brain the sequences were isolated - F=frontal lobe, G=basal ganglia, M=medial temporal lobe, and T=non-medial temporal lobe, the next number is the clone identifier, and G represents the genetic target the sequences were amplified from - D=DNA and R=RNA.


Input File Format

The input file has the same format as input files for the phylogenetic program PHYLIP.  The input file should be saved as a text file, and the file name must end in .paq for the program to recognize the input file.  The format requires that the first line of the file consists of two integers; the first integer representing the number of sequences in the file and the second integer representing the number of characters in each sequence. Following the first line are the sequences. 

For each sequence, the first 11 characters are reserved for the sequence name, and the sequence data follows.  The sequences should be non-interleved (not segmented but the complete sequence).  Carriage returns should only be used at the end of each sequence.   

Data restrictions: 

The program will accept nucleotide or amino acid sequences.  The program is case-sensitive, but there should be no problems if you are consistent with your case usage.  Each sequence must have the same number of characters, and it is assumed that the sequences are aligned.  Ambiguous characters should be indicated using 'n' or 'N'.  Gaps should be indicated by '-' or '.'. Comparisons between sequences do not consider gaps or ambiguous characters when calculating distance between sequences. 


Running the Program

UNIX 
	The program file should reside in the same directory as the data file. After compiling the program, execute the program.  The first prompt is for the file name of the input file.  Next, the user can decide whether to consider or ignore gap positions in the calculation of distance between sequences.  Once that is entered, some information about the input file is displayed, including: the minimum number of genetic changes (all genetic changes, nucleotide or amino acid are referred to as nucleotide changes in the program), the maximum number of genetic changes, and the base frequencies in the data set.  The program also indicates a region where no sequences differ by a number of genetic changes. This information may be very useful for a starting radius value. 

Windows 
	The executable file will not operate properly if it is on the Desktop, so it should be on the local hard drive.  Once the program has been started, you must open the input file (with .paq ending to the file name).  You can choose to consider or ignore gap positions in the calculation of distance between sequences.  The genetic data can be seen if the “Data Matrix” button is clicked, while the “Genet dist.” button will show the number of difference among the sequences.  Clustering can b examined by inputting a range of starting and ending radius values and clicking the “Group” button.  Once the “Group” button has been clicked, the main window will indicate the radius, the number of clusters found, and for each cluster, the center sequence and its neighbors are shown.  Any variants not contained in any of the clusters are also indicated.  The right and left arrows by the radius value allow the user to change the radius value to see the changes in clustering that result.  The order of the sequences can be changed to correspond to the clustering found by the various radius values by clicking on the “Nt Matrix new” button.  The scale of the window can be changed to help view the genetic distance matrix.


Program Menu

Here are the options for the program: 
1. Display the genetic distance matrix 
2. Display the clusters for a range of radius values
3. Display the clusters for a single radius value 
4. Search a cluster for sub-clusters 
5. Display a rearranged genetic distance matrix 
6. Display all potential clusters for a single radius value 
7. Display the average distance from all centers in a range of radius values 
8. Display all genetic distance from a single sequence  
9. Display the genetic distance between two sequences 

Option 1:  This option displays a matrix of genetic distances between all sequences.  The matrix is symmetrical with zeroes along the diagonal. 

Option 2:  For a range of radius values, the program automatically finds distinct clusters (not identical or subsets of other clusters). 

Option 3:  Same as option 2, only for a single radius value. 

Option 4:  First you enter a radius value, and the program automatically finds the distinct clusters.  You can choose which cluster you want to investigate further.  Then you can choose a smaller radius value, and the program automatically finds the distinct clusters from the sequences which were contained in the previous cluster you picked.  This allows you to search for possible sub-clusters within larger clusters. 

Option 5:  The genetic distance matrix (from option 1) contains a lot of data which may not be detectable to the eye.  Here, you choose a radius value, the program automatically finds the distinct clusters, and re-arranges the genetic distance matrix so that sequences in the same cluster are next to each other in the distance matrix.  The result is a matrix which contains blocks of small genetic distances along the matrix diagonal and larger distances elsewhere. This option allows the user to infer possible clusters and sub-clusters. 

Option 6:  This option displays all the potential clusters.  Each sequences is used as the center of a cluster, and the number and identity of other sequences within the cluster (neighbors) are shown.  The average distance from the center and all the neighbors is also shown.  This option is useful when the data clusters overlap. 

Option 7:  For a range of radius values, the average distance from all centers is displayed.  This output is helpful in determining the best radius values for each cluster.  **Note: in the "optimal" clustering of sequences, the clusters may each have different radius values. 

Option 8:  For a specific sequence, the distance from that sequence to all other sequences is displayed. 

Option 9:  For two specified sequences, the distance between the two sequences is displayed.


 
Literature Cited

Charleston, M.A. (1998) Spectrum: spectral analysis of phylogenetic data. Bioinformatics, 14,98-99. 

Dopazo, J., A. Dress, and A. von Haeseler. (1993) Split decomposition: a technique to analyze viral evolution. Proc Natl Acad Sci USA, 90,10320-10324. 

Eigen,M. (1993) Viral Quasispecies. Sci. Am., 269, 42-49. 

Felsenstein,J. (1993) PHYLIP: Phylogeny Inference Package, Version 3.5. University of Washington, Seattle, WA. 

Gower, J.C. (1966) Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika, 53,325-333. 

Guenoche, A. (1986) Graphical representation of a boolean array. Computers and the humanities, 20,277-281. 

Higgins,D.G. (1992) Sequence ordinations: a multivariate analysis approach to analysing large sequence data sets. Comput. Applic. Biosci., 8, 15-22. 

Mardia, K.V., J.T. Kent, and J.M. Bibby. (1979) Multivariate analysis. Academic Press, New York. 

Shapshak,P., Segal,D.M., Crandall,K.A., Fujimura,R.K., Zhang,B.T., Xin,K.Q., Okuda,K. Petito,C.K., Eisdorfer,C. and Goodkin K. (1999) Independent evolution of HIV type 1 in different brain regions. AIDS Res. Hum. Retroviruses, 15, 811- 820. 

Swofford,D.L. (1999) PAUP*. Phylogenetic Analysis Using Parsimony (*and Other Methods). Version 4. Sinauer Associates, Sunderland, Massachusetts. 

van Heel,M. (1991) A new family of powerful multivariate statistical sequence analysis techniques. J. Mol. Biol., 220, 877-887.
