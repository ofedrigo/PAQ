// PAQ.cpp    create a matrix of number of n.t. differences between 
//	      sequences.  inputs the radius for clustering, and lists
//	      the sequences that are w/in the radius, etc.
//	      Note: the average distance is now calculated as the sum
//	      of the distance squared, normalized by the number of nbrs
// edited: 8/9/2000

#include <iostream.h>
#include <iomanip.h>		// for setprecision()
#include <fstream.h>
#include <ctype.h>              // for isspace()
#include <math.h>

const int inputlth = 20;	// length of input file name
const int sequences = 200;	// approx. # of sequences
const int maxradius = 300;	// approx. maximum radius
const int sites = 1500;		// approx. # of sites/sequence	
const int num_char = 128;       // number of chars
const char min_char = ' ';      // first printable char
const char max_char = 'Z';      // last printable char
const int skip = 10;		// length of sequence name (ignored)

typedef int intmatrix[sequences] [sequences];	// matrix of nt changes
typedef char charmatrix[sequences] [sites];	// data (nt) matrix
typedef int infomatrix[maxradius][sequences][sequences];
typedef float avdistmatrix[maxradius][sequences];  // matrix of av. dist.
typedef int clusmat[maxradius][sequences][sequences]; // matrix of clusters
typedef float clusdist[maxradius][sequences];	// av. dist. per cluster

/**************************************************************************
infomatrix[radius][center][i]:  when i=0, this indicates the number of
variants contained in the cluster with the given center and radius.
the components where i>0 contains the identity of the variants within
the cluster.
avdistmatrix[radius][center]:  contains the average distance from all
variants within the cluster to the given center and radius.
clus[radius][i][j]:  when i=j=0, this indicates the number of clusters
present with the given radius.  when i=1, j=2,3,... indicates all the
variants that belong to any cluster (and thus all variants that do not
belong to any cluster).  when i>1 and j>1, this identifies the clusters
(i=2 means the 1st cluster) and the variants within the clusters (j).
clusdist[radius][cluster#]:  this contains the average distance for the
cluster, given the radius.
*************************************************************************/

void changeMatrix(int& column, int& numsites, int& maxchange, int& row,
		int& numseqs, int& j, char datamatrix[sequences][sites],
		int& numdiff, int matrix[sequences][sequences]);
void GAPchangeMatrix(int& column, int& numsites, int& maxchange, int& row,
                int& numseqs, int& j, char datamatrix[sequences][sites],
                int& numdiff, int matrix[sequences][sequences]); 
void printMatrix(int matrix[sequences][sequences], int& row, 
		int& numseqs, int& column, int frequcount[]);
void AvDistanceInfo(int& maxchange, int& iRange1, int& iRange2,
		int& row, int& count,int& dist, int& column,
		int info[maxradius][sequences][sequences], 
		int matrix[sequences][sequences], 
		float avdist[maxradius][sequences], int& numseqs);
void Sort(int iRange1, int iRange2, int& numseqs, 
		int info[maxradius][sequences][sequences], 
		float avdist[maxradius][sequences],
		int clus[maxradius][sequences][sequences]);
void Sortsmall(int iRange1, int iRange2, int& numseqs,
                int info[maxradius][sequences][sequences],
                float avdist[maxradius][sequences],
                int clus[maxradius][sequences][sequences],
		int origclus[sequences]);
void SubCluster(int& numseqs, int& numsites, int& iSubCluster, int& j, 
		int& row, int& column, int& numdiff, int& x, int& iRange1,
		int& newseqs, char smalldata[sequences][sites],
		char datamatrix[sequences][sites],
		int origclus[sequences], 
		int clus[maxradius][sequences][sequences],
		int smallchange[sequences][sequences]);
void Resort(int clus[maxradius][sequences][sequences],
		char datamatrix[sequences][sites],
		int& numsites, int& numseqs, int& iRange1,
		char sorteddata[sequences][sites], 
		int seqorder[sequences],
		int sortedmat[sequences][sequences]);
void CompareClust(float avdist[maxradius][sequences], int& radius, int& count,
		int findclus[sequences][sequences], int& maxvar,
		int info[maxradius][sequences][sequences],int& maxcenter,
		float clusad[maxradius][sequences], int tempvar[maxradius],
		int& numseqs); 
void AllClusters(int radius, int& numseqs,
                int info[maxradius][sequences][sequences],
                int matrix[sequences][sequences]);
void AvDistance(int radius, int& numseqs,
		int info[maxradius][sequences][sequences],
		int matrix[sequences][sequences]);
void SequenceNT(int sequence, int matrix[sequences][sequences],
                int numseqs);
void CompareTwo(int sequence, int sequence2,
                int matrix[sequences][sequences]);
void Menu();


void zeromatrix( intmatrix );	// zero out components of intmatrix
void zerofreq(int []);		// zero out freq. of each nt
void zeroinfo( infomatrix );	// zero out entries of infomatrix
void zeroavdist( avdistmatrix );  // zero out entries of avdistmatrix
void Print(const int[]);	// print out freq. of each nt

int main()
{
intmatrix	matrix;			// our nt change matrix
intmatrix	findclus;		// to define clusters
intmatrix	sortedmat;		// sorted nt change matrix
intmatrix	smallchange;		// 
charmatrix	datamatrix;		// matrix of nt data
charmatrix	sorteddata;		// sorted data matrix
charmatrix	smalldata;		// subset of datamatrix
char		amino;			// amino acid read
clusmat		clus;			// cluster data
int		row;			// row counter
int		column;			// column counter
int		seqorder[sequences];	// order of seq in resorted mat
int		i;			// counters
int		j;
int		k;
int		l;
int		x;
int		y;
int		freqcount[num_char];    // freq. counts for each nt
ifstream	infile;			// input file
int		numdiff;		// # of diff. between 2 sequences
char		file[inputlth];		// length of input file name
int		numseqs;		// # of sequences
int		newseqs;		// # seq. in smalldata matrix
int		numsites;		// # of sites
int		radius;			// radius for clustering
int		count;			// count of neighbors
int		dist;			// distance from center
int		minchange;		// min # of nt changes in data
int		maxchange;		// max # of nt changes in data
int		novar;			// # of var. in cluster
int		tempvar[maxradius];	// temporary array
int		origclus[sequences];	// order of orig. sequences
int		maxvar;			// max # of variants
int		maxcenter;		// center with max # nbrs
int		cluster;		// # of clusters
int		different;		// boolean variable
int		numclus;		// # of clusters present
int		notclus;		// # var. not in any clusters
char		UserInput;
int             iRange1;
int             iRange2;
bool		test = true;
infomatrix	info;
avdistmatrix	avdist;
clusdist	clusad;
int		iSubCluster;
int		sequence;
int		sequence2;

// zero out all matrices
zeromatrix( matrix );
zerofreq(freqcount);
zeroinfo(info);


// request the input file
cout << "What is the input file name?   : ";
cin >> file;
infile.open(file);
if (!infile)
{
  cout << "*** Cannot open the input file! ***" << endl;
  return 1;
}
infile >> numseqs;		// read in # of sequences
infile >> numsites;		// read in # of sites
infile.ignore(100,'\n');	// skip to next line of input

// read data from input file into datamatrix
for (row = 0; row <=(numseqs-1); row++)  // loop for # of sequences
{
infile.ignore(skip,'\n');          //ignore the sequence name
for (column= 0; column <=(numsites-1); column++)	
				// loop for length of sequence
// this loop counts the number of times a nt or amino acid occurs
// this is used later to calculate nt or aa frequencies
{
infile.get(amino);
if (isprint(amino))
freqcount[amino]++;		// increment count for appropriate nt
datamatrix[row][column] = amino;	// copy nt to data matrix
}
infile.ignore((numsites+400),'\n');	// jump to next row
}

char cGapDec;
cout << "Would you like to include gaps in the change matrix?\n";
cout << "Type 'y' for yes or 'n' for no.\n";
cin >> cGapDec;
if(cGapDec =='n' || cGapDec == 'N') {
changeMatrix(column, numsites, maxchange, row, numseqs, j,
             datamatrix, numdiff, matrix);
}
else {
GAPchangeMatrix(column, numsites, maxchange, row, numseqs, j,
             datamatrix, numdiff, matrix);
}
Print(freqcount);

	while(test) {
	Menu();	
	cin >> UserInput;
	
	switch(UserInput) {
	case '1':
	  printMatrix(matrix, row, numseqs, column, freqcount);  
	break;
	case '2':
	cout << "Enter the starting radius value:  ";
	cin >> iRange1;
	cout << endl;
	cout << "Enter the ending radius value:    ";
	cin >> iRange2;
	AvDistanceInfo(maxchange, iRange1, iRange2, row, count, dist,
		column,info, matrix, avdist, numseqs);
	Sort(iRange1, iRange2, numseqs, info, avdist, clus);
	break;
	case '3':
	cout << "Enter the  radius value:  ";
	cin >> iRange1;
	iRange2 = iRange1;
	AvDistanceInfo(maxchange, iRange1, iRange2, row, count, dist,
		column,info, matrix, avdist, numseqs);
	Sort(iRange1, iRange2, numseqs, info, avdist, clus);
	break;
	case '4': 
	cout << "\nEnter the radius to identify clusters:  ";
        cin >> iRange1;
	iRange2 = iRange1;
	AvDistanceInfo(maxchange, iRange1, iRange2, row, count, dist,
		column,info, matrix, avdist, numseqs);
	Sort(iRange1, iRange2, numseqs, info, avdist, clus);
	cout << "\n\n\nEnter the cluster number you would like to ";
	cout << "further investigate:  ";
	cin >> iSubCluster;
	SubCluster(numseqs, numsites, iSubCluster, j, row, column, 
		numdiff, x, iRange1, newseqs, smalldata, datamatrix, 
		origclus, clus, smallchange);
	zeromatrix ( smallchange );
	if(cGapDec == 'n' || cGapDec == 'N') {
	changeMatrix(column, numsites, maxchange, row, newseqs, j,
		smalldata, numdiff, smallchange);
	}
	else {
	GAPchangeMatrix(column, numsites, maxchange, row, newseqs, j,
        smalldata, numdiff, smallchange);
	}
	cout << "\nNow enter a smaller radius value:  ";
	cin >> i;
	cout << "For radius = " << iRange1 << ", cluster #";
	cout << iSubCluster << " had center ";
	cout << clus[iRange1][iSubCluster][0];
	cout << " and contained " << clus[iRange1][iSubCluster][1];
	cout << " other variants," << endl;
	cout << "including:  " ;
	for (int j=2; j<=clus[iRange1][iSubCluster][1]+2; j++)
	  {
	    cout << clus[iRange1][iSubCluster][j] << " ";
	  }
	cout << endl;
	cout << "\nFurther analysis of cluster #" << iSubCluster;
	cout << " found:" << endl;
	iRange1 = i;
	iRange2 = iRange1;
	AvDistanceInfo(maxchange, iRange1, iRange2, row, count, dist,
                column,info, smallchange, avdist, newseqs);
        Sortsmall(iRange1, iRange2, newseqs, info, avdist, clus,
		origclus);
	break;
	case '5':
	cout << "\nEnter a radius value:  ";
        cin >> iRange1;
	iRange2 = iRange1;
        AvDistanceInfo(maxchange, iRange1, iRange2, row, count, dist,
		column,info, matrix, avdist, numseqs);
        Sort(iRange1, iRange2, numseqs, info, avdist, clus);
	Resort(clus, datamatrix, numsites, numseqs, iRange1, 
		sorteddata, seqorder, sortedmat);
	printMatrix(sortedmat, row, numseqs, column, freqcount);
        break;
        case '6':
        cout << "\nEnter the radius value:  ";
        cin >> iRange1;
        iRange2 = iRange1;
        AvDistanceInfo(maxchange, iRange1, iRange2, row, count, dist,
                column,info, matrix, avdist, numseqs);
        AllClusters(iRange1, numseqs, info, matrix);
        break;
        case '7':
        cout << "\nEnter the starting radius value:  ";
        cin >> iRange1;
	cout << endl;
        cout << "Enter the ending radius value:    "; 
        cin >>  iRange2;
	for(int z=iRange1; iRange1<=iRange2; iRange1++) {
        cout <<"\n\nFor the radius "<<iRange1<< ":\n\n";
        AvDistanceInfo(maxchange, iRange1, iRange1, row, count, dist,
                column,info, matrix, avdist, numseqs);
        AvDistance(iRange1, numseqs, info, matrix);
        }
        break;
        case '8':
        cout << "\nEnter the sequence number to be analyzed:  ";
        cin >> sequence;
	cout << endl;
        SequenceNT(sequence, matrix,numseqs);
        break;
        case '9':
        cout << "\nEnter the first sequence number:   ";
        cin >> sequence;
	cout << endl;
        cout << "Enter the second sequence number:  ";
        cin >> sequence2;
	cout << endl;
        CompareTwo(sequence, sequence2,matrix);
        break;
        case '0':
        test = false;
        break;
        default:
        cout << "NOT AN OPTION!\n";
        cout << "\n\nPlease specify another option.\n";
        cin >> UserInput;
	break;
	}
	}

return 0;
}

void Menu() 
{
cout <<"\n\n\n\n\n\n";
cout << "	<<<<<<<<<	MENU	>>>>>>>>>>>>>>\n\n\n\n";
cout << "\n	1. Display the genetic distance matrix.\n\n";
cout << "	2. Display the groups for a range of radius values.\n\n";
cout << "	3. Display the groups for a single radius value.\n\n";
cout << "	4. Search a group for sub-groups.\n\n";
cout << "	5. Display a rearranged genetic distance matrix.\n\n";
cout << "	6. Display all potential groups for a single\n";
cout << "	   radius value.\n\n";
cout << "	7. Display the average distance from all centers\n";
cout << "	   in a range of radius values.\n\n";
cout << "	8. Display all the genetic distance from a\n";
cout << "	   single sequence.\n\n";
cout << "	9. Display the genetic distance between two sequences.\n\n";
cout << "	0. Quit\n\n\n";
cout << "	Choose an option from the menu:   "; 
} 


void GAPchangeMatrix(int& column, int& numsites, int& maxchange, int& row,
                int& numseqs, int& j, char datamatrix[sequences][sites],
                int& numdiff, int matrix[sequences][sequences])
{
int maxgap;
int change[1000];
int columns;
int gaplth;
int x;
int minchange;
// create a matrix, where cells have # of differences between sequences


for (int i=0; i<=999; i++)
    change[i]=0;

// create a matrix, where cells have # of differences between sequences
maxchange = 0;     // records the largest # of nt changes
minchange = 1000;  // records the smallest # of nt changes
for (row= 0; row <=(numseqs-2); row++)  //loop for all rows (-1)
{
for (j=row+1; j<=(numseqs-1); j++)   //loop for 2nd comparison row
 { numdiff=0;
   for (column=0; column<=(numsites-1); column++)   //loop for columns
    {
            if (datamatrix[row][column] != datamatrix[j][column])
                {  numdiff++;
                }
         }
	
  matrix[row][j]=numdiff;
  matrix[j][row]=numdiff;
  change[numdiff]++;
  if (numdiff >0)
    {
      if (numdiff < minchange)
      minchange = numdiff;
    }
  if (numdiff > maxchange)
        maxchange = numdiff;
 }
}
maxgap = 0;
gaplth = 0;
x = 0;
for (int i=minchange; i<=maxchange; i++)
  {
    if (change[i] == 0)
        x++;
    if (x>gaplth)
      {
        gaplth = x;
        maxgap = i;
      }
    if (change[i] > 0)
        x=0;
  }

cout << "The min number of nt changes was " << minchange << endl;
cout << "The max number of nt changes was " << maxchange << endl;
cout << "The largest region with no nt change was of length " << gaplth;
if (gaplth > 0)
  {
    cout << " and occurred between " << (maxgap-gaplth+1);
    cout << " and " << maxgap << " (inclusive)" << endl;
  }
else
  cout << endl;
}







void changeMatrix(int& column, int& numsites, int& maxchange, int& row,
		int& numseqs, int& j, char datamatrix[sequences][sites],
		int& numdiff, int matrix[sequences][sequences])
{	
int maxgap;
int change[1000];
int columns;
int gaplth;
int x;
int minchange;
// create a matrix, where cells have # of differences between sequences


for (int i=0; i<=999; i++)
    change[i]=0;

// create a matrix, where cells have # of differences between sequences
maxchange = 0;     // records the largest # of nt changes
minchange = 1000;  // records the smallest # of nt changes
for (row= 0; row <=(numseqs-2); row++)  //loop for all rows (-1)
{
for (j=row+1; j<=(numseqs-1); j++)   //loop for 2nd comparison row
 { numdiff=0;
   for (column=0; column<=(numsites-1); column++)   //loop for columns
    { 
	if (datamatrix[row][column] != '*' && datamatrix[j][column] !='*')
	  {
           if (datamatrix[row][column] != '-' && datamatrix[j][column]!='-')
           {

            if (datamatrix[row][column] != datamatrix[j][column])
                {  numdiff++;
                }
           }
          }
        }

  matrix[row][j]=numdiff;
  matrix[j][row]=numdiff;
  change[numdiff]++;
  if (numdiff >0)
    {
      if (numdiff < minchange)
      minchange = numdiff;
    }
  if (numdiff > maxchange)
        maxchange = numdiff;
 }
}
maxgap = 0;
gaplth = 0;
x = 0;
for (int i=minchange; i<=maxchange; i++)
  {
    if (change[i] == 0)
        x++;
    if (x>gaplth)
      {
        gaplth = x;
        maxgap = i;
      }
    if (change[i] > 0)
        x=0;
  }

cout << "The min number of nt changes was " << minchange << endl;
cout << "The max number of nt changes was " << maxchange << endl;
cout << "The largest region with no nt change was of length " << gaplth;
if (gaplth > 0)
  {
    cout << " and occurred between " << (maxgap-gaplth+1);
    cout << " and " << maxgap << " (inclusive)" << endl;
  }
else
  cout << endl;
}

void SequenceNT(int sequence, int matrix[sequences][sequences],
                int numseqs)
{
for(int mycounter = 0; mycounter <= numseqs-1; mycounter++) {
        cout <<"\nNumber of changes from sequence number ";
	cout << setw(4) << mycounter+1 << " :  ";
        cout << matrix[mycounter][sequence-1];
        }

}

void CompareTwo(int sequence, int sequence2, int
                matrix[sequences][sequences])
{
cout << "\nThe number of changes between sequence number " << sequence
        << " \nand sequence number "<< sequence2 <<" : "
        << matrix[sequence-1][sequence2-1] << "\n";
}

void printMatrix(int matrix[sequences][sequences], int& row,
                int& numseqs, int& column, int freqcount[])
{
// print out the change matrix

cout << endl;
cout << "nt differences between sequences:" << endl;
for (row=0; row<=(numseqs-1); row++)
{
for (column=0; column<=(numseqs-1); column++)
{
        cout << setw(4) << matrix [row] [column];
}
cout << endl;
}

cout << endl;
}

void AvDistanceInfo(int& maxchange, int& iRange1, int& iRange2,
                int& row, int& count,int& dist, int& column,
                int info[maxradius][sequences][sequences],
                int matrix[sequences][sequences],
                float avdist[maxradius][sequences], int& numseqs)
{
int j;
int k;

// a loop to fill up the info and avdist matrices with a range of radius
for (int radius =iRange1; radius<=iRange2; radius++)
{
 for (row = 0; row <= (numseqs-1); row++)
  {
        count = 0;
        dist = 0;
        for (column = 0; column<= (numseqs-1); column++)
         {
           if (matrix[row][column] <= radius)
                {
                  count++;
		  info[radius][row][count]=(column+1);
                  dist = dist + (matrix[row][column]*matrix[row][column]);
                }
         }
        if (count != 1)
         {
          info[radius][row][0]=count-1;  // # of nbrs in cluster
          avdist[radius][row] = (float(dist)/float(count-1));
         }
  }
}   // end for loop to fill the info and avdist matrices
}

void AllClusters(int radius, int& numseqs,
                 int info[maxradius][sequences][sequences],
                int matrix[sequences][sequences])
{
// a loop to test neighbors w/in a radius
        if (radius != 0)
        {
        for (int row = 0; row <= (numseqs-1); row++)
        {
        cout << endl;
        cout <<"Center #"<< (row+1) <<" has these neighbors: "<< endl;
        int count = 0;
        int dist = 0;
        for (int column = 0; column<= (numseqs-1); column++)
         {
           if (matrix[row][column] <= radius)
                { cout << (column+1) << "  ";
                  count++;
                  dist = dist + (matrix[row][column]*matrix[row][column]);
                }
         }
        cout << endl;
        cout << "Center #" << (row+1) << " has " << (count-1);
        cout << " neighbors";
        if (count != 1)
         {
	  cout << ", average distance from center = ";
          cout << float(dist)/float(count-1) << endl;
         }
        else
          cout << endl;
        }
        }
}

void AvDistance(int radius, int& numseqs,
                 int info[maxradius][sequences][sequences],
                int matrix[sequences][sequences])
{
// a loop to test neighbors w/in a radius
        if (radius != 0)
        {
        for (int row = 0; row <= (numseqs-1); row++)
        {
        int count = 0;
        int dist = 0;
        for (int column = 0; column<= (numseqs-1); column++)
         {
           if (matrix[row][column] <= radius)
                { 
                  count++;
                  dist = dist + (matrix[row][column]*matrix[row][column]);
                }
         }
        cout << "For center #" << setw(3) << (row+1);
        cout << ", average distance from center = ";
        if (count != 1)
         {
          cout << float(dist)/float(count-1) << endl;
         }
        else
          cout << "0" << endl;
        }
        }
}

void Sort(int iRange1, int iRange2, int& numseqs, 
                int info[maxradius][sequences][sequences],
                float avdist[maxradius][sequences],
		int clus[maxradius][sequences][sequences])
{
int x;
int y;
int j;
int k;
int findclus[sequences][sequences];
int tempvar[sequences];
float clusad[maxradius][sequences];
int maxcenter;
int maxvar;
int cluster;
int count;
int radius;
int notclus;
int overlap;
//  sorts the # of variants in clusters in descending order.  finds
//  the clusters for the radius values

for (radius = iRange1; radius<=iRange2; radius ++)
{
  for (int i=0; i<=(numseqs-1); i++)
    {
        tempvar[i] = info[radius][i][0];
    }

// tempvar[i] holds the # of variants w/in the cluster w/ center I

  zeromatrix (findclus);
  maxvar = 99;
  cluster = 0;
  count = 0;
  while (maxvar >0)
    {
        x=0;
        y=0;
        maxvar = 0;
        maxcenter = 0;
        for (int center = 0; center <= (numseqs-1); center++)
          {
            if (tempvar[center]!=0)
              {
               if (tempvar[center] > maxvar)
                {
                  maxvar = tempvar[center];
                  maxcenter = center;
                }
            if (tempvar[center] == maxvar)
		{
                if (avdist[radius][center] < avdist[radius][maxcenter])
                   maxcenter = center;
              }
          }  // end the for loop
    }  // found the cluster w/ max. # of variants inside

      CompareClust(avdist, radius, count,findclus, maxvar,
                info, maxcenter,clusad, tempvar, numseqs);
//I took out two brackets that let it compile here
//Rob, 11/12/99
}
clus[radius][0][0] = (count-1);  // # of distinct clusters

cout << endl;
cout << "For radius = " << radius << ", there were ";
cout << (count-1) << " distinct groups present." << endl;

for (int i=1; i<=(count-1); i++)
  {
    clus[radius][i][0] = findclus[i][0];  // center of cluster
    clus[radius][i][1] = findclus[i][1];  // # of nbrs in cluster
    for (int j=2; j<=findclus[i][1]+2; j++)
{
        clus[radius][i][j] = info[radius][findclus[i][0]-1][j-1];
      }
  }

for (int i=1; i<=(count-1); i++)
  {
    cout << "Group #" <<  i << " has center = ";
    cout << clus[radius][i][0];
    if (clus[radius][i][1]>0)
      {
        cout << ", " << clus[radius][i][1];
        cout << " neighbors," << endl;
	cout << "and average distance from center = ";
        cout << avdist[radius][clus[radius][i][0]-1] << endl;
      }
    else
      {
	cout << ", and " << clus[radius][i][1];
        cout << " neighbors." << endl;
      }
    cout << "This group contains:  ";
    for (int j=2; j<=findclus[i][1]+2; j++)
      {
        cout << clus[radius][i][j] << " ";
      }
    cout << endl;
  }

notclus = 0;
for (int i=1; i<= clus[radius][0][0]; i++)
 {
for (int j=2; j<= numseqs+1; j++)
      {
//      cout << findclus[i][j] << " ";
        if (findclus[i][j]>0)
          {
            findclus[0][j]++;
            notclus++;  // indicates whether any var. not in clus.
          }
      }
// cout << endl;
  }

notclus=0;  // # of var. not in any clusters
overlap=0;  // indicates if any overlap
for (j=2; j<= numseqs+1; j++)
  {
    if (findclus[0][j]==0)
      notclus++;
    if (findclus[0][j]>1)
      overlap++;
  }
if (overlap>0)
  {
    cout << "Of these variants, the following belonged to more ";
    cout << "than one group: " << endl;
    for (j=2; j<= numseqs+1; j++)
      {
        if (findclus[0][j]>1)
          cout << j-1 << " ";
      }
    cout << endl;
  }
cout << "There were " << notclus ;
cout << " variants (" << float(100*notclus/numseqs);
cout << "%) that were not contained in any groups:" << endl;
if (notclus>0)
  {
    k = 2;
    for (j=2; j<= numseqs+1; j++)
      {
        if (findclus[0][j]==0)
          {
            cout << j-1 << " ";
            clus[radius][clus[radius][0][0]+1][k] = j-1;
            k++;
          }
      }
  }
cout << endl;
clus[radius][clus[radius][0][0]+1][1] = notclus;

}  // end the for loop searching for clusters
}

void Sortsmall(int iRange1, int iRange2, int& numseqs,
                int info[maxradius][sequences][sequences],
                float avdist[maxradius][sequences],
                int clus[maxradius][sequences][sequences],
		int origclus[sequences])
{

int x;
int y;
int j;
int k;
int findclus[sequences][sequences];
int tempvar[sequences];
float clusad[maxradius][sequences];
int maxcenter;
int maxvar;
int cluster;
int count;
int radius;
int notclus;
int overlap;

//  sorts the # of variants in clusters in descending order.  finds
//  the clusters for the radius values

for (radius = iRange1; radius<=iRange2; radius ++)
{
  for (int i=0; i<=(numseqs-1); i++)
    {
        tempvar[i] = info[radius][i][0];
    }

// tempvar[i] holds the # of variants w/in the cluster w/ center I

  zeromatrix (findclus);
  maxvar = 99;
  cluster = 0;
  count = 0;
  while (maxvar >0)
    {
        x=0;
        y=0;
        maxvar = 0;
        maxcenter = 0;
        for (int center = 0; center <= (numseqs-1); center++)
          {
            if (tempvar[center]!=0)
              {
               if (tempvar[center] > maxvar)
                {
                  maxvar = tempvar[center];
                  maxcenter = center;
                }
               if (tempvar[center] == maxvar)
                {
                 if (avdist[radius][center] < avdist[radius][maxcenter])
                   maxcenter = center;
                }
               }  // end the for loop
          }  // found the cluster w/ max. # of variants inside

        CompareClust(avdist, radius, count,findclus, maxvar,
                info, maxcenter,clusad, tempvar, numseqs);
//I took out two brackets that let it compile here
//Rob, 11/12/99
    }
  clus[radius][0][0] = (count-1);  // # of distinct clusters

  cout << endl;
  cout << "For radius = " << radius << ", there were ";
  cout << (count-1) << " distinct sub-groups present." << endl;

  for (int i=1; i<=(count-1); i++)
    {
      clus[radius][i][0] = findclus[i][0];  // center of cluster
      clus[radius][i][1] = findclus[i][1];  // # of nbrs in cluster
      for (int j=2; j<=findclus[i][1]+2; j++)
	{
          clus[radius][i][j] = info[radius][findclus[i][0]-1][j-1];
        }
    }

  for (int i=1; i<=(count-1); i++)
    {
      cout << "Sub-group #" <<  i << " has center ";
      cout << origclus[clus[radius][i][0]-1] << ", with ";
      cout << clus[radius][i][1];
      cout << " neighbors, including: ";
      for (int j=2; j<=findclus[i][1]+2; j++)
        {
	  cout << origclus[clus[radius][i][j]-1] << " ";
        }
      cout << endl;
    }

  notclus = 0;
  for (int i=1; i<= clus[radius][0][0]; i++)
    {
      for (int j=2; j<= numseqs+1; j++)
        {
//        cout << findclus[i][j] << " ";
          if (findclus[i][j]>0)
            {
              findclus[0][j]++;
              notclus++;  // indicates whether any var. not in clus.
            }
        }
//    cout << endl;
     }

  if (notclus >0)
    {
      notclus=0;  // # of var. not in any clusters
      overlap=0;  // indicates if any overlap
      for (j=2; j<= numseqs+1; j++)
        {
          if (findclus[0][j]==0)
            notclus++;
          if (findclus[0][j]>1)
            overlap++;
        }
      if (overlap>0)
        {
          cout << "Of these variants, the following belonged to more ";
          cout << "than one sub-group: " << endl;
          for (j=2; j<= numseqs+1; j++)
            {
              if (findclus[0][j]>1)
                cout << origclus[j-2] << " ";
            }
          cout << endl;
        }
      cout<< "There were " << notclus ;
      cout<< " variants (" << float(100*notclus/numseqs);
      cout<< "%) that were not contained in any sub-groups:"<<endl;
      for (j=2; j<= numseqs+1; j++)
        {
          if (findclus[0][j]==0)
            cout << origclus[j-2] << " ";
	}
      cout << endl;
    }

}  // end the for loop searching for clusters
}

void Resort(int clus[maxradius][sequences][sequences],
                char datamatrix[sequences][sites],
                int& numsites, int& numseqs, int& iRange1,
                char sorteddata[sequences][sites],
		int seqorder[sequences],
		int sortedmat[sequences][sequences])
{
int i;
int j;
int k;
int l;
int different;
int radius;
int row;
int column;
int numdiff;
int x;
char key;

radius = iRange1;
k=0;
cout << "\nThe new order of the sequences for the rearranged nucleotide";
cout << endl << "change matrix is:   ";
i=1;
  for (j=2; j<=clus[radius][i][1]+2; j++)
    {
      for (x=0; x<=numsites-1; x++)
        {
          sorteddata[k][x]=datamatrix[clus[radius][i][j]-1][x];
          seqorder[k]=clus[radius][i][j];
        }
      k++;
      cout << clus[radius][i][j] << " ";
    }
// cout << endl;

for (i=2; i<= clus[radius][0][0]+1; i++)
  {
    for (j=2; j<=clus[radius][i][1]+2; j++)
      {
        different = 0;
        for (l=0; l<=k; l++)
          {
            if (clus[radius][i][j] == seqorder[l])
              {
                l=k;
                different = 0;
              }
            else
              {
                different = 1;
              }
          }
        if (different == 1)
          {
            for (x=0; x<=numsites-1; x++)
              {
                sorteddata[k][x]=datamatrix[clus[radius][i][j]-1][x];
                seqorder[k]=clus[radius][i][j];
              }
            k++;
            cout << clus[radius][i][j] << " ";
          }
      }
//  cout << endl;
  }

for (row= 0; row <=(numseqs-2); row++)  //loop for all rows (-1)
{
 for (j=row+1; j<=(numseqs-1); j++)   //loop for 2nd comparison row
  {
   numdiff=0;
   for (column=0; column<=(numsites-1); column++)   // columns
   {if (sorteddata[row][column]!='*' && sorteddata[j][column]!='*')
    {
     if (sorteddata[row][column]!='-' && sorteddata[j][column]!='-')
      {
       if (sorteddata[row][column] != sorteddata[j][column])
        {  numdiff++;
        }
      }
    }
  }
 sortedmat[row][j]=numdiff;
 sortedmat[j][row]=numdiff;
 }
}

      cout << "\n\nPress any key, then 'return' to see the ";
      cout << "rearranged nt change matrix  ";
      cin >> key;
}

void SubCluster(int& numseqs, int& numsites, int& iSubCluster, int& j, 
		int& row, int& column, int& numdiff, int& x, int& iRange1,
                int& newseqs, char smalldata[sequences][sites],
                char datamatrix[sequences][sites],
                int origclus[sequences],
                int clus[maxradius][sequences][sequences],
                int smallchange[sequences][sequences])
{
int k;
int radius;

radius = iRange1;
k=0;
for (j=2; j<=clus[radius][iSubCluster][1]+2; j++)
  {
    for (x=0; x<=numsites-1; x++)
      {
        smalldata[k][x]=datamatrix[clus[radius][iSubCluster][j]-1][x];
      }
    k++;
  }

newseqs=clus[radius][iSubCluster][1]+1;
for (j=0; j<=clus[radius][iSubCluster][1]+1; j++)
  {
    origclus[j]=clus[radius][iSubCluster][j+2];
  }
/**
cout << "iSubCluster = " << iSubCluster << endl;
cout << "radius = " << radius << endl;
cout << "# seqs = " << newseqs ;
cout << endl;
**/
}



void CompareClust(float avdist[maxradius][sequences], int& radius, 
		int& count, int findclus[sequences][sequences], 
		int& maxvar, int info[maxradius][sequences][sequences],
		int& maxcenter, float clusad[maxradius][sequences],
		int tempvar[maxradius], int& numseqs)
{
int i;
int x;
int y;
int k;

clusad[radius][count+1] = avdist[radius][maxcenter];
  if (maxvar > 0)
    {
        findclus[count+1][0] = maxcenter+1;  // center of cluster
        findclus[count+1][1] = maxvar;       // # of nbrs in cluster

// this loop places a 1 in the component correlating to the variant
// which is present in the cluster
        for (i=1; i<=(info[radius][maxcenter][0]+1); i++)
          {
            findclus[count+1][info[radius][maxcenter][i]+1]=1;
          }

// for the 1st cluster, just copy info, no need to compare clusters,
// but for 2nd cluster, need to compare variants w/in clusters

        k=0;
if (count > 0)
          {
            while (k<count)
              {
                // cout << "k1= " << k << " ";
                x=0;
                y=0;
                for (i=2; i<=numseqs+1; i++)
                  {
                    if (findclus[k+1][i]>findclus[count+1][i])
                        x++;
                    if (findclus[k+1][i]<findclus[count+1][i])
                        y++;
                  }
                if (y==0 && x>0)  // prev clus contains the present clus
                  {
                    k=count;
                    for (i=0; i<=numseqs+1; i++)
                      {
                        findclus[count+1][i]=0;
		      }
                    count--;
                  }
                if (x==0 && y==0)  // clusters are same, but we already
                                   // know prev. clus smaller av. dist
                  {
                    k=count;
                    //  cout << "SAME ";
                    for (i=0; i<=numseqs+1; i++)
                      {
                        findclus[count+1][i]=0;
                      }
                    count--;
                  }
                //  cout << "k2= " << k << endl;
                k++;
          }
    }

  }
tempvar[maxcenter]=0; // removes the max for the next loop
count++;
   // end of loop to find clusters for a radius

}

void zerofreq(/* out */ int freqcount[])
// zeroes out positions min_char through max_char
{
        char index;     // loop control and index variable
        for (index = min_char; index <= max_char; index++)
                freqcount[index] = 0;
}

void Print(/* in */ const int freqcount[])  // list of char counts
{
        char index;     // loop control and index variable
        int  tot = 0;	// total number of changes
	int  aa;	// change of a.a. char to integer

        for (index = min_char + 1; index <= max_char; index++)
                tot = tot + freqcount[index];

        for (index = min_char; index <= max_char; index++)
//              if (index == min_char)
//                      cout << "There were " << freqcount[index]
//                           << " gaps" << endl;
//              else if (freqcount[index] > 0)
		if (freqcount[index] > 0)    
                        cout << index << " occurred "
                        << setw(3) << freqcount[index]
                        << " time(s)"
                        << " at a frequency of "
                        << setw(5)
			<< float(freqcount[index])/float(tot)*100.00
                        << "%" << endl;
}

void zeromatrix(/* out */ intmatrix matrix )
{
int rows;
int columns;

for (rows=0; rows<70; rows++)
	for (columns=0; columns<=70; columns++)
		matrix [rows] [columns] =0;
}

void zeroinfo(/* out */ infomatrix matrix )
{
int rows;
int columns;
int j;

for (rows=0; rows<200; rows++)
        for (columns=0; columns<100; columns++)
		for (j=0; j<100; j++)
                  matrix [rows] [columns] [j] = 0;
}

void zeroavdist(/* out */ avdistmatrix matrix )
{
int rows;
int columns;

for (rows=0; rows<200; rows++)
        for (columns=0; columns<=100; columns++)
                matrix [rows] [columns] =0;
}


