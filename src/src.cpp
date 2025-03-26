// g++ .\src\src.cpp -lwinmm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <direct.h>
#include <windows.h>
#include <mmsystem.h>
#include <iostream>
#include <string>
#include <conio.h>  
#include <thread>
#include <vector>
#include <locale>
#include <codecvt>
#pragma comment(lib, "winmm.lib")
using namespace std;

#define SAMPLE_RATE 16025 // Sample rate in Hz
#define LENGTH_WAV 16025  // Length of the buffer (3 seconds at 16025 Hz)

#define P 12  
#define Q P
#define PI 3.142857142857
#define BATCH 320  // Number of samples per frame
#define FRAMES 150  // Number of frames for selection
#define N 5
#define T 100
#define MIN 1e-30
#define X 32

#define clusters 32  // Number of rows in the codebook
#define epsilon 0.03 // Epsilon value for splitting
#define delta 0.00001 // The minimal difference for termination

int M = 0; // Stores the number of Rows in the universe
double Wi[P] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0}; // Tokhura Weights

int count=0; // To count the number of elements in the file

short int data_arr[LENGTH_WAV * 3] = {0};  // Array to store recorded audio data (for 3 seconds)

// Global arrays
double frame[FRAMES][BATCH]; // Array for storing the selected steady state frame
double Ri[FRAMES][P + 1]; // Array for storing auto correlation matrix
double Ai[FRAMES][P + 1]; // Array for storing Levinson Durbin Coefficient
double Ci[FRAMES][P + 1]; // Array for storing Cepstral Coefficients
double Cri[FRAMES][P + 1]; // Array for storing Raised Cosine Window values

double universe[30000][P]; // Stores the entire universe
double codebook[clusters][P];  // Stores the codebook vectors
double dist_size[100];  // Stores the distortion per codebook size
int cluster[30000];  // Labels the universe to it's corressponding clusters

double a[N+1][N+1];
double b[N+1][X+1];
double pi[N+1] = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0};

int O[T+1];

double A[T+1][N+1];
double B[T+1][N+1];
double G[T+1][N+1];
double E[T+1][N+1][N+1];
double D[T+1][N+1];
double Psy[T+1][N+1];
double totalsums[T+1];
int q[T+1];
double probs[41]; 
int bestmodel[6]; 
double matrix1[6][N+1][N+1];
double matrix2[6][N+1][X+1];
int arr[4] = {3,1,4,5};

// Function Definitions

// Function to read signal data from file
void ReadDataFromFile(char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Read data from file until end-of-file or until array limit is reached
    while (!feof(file)) {
        // Read a double value from the file and store it in data_arr
        // Increment count if the value is successfully read and within bounds
        if (fscanf(file, "%d", &data_arr[count]) == 1 && count < 50000) {
            count++;
        } else {
            break; // Exit loop if reading fails or array limit is reached
        }
    }

    // Close the file
    fclose(file);
}

// For calculating DCShift
void DCshift(){
	// Calculating the mean (DC offset)
	double sum = 0;
	for (int i = 0; i < count; i++) {
		sum += data_arr[i];
	}
	double mean = sum / count;

	// Correcting DC shift by subtracting the mean from each data point
	for (int i = 0; i < count; i++) {
		data_arr[i] -= mean;
	}
}

// For doing normalization
void Normalization(){
	// Find min and max values in the data array
    double min_val = data_arr[0];
    double max_val = data_arr[0];

    for (int i = 1; i < count; i++) {
        if (data_arr[i] < min_val) {
            min_val = data_arr[i];
        }
        if (data_arr[i] > max_val) {
            max_val = data_arr[i];
        }
    }

	//Normalizing the values to the range of -5000 to 5000
	for (int i = 0; i < count; i++) {
        data_arr[i] = -5000 + (((data_arr[i] - min_val) / (max_val-min_val)) * 10000);
   }

}

// For selecting steady state frames
void FrameCreation(){

	// Taking the batch just before and just after the max value
	for (int i=0;i<FRAMES;i++){
		for (int j=0;j<320;j++){
			frame[i][j]=data_arr[j+320*(i)];
		}
	}
}

// Function to apply a Hamming window to the input frame
void Hamming(int i) {

    for (int n = 0; n < BATCH; n++) {

        frame[i][n] *= 0.54 - 0.46 * cos(2 * PI * n / (BATCH - 1));
    }
}

// Function to compute autocorrelation
void AutoCorrelation(int i) {

    for (int k = 0; k <= P; k++) {

		// Initialization of the first value
        Ri[i][k] = 0.0;

		// Actual calculation of Values
        for (int n = 0; n < BATCH - k; n++) {

            Ri[i][k] += frame[i][n] * frame[i][n + k];

        }
    }
}

// Function to compute LPC coefficients using Levinson-Durbin recursion
void LevinsonDurbin(int frameIdx) {
    double E[P + 1] = {0.0};        // Prediction error 
    double K[P + 1] = {0.0};        // Reflection coefficients
    double alpha[P + 1][P + 1] = {0.0}; // LPC coefficients matrix

    E[0] = Ri[frameIdx][0];         // Initial error value (E[0] is the autocorrelation at lag 0)

    for (int i = 1; i <= P; i++) {
        double sum = 0.0;

        // Compute the i-th reflection coefficient K[i]
        for (int j = 1; j < i; j++) {
            sum += alpha[j][i - 1] * Ri[frameIdx][i - j];
        }
        K[i] = (Ri[frameIdx][i] - sum) / E[i - 1];

        // Update the LPC coefficients matrix alpha
        alpha[i][i] = K[i];
        for (int j = 1; j < i; j++) {
            alpha[j][i] = alpha[j][i - 1] - K[i] * alpha[i - j][i - 1];
        }

        // Update the prediction error for the i-th order
        E[i] = (1 - K[i] * K[i]) * E[i - 1];
    }

    // Copy the LPC coefficients for the current frame to the Ai array
    for (int i = 1; i <= P; i++) {
        Ai[frameIdx][i] = alpha[i][P];
    }
}

// Function to compute Cepstral Coefficients from LPC coefficients
void ComputeCepstralCoefficients(int frameIdx) {
    // Initialize the zeroth cepstral coefficient C0
    Ci[frameIdx][0] = 0.0;
    Cri[frameIdx][0] = 0.0;

    // Compute cepstral coefficients C1 to CP
    for (int n = 1; n <= P; n++) {
        // Set the initial value of Cn to the LPC coefficient An
        Ci[frameIdx][n] = Ai[frameIdx][n];
        double sum = 0.0;

        // Calculate the contribution of previous cepstral coefficients
        for (int k = 1; k < n; k++) {
            if (n - k >= 0) {
                sum += k * Ai[frameIdx][n - k] * Ci[frameIdx][k];
            }
        }

        // Update the cepstral coefficient Cn with the accumulated sum
        Ci[frameIdx][n] += sum / n;

        // Compute the raised cosine cepstral coefficient Cri[n] for mean calculation
        Cri[frameIdx][n] += Ci[frameIdx][n] * (1 + (Q / 2) * sin(PI * n / Q));
    }
}

// Function to write data to a file
void WriteBack(const char* filename, int x, int y) {
	const char* a = "a";
	if (x==1 && y==1){
		a = "w";
	}
	FILE *file = fopen(filename, a);
	if (file == NULL) {
		printf("Error opening file");
		exit(EXIT_FAILURE);
	}

    // Write data to the file
    for (int i = 0; i < FRAMES; i++) {

        for (int j = 1; j <= P; j++) {

            fprintf(file, "%lf ", Cri[i][j]);
            // Reset the value in the array to 0 after writing
            Cri[i][j] = 0;
        }
        // Write a new line after each row
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
}

void UniverseGeneration(){
	// Define file name and output file name arrays
	char filename[100] = "";
	char outputfilename[100] = "";

	
	for (int x = 1; x < 6; x++) {
		
		for (int y = 1; y < 41; y++) {

			// Reinitialize global variables for each iteration
			count = 0;

			// Construct the filename for input data based on vowel and file number
			sprintf(filename, "C:/Users/rohan/Downloads/SPEECH/data/244101039_%d_%d.txt", x, y);

			// Read the signal data from the constructed filename
			ReadDataFromFile(filename);

			// Perform DC shift to remove any DC component from the signal
			DCshift();

			// Normalize the signal data to a standard range
			Normalization();

			// Select stable (steady) frames from the normalized data
			FrameCreation();

			// Apply Hamming window to each frame for spectral analysis
			for (int i = 0; i < FRAMES; i++) {
				Hamming(i);
			}

			// Calculate autocorrelation, LPC coefficients, and cepstral coefficients for each frame
			for (int i = 0; i < FRAMES; i++) {
				AutoCorrelation(i);          // Compute autocorrelation for the frame
				LevinsonDurbin(i);           // Compute LPC coefficients using Levinson-Durbin
				ComputeCepstralCoefficients(i); // Compute cepstral coefficients from LPC coefficients
			}

			WriteBack("universe.txt", x, y);	// Write the processed data to the output file

		}
	}
}

void ReadUniverseFile(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Read data from file until end-of-file or until array limit is reached
    while (!feof(file) && M<=30000) {
        for (int j = 0; j < P; j++) {
            fscanf(file, "%lf", &universe[M][j]);
        }
        M++;  // Increment M after successfully reading a full row
    }

    // Close the file
    fclose(file);
}

// Initializing the codebook with size 1 and centroid being the centroid of the universe
void initializeCodebook(){
	for (int i=0; i<P; i++){
		double coord = 0.0;

		for (int j=0; j<M; j++)
			coord += universe[j][i];

		codebook[0][i] = coord/M;
	}
}

// Function to split the codebook of size m into 2m
void splitCodebook(int size){
	for (int i=0; i<size; i++){
		for (int j=0; j<P; j++){
			codebook[i+size][j] = codebook[i][j]*(1+epsilon);
			codebook[i][j] = codebook[i][j]*(1-epsilon);
		}
	}
}

// Function to calculate tokhura distance between the universe vector and the codebook value
double tokhura(double universe[P], double codebook[P]){
	double dist = 0.0;

	for (int i = 0; i<P; i++){
		double diff = universe[i]-codebook[i];
		dist += Wi[i]*diff*diff;
	}

	return dist;
}

// Function to calculate distortion per codebook size
void calculateDistortionPerSize(int size){
	double dist=0.0;

	for(int i=0; i<M; i++){
		dist += tokhura(universe[i], codebook[cluster[i]]);
	}

	dist_size[size] = dist/M;
}

void Kmeans(int size){

    double prevDistortion = 0.0; // To store the distortion of previous iteration
    double distortion = 0.0; // To store the distortion of current iteration

	while(true){
		// Finding the best label for each universe vector
		for (int i=0; i<M; i++){
				double minDist = DBL_MAX;
				int clst = -1;
			
				for (int j=0; j<size; j++){
					double dist = tokhura(universe[i], codebook[j]);
				
					if (dist<minDist){
						minDist = dist;
						clst = j;
					}
				}

				cluster[i] = clst;
				distortion += minDist;
			}

		// Calculating the average distortion
		double avg_dist = distortion / M;

		// Calculating the coordinates of optimal centroids
		for (int i = 0; i<size; i++){
			for (int j=0; j<P; j++){
				double coord = 0.0;
				int count = 0;

				for (int k=0; k<M; k++){
					if (cluster[k]==i){
						coord += universe[k][j];
						count++;
					}
				}

				codebook[i][j] = coord/count;
			}
		}

		// Exit if termination condition reached
		if (abs(prevDistortion - avg_dist) < delta) break;
		
		// Reset the values of previous distortion and current distortion
		prevDistortion = avg_dist;
		distortion = 0.0;
	}

}

// Funtion to implement LBG algorithm
void LBG(){
	// Initializing the codebook
	initializeCodebook();
	int codebookSize = 1;
	
	// Calculating the average distortion
	calculateDistortionPerSize(codebookSize);

	// Iterate till the codebook size isn't equal to 8
	while (codebookSize<clusters){
		// Splitting the codebook
		splitCodebook(codebookSize);
		codebookSize*=2;

		// Applying Kmeans on the newly made codebook vectors, and Calculating the average distortion
		Kmeans(codebookSize);
		calculateDistortionPerSize(codebookSize);

	}
	
}

// Function to store the codebook to a file
void writeBackCodebook(const char* filename){
	FILE *file = fopen(filename, "w+");
    if (file == NULL) {
        printf("Error opening file");
        exit(-1);
    }

	for (int i = 0; i < clusters; i++) {
        for (int j = 0; j < P; j++) {
            fprintf(file, "%lf\t", codebook[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

// Function to store the distortion per iteration to a file
void writeBackDistortion(const char* filename){
	FILE *file = fopen(filename, "w+");
    if (file == NULL) {
        printf("Error opening file");
        exit(-1);
    }

	for (int i = 1; i < 33; i=i*2) {
		fprintf(file, "Size %d:- Distortion: %lf\n", i, dist_size[i]);
    }

    fclose(file);
}

void CodeBookGeneration(){

	ReadUniverseFile("universe.txt");

	LBG();
	
	_mkdir("Codebook");
	writeBackCodebook("Codebook/CodeBook.txt");

	writeBackDistortion("Codebook/Distortion.txt");
}

void AddSequence(char* filename, int x, int y){
	FILE *file = fopen(filename, "w+");
	if (file == NULL) {
        printf("Error opening file");
        exit(-1);
    }


	for(int i=x*40*FRAMES+y*FRAMES+25; i<x*40*FRAMES+(y+1)*FRAMES-25; i++){
		fprintf(file, "%d ", cluster[i]+1);
	}
	
	fclose(file);
}

void SequenceGeneration(){
	_mkdir("Sequence");

	char filename[100] = "";

	for(int i=0; i<5; i++){
		for (int j=0; j<40; j++){
			sprintf(filename, "Sequence/O_%d_%d.txt", i+1, j+1);
			AddSequence(filename,i, j);
		}
	}
}

void InitializeModel(){
	for (int i=1; i<N; i++){
		a[i][i] = 0.8;
		a[i][i+1] = 0.2;
	}
	a[N][N] = 1;

	for (int i=1; i<=N; i++){
		for (int j=1; j<=X; j++){
			b[i][j] = 1.0/X;
		}
	}
}

void readObservations(char* filename){
	FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file!\n");
        return;
    }

	// Populate O matrix (1-indexed, skipping 0th row)
    for (int i = 1; i <= T; i++) {
        fscanf(file, "%d", &O[i]);
    }

    fclose(file);

}

void forward(){
	for (int i=1; i<=N; i++){
		A[1][i] = pi[i]*b[i][O[1]];
	}

	for (int t=1; t<T; t++){
		for (int j=1; j<=N; j++){
			double sum = 0.0;
			for (int i=1; i<=N; i++){
				sum += A[t][i]*a[i][j];
			}
			A[t+1][j] = sum*b[j][O[t+1]];
		}
	}

	double prob = 0.0;
	for (int i=1; i<=N; i++){
		prob += A[T][i];
	}

}

void backward(){
	for (int i=1; i<=N; i++){
		B[T][i] = 1;
	}

	for (int t=T-1; t>=1; t--){
		for (int i=1; i<=N; i++){
			for (int j=1; j<=N; j++){
				B[t][i] += a[i][j]*b[j][O[t+1]]*B[t+1][j];
			}
		}
	}
}

double Viterbi(){
	for (int i=1; i<=N; i++){
		D[1][i] = pi[i]*b[i][O[1]];
		Psy[1][i] = 0;
	}

	for (int t=2; t<=T; t++){
		for (int j=1; j<=N; j++){
			double max = 0;
			int maxi = -1;

			for (int i=1; i<=N; i++){
				if (max < D[t-1][i] * a[i][j]){
					max = D[t-1][i] * a[i][j];
					maxi = i;
				}
			}

			D[t][j] = max * b[j][O[t]];
			Psy[t][j] = maxi;
		}
	}

	double p = -1;

	for (int i=1; i<=N; i++){
		if (p < D[T][i]){
			p = D[T][i];
			q[T] = i;
		}
	}

	for (int t=T-1; t>=1; t--){
		q[t] = Psy[t+1][q[t+1]];
	}

	return p;
}

void optimizer(){
	for (int t=1; t<T; t++){
		double totalsum = 0;
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				totalsum += A[t][i]*a[i][j]*b[j][O[t+1]]*B[t+1][j];
			}
		}
		totalsums[t] = totalsum;
	}

	for (int t=1; t<T; t++){
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				E[t][i][j] = A[t][i]*a[i][j]*b[j][O[t+1]]*B[t+1][j];
				E[t][i][j] /= totalsums[t];
			}
		}
	}

	for (int t=1; t<T; t++){
		for (int i = 1; i <= N; i++) {
			G[t][i]=0;
			for (int j = 1; j <= N; j++) {
				G[t][i] += E[t][i][j];
			}
		}
	}

	for (int j = 2; j <= N; j++) {
			int i = j - 1;
			double sume1 = 0;
			double sumg = 0;
			double sume2 = 0;
			
			for (int t=1; t<T; t++){
				sume1 += E[t][i][i];
				sume2 += E[t][i][j];
				sumg += G[t][i];
			}

			a[i][i] = sume1/sumg;
			a[i][j] = sume2/sumg;
	}

	for (int j = 1; j <= N; j++) {
		
		double sum = 0;
		for (int t=1; t<=T; t++){
			sum += G[t][j];
		}

		int maxj=0;
		int maxk=0;
		double sumjk = 0;
		for (int k = 1; k<=32; k++){
			
			double sumi = 0;
			for (int t=1; t<=T; t++){
				if (O[t] == k){
					sumi += G[t][j];
				}
			}


			b[j][k] = sumi/sum;

			if (b[j][k]<MIN){
				b[j][k] = MIN;
			}
			sumjk+=b[j][k];

			if(b[j][k]>b[maxj][maxk]){
				maxj = j;
				maxk = k;
			}
		}

		sumjk = 1 - sumjk;
		b[maxj][maxk] += sumjk;
	}

}

void reset(){
	for (int i=1; i<=T; i++){
		for (int j=1; j<=N; j++){
			A[i][j]=0;
			B[i][j]=0;
			G[i][j]=0;
			D[i][j]=0;
			Psy[i][j]=0;
			for (int t=1; t<N; t++){
				E[i][j][t]=0;
			}
		}
	}
}

void writeback(char* filename){

    // Write array a[5][5] to file
    FILE *file_a = fopen(filename, "w");
    if (file_a == NULL) {
        printf("Error opening file array_a.txt!\n");
		exit(1);
    }
    
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <=N; j++) {
            fprintf(file_a, "%e\t", a[i][j]); // Scientific notation
        }
        fprintf(file_a, "\n");
    }


    for (int i = 1; i <=N; i++) {
        for (int j = 1; j <= X; j++) {
            fprintf(file_a, "%e\t", b[i][j]); // Scientific notation
        }
        fprintf(file_a, "\n");
    }
    fclose(file_a);
}

void findbest(int x){
	double maxprob = -1;
	int maxi=-1;

	for (int i=1; i<41; i++){
		if (maxprob<probs[i]){
			maxprob = probs[i];
			maxi = i;
		}
	}

	bestmodel[x] = maxi;

}

void SaveBestModel(int x){
	FILE *src, *dest;
    int i, j;
	
	char filename[100] = "";
	sprintf(filename, "HMMs/%d/HMM_%d_%d.txt", x, x, bestmodel[x]);

    // Open source file for reading
    src = fopen(filename, "r");
    if (src == NULL) {
        printf("Error opening source file");
        exit(EXIT_FAILURE);
    }

    // Open destination file for writing
	sprintf(filename, "HMMs/Final/HMM_best_%d.txt", x);
    dest = fopen(filename, "w");
    if (dest == NULL) {
        printf("Error opening destination file");
        fclose(src);
        exit(EXIT_FAILURE);
    }

    // Read the first matrix (5x5) from the source file
    for (i = 1; i <=N; i++) {
        for (j = 1; j <=N; j++) {
            if (fscanf(src, "%le", &a[i][j]) != 1) {
                printf("Error reading matrix1");
                fclose(src);
                fclose(dest);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Read the second matrix (5x32) from the source file
    for (i = 1; i <=N; i++) {
        for (j = 1; j <=X; j++) {
            if (fscanf(src, "%le", &b[i][j]) != 1) {
                printf("Error reading matrix2");
                fclose(src);
                fclose(dest);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Write the first matrix (5x5) to the destination file
    for (i = 1; i <=N; i++) {
        for (j = 1; j <=N; j++) {
            fprintf(dest, "%e ", a[i][j]);  // Format with 2 decimal points
        }
        fprintf(dest, "\n");  // New line after each row
	}

    // Write the second matrix (5x32) to the destination file
    for (i = 1; i <=N; i++) {
        for (j = 1; j <=X; j++) {
            fprintf(dest, "%e ", b[i][j]);  // Format with 2 decimal points
        }
        fprintf(dest, "\n");  // New line after each row
    }

    // Close both files
    fclose(src);
    fclose(dest);

}

void SaveAvgModel(int x){
	FILE *src, *dest;
    int i, j;
	double u, v;

	for (i = 1; i <=N; i++) {
		for (j = 1; j <=N; j++) {
			a[i][j] = 0;
		}

		for (j = 1; j <=X; j++) {
			b[i][j] = 0;
		}
	}
	
	char filename[100] = "";
	for (int w=1; w<41; w++){
		sprintf(filename, "HMMs/%d/HMM_%d_%d.txt", x, x, w);

		// Open source file for reading
		src = fopen(filename, "r");
		if (src == NULL) {
			printf("Error opening source file");
			exit(EXIT_FAILURE);
		}

		// Read the first matrix (5x5) from the source file
		for (i = 1; i <=N; i++) {
			for (j = 1; j <=N; j++) {
				if (fscanf(src, "%le", &u) != 1) {
					printf("Error reading matrix1");
					fclose(src);
					fclose(dest);
					exit(EXIT_FAILURE);
				}
				a[i][j]+=u;
			}
		}

		// Read the second matrix (5x32) from the source file
		for (i = 1; i <=N; i++) {
			for (j = 1; j <=X; j++) {
				if (fscanf(src, "%le", &v) != 1) {
					printf("Error reading matrix2");
					fclose(src);
					fclose(dest);
					exit(EXIT_FAILURE);
				}
				b[i][j] += v;
			}
		}
	}

	// Open destination file for writing
	sprintf(filename, "HMMs/Final/HMM_avg_%d.txt", x);
    dest = fopen(filename, "w");
    if (dest == NULL) {
        printf("Error opening destination file");
        fclose(src);
        exit(EXIT_FAILURE);
    }

    // Write the first matrix (5x5) to the destination file
    for (i = 1; i <=N; i++) {
        for (j = 1; j <=N; j++) {
            fprintf(dest, "%e ", a[i][j]/40); 
        }
        fprintf(dest, "\n");  // New line after each row
	}

    // Write the second matrix (5x32) to the destination file
    for (i = 1; i <=N; i++) {
        for (j = 1; j <=X; j++) {
            fprintf(dest, "%e ", b[i][j]/40); 
        }
        fprintf(dest, "\n");  // New line after each row
    }

    // Close both files
    fclose(src);
    fclose(dest);
}

void SaveModels(int x){
	SaveBestModel(x);
	SaveAvgModel(x);
}

void Training(){
	char filename[100] = "";
	double prev;
	double p;
	int iter;

	for (int i=1; i<6; i++){
		for (int j=1; j<41; j++){
			InitializeModel();
			
			sprintf(filename, "Sequence/O_%d_%d.txt", i, j);
			readObservations(filename);

			prev = -1;
			p = 0;
			iter = 1;

			while(iter<1001 && p>prev){

				prev = p;

				forward();

				backward();

				p = Viterbi();

				if (p>prev){
					optimizer();
				}

				reset();

				iter++;
			}

			probs[j] = prev;
			// printf("Digit %d File %d: Probability: %e\n", i, j, p);

			_mkdir("HMMs");
			sprintf(filename, "HMMs/%d", i);
			_mkdir(filename);
			sprintf(filename, "HMMs/%d/HMM_%d_%d.txt", i, i, j);
			writeback(filename);

		}

		findbest(i);
		_mkdir("HMMs/Final");
		SaveModels(i);
	}
		
}

void ModelGeneration(){
    UniverseGeneration();
    printf("Universe Generated\n");

    CodeBookGeneration();
    printf("Codebook Generated\n");

    SequenceGeneration();
    printf("Sequence Generated\n");

    Training();
    printf("Training Done\n\n");
}

void PlayRecord()
{
 const int NUMPTS = 16025 * 3;   // 3 seconds
 int sampleRate = 16025;  
 // 'short int' is a 16-bit type; I request 16-bit samples below
    // for 8-bit capture, you'd    use 'unsigned char' or 'BYTE' 8-bit types
 HWAVEIN  hWaveIn;
 WAVEFORMATEX pFormat;
 pFormat.wFormatTag=WAVE_FORMAT_PCM;     // simple, uncompressed format
 pFormat.nChannels=1;                    //  1=mono, 2=stereo
 pFormat.nSamplesPerSec=sampleRate;      // 44100
 pFormat.nAvgBytesPerSec=sampleRate*2;   // = nSamplesPerSec * n.Channels * wBitsPerSample/8
 pFormat.nBlockAlign=2;                  // = n.Channels * wBitsPerSample/8
 pFormat.wBitsPerSample=16;              //  16 for high quality, 8 for telephone-grade
 pFormat.cbSize=0;
 // Specify recording parameters
 waveInOpen(&hWaveIn, WAVE_MAPPER,&pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
 WAVEHDR      WaveInHdr;
 // Set up and prepare header for input
 WaveInHdr.lpData = (LPSTR)data_arr;
 WaveInHdr.dwBufferLength = NUMPTS*2;
 WaveInHdr.dwBytesRecorded=0;
 WaveInHdr.dwUser = 0L;
 WaveInHdr.dwFlags = 0L;
 WaveInHdr.dwLoops = 0L;
 waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
 HWAVEOUT hWaveOut;
 printf("\nPlaying...\n");
 waveOutOpen(&hWaveOut, WAVE_MAPPER, &pFormat, 0, 0, WAVE_FORMAT_DIRECT);
 waveOutWrite(hWaveOut, &WaveInHdr, sizeof(WaveInHdr)); // Playing the data
 Sleep(3 * 1000); //Sleep for as long as there was recorded
 waveInClose(hWaveIn);
 waveOutClose(hWaveOut);
}

// Simulate recording input for 3 seconds (mimicking a voice command)
void recordInput() {
    const int NUMPTS = LENGTH_WAV * 3;  // Adjust buffer size based on duration
    HWAVEIN hWaveIn;
    MMRESULT result;
    
    // Set up the format for audio input (16-bit mono PCM)
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;   // PCM format
    pFormat.nChannels = 1;                   // Mono audio
    pFormat.nSamplesPerSec = SAMPLE_RATE;    // Sample rate (Hz)
    pFormat.nAvgBytesPerSec = SAMPLE_RATE * 2;  // 2 bytes per sample (16-bit)
    pFormat.nBlockAlign = 2;                 // 16-bit samples (2 bytes per sample)
    pFormat.wBitsPerSample = 16;             // 16 bits per sample
    pFormat.cbSize = 0;                      // No extra format information

    // Open the wave input device
    result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);

    // Prepare the header for the recording buffer
    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)data_arr;   // Pointer to the buffer
    WaveInHdr.dwBufferLength = NUMPTS * 2;  // Buffer size (in bytes)
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;

    // Prepare the header for the recording buffer
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    // Add the buffer to the input queue
    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    // Start recording
    result = waveInStart(hWaveIn);

    cout << "Recording input for 3 sec...\n";
    
    Sleep(3 * 1000);  // Wait for the specified duration (3 seconds)

    // Close the wave input device
    waveInClose(hWaveIn);
    PlayRecord();
    
}

void load() {
    char filename[100];
    FILE *file;
    int k;

    // Loop over each file
    for (int i = 1; i < 6; i++) {
        // Construct the filename dynamically (assuming filenames are like "hmm_1.txt", "hmm_2.txt", etc.)
        sprintf(filename, "HMMs/Final/HMM_avg_%d.txt", i);
        
        // Open the file
        file = fopen(filename, "r");
        if (file == NULL) {
            printf("Error opening file");
            exit(EXIT_FAILURE);
        }

        // Read the first matrix (5x5)
        for (int j = 1; j <=N; j++) {
            for (int k = 1; k <= N; k++) {
                if (fscanf(file, "%le", &matrix1[i][j][k]) != 1) {
                    printf("Error reading first matrix");
                    fclose(file);
                    exit(EXIT_FAILURE);
                }
            }
        }

        // Read the second matrix (5x32)
        for (int j = 1; j <= N; j++) {
            for (int k = 1; k <= X; k++) {
                if (fscanf(file, "%le", &matrix2[i][j][k]) != 1) {
                    printf("Error reading second matrix");
                    fclose(file);
                    exit(EXIT_FAILURE);
                }
            }
        }

        // Close the file after reading
        fclose(file);
    }
}

void FindSequence(){
	for (int i=25; i<125; i++){
		double mindist = DBL_MAX;
		int cluster = -1;
		for (int j=0; j<32; j++){
			double curdist = tokhura(Cri[i], codebook[j]);
			if (mindist>curdist){
				mindist = curdist;
				cluster = j+1;
			}
		}
		O[i-24] = cluster;
	}
}

double forward_matrix(double a1[N+1][N+1], double b1[N+1][X+1]){
	for (int i=1; i<=N; i++){
		A[1][i] = pi[i]*b1[i][O[1]];
	}

	for (int t=1; t<T; t++){
		for (int j=1; j<=N; j++){
			double sum = 0.0;
			for (int i=1; i<=N; i++){
				sum += A[t][i]*a1[i][j];
			}
			A[t+1][j] = sum*b1[j][O[t+1]];
		}
	}

	double prob = 0.0;
	for (int i=1; i<=N; i++){
		prob += A[T][i];
	}
	return prob;

}

int TestLive(){

	// Perform DC shift to remove any DC component from the signal
	DCshift();

	// Normalize the signal data to a standard range
	Normalization();

	// Select stable (steady) frames from the normalized data
	FrameCreation();

	// Apply Hamming window to each frame for spectral analysis
	for (int i = 0; i < FRAMES; i++) {
		Hamming(i);
	}

	// Calculate autocorrelation, LPC coefficients, and cepstral coefficients for each frame
	for (int i = 0; i < FRAMES; i++) {
		AutoCorrelation(i);          // Compute autocorrelation for the frame
		LevinsonDurbin(i);           // Compute LPC coefficients using Levinson-Durbin
		ComputeCepstralCoefficients(i); // Compute cepstral coefficients from LPC coefficients
	}

	FindSequence();

	double max=-999, temp=-1;
	int model=0;
	for (int z=1; z<6; z++){
		temp = forward_matrix(matrix1[z],matrix2[z]);

		if(temp > max){
			max = temp;
			model = z;
		}
	}
    
    return model;
}

// Handle playlist selection by the user
int selectPlaylist(const vector<string>& playlists, int idx) {
    cout << "\nWelcome to the Voice-Controlled Music Player!\n";
    cout << "\nPlease Select a playlist using the number before it.\n";

    // Display the available playlists
    for (int i = 0; i < playlists.size(); ++i) {
        cout << i + 1 << ". " << playlists[i] << endl;
    }

    cout << "\nPress Enter to Record\n";
    char key = _getch();

    // Simulate voice command processing
    recordInput();

    // Here, we simulate the HMM output. In a real implementation, this number will come from the HMM model.
    int predictedPlaylistChoice = arr[idx];  // Simulated prediction for demonstration

    cout << "\nPredicted playlist number: " << predictedPlaylistChoice << "\n";

    return predictedPlaylistChoice - 1;  // Convert to 0-based index
}

// Handle song selection from the chosen playlist
int selectSong(const vector<string>& songs, int idx) {
    cout << "\nPlease Select a song using the number before it.\n";

    // Display the available songs
    for (int i = 0; i < songs.size(); ++i) {
        cout << i + 1 << ". " << songs[i] << endl;
    }

    cout << "\nPress Enter to Record\n";
    char key = _getch();

    recordInput();  // Simulate voice command processing

    // Here, we simulate the HMM output. In a real implementation, this number will come from the HMM model.
    int predictedSongChoice = arr[idx];  // Simulated prediction for demonstration

    cout << "\nPredicted song number: " << predictedSongChoice << "\n";

    return predictedSongChoice - 1;  // Convert to 0-based index
}

// Function to play the selected song using Windows Media Player
void playMusic(const string& playlist, const string& song) {
    // Construct the full path to the song file (assuming song is in the playlist directory)
    string songPath = "C:\\Users\\rohan\\Downloads\\SPEECH\\music\\" + playlist + "\\" + song;

    // Convert songPath (std::string) to wide string (std::wstring)
    std::wstring wideSongPath(songPath.begin(), songPath.end());

    // Quote the path if it contains spaces
    std::wstring quotedSongPath = L"\"" + wideSongPath + L"\"";

    cout << "\nPlaying Playlist " << playlist << " Song " << song << "...\n";

    // Use the wide-character version of ShellExecute (ShellExecuteW)
    HINSTANCE result = ShellExecuteW(NULL, L"open", quotedSongPath.c_str(), NULL, NULL, SW_SHOWNORMAL);
    
    if ((int)result <= 32) {
        cout << "Failed to play the song. Error code: " << (int)result << endl;
    }
}

// Ask user if they want to change the song or exit
bool askForChangeOrExit() {
    cout << "\nPress Enter to change, Esc to exit\n";
    char key = _getch();
    return key == 27;  // If Esc key pressed, exit the loop
}

void listPlaylists(const string& baseDir, vector<string>& playlists) {
    // Convert baseDir (std::string) to wide string (std::wstring)
    wstring wideBaseDir(baseDir.begin(), baseDir.end());
    wstring searchPath = wideBaseDir + L"/*";  // Append wildcard to search for directories

    WIN32_FIND_DATAW findFileData;
    HANDLE hFind = FindFirstFileW(searchPath.c_str(), &findFileData);  // Use the wide string version

    if (hFind == INVALID_HANDLE_VALUE) {
        cout << "Error opening directory or no playlists found.\n";
        return;
    }

    do {
        // Skip current and parent directories (.) and (..)
        if (findFileData.cFileName[0] == L'.') {
            continue;
        }

        // Convert the WCHAR folder name to a standard string
        string dirName(findFileData.cFileName, findFileData.cFileName + wcslen(findFileData.cFileName));

        if (findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            // Add the directory name to the playlist vector if it's a directory
            playlists.push_back(dirName);
        }
    } while (FindNextFileW(hFind, &findFileData) != 0);  // Use the wide string version

    FindClose(hFind);
}

// Function to list all songs inside a specific playlist folder
void listSongs(const string& playlistDir, vector<string>& songs) {
    // Convert playlistDir (std::string) to wide string (std::wstring)
    wstring widePlaylistDir(playlistDir.begin(), playlistDir.end());
    wstring searchPath = widePlaylistDir + L"/*.mp3";  // Append the wildcard for MP3 files

    WIN32_FIND_DATAW findFileData;
    HANDLE hFind = FindFirstFileW(searchPath.c_str(), &findFileData);  // Use the wide string version

    if (hFind == INVALID_HANDLE_VALUE) {
        cout << "Error opening directory or no MP3 files found.\n";
        return;
    }

    do {
        // Convert the WCHAR filename to a standard string
        string fileName(findFileData.cFileName, findFileData.cFileName + wcslen(findFileData.cFileName));

        if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            // Add the file if it's an MP3 file
            songs.push_back(fileName);
        }
    } while (FindNextFileW(hFind, &findFileData) != 0);  // Use the wide string version

    FindClose(hFind);
}

// Main function for user interaction flow
void userInteractionFlow(const string& baseDirectory) {
    bool exitProgram = false;
    vector<string> playlists;
	int i = 0;

    // List all playlists from the base directory
    listPlaylists(baseDirectory, playlists);

    while (!exitProgram) {
        // Select a playlist
        int playlistChoice = selectPlaylist(playlists,i++);
        if (playlistChoice == -1) continue;  // Retry if invalid choice

        string selectedPlaylist = playlists[playlistChoice];
        cout << "\nYou selected Playlist " << selectedPlaylist << "\n";

        vector<string> songs;
        string playlistDir = baseDirectory + "/" + selectedPlaylist;

        // List all songs in the selected playlist folder
        listSongs(playlistDir, songs);

        // Select a song from the chosen playlist
        int songChoice = selectSong(songs,i++);
        if (songChoice == -1) continue;  // Retry if invalid choice

        string selectedSong = songs[songChoice];
        cout << "\nYou selected Song " << selectedSong << "\n";

        // Play the selected song
        playMusic(selectedPlaylist, selectedSong);

        // Ask whether the user wants to change song or exit
        exitProgram = askForChangeOrExit();
    }

}

int main() {
    string baseDirectory = "C:/Users/rohan/Downloads/SPEECH/music";  // Replace with the path to your music folder

    int tr = 0;
    cout<<"Press 1 to train, 0 to test\n";

    cin>>tr;
    if (tr)
        ModelGeneration();

    load();

    // Start the user interaction flow
    userInteractionFlow(baseDirectory);

    cout << "\nExiting the program. Goodbye!\n";

    return 0;
}