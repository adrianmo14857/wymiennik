#include "exch.h"

unsigned  seed = chrono::system_clock::now().time_since_epoch().count();
mt19937 gen(seed);
uniform_real_distribution<double> dist(0.0, 1.0);
double myRand() { return dist(gen); }


vec Markov(int num, const matrix& H, const vec& b)
{
	const double SMALL = 1.0e-4;     // Termination criterion for random walk

	int N = H.size();
	vec result(N);

	// Compute row sums, transition matrix and cumulative probabilities
	vec rowSum(N);
	matrix p(N, vec(N));
	matrix sump(N, vec(N, 0));

	for (int i = 0; i < N; i++)
	{
		rowSum[i] = 0;
		for (int j = 0; j < N; j++) rowSum[i] += abs(H[i][j]);
		for (int j = 0; j < N; j++)
		{
			p[i][j] = abs(H[i][j]) / rowSum[i];
			if (j == 0) sump[i][j] = p[i][j];
			else          sump[i][j] = sump[i][j - 1] + p[i][j];
		}
	}

	// Now solve by random walk for each coordinate
	for (int i = 0; i < N; i++)
	{
		double sumX = 0.0;
		for (int trial = 1; trial <= num; trial++)
		{
			int index = i;
			double W = 1, X = b[i];                 // W = multiplied weights; X = displacement on random walk
			while (abs(W) > SMALL)
			{
				double r = myRand();                 // Random number on [0,1)
				int j = 0;
				while (r > sump[index][j]) j++;    // First column where cumulative probability exceeds r

				W *= rowSum[index];
				if (H[index][j] < 0) W = -W;
				X += W * b[j];
				index = j;
			}
			sumX += X;
		}
		result[i] = sumX / num;
	}

	return result;
}


vec calculateMarkov(int num,matrix A,vec b) {

	int N = A.size();
	double Amax = 0.0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) if (abs(A[i][j]) > Amax) Amax = abs(A[i][j]);
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++) A[i][j] /= Amax;
		b[i] /= Amax;
	}

	matrix H = A;
	for (int i = 0; i < N; i++) H[i][i]++;

	vec X = Markov(num, H, b);

	return X;
}