#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <math.h>
//#include <fstream>
#include <thread>
#include <future>

using namespace std;

double MIN; // for the recursive function
vector<double> solution;
constexpr auto pi = 3.14;
constexpr auto ITERATIONS = 30;
typedef chrono::high_resolution_clock Clock;
//ofstream output("date.txt");

// 1. Rastrigin's function (f : [-5.12, 5.12])
double Rastrigin(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += x[i] * x[i] - 10 * cos(2 * pi * x[i]);
	}
	res += (double)10 * x.size();
	return res;
}

// 2. Michalewicz's function (f : [0, pi])
double Michalewicz(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += sin(x[i]) * pow(sin(((i + 1) * x[i] * x[i]) / pi), 20);
	}
	return -res;
}

// 3. Sphere function (f : [-5.12, 5.12])
double Sphere(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += x[i] * x[i];
	}
	return res;
}

// 4. Rosenbrock's function (f : [-5, 10])
double Rosenbrock(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size() - 1; ++i) {
		res += 100 * pow(x[i] * x[i] - x[i + 1], 2) + (x[i] - 1) * (x[i] - 1);
	}
	return res;
}




void recursiveCalculation(double f(vector<double>), int recLevel, vector<double> &x, double lowerBound, double upperBound, int dim) {
	if (recLevel != dim) {
		for (double i = lowerBound; i <= upperBound; i += 0.02) {
			x.push_back(i);
			recursiveCalculation(f, recLevel + 1, x, lowerBound, upperBound, dim);
			x.erase(x.begin() + recLevel);
		}
	}
	else {
		double res = f(x);
		if (res < MIN) {
			MIN = res;
			solution = x;
		}
	}
}

void deterministicApproach(double f(vector<double>), double lowerBound, double upperBound, int dim) {
	vector<double> x;
	MIN = DBL_MAX;

	auto start = Clock::now();
	recursiveCalculation(f, 0, x, lowerBound, upperBound, dim);
	auto finish = Clock::now();

	int time = chrono::duration_cast<chrono::milliseconds>(finish - start).count();
	
	cout << "Min: " << MIN << " At: x=(";
	for (int i = 0; i < solution.size(); ++i) {
		cout << solution[i] << ",";
	}
	cout << ") Time: " << time << "ms" << endl;
}




pair<double, vector<double>> lazyMinCalculation(double f(vector<double>), double lowerBound, double upperBound, int dim) {
	double min = DBL_MAX;
	vector<double> sol;
	for (double i = lowerBound; i <= upperBound; i += 0.002) {
		vector<double> x;
		for (int j = 0; j < dim; ++j) {
			x.push_back(i);
		}
		double res = f(x);
		if (res < min) {
			min = res;
			sol = x;
		}
	}
	return { min, sol };
}

void lazyDeterministicApproach(double f(vector<double>), double lowerBound, double upperBound, int dim) {
	auto start = Clock::now();
	pair<double, vector<double>> res = lazyMinCalculation(f, lowerBound, upperBound, dim);
	auto finish = Clock::now();

	int time = chrono::duration_cast<chrono::milliseconds>(finish - start).count();

	cout << "Min: " << res.first << " At: x=(";
	for (int i = 0; i < res.second.size(); ++i) {
		cout << res.second[i] << ",";
	}
	cout << ") Time: " << time << "ms" << endl;
}



double generateRandomDouble(double lowerBound, double upperBound) {
	uniform_real_distribution<double> unif(lowerBound, upperBound);
	random_device rd;
	return unif(rd);
}

pair<double, vector<double>> randomSearch(double f(vector<double>), double lowerBound, double upperBound, int dim) {
	double min = DBL_MAX;
	double res;
	vector<double> sol;
	for (int i = 0; i < 500000; ++i) {
		auto start = Clock::now();
		vector<double> x;
		for (int j = 0; j < dim; ++j) {
			x.push_back(generateRandomDouble(lowerBound, upperBound));
		}
		res = f(x);
		if (res < min) {
			min = res;
			sol = x;
		}
		auto finish = Clock::now();
	}
	return { min, sol };
}

void nedeterministicApproach(double f(vector<double>), double lowerBound, double upperBound, int dim) {
	double minSol = DBL_MAX;
	vector<double> minSolValues;
	double maxSol = -DBL_MAX;
	vector<double> maxSolValues;
	double avgSol = 0;

	int minTime = INT_MAX;
	int maxTime = INT_MIN;
	int avgTime = 0;

	for (int i = 0; i < ITERATIONS; ++i) {
		auto start = Clock::now();

		auto future1 = async(randomSearch, f, lowerBound, upperBound, dim);
		auto future2 = async(randomSearch, f, lowerBound, upperBound, dim);

		pair<double, vector<double>> sol1 = future1.get();
		pair<double, vector<double>> sol2 = future2.get();

		double sol;
		vector<double> solValues;

		if (sol1.first < sol2.first) {
			sol = sol1.first;
			solValues = sol1.second;
		}
		else {
			sol = sol2.first;
			solValues = sol2.second;
		}

		if (sol < minSol) {
			minSol = sol;
			minSolValues = solValues;
		}
		if (sol > maxSol) {
			maxSol = sol;
			maxSolValues = solValues;
		}
		avgSol += sol;

		auto finish = Clock::now();

		int time = chrono::duration_cast<chrono::milliseconds>(finish - start).count();
		if (time < minTime) minTime = time;
		if (time > maxTime) maxTime = time;
		avgTime += time;

		cout << i + 1 << ". Min: " << sol << " At: x=(";
		for (int i = 0; i < solValues.size(); ++i) {
			cout << solValues[i] << ",";
		}
		cout << ") Time: " << time << "ms" << endl;
	}

	cout << endl << "Stats:" << endl;
	cout << "Min solution: " << minSol << " At: x=(";
	for (int i = 0; i < minSolValues.size(); ++i) {
		cout << minSolValues[i] << ",";
	}
	cout << ")" << endl;
	cout << "Max solution: " << maxSol << " At: x=(";
	for (int i = 0; i < maxSolValues.size(); ++i) {
		cout << maxSolValues[i] << ",";
	}
	cout << ")" << endl;
	cout << "Avg solution: " << avgSol / ITERATIONS << endl << endl;

	cout << "Min time: " << minTime << "ms" << endl
		<< "Max time: " << maxTime << "ms" << endl
		<< "Avg time: " << avgTime / ITERATIONS << "ms" << endl;
}

int main()
{
	int n;
	cout << "Enter the number of dimensions: "; cin >> n;

	
	cout << "RASTRIGIN RESULTS FOR " << n << " DIMENSIONS..." << endl << endl;
	cout << "=============LAZY============" << endl;
	lazyDeterministicApproach(Rastrigin, -5.12, 5.12, n);
	cout << "=======NEDETERMINISTIC=======" << endl;
	nedeterministicApproach(Rastrigin, -5.12, 5.12, n);
	cout << "========DETERMINISTIC========" << endl;
	deterministicApproach(Rastrigin, -5.12, 5.12, n);
	
	
	cout << "MICHALEWICZ RESULTS FOR " << n << " DIMENSIONS..." << endl << endl;
	cout << "=============LAZY============" << endl;
	lazyDeterministicApproach(Michalewicz, 0, pi, n);
	cout << "=======NEDETERMINISTIC=======" << endl;
	nedeterministicApproach(Michalewicz, 0, pi, n);
	cout << "========DETERMINISTIC========" << endl;
	deterministicApproach(Michalewicz, 0, pi, n);
	

	cout << "SPHERE RESULTS FOR " << n << " DIMENSIONS..." << endl << endl;
	cout << "=============LAZY============" << endl;
	lazyDeterministicApproach(Sphere, -5.12, 5.12, n);
	cout << "=======NEDETERMINISTIC=======" << endl;
	nedeterministicApproach(Sphere, -5.12, 5.12, n);
	cout << "========DETERMINISTIC========" << endl;
	deterministicApproach(Sphere, -5.12, 5.12, n);
	
	
	cout << "ROSENBROCK RESULTS FOR " << n << " DIMENSIONS..." << endl << endl;
	cout << "=============LAZY============" << endl;
	lazyDeterministicApproach(Rosenbrock, -5, 10, n);
	cout << "=======NEDETERMINISTIC=======" << endl;
	nedeterministicApproach(Rosenbrock, -5, 10, n);
	cout << "========DETERMINISTIC========" << endl;
	deterministicApproach(Rosenbrock, -5, 10, n);
	
}