#include<iostream>
#include<vector>
#include<chrono>
#include <thread>
using namespace std;
using namespace std::chrono;
class arrayData { // used for storing indexs and sum of sub arrays

public:
	arrayData();
	arrayData(const int&, const int&, const int&);
	arrayData(const arrayData &);
	~arrayData();
	void printData(); // prints the indexes and sum for testing
	int lowIndex;
	int highIndex;
	int maxSum;

};
const int bounds = 300; // will make random numbers between -bounds-1 and bounds 
const int m = 10; // number of trials for each run.
void vprint(const vector<int>&); //  prints vector for testing
vector<int> makeranv(const int &); // takes in size and makes a random vecotr of ints from +/- bounds of that size, returns that vector
arrayData BRUTE_FORCE(const vector<int> &); // calculates max sub array indexes and sum using brute force
arrayData FIND_MAX_SUBARRAY(const vector<int> &, const int &low, const int &high); // calculates max sub array indexes and sum using divide-and-conquer
arrayData FIND_MAX_CROSSING_SUBARRAY(const vector<int> &A, const int &, const int &, const int &); // used to calculate the max sub array across the subaraay in linear time
void main() {

	vector<vector<int>> A, run1, run2; //run1=Bru, etc...
	vector<vector<double>>  run1T, run2T;
	vector<double>run1Tavg, run2Tavg;
	vector<double> tempV;
	duration<double, std::milli> time; // total time in milliseconds
	high_resolution_clock::time_point t1; // start time
	high_resolution_clock::time_point t2; // end time

	for (size_t i = 5000; i <= 100000; i += 5000) { // initializing all the test cases
		A.push_back(makeranv(i));
	}
	cout << "ALG1 BEGIN" << endl;
	for (size_t i = 0; i < A.size(); i++) {// run2/BRUTE-FORCE
		cout << "run " << i + 1 << " " << endl;
		for (size_t j = 0; j < m; j++) {

			t1 = high_resolution_clock::now();
			BRUTE_FORCE(A[i]).printData();
			t2 = high_resolution_clock::now();
			time = t2 - t1;
			tempV.push_back(time.count());
		}
		run1T.push_back(tempV);
		tempV.clear();
	}
	cout << "ALG1 END, ALG2 BEGIN" << endl;
	for (size_t i = 0; i < A.size(); i++) {// run2/DIVIDE-AND-CONQUER
		cout << "run " << i + 1 << " " << endl;
		for (size_t j = 0; j < m; j++) {

			t1 = high_resolution_clock::now();
			FIND_MAX_SUBARRAY(A[i], 0, A[i].size() - 1).printData();
			t2 = high_resolution_clock::now();
			time = t2 - t1;
			tempV.push_back(time.count());
		}
		run2T.push_back(tempV);
		tempV.clear();
	}

	// setting all to zero for later  average calculations
	for (size_t i = 0; i < A.size(); i++) {
		run1Tavg.push_back(0);
		run2Tavg.push_back(0);
	}

	for (size_t i = 0; i < A.size(); i++) { // calulating averge RT for run1/BRUTE-FORCE
		for (size_t j = 0; j < m; j++)
		{
			run1Tavg[i] += run1T[i][j];
		}
		run1Tavg[i] /= m;
	}
	for (size_t i = 0; i < A.size(); i++) {// DIVIDE-AND-CONQUER average rt
		for (size_t j = 0; j < m; j++) {
			run2Tavg[i] += run2T[i][j];
		}
		run2Tavg[i] /= m;
	}
	cout << "Run times for BRUTE-FORCE:" << endl;
	for (size_t i = 0; i < A.size(); i++) {// output
		cout << run1Tavg[i] << " ";
	}cout << endl;
	cout << "Run times for DIVIDE-AND-CONQUER:" << endl;
	for (size_t i = 0; i < A.size(); i++) {// output
		cout << run2Tavg[i] << " ";
	}cout << endl;

	
}
arrayData BRUTE_FORCE(const vector<int> &A) { // RT=O(n^2) 
	int max = 0;
	int maxI = 0;
	int maxJ = -1;
	for (size_t i = 0; i < A.size(); i++) {
		int sum = 0;
		for (size_t j = i; j < A.size(); j++) {
			sum += A[j]; //sum is that of A[i..j]
			if (sum>max){ // stores the indexes of the max sub array so far
				max = sum;
				maxI = i;
				maxJ = j;
			}
		}
	}
	arrayData temp(maxI, maxJ, max);
	return temp;
}
arrayData FIND_MAX_SUBARRAY(const vector<int> &A, const int &low, const int &high) { // RT = O(nlgn)

	if (high == low) {
		arrayData temp(low, high, A[low]);
		return temp;
	}
	else {
		int mid = (low + high) / 2;

		arrayData left = FIND_MAX_SUBARRAY(A, low, mid);
		arrayData right = FIND_MAX_SUBARRAY(A, mid + 1, high);
		arrayData cross = FIND_MAX_CROSSING_SUBARRAY(A, low, mid, high);

		if (left.maxSum >= right.maxSum && left.maxSum >= cross.maxSum) { // returns sub array left data if has greatest sum
			return left;
		}
		else if (right.maxSum >= left.maxSum && right.maxSum >= cross.maxSum) { // returns sub array right data if has greatest sum
			return right;
		}
		else { // returns sub array cross data if has greatest sum
			return cross;
		}

	}
}
arrayData FIND_MAX_CROSSING_SUBARRAY(const vector<int> &A, const int &low, const int &mid, const int &high) { //RT= O(n) 
	int leftSum = -INT_MAX;
	int rightSum = -INT_MAX;
	int sum = 0;
	int indexL, indexR;

	for (int i = mid; i >= low; i--) {
		sum += A[i];
		if (sum > leftSum) {
			leftSum = sum;
			indexL = i;
		}
	}
	sum = 0;

	for (int j = mid + 1; j <= high; j++) {
		sum += A[j];
		if (sum > rightSum) {
			rightSum = sum;
			indexR = j;
		}
	}
	arrayData temp(indexL, indexR, leftSum + rightSum); // set of data (indexes and sum)
	return temp;
}
arrayData::arrayData() { //constructor
	lowIndex = highIndex = maxSum = NULL;
}
arrayData::arrayData(const int &low, const int &high, const int &sum) { // constructor
	lowIndex = low;
	highIndex = high;
	maxSum = sum;
}
arrayData::arrayData(const arrayData &data) { // copy constructor
	lowIndex = data.lowIndex;
	highIndex = data.highIndex;
	maxSum = data.maxSum;
}
arrayData::~arrayData() { //destructor
	lowIndex = highIndex = maxSum = NULL;
}
void arrayData::printData() { // for testing
	cout << "lowIndex: " << lowIndex
		<< " highIndex: " << highIndex
		<< " maxSum: " << maxSum << endl;
}
void vprint(const vector<int>&A) { // prints vector
	for (size_t i = 0; i < A.size(); i++) {
		cout << A[i] << " ";
	}
	cout << endl;
}
vector<int> makeranv(const int &size) { // takes in size and makes a random vector of ints from -bound-1 and  bounds of that size, returns vector
	vector<int> V;
	for (int i = 0; i < size; i++) {
		V.push_back((rand() % (2 * bounds) + 1) - bounds);
	}
	return V;
}
