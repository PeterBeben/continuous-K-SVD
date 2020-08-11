//     Copyright (C) 2019 Piotr (Peter) Beben <pdbcas@gmail.com>
//     See LICENSE included.

#include "MatchingPursuit.h"
#include "OrthogonalPursuit.h"
#include "ksvd.h"
#include "ksvd_dct2D.h"
#include "cosine_transform.h"

#include <Eigen/Dense>
#include <QCoreApplication>
#include <iostream>
#include <ctime>
#include <omp.h>

using std::cout;
using std::endl;
using std::vector;
using Eigen::MatrixXf;
using Eigen::Matrix;
using Eigen::VectorXf;
using Eigen::Index;
using Eigen::Map;
using Eigen::Dynamic;

using vectorf = vector<float, Eigen::aligned_allocator<float>>;


void test_ksvd();
void test_ksvd_dct2D();

//-----------------------------------------------------------

int main(int argc, char *argv[])
{

	QCoreApplication a(argc, argv);
	Eigen::initParallel();

	test_ksvd();
	test_ksvd_dct2D();

	return a.exec();
}


//-----------------------------------------------------------


void test_ksvd()
{

	Index ndim = 25;             // Dimension of atoms and signals.
	Index nsig = 500;            // No. of signal.
	Index natm = 100;            // No. of atoms.
	Index latm = 8;              // Sparsity constraint.
	int maxiters = 8;            // No. of K-SVD iterations.
	bool useOpenMP = true;



	MatchingPursuit mp;
	OrthogonalPursuit op;

	srand(2);
	//std::srand(std::time(nullptr));

	MatrixXf Y = MatrixXf::Random(ndim,nsig);
	MatrixXf D = MatrixXf::Random(ndim,natm);
	MatrixXf X = MatrixXf::Random(natm,nsig);

	D.colwise().normalize();

	MatrixXf Dorig = D;
	MatrixXf Xorig = X;

	clock_t c_start;
	std::clock_t c_end;
	double time_elapsed_ms;


	cout << "\n\n** Testing ksvd **" ;
	cout << "\n<Error measured as average coordinate difference>\n";

	cout << "\n* Using matching pursuit *" ;
	cout << "\nInitial error: " ;
	cout << (Y-(D*X)).cwiseAbs().sum()/(ndim*nsig);

	c_start = std::clock();
	ksvd(useOpenMP, Y, latm, maxiters, 0.0, 1, mp, D, X);
	c_end = std::clock();

	cout << "\nError after " << maxiters <<" iterations: " ;
	cout << (Y-(D*X)).cwiseAbs().sum()/(ndim*nsig);
	time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "\nCPU time used: " << time_elapsed_ms << " ms\n";


	D = Dorig;
	X = Xorig;

	cout << "\n* Using orthogonal pursuit *" ;
	cout << "\nInitial error: " ;
	cout << (Y-(D*X)).cwiseAbs().sum()/(ndim*nsig);

	c_start = std::clock();
	ksvd(useOpenMP, Y, latm, maxiters, 0.0, 1, op, D, X);
	c_end = std::clock();

	cout << "\nError after " << maxiters <<" iterations: " ;
	cout << (Y-(D*X)).cwiseAbs().sum()/(ndim*nsig);
	time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "\nCPU time used: " << time_elapsed_ms << " ms\n";



}


//-----------------------------------------------------------


void test_ksvd_dct2D()
{

	Index nfreq = 5;             // No. of frequencies in cosine transform
	Index ndim = nfreq*nfreq;    // Dimension of atoms.
	Index nsig = 500;            // No. of signal.
	Index natm = 100;            // No. of atoms.
	Index latm = 6;              // Sparsity constraint.
	Index maxSamples = ndim;     // Max. signal samples.
	int maxiters = 18;           // No. of K-SVD iterations.
	bool useOpenMP = true;


	MatchingPursuit mp;
	OrthogonalPursuit op;

	srand(2);
	//std::srand(std::time(nullptr));

	vector<Index> ndimrand(nsig);
	for(Index isig=0; isig < nsig; ++isig){
		ndimrand[isig] = 1 + std::rand()/((RAND_MAX + 1u)/maxSamples);
	}

	vector<VectorXf> Ya(nsig);
	vector<VectorXf> Ua(nsig);
	vector<VectorXf> Va(nsig);
	MatrixXf D = MatrixXf::Random(ndim,natm);
	MatrixXf X = MatrixXf::Random(natm,nsig);

	D.colwise().normalize();

	for(Index isig=0; isig < nsig; ++isig){
		Index ndimsig = ndimrand[isig];
		Ya[isig] = VectorXf::Random(ndimsig);
		Ua[isig] = VectorXf::Random(ndimsig);
		Va[isig] = VectorXf::Random(ndimsig);
	}

	MatrixXf Dorig = D;
	MatrixXf Xorig = X;

	clock_t c_start;
	std::clock_t c_end;
	double time_elapsed_ms;


	cout << "\n\n** Testing ksvd_dct2D **" ;
	cout << "\n<Error measured as average (coord. diff., cosine angle, length vect. diff.)>\n" ;

	cout << "\n* Using matching pursuit *" ;
	cout << "\nInitial error: " ;
	print_error_dct2D(Ya, Ua, Va, D, X, nfreq);

	c_start = std::clock();
    ksvd_dct2D(useOpenMP, Ya, Ua, Va, nfreq, latm, maxiters, 0.0, mp, D, X);
	c_end = std::clock();

	cout << "\nError after " << maxiters <<" iterations: " ;
	print_error_dct2D(Ya, Ua, Va, D, X, nfreq);
	time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "\nCPU time used: " << time_elapsed_ms << " ms\n";


	D = Dorig;
	X = Xorig;

	cout << "\n* Using orthogonal pursuit *" ;
	cout << "\nInitial error: " ;
	print_error_dct2D(Ya, Ua, Va, D, X, nfreq);

	c_start = std::clock();
    ksvd_dct2D(useOpenMP, Ya, Ua, Va, nfreq, latm, maxiters, 0.0, op, D, X);
	c_end = std::clock();

	cout << "\nError after " << maxiters <<" iterations: " ;
	print_error_dct2D(Ya, Ua, Va, D, X, nfreq);
	time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "\nCPU time used: " << time_elapsed_ms << " ms\n";



}
