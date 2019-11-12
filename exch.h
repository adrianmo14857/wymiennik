#ifndef EXCH_H
#define EXCH_H

#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <vector>

using namespace std;
using vec = vector<double>;
using matrix = vector<vec>;

const double Mm = 3000;	//Masa wymiennika po stronie pierwotnej
const double Mco = 3000;	//Masa wymiennika po stronie wtornej
const double cwym = 2700;	//zastepcze cieplo wlasciwe wymiennika
const double ro = 1000;	//gestosc wody
const double cw = 4200;	//cieplo wlasciwe wody
const double kw = 250000;	//wspolczynnik przenikania ciepla

double myRand();
vec Markov(int num, const matrix& H, const vec& b);
vec calculateMarkov(int num, matrix A, vec b);


#endif //EXCH_H
