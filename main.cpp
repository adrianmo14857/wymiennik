#include "exch.h"

int main() {
	
	double Fzm = 22;
	double To = -20;
	double Tzm = 70 - 2.5 * (To - 6);
	double Fzco = 41;
	double Tpco = 45;

	double A1 = (-Fzm * ro * cw - kw) / (Mm * cwym);
	double B1 = -kw / (Mm * cwym);
	double WW1 = ((Fzm * ro * cw) / (Mm * cwym)) * Tzm;
	double A2 = kw / (Mco * cwym);
	double B2 = -(Fzco * ro * cw + kw) / (Mco * cwym);
	double WW2 = ((Fzco * ro * cw) / (Mco * cwym)) * Tpco;

	
	matrix A = { { A1, A2},
				{ B1, B2},
	};
	vec b = { WW1, WW2 };

	vec X = calculateMarkov(1000, A, b);

	for (double x : X) cout << x << '\t';

	cout << endl;
		
	return 0;
}