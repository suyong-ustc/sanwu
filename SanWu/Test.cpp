#include "Test.h"
#include <armadillo>
#include "DIC/DICAlgorithm.h"

using namespace arma;

bool TestNormalizeVectorize()
{
	mat a = randu(3, 3);

	mat w = zeros(3, 3);
	for (int r = 1; r < 3; ++r)
		for (int c = 1; c < 3; ++c)
			w(r, c) = 1;

	const double zncc = ZNCC(a, a, w);
	std::cout << "the correaltion function is " << zncc << std::endl;

	return true;
}
