/**
* @file test-0038.c
*
* @brief Plus of the constatnts/variables.
*/

#include <limits.h>
#include <math.h>

int main(int argc, const char *argv[]) {
	// Plus of constatnts.
	int a = 4 + 5;

	// Plus of a variable and constant.
	int b = a + 1;

	// Plus of variables.
	int c = a + b;

	// Overflows.
	int d = INT_MAX + 1;
	unsigned e = UINT_MAX + 1;
	unsigned char f = 255 + 255;

	// Special float values;
	float nan = NAN;
	float inf = INFINITY;

	//  INF + INF = INF
	float g = inf + inf;

	//  INF + (-INF) = -NAN
	float h = inf + -inf;

	//  INF + NAN = NAN
	float i = inf + nan;

	//  INF + number = INF
	float j = inf + 10;

	// -INF + INF = -NAN
	float k = -inf + inf;

	// -INF + (-INF) = -INF
	float l = -inf + -inf;

	// -INF + NAN = NAN
	float m = -inf + -nan;

	// -INF + number = -INF
	float n = -inf + 158;

	//  NAN + INF = NAN
	float o = nan + inf;

	//  NAN + (-INF) = NAN
	float p = nan + -inf;

	//  NAN + NAN = NAN
	float q = nan + nan;

	//  NAN + number = NAN
	float r = nan + 14789;

	//  number + INF = INF
	float s = 1589 + inf;

	//  number + (-INF) = -INF
	float t = 148997 + -inf;

	//  number + NAN = NAN
	float u = 148997 + nan;

	return 0;
}

