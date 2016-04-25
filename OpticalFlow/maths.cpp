//******************************************************************************
// HSC - Histogram of Shearlet Coefficients.
// Copyright (C) 2010-2011 by William Robson Schwartz.
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/lgpl.html>
//******************************************************************************

#include "headers.h"
#include <xmmintrin.h>
#include <emmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <emmintrin.h>
#include "maths.h"

namespace HSCNameSpace{


	float Maths::DotProductSSENotMultof4(float *v, float *v2, int n) {
		__m128 a, b, c, d;
		float out_sse[4];
		int i, n2;

		d = _mm_setzero_ps();
		if ((n & 3) == 0) { // multiple of 4
			for (i = 0; i < n; i += 4) {
				a = _mm_loadu_ps(v + i);
				b = _mm_loadu_ps(v2 + i);
				c = _mm_mul_ps(a, b);
				d = _mm_add_ps(c, d);
			}
			_mm_storeu_ps(out_sse, d);
		}
		else {  // n not multiple of 4
			n2 = n - 4;
			for (i = 0; i < n2; i += 4) {
				a = _mm_loadu_ps(v + i);
				b = _mm_loadu_ps(v2 + i);
				c = _mm_mul_ps(a, b);
				d = _mm_add_ps(c, d);
			}
			_mm_storeu_ps(out_sse, d);
			n2 = n - (n & 0x3);
			// do the remaining elements
			for (i = n2; i < n; i++) {
				out_sse[0] += v[i] * v2[i];
			}
		}

		return out_sse[0] + out_sse[1] + out_sse[2] + out_sse[3];
	}




	void Maths::DivideVectorsSSE(float *v1, float value, float *outputvect, int n) {
		__m128 a, b, c;
		int i, n2;

		if ((n & 3) == 0) { // multiple of 4
			b = _mm_load1_ps(&value);
			for (i = 0; i < n; i += 4) {
				a = _mm_loadu_ps(v1 + i);
				c = _mm_div_ps(a, b);
				_mm_storeu_ps(outputvect + i, c);
			}
		}
		else {
			n2 = n - 4;
			b = _mm_load1_ps(&value);
			for (i = 0; i < n2; i += 4) {
				a = _mm_loadu_ps(v1 + i);
				c = _mm_div_ps(a, b);
				_mm_storeu_ps(outputvect + i, c);
			}

			n2 = n - (n & 0x3);
			for (i = n2; i < n; i++) {
				outputvect[i] = v1[i] / value;
			}
		}
	}
}