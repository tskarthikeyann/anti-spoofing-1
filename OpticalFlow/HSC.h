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

#ifndef HSC_H
#define HSC_H
#include "misc.h"
#include "maths.h"

namespace HSCNameSpace{
	class HSC {
		int nlevels;
		int norient;
		int maskSize;
		bool useCells;
		string LapFilter;	// filter used for Laplacian decomposition
		string normType;	// normalization type (global, local, local+global, none)
		Vector<float> *g;
		Vector<float> *h;
		vector<Matrix<float> *> matrixShear;
		vector<vector<Matrix<float> *> > dFeatures;	// keep computed features for the current image

		// maths for speed-up (for normalization)
		Maths mat;

		// sum a region of an matrix
		float SumRegion(Matrix<float> *m, int x0, int y0, int x1, int y1);

		// extract features withour selection
		void ExtractFeatures(Position *pos, float *feat);

		// reset some structures after extraction features from an image
		void ResetStructures();

	public:
		HSC();

		// add new image for this extraction method
		void AddNewImage(IplImage *img);

		// initilize extraction method
		void InitializeExtractionMethod();

		// get number of features per block
		int GetNFeatures();

		// set parameters
		void SetParameter(string name, string value);

		// function to extract feature vectors
		Vector<float> *ExtractFeatures(int blockW, int blockH, int strideX, int strideY);

		// function to extract feature vectors considering image as a whole (not dividing in blocks)
		Vector<float> *ExtractFeatures();
	};
}
#endif