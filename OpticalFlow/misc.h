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

#ifndef MISC_H
#define MISC_H
#include "maths.h"

namespace HSCNameSpace{
	// get pixel from 1-channel image
#define GetPixel1Channel(img, x, y) ((int) ((uchar *)((img)->imageData + (img)->widthStep * (y)))[x])

	// crop image into a matrix
	Matrix<float> *CropImgMatrix(IplImage *img, int x0, int y0, int x1, int y1);

	// float vector to a matrix (assuming that vector was built columnwise)
	Matrix<float> *FloattoMatrix(float *v, int nrows, int ncols);

	/**********************************/
	/* Laplacian pyramid decompostion */
	/**********************************/
	// generate filters for the Laplacian pyramid
	// h filter
	Vector<float> *LapPyrpfiltersH(string fname);

	// g filter
	Vector<float> *LapPyrpfiltersG(string fname);

	// Laplacian Pyramid Decomposition
	void LapPyrlpdec(Matrix<float> *x, Vector<float> *h, Vector<float> *g, Matrix<float> **xlo, Matrix<float> **xhi);

	// 2D seperable filtering with extension handling
	Matrix<float> *LapPyrsefilter2(Matrix<float> *x, Vector<float> *f1, Vector<float> *f2, vector<int> &shift);

	// 2D extension
	Matrix<float> *LapPyrextend2(Matrix<float> *x, int ru, int rd, int cl, int cr);

	// 
	vector<int> *LapPyrgetPerIndices(int lx, int lb, int le);

	// Two dimensional matrix convolution (separable filters).
	Matrix<float> *LapPyrconv2Separable(Vector<float> *h1, Vector<float> *h2, Matrix<float> *A, string shape);

	// Two dimensional matrix convolution
	Matrix<float> *LapPyrconv2(Matrix<float> *A, Matrix<float> *B, string shape);

	// execute shearlets given an input image
	void Shearletsapply(Matrix<float> *x, int ndirections, Matrix<float> **xlo, vector<Matrix<float> *> &d,
		vector<Matrix<float> *> &shearMat, Vector<float> *g, Vector<float> *h, int nlevels);

}
#endif