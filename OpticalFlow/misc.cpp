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
#include "maths.h"
#include "misc.h"
#if defined(_WIN32)
#include <windows.h> 
#include <ctime>
#endif
#if defined(linux)
#include <stddef.h>
#include <stdio.h> 
#include <sys/types.h>
#include <dirent.h>
#endif

namespace HSCNameSpace{


	Matrix<float> *CropImgMatrix(IplImage *img, int x0, int y0, int x1, int y1) {
		Matrix<float> *m;
		int x, y;

		CV_FUNCNAME("CropImgMatrix");

		if (img->nChannels > 1) {
			DET_ERROR("image must have a single channel");
		}

		m = new Matrix<float>(y1 - y0 + 1, x1 - x0 + 1);
		for (y = y0; y <= y1; y++) {
			for (x = x0; x <= x1; x++) {
				//m->GetElement(x, y)
				m->SetValue(x - x0, y - y0, (float)GetPixel1Channel(img, x, y));
			}
		}
		cv::Mat mattemp = cv::Mat(m->GetNCols(), m->GetNRows(), CV_32FC1, m->GetData());
		cv::imwrite("firstFeaturewithin.jpg", mattemp);

		cvSaveImage("img.jpg", img);

		return m;
	}


	Matrix<float> *FloattoMatrix(float *v, int nrows, int ncols) {
		Matrix<float> *m;
		int y, x, idx;

		m = new Matrix<float>(nrows, ncols);

		idx = 0;
		for (x = 0; x < ncols; x++) {
			for (y = 0; y < nrows; y++) {
				m->SetValue(x, y, v[idx++]);
			}
		}

		return m;
	}


	/**********************************/
	/* Laplacian pyramid decompostion */
	/**********************************/
	Vector<float> *LapPyrpfiltersH(string fname) {
		Vector<float> *v = NULL;

		CV_FUNCNAME("LapPyrpfiltersH");

		if (fname.compare("9/7") == 0 || fname.compare("9-7") == 0) {
			v = new Vector<float>(9);
			v->SetElement(0, 0.037828455506995f);
			v->SetElement(1, -0.02384946501938f);
			v->SetElement(2, -0.11062440441842f);
			v->SetElement(3, 0.37740285561265f);
			v->SetElement(4, 0.8526986790094f);
			v->SetElement(5, 0.37740285561265f);
			v->SetElement(6, -0.11062440441842f);
			v->SetElement(7, -0.02384946501938f);
			v->SetElement(8, 0.037828455506995f);
		}

		else if (fname.compare("5-3") == 0 || fname.compare("5/3") == 0) {
			v = new Vector<float>(5);
			v->SetElement(0, -0.176776695296637f);
			v->SetElement(1, 0.353553390593274f);
			v->SetElement(2, 1.06066017177982f);
			v->SetElement(3, 0.353553390593274f);
			v->SetElement(4, -0.176776695296637f);
		}

		else if (fname.compare("Burt") == 0){
			v = new Vector<float>(5);
			v->SetElement(0, -0.0707106781186548f);
			v->SetElement(1, 0.353553390593274f);
			v->SetElement(2, 0.848528137423857f);
			v->SetElement(3, 0.353553390593274f);
			v->SetElement(4, -0.0707106781186548f);
		}

		else {
			DET_ERROR("Unknown parameter for pfilters");
		}

		return v;
	}


	Vector<float> *LapPyrpfiltersG(string fname) {
		Vector<float> *v = NULL;

		CV_FUNCNAME("LapPyrpfiltersG");

		if (fname.compare("9/7") == 0 || fname.compare("9-7") == 0) {
			v = new Vector<float>(7);
			v->SetElement(0, -0.064538882628938f);
			v->SetElement(1, -0.040689417609558f);
			v->SetElement(2, 0.41809227322221f);
			v->SetElement(3, 0.78848561640566f);
			v->SetElement(4, 0.41809227322221f);
			v->SetElement(5, -0.040689417609558f);
			v->SetElement(6, -0.064538882628938f);
		}

		else if (fname.compare("5-3") == 0 || fname.compare("5/3") == 0) {
			v = new Vector<float>(3);
			v->SetElement(0, 0.353553390593274f);
			v->SetElement(1, 0.707106781186547f);
			v->SetElement(2, 0.353553390593274f);
		}

		else if (fname.compare("Burt") == 0){
			v = new Vector<float>(7);
			v->SetElement(0, -0.0151522881682832f);
			v->SetElement(1, -0.0757614408414158f);
			v->SetElement(2, 0.368705678761557f);
			v->SetElement(3, 0.858629662869379f);
			v->SetElement(4, 0.368705678761557f);
			v->SetElement(5, -0.0757614408414158f);
			v->SetElement(6, -0.0151522881682832f);
		}

		else {
			DET_ERROR("Unknown parameter for pfilters");
		}

		return v;
	}


	// Laplacian Pyramid Decomposition
	void LapPyrlpdec(Matrix<float> *x, Vector<float> *h, Vector<float> *g, Matrix<float> **c, Matrix<float> **d) {
		int i, xa, ya, xx, yy;
		vector<int> selectedCols, selectedRows;
		int adjust;
		Matrix<float> *xlo, *maux;
		vector<int> shift(2);

		shift[0] = 0;
		shift[1] = 0;

		xlo = LapPyrsefilter2(x, h, h, shift);

		// select rows and columns
		for (i = 0; i < xlo->GetNRows(); i += 2) {
			selectedRows.push_back(i);
		}
		for (i = 0; i < xlo->GetNCols(); i += 2) {
			selectedCols.push_back(i);
		}
		(*c) = xlo->GetSelectedColsRows(selectedCols, selectedRows);

		adjust = (g->GetNElements() + 1) % 2;


		delete xlo;
		xlo = new Matrix<float>(x->GetNRows(), x->GetNCols());
		xlo->SetAllElements(0.0f);


		ya = 0;
		for (yy = 0; yy < (*c)->GetNRows(); yy++) {
			xa = 0;
			for (xx = 0; xx < (*c)->GetNCols(); xx++) {
				xlo->SetValue(xa, ya, (*c)->GetElement(xx, yy));
				xa += 2;
			}
			ya += 2;
		}


		shift[0] = adjust;
		shift[1] = adjust;

		maux = LapPyrsefilter2(xlo, g, g, shift);

		(*d) = new Matrix<float>(x->GetNRows(), x->GetNCols());

		for (yy = 0; yy < x->GetNRows(); yy++) {
			for (xx = 0; xx < x->GetNCols(); xx++) {
				(*d)->SetValue(xx, yy, x->GetElement(xx, yy) - maux->GetElement(xx, yy));
			}
		}

		delete xlo;
		delete maux;
	}



	// 2D seperable filtering with extension handling
	Matrix<float> *LapPyrsefilter2(Matrix<float> *x, Vector<float> *f1, Vector<float> *f2, vector<int> &shift) {
		float lf1, lf2;
		Matrix<float> *y, *mret;

		lf1 = ((float)f1->GetNElements() - 1.0f) / 2.0f;

		lf2 = ((float)f2->GetNElements() - 1.0f) / 2.0f;

		y = LapPyrextend2(x, (int)floor(lf1) + shift[0], (int)ceil(lf1) - shift[0],
			(int)floor(lf2) + shift[1], (int)ceil(lf2) - shift[1]);

		mret = LapPyrconv2Separable(f1, f2, y, "valid");
		delete y;

		return mret;
	}



	Matrix<float> *LapPyrextend2(Matrix<float> *x, int ru, int rd, int cl, int cr) {
		vector<int> *I;
		Matrix<float> *maux, *mret;
		int rx, cx;

		rx = x->GetNRows();
		cx = x->GetNCols();

		// for rows
		I = LapPyrgetPerIndices(rx, ru, rd);
		maux = x->GetSelectedRows(I);
		delete I;

		// for columns
		I = LapPyrgetPerIndices(cx, cl, cr);
		mret = maux->GetSelectedCols(I);
		delete I;
		delete maux;

		return mret;
	}



	vector<int> *LapPyrgetPerIndices(int lx, int lb, int le) {
		vector<int> *I;
		int i;

		CV_FUNCNAME("LapPyrgetPerIndices");

		I = new vector<int>;

		for (i = 0; i < lb; i++) {
			I->push_back(lx - lb + i);
		}
		for (i = 0; i < lx; i++) {
			I->push_back(i);
		}
		for (i = 0; i < le; i++) {
			I->push_back(i);
		}

		if (lx < lb || lx < le) {
			DET_ERROR("Conditions not checked");
		}

		return I;
	}


	// Two dimensional matrix convolution (separable filters)
	Matrix<float> *LapPyrconv2Separable(Vector<float> *h1, Vector<float> *h2, Matrix<float> *A, string shape) {
		int x, y, j;
		int mc, nc;
		Matrix<float> *C, *C2, *maux;
		float value;
		vector<int> selectCols, selectRows;

		CV_FUNCNAME("LapPyrconv2Separable");

		mc = max(max(A->GetNRows() + h1->GetNElements() - 1, A->GetNRows()), h1->GetNElements());
		nc = max(max(A->GetNCols() + h2->GetNElements() - 1, A->GetNCols()), h2->GetNElements());

		C = new Matrix<float>(A->GetNRows(), nc);
		C->SetAllElements(0.0f);

		// convolve rows
		for (y = 0; y < C->GetNRows(); y++) {
			for (x = 0; x < C->GetNCols(); x++) {
				value = 0;
				for (j = 0; j < h1->GetNElements(); j++) {
					if (x - j >= 0 && x - j < A->GetNCols())
						value += A->GetElement(x - j, y) * h1->GetElement(j);
				}
				C->SetValue(x, y, value);
			}
		}
		C2 = new Matrix<float>(mc, C->GetNCols());
		C2->SetAllElements(0.0f);

		// convolve cols
		for (y = 0; y < C2->GetNRows(); y++) {
			for (x = 0; x < C2->GetNCols(); x++) {
				value = 0;
				for (j = 0; j < h2->GetNElements(); j++) {
					if (y - j >= 0 && y - j < C->GetNRows())
						value += C->GetElement(x, y - j) * h2->GetElement(j);
				}
				C2->SetValue(x, y, value);
			}
		}

		if (shape.compare("valid") == 0) {
			for (j = h1->GetNElements() - 1; j < C2->GetNCols() - h1->GetNElements() + 1; j++) {
				selectCols.push_back(j);
			}
			for (j = h2->GetNElements() - 1; j < C2->GetNRows() - h2->GetNElements() + 1; j++) {
				selectRows.push_back(j);
			}

			maux = C2->GetSelectedColsRows(selectCols, selectRows);
			delete C2;
			C2 = maux;
		}

		else if (shape.compare("same") == 0) {
			for (j = (h1->GetNElements() / 2); j < C2->GetNCols() - (h1->GetNElements() / 2); j++) {
				selectCols.push_back(j);
			}
			for (j = (h2->GetNElements() / 2); j < C2->GetNRows() - (h2->GetNElements() / 2); j++) {
				selectRows.push_back(j);
			}

			maux = C2->GetSelectedColsRows(selectCols, selectRows);
			delete C2;
			C2 = maux;
		}

		else if (shape.compare("full") == 0) {
			;
		}

		else {
			DET_ERROR("invalid shape type for convolution");
		}

		delete C;
		return C2;
	}




	// Two dimensional matrix convolution 
	Matrix<float> *LapPyrconv2(Matrix<float> *A, Matrix<float> *B, string shape) {
		int x, y, i, j;
		int mc, nc;
		Matrix<float> *C, *maux;
		float value;
		vector<int> selectCols, selectRows;

		CV_FUNCNAME("LapPyrconv2");

		mc = max(max(A->GetNRows() + B->GetNRows() - 1, A->GetNRows()), B->GetNRows());
		nc = max(max(A->GetNCols() + B->GetNCols() - 1, A->GetNCols()), B->GetNCols());

		C = new Matrix<float>(mc, nc);
		C->SetAllElements(0.0f);

		// apply convolution
		for (y = 0; y < C->GetNRows(); y++) {
			for (x = 0; x < C->GetNCols(); x++) {
				value = 0;
				for (j = 0; j < B->GetNRows(); j++) {
					for (i = 0; i < B->GetNCols(); i++) {
						if (x - i >= 0 && x - i < A->GetNCols() && y - j >= 0 && y - j < A->GetNRows())
							value += A->GetElement(x - i, y - j) * B->GetElement(i, j);
					}
				}
				C->SetValue(x, y, value);
			}
		}

		if (shape.compare("same") == 0) {
			for (j = (B->GetNCols() / 2); j < C->GetNCols() - (B->GetNCols() / 2) + 1; j++) {
				selectCols.push_back(j);
			}
			for (j = (B->GetNRows() / 2); j < C->GetNRows() - (B->GetNRows() / 2) + 1; j++) {
				selectRows.push_back(j);
			}

			maux = C->GetSelectedColsRows(selectCols, selectRows);
			delete C;
			C = maux;
		}

		else if (shape.compare("full") == 0) {
			;
		}

		else {
			DET_ERROR("invalid shape type for convolution");
		}

		return C;
	}





	void Shearletsapply(Matrix<float> *x, int ndirections, Matrix<float> **xlo, vector<Matrix<float> *> &d,
		vector<Matrix<float> *> &shearMat, Vector<float> *g, Vector<float> *h, int nlevels) {
		Matrix<float> *xhi, *maux;
		int i;



		// execute pyramid decomposition
		LapPyrlpdec(x, h, g, xlo, &xhi);

		//todo for testing
		cv::Mat mattempshearlet = cv::Mat(xhi->GetNCols(), xhi->GetNRows(), CV_32FC1, xhi->GetData());
		cv::imwrite("ConvolutionOriginal.jpg", mattempshearlet);

		// apply convolution for each matrix
		for (i = 0; i < (int)shearMat.size(); i++) {
			maux = LapPyrconv2(xhi, shearMat.at(i), "same");
			d.push_back(maux);

			//todo: for testing 
			char shearletfilename[40];
			sprintf(shearletfilename, "Level%dConvolution%d.jpg", nlevels, i);
			cv::Mat mattempshearlet = cv::Mat(maux->GetNCols(), maux->GetNRows(), CV_32FC1, maux->GetData());
			cv::imwrite(shearletfilename, mattempshearlet);

		}

		delete xhi;
	}
}