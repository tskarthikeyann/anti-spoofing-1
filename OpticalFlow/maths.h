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

#ifndef MATHS_H
#define MATHS_H
#include "headers.h"

#define EPS 0.00000001

namespace HSCNameSpace{
	template <class T> class Vector;



	//////////////////
	// class Vector //
	//////////////////
	template <class T>
	class Vector {
		T *realp; // real pointer before align
		T *v;     // data poiter
		int n;        // number of elements

		// allocate memory for the vector
		void AllocateMemory(int n) {
			CV_FUNCNAME("Vector::AllocateMemory");

			realp = (T *)malloc(sizeof(T)* n);

			if (realp == NULL)
				DET_ERROR("No available memory to allocate vector");

			memset(realp, 0, sizeof(T)* n);
			v = realp;
			this->n = n;
			// if change here, check setData function
		}


		void LoadVector(char *filename) {
			FILE *f;
			int c, r;

			f = fopen(filename, "rb");
			fread(&r, sizeof(int), 1, f);
			fread(&c, sizeof(int), 1, f);
			AllocateMemory(max(r, c));

			fread(v, sizeof(T), max(r, c), f);

			fclose(f);
		}


	public:
		// Initialize
		Vector(int n) { AllocateMemory(n); }

		Vector(bool onlyHeader) {
			realp = NULL;
			v = NULL;
			n = 0;
		}

		Vector(string filename) { LoadVector((char *)filename.c_str()); }

		// load a vector
		Vector(char *filename) { LoadVector(filename); }

		// release structures
		~Vector() {
			free(realp);
		}

		// return number of bytes used for this vector
		int GetNumBytes() { return n * sizeof(T); }


		// retrieve element at row y and column x
		inline T GetElement(int x) { return v[x]; }
		inline T operator[](int x) { return v[x]; }


		// set data to a vector - at your risk!!
		void SetData(T *pointer) {
			free(realp);
			realp = pointer;
			v = pointer;
		}

		void SetDataCopy(T *data) {
			memcpy(realp, data, sizeof(T)* n);
		}


		// set data to a vector - at your risk - DONOT FREE PREVIOUS DATA!!
		void SetData(T *pointer, int nelements) {
			realp = pointer;
			v = pointer;
			n = nelements;
		}


		// Concatenate data of b in this vector
		void Concatenate(Vector<T> *b) {
			T *tmp;

			tmp = (T *)malloc(sizeof(T)* (this->n + b->GetNElements()));
			memcpy(tmp, realp, sizeof(T)* this->n);
			memcpy(tmp + this->n, b->GetData(), sizeof(T)* b->GetNElements());
			delete realp;
			this->n = this->n + b->GetNElements();
			this->realp = tmp;
			this->v = tmp;
		}


		// set value
		inline void SetElement(int x, T value) { v[x] = value; }

		// set all value
		void SetAllElements(T value) {
			int i;

			for (i = 0; i < n; i++)
				SetElement(i, value);
		}

		// set n elements from p with value 
		void SetRangeElements(int p, int n, T value) {
			int i;

			for (i = p; i < n + p; i++)
				SetElement(i, value);
		}

		// retrieve matrix pointer
		inline T *GetData() { return v; }

		// retrive number of elements
		inline int GetNElements() { return n; }

		// reset elements of vector
		inline void ResetVector() { memset(v, 0, sizeof(T)* n); }

		// select elements from the vector
		inline Vector<T> *SelectElements(vector<int> *selectedElems) {
			Vector<T> *sv;
			int i;

			sv = new Vector<T>((int)selectedElems->size());
			for (i = 0; i < sv->GetNElements(); i++) {
				sv->SetElement(i, v[selectedElems->at(i)]);
			}

			return sv;
		}

		// print vector
		void Print() {
			int i;

			for (i = 0; i < this->n; i++) {
				printf("%5.3f ", GetElement(i));
			}
			printf("\n");
		}

		// write vector
		void Write(char *filename) {
			int i;
			FILE *f;
			T value;

			f = fopen(filename, "wb");
			// header
			i = 1;
			fwrite(&i, 1, sizeof(int), f);
			fwrite(&n, 1, sizeof(int), f);

			// data
			for (i = 0; i < n; i++) {
				value = GetElement(i);
				fwrite(&value, 1, sizeof(T), f);
			}

			fclose(f);
		}

		// copy this vector to a new one
		Vector *Copy() {
			Vector<T> *v1;
			T *data;
			v1 = new Vector<T>(n);
			data = v1->GetData();
			memcpy(data, this->GetData(), sizeof(T)* n);
			return v1;
		}


		// resize vector
		void Resize(int n) {
			v = realloc(v, sizeof(T)* n);
			realp = v;
			this->n = n;
		}


		// compute mean
		T Mean() {
			int i;
			T value = 0;

			for (i = 0; i < this->n; i++) {
				value += this->GetElement(i);
			}

			return value / (float) this->n;
		}

		// compute var
		T Std() {
			T mean;
			T sum2, sumc;
			T var;
			int i;

			mean = this->Mean();
			sum2 = 0;
			sumc = 0;
			for (i = 0; i < this->n; i++) {
				sum2 += (this->GetElement(i) - mean) * (this->GetElement(i) - mean);
				sumc += (this->GetElement(i) - mean);
			}

			var = (sum2 - ((sumc * sumc) / (float)this->n)) / ((float)this->n - 1);
			return sqrt(var);
		}


		// compute L2 norm of a vector
		T L2Norm() {
			int i;
			T norm = 0;

			for (i = 0; i < this->n; i++) {
				norm += this->GetElement(i) * this->GetElement(i);
			}
			return sqrt(norm);
		}


		// add all values of the vector
		T Sum() {
			int i;
			T svalue = 0;

			for (i = 0; i < this->n; i++) {
				svalue += this->GetElement(i);
			}

			return svalue;
		}

		void Write(string filename) { Write((char *)filename.c_str()); }
	};








	//////////////////
	// class Matrix //
	//////////////////
	template <class T>
	class Matrix { // matrix is column based
		T *realp; // real pointer before align
		T *m;     // data pointer
		int r, c;     // number of rows and columns
		int idx;      // index to get elements sequentially

		// allocate memory for the matrix
		inline void AllocateMemory(int r, int c) {

			CV_FUNCNAME("Matrix::AllocateMemory");

			realp = new T[r * c];

			if (realp == NULL)
				DET_ERROR("No available memory to allocate matrix");


			//	realp = (T *) malloc(sizeof(T) * r * c);
			memset(realp, 0, sizeof(T)* (r * c));
			m = realp;
			this->r = r;
			this->c = c;
			this->idx = 0;
		}

		void LoadMatrix(char *filename) {
			FILE *f;
			int x, y, r, c;
			T var;

			f = fopen(filename, "rb");
			fread(&r, sizeof(int), 1, f);
			fread(&c, sizeof(int), 1, f);
			AllocateMemory(r, c);

			// assuming file stores rowise data
			for (y = 0; y < r; y++) {
				for (x = 0; x < c; x++) {
					fread(&var, sizeof(T), 1, f);
					SetValue(x, y, var);
				}
			}

			fclose(f);
		}


		// concatenate two matrices in rows
		void LoadMatrix(char *filename, char *filename2) {
			FILE *f, *f2;
			int x, y, r, c, r2, c2;
			T var;

			f = fopen(filename, "rb");
			fread(&r, sizeof(int), 1, f);
			fread(&c, sizeof(int), 1, f);

			f2 = fopen(filename2, "rb");
			fread(&r2, sizeof(int), 1, f2);
			fread(&c2, sizeof(int), 1, f2);

			// check dimensions
			if (c2 != c) {
				printf("Error: matrices need to have same number of columns!\n");
				exit(2);
			}

			AllocateMemory(r + r2, c);

			// assuming file stores rowise data
			for (y = 0; y < r; y++) {
				for (x = 0; x < c; x++) {
					fread(&var, sizeof(T), 1, f);
					SetValue(x, y, var);
				}
			}
			fclose(f);

			for (y = r; y < r + r2; y++) {
				for (x = 0; x < c; x++) {
					fread(&var, sizeof(T), 1, f2);
					SetValue(x, y, var);
				}
			}


			fclose(f2);
		}


	public:

		// clear structures
		~Matrix() {
			//printf("cleaning matrix\n");
			delete this->realp;
			//fflush(stdout);
			this->realp = NULL;
			this->m = NULL;
			this->r = -1;
			this->c = -1;
		}

		// Initialize
		Matrix(int r, int c) { AllocateMemory(r, c); }

		Matrix(string filename) { LoadMatrix((char *)filename.c_str()); }

		// load a matrix
		Matrix(char *filename) { LoadMatrix(filename); }

		// load matrix from a stream
		//Matrix(void *stream, int nbytes) { LoadMatrix(stream, nbytes); }

		// load a matrices and concatenate them
		Matrix(char *filename, char *filename2) { LoadMatrix(filename, filename2); }
		Matrix(string filename, string filename2) { LoadMatrix((char *)filename.c_str(), (char *)filename2.c_str()); }

		// return number of bytes used for this matrix
		int GetNumBytes() { return r * c * sizeof(T); }

		// retrieve element at row y and column x
		inline T GetElement(int x, int y) { return m[x * r + y]; }

		// set value
		inline void SetValue(int x, int y, T value) { m[x * r + y] = value; }

		// retrieve matrix pointer
		inline T *GetData() { return m; }

		// retrieve pointer for a column (column starts on 0)
		inline T *GetColumn(int i) { return m + r*i; }

		// retrieve a row of the matrix in a new vector
		inline Vector<T> *GetRow(int i) {
			Vector<T> *v;
			int j;

			v = new Vector<T>(this->c);
			for (j = 0; j < this->c; j++) {
				v->SetElement(j, this->GetElement(j, i));
			}

			return v;
		}

		// set position to get elements sequentially
		inline void SetPosition(int x, int y) { idx = x * r + y; }

		// set all value
		void SetAllElements(T value) {
			int y, x;

			for (x = 0; x < c; x++)
			for (y = 0; y < r; y++)
				SetValue(x, y, value);
		}


		// abs all values in the matrix
		void AbsAllElements() {
			int nelems;
			int i;

			nelems = r * c;
			for (i = 0; i < nelems; i++) {
				m[i] = abs(m[i]);
			}
		}


		// get current element, increment postion
		inline float GetCurrentElement() { return m[idx++]; }

		// set value at current position, increment it
		inline void SetCurrentElement(T value) { m[idx++] = value; }

		// set all values to 0
		inline void ZeroAll() { memset(m, 0, sizeof(T)* (r * c)); }

		// get number of rows
		inline int GetNRows() { return r; }

		// get number of columns
		inline int GetNCols() { return c; }

		Matrix<T> *GetSelectedRows(Vector<int> *selectedrows) {
			Matrix<T> *m;
			int i, j, k;

			m = new Matrix<T>(selectedrows->GetNElements(), this->GetNCols());
			for (i = 0; i < selectedrows->GetNElements(); i++) {
				j = selectedrows->GetElement(i);
				for (k = 0; k < this->GetNCols(); k++) {
					m->SetValue(k, i, this->GetElement(k, j));
				}
			}

			return m;
		}


		Matrix<T> *GetSelectedRows(vector<int> *selectedrows) {
			Matrix<T> *m;
			int i, j, k;

			m = new Matrix<T>((int)selectedrows->size(), this->GetNCols());
			for (i = 0; i < (int)selectedrows->size(); i++) {
				j = selectedrows->at(i);
				for (k = 0; k < this->GetNCols(); k++) {
					m->SetValue(k, i, this->GetElement(k, j));
				}
			}

			return m;
		}


		Matrix<T> *GetSelectedCols(vector<int> *selectedcols) {
			Matrix<T> *m;
			int i, j, k;

			m = new Matrix<T>(this->GetNRows(), (int)selectedcols->size());
			for (i = 0; i < (int)selectedcols->size(); i++) {
				j = selectedcols->at(i);
				for (k = 0; k < this->GetNRows(); k++) {
					m->SetValue(i, k, this->GetElement(j, k));
				}
			}

			return m;
		}

		// select columns and row of the matrix
		Matrix<T> *GetSelectedColsRows(vector<int> &selectedcols, vector<int> &selectedrows) {
			Matrix<T> *m, *maux;

			maux = this->GetSelectedRows(&selectedrows);
			m = maux->GetSelectedCols(&selectedcols);
			delete maux;

			return m;
		}

		// select columns in the range [lowIdx, highIdx] (closed interval!)
		Matrix<T> *GetSelectedCols(int lowIdx, int highIdx) {
			vector<int> selectedcols;
			int i;

			for (i = lowIdx; i <= highIdx; i++) {
				selectedcols.push_back(i);
			}

			return GetSelectedCols(&selectedcols);
		}




		// contatenate rows of matrices, first m1 then m2
		Matrix<T> *ConcatenateMatricesRows(Matrix<T> *m1, Matrix<T> *m2) {
			Matrix<T> *m;
			int x, y, idxrow;

			// only copy the second matrix
			if (m1 == NULL) {
				m = m2->Copy();
				return m;
			}

			// test compatibility
			if (m1->GetNCols() != m2->GetNCols()) {
				printf("Incompatible number of columns!\n");
				exit(2);
			}

			m = new Matrix<T>(m1->GetNRows() + m2->GetNRows(), m1->GetNCols());

			idxrow = 0;
			// m1
			for (y = 0; y < m1->GetNRows(); y++) {
				for (x = 0; x < m1->GetNCols(); x++) {
					m->SetValue(x, idxrow, m1->GetElement(x, y));
				}
				idxrow++;
			}
			// m2
			for (y = 0; y < m2->GetNRows(); y++) {
				for (x = 0; x < m2->GetNCols(); x++) {
					m->SetValue(x, idxrow, m2->GetElement(x, y));
				}
				idxrow++;
			}

			return m;
		}


		// contatenate cols of matrices, first m1 then m2
		Matrix<T> *ConcatenateMatricesCols(Matrix<T> *m1, Matrix<T> *m2) {
			Matrix<T> *m;
			int x, y, idxcol;

			// only copy the second matrix
			if (m1 == NULL) {
				m = m2->Copy();
				return m;
			}

			// test compatibility
			if (m1->GetNRows() != m2->GetNRows()) {
				printf("Incompatible number of rows!\n");
				exit(2);
			}

			m = new Matrix<T>(m1->GetNRows(), m1->GetNCols() + m2->GetNCols());

			idxcol = 0;
			// m1
			for (x = 0; x < m1->GetNCols(); x++) {
				for (y = 0; y < m1->GetNRows(); y++) {
					m->SetValue(idxcol, y, m1->GetElement(x, y));
				}
				idxcol++;
			}
			// m2
			for (x = 0; x < m2->GetNCols(); x++) {
				for (y = 0; y < m2->GetNRows(); y++) {
					m->SetValue(idxcol, y, m2->GetElement(x, y));
				}
				idxcol++;
			}

			return m;
		}


		// sum elements
		inline T Sum() {
			int i, j;
			T sum = 0;
			for (i = 0; i < this->r; i++) {
				for (j = 0; j < this->c; j++) {
					sum += GetElement(j, i);
				}
			}
			return sum;
		}



		// print matrix
		void Print() {
			int i, j;

			for (i = 0; i < this->r; i++) {
				for (j = 0; j < this->c; j++) {
					printf("%5.3f ", GetElement(j, i));
				}
				printf("\n");
			}
		}

		// write to a file - write row-wise
		void Write(char *filename) {
			int i, j;
			FILE *f;
			float value;

			f = fopen(filename, "wb");
			if (f == NULL)  {
				printf("Couldn't open '%s' for writing\n", filename);
				exit(2);
			}

			// header
			fwrite(&r, 1, sizeof(int), f);
			fwrite(&c, 1, sizeof(int), f);

			// data
			for (i = 0; i < r; i++) {
				for (j = 0; j < c; j++) {
					value = (float)GetElement(j, i);
					fwrite(&value, 1, sizeof(float), f);
				}
			}

			fclose(f);
		}



		// write to a file - write row-wise
		void Write(string filename, int nrows) {
			int i, j;
			FILE *f;
			float value;

			f = fopen((char *)filename.c_str(), "wb");
			// header
			fwrite(&nrows, 1, sizeof(int), f);
			fwrite(&c, 1, sizeof(int), f);

			// data
			for (i = 0; i < nrows; i++) {
				for (j = 0; j < c; j++) {
					value = (float)GetElement(j, i);
					fwrite(&value, 1, sizeof(float), f);
				}
			}

			fclose(f);
		}

		void Write(string filename) { Write((char *)filename.c_str()); }

		// copy this matrix to a new one
		Matrix *Copy() {
			Matrix *m1;
			T *data;
			m1 = new Matrix<T>(r, c);
			data = m1->GetData();
			memcpy(data, this->GetData(), sizeof(T)* r * c);
			return m1;
		}

		Vector<T> *RowsMean() {
			int i, j;
			Vector<T> *v;

			v = new Vector<T>(c);
			v->SetAllElements((T)0);

			for (i = 0; i < r; i++) {
				for (j = 0; j < c; j++) {
					v->SetElement(j, v->GetElement(j) + this->GetElement(j, i));
				}
			}

			for (j = 0; j < c; j++) {
				v->SetElement(j, v->GetElement(j) / (T)r);
			}

			return v;
		}


		void SetRow(Vector<T> *InVector, int r) {
			int i;
			// check sizes
			if (InVector->GetNElements() != this->GetNCols()) {
				printf("error diff n columns\n");
				exit(2);
			}

			for (i = 0; i < InVector->GetNElements(); i++) {
				this->SetValue(i, r, InVector->GetElement(i));
			}
		}

		// set c-th column with data
		void SetColumn(T *data, int c) {
			T *colptr;

			colptr = this->GetColumn(c);
			memcpy(colptr, data, this->GetNRows() * sizeof(T));
		}


		// set selected columns in the range [lowIdx, highIdx] (closed interval!)
		void SetSelectedCols(Matrix<T> *mdata, int lowIdx, int highIdx) {
			vector<int> selectedcols;
			int i, idx;

			CV_FUNCNAME("Matrix::SetSelectedCols");

			if ((mdata->GetNCols() != highIdx - lowIdx + 1) || (mdata->GetNRows() != this->GetNRows()))
				DET_ERROR("Invalid number of columns/rows");

			for (idx = 0, i = lowIdx; i <= highIdx; i++) {
				this->SetColumn(mdata->GetColumn(idx++), i);
			}
		}


		// retrieve a new matrix with selected rows
		Matrix<T> * SetSelectedRows(vector<int> &rowIdx) {
			Vector<T> *v;
			Matrix<T> *m;
			int i;

			// new matrix
			m = new Matrix<T>((int)rowIdx.size(), this->GetNCols());

			for (i = 0; i < (int)rowIdx.size(); i++) {
				v = this->GetRow(rowIdx[i]);
				m->SetRow(v, i);
				delete v;
			}

			return m;
		}
	};







	class Maths {

	public:
		// project vector v according to matrix m and result in res
		//void Project(Vector *v, Matrix *m, Vector *res);
		//void Project2(Vector *v, Matrix *m, Vector *res);
		float DotProduct(float *v, float *v2, int n);
		float StepDotProduct(float *v, float *v2, int r0, int r1);

		float DotProductSSEMultof4(float *v, float *v2, int n);

		// v1 ./ value 
		void DivideVectorsSSE(float *v1, float value, float *outputvect, int n);

		// dot product when n is not multiple of 4
		float DotProductSSENotMultof4(float *v, float *v2, int n);
		float dot_product_sse(float * vec1, float * vec2, int N);
		//Vector *Project(Vector *v, Matrix *m);

		void ZscoreSSE(float *data, float *mean, float *std, float *outputvect, int n);

		// v1 - v2 (WARNING: this operation add values outside the range (assume extra space was allocated), when n not multiple of 4)
		void SubtractVectorsSSE(float *v1, float *v2, float *outputvect, int n);



		// not aligned operations
		void SubtractVectorsSSENotAligned(float *v1, float *v2, float *outputvect, int n);
		void SubtractVectorsSSENotAligned(short *v1, short *v2, short *outputvect, int n);
		void AddVectorsNotAligned(float *v1, float *v2, float *outputvect, int n);


		inline void AddVectorsNotAligned(short *v1, short *v2, short *outputvect, int n);

		// keep the max between values of a vector and a constant
		// v1 = max(v1, value)
		void KeepMaxSSE(float *v1, float value, float *outputvect, int n);


		float sse3_inner(const float* a, const float* b, unsigned int size);
		// add vector values of a vector using SSE, needs a auxiliar vector 


		// not aligned but n multiple of 4
		void SubtractVectorsSSENotAlignedMult4(float *v1, float *v2, float *outputvect, int n);
		void AddVectorsNotAlignedMult4(float *v1, float *v2, float *outputvect, int n);

		// operations to find the correct bin for co-occurrence matrix
		void FindCoocBins(unsigned char *v1, unsigned char *v2, unsigned char *ret1, unsigned char *ret2, int n);

		void Test();
		void TestFunctions();
	};

}

#endif