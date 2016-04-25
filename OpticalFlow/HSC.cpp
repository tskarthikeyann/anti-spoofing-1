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
#include "misc.h"
#include "HSCmatrices.h"
#include "HSC.h"
namespace HSCNameSpace{

	HSC::HSC() {

		norient = 8;
		maskSize = 8;
		LapFilter = "Burt";
		normType = "global";
		nlevels = 2;
		g = NULL;
		h = NULL;
		useCells = true;
	}


	void HSC::ResetStructures() {
		vector<Matrix<float> *> d;
		int s, i;

		for (s = 0; s < (int)dFeatures.size(); s++) {
			d = dFeatures.at(s);

			for (i = 0; i < (int)d.size(); i++) {
				delete d.at(i);
			}
		}

		dFeatures.clear();
	}



	int HSC::GetNFeatures() {

		if (useCells == true)
			return 4 * nlevels * norient;
		else
			return nlevels * norient;
	}


	float HSC::SumRegion(Matrix<float> *m, int x0, int y0, int x1, int y1) {
		float value = 0;
		int x, y;

		for (y = y0; y <= y1; y++) {
			for (x = x0; x <= x1; x++) {
				value += fabs(m->GetElement(x, y));
			}
		}

		return value;
	}


	void HSC::InitializeExtractionMethod() {
		char str[1024];

		// load filters
		h = LapPyrpfiltersH(LapFilter);
		g = LapPyrpfiltersG(LapFilter);

		sprintf(str, "M_%d_%dx%d", norient, maskSize, maskSize);
		matrixShear = GetHSCMatrix(str);

		// debug parameters
		printf("[HSC] [nlevels: %d] [norient: %d] [Filter: %s] [norm: %s] [Cells: %d]\n", nlevels, norient,
			LapFilter.c_str(), normType.c_str(), useCells);
	}



	void HSC::AddNewImage(IplImage *img) {
		vector<Matrix<float> *> d;
		Matrix<float> *m, *xlo = NULL;
		IplImage *dst;
		int i;

		// reset structures first
		ResetStructures();

		dst = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);

		// convert to a gray scale image
		cvCvtColor(img, dst, CV_RGB2GRAY);

		//todo:for update
		//cvEqualizeHist(dst, dst);

		// compute features for the entire image
		m = CropImgMatrix(dst, 0, 0, img->width - 1, img->height - 1);

		//todo:for testing,image,m is a gray pic and 
		cv::Mat mattemp = cv::Mat(m->GetNRows(), m->GetNCols(), CV_32FC1, m->GetData());
		cv::imwrite("firstFeature.jpg", mattemp);

		// release image
		cvReleaseImage(&dst);

		// apply shearlet transform
		for (i = 0; i < nlevels; i++) {

			//m:image norient:orientation,xlo:different levels ,d is convolutioned feature,h,g:lap filter
			Shearletsapply(m, norient, &xlo, d, matrixShear, g, h, i);


			//todo: for testing 
			char filename[40];
			sprintf(filename, "Scale%d.jpg", i);
			cv::Mat mattemp = cv::Mat(xlo->GetNCols(), xlo->GetNRows(), CV_32FC1, xlo->GetData());
			cv::imwrite(filename, mattemp);


			//todo:dFeatures stores for features in 8 direction?
			// store results to be used later
			dFeatures.push_back(d);


			// clear structures and set for the next level
			delete m;
			d.clear();

			// set downsampled image
			m = xlo;
		}

		// release image matrix
		delete xlo;

	}




	void HSC::ExtractFeatures(Position *pos, float *feat) {
		int i, s, x0, y0, hh, ww, x1, y1;
		vector<Matrix<float> *> d;
		Matrix<float> *mtmp;
		float value;
		float values[4];
		int idx = 0;
		float valueSum = 0;
		float L2;
		int nfeatLevel;	// number of features in one level
		int xhalf, yhalf;

		CV_FUNCNAME("HSC::ExtractFeatures");

		// define number of features ina single level
		if (useCells == true)
			nfeatLevel = 4 * norient;
		else
			nfeatLevel = norient;

		// set initial coordinates to extract features
		x0 = pos->x0;
		y0 = pos->y0;
		ww = pos->w;
		hh = pos->h;

		// retrieve computed features for each scale
		for (s = 0; s < (int)dFeatures.size(); s++) {

			// retrieve features
			d = dFeatures.at(s);

			// set block location
			x0 = x0 >> s;
			y0 = y0 >> s;
			hh = hh >> s;
			ww = ww >> s;
			x1 = x0 + ww - 1;
			y1 = y0 + hh - 1;

			// compute the histogram for each orientation
			for (i = 0; i < (int)d.size(); i++) {

				mtmp = d.at(i);

				if (useCells == false) {
					value = SumRegion(mtmp, x0, y0, x1, y1);
					feat[idx + i] = value;
					valueSum += value;
				}
				else {
					xhalf = x0 + ((x1 - x0 + 1) / 2);
					yhalf = y0 + ((y1 - y0 + 1) / 2);
					values[0] = SumRegion(mtmp, x0, y0, xhalf - 1, yhalf - 1);
					values[1] = SumRegion(mtmp, xhalf, y0, x1, yhalf - 1);
					values[2] = SumRegion(mtmp, x0, yhalf, xhalf - 1, y1);
					values[3] = SumRegion(mtmp, xhalf, yhalf, x1, y1);
					valueSum += values[0];

					// set them into the feature vector
					feat[idx + i] = values[0];
					feat[idx + i + norient] = values[1];
					feat[idx + i + (2 * norient)] = values[2];
					feat[idx + i + (3 * norient)] = values[3];
				}
			}
			idx += nfeatLevel;
		}


		if (valueSum > 0) {
			// L2
			if (normType.compare("local") == 0 || normType.compare("local+global") == 0) { // normalization for each level
				idx = 0;
				for (i = 0; i < nlevels; i++) {  // local normalization
					L2 = mat.DotProductSSENotMultof4(feat + idx, feat + idx, nfeatLevel);
					L2 = sqrt(L2 + (float)EPS);

					// compute L2-norm
					mat.DivideVectorsSSE(feat + idx, L2, feat + idx, nfeatLevel);
					idx += nfeatLevel;
				}
			}
			if (normType.compare("global") == 0 || normType.compare("local+global") == 0) {
				L2 = mat.DotProductSSENotMultof4(feat, feat, this->GetNFeatures());
				L2 = sqrt(L2 + (float)EPS);

				// compute L2-norm
				mat.DivideVectorsSSE(feat, L2, feat, this->GetNFeatures());
			}
			if (normType.compare("none") == 0) {
				;
			}
		}
	}



	Vector<float> *HSC::ExtractFeatures(int blockW, int blockH, int strideX, int strideY) {
		Vector<float> *features;
		vector<Position> blocks;
		int x, y, idx, w, h, i;
		Position pos;

		CV_FUNCNAME("HSC::ExtractFeatures");

		if (dFeatures.size() == 0)
			DET_ERROR("call AddNewImage() first!");

		// create list of blocks
		h = dFeatures.at(0).at(0)->GetNRows();
		w = dFeatures.at(0).at(0)->GetNCols();

		for (y = 0; y < h - blockH + 1; y += strideY) {
			for (x = 0; x < w - blockW + 1; x += strideX) {
				pos.x0 = x;
				pos.y0 = y;
				pos.y1 = y + blockH - 1;
				pos.x1 = x + blockW - 1;
				pos.h = blockH;
				pos.w = blockW;
				blocks.push_back(pos);
			}
		}

		// allocate feature vector
		features = new Vector<float>((int)blocks.size() * GetNFeatures());

		// extract features
		idx = 0;
		for (i = 0; i < (int)blocks.size(); i++) {
			pos = blocks.at(i);
			ExtractFeatures(&pos, features->GetData() + idx);
			idx += GetNFeatures();
		}

		return features;
	}




	Vector<float> *HSC::ExtractFeatures() {
		Vector<float> *features;
		vector<Position> blocks;
		int w, h;
		Position pos;

		CV_FUNCNAME("HSC::ExtractFeatures");

		if (dFeatures.size() == 0)
			DET_ERROR("call AddNewImage() first!");

		// allocate feature vector
		features = new Vector<float>(GetNFeatures());

		// set block size as the image size
		h = dFeatures.at(0).at(0)->GetNRows();
		w = dFeatures.at(0).at(0)->GetNCols();
		pos.x0 = 0;
		pos.y0 = 0;
		pos.x1 = w - 1;
		pos.y1 = h - 1;
		pos.h = h;
		pos.w = w;

		// extract features
		ExtractFeatures(&pos, features->GetData());

		return features;
	}




	void HSC::SetParameter(string name, string value) {

		CV_FUNCNAME("HSC::SetParameter");

		// check parameters
		if (name.compare("LapFilter") == 0) {
			LapFilter = value;
		}

		else if (name.compare("normalization") == 0) {
			normType = value;
		}

		else if (name.compare("nlevels") == 0) {
			nlevels = atoi(value.c_str());
		}

		else if (name.compare("norient") == 0) {
			norient = atoi(value.c_str());
		}

		else if (name.compare("maskSize") == 0) {
			maskSize = atoi(value.c_str());
		}

		else if (name.compare("useCell") == 0) {
			if (value.compare("true") == 0)
				useCells = true;
			else
				useCells = false;
		}

		else {
			printf("Parameter '%s' invalid\n", name.c_str());
			exit(2);
		}
	}
}