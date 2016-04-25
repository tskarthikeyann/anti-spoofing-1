#include <iostream>  
#include "highgui.h"
#include "opencv2/opencv.hpp"  
#include"string"
#include <stdio.h>
#include "headers.h"
#include "HSC.h"

//#include"FreeImage.h"

using namespace cv;
using namespace std;
using namespace HSCNameSpace;

#define UNKNOWN_FLOW_THRESH 1e9 
const char* face_cascade_name = "d:\\opencv\\sources\\data\\haarcascades\\haarcascade_frontalface_alt2.xml";
CascadeClassifier face_cascade;
string window_name = "人脸识别";

//todo for testing
int frameNIndex = 0;

void resizePic(Mat & image, CvSize size){

	IplImage faceImage = IplImage(image);
	IplImage * resizedImage;
	resizedImage = cvCreateImage(size, IPL_DEPTH_8U, 3);
	cvResize(&faceImage, resizedImage, CV_INTER_LINEAR);
	image = Mat(resizedImage);
}

void resizePic_Gray(Mat & image, CvSize size){

	IplImage faceImage = IplImage(image);
	IplImage * resizedImage;
	resizedImage = cvCreateImage(size, IPL_DEPTH_8U, 1);
	cvResize(&faceImage, resizedImage, CV_INTER_LINEAR);
	image = Mat(resizedImage);
}

bool detectFace(Mat frame, Mat & dectedFace, CvSize size){
	std::vector<Rect>facesRect;
	std::vector<Mat> facesMat;
	bool isFaceGet = 0;

	Mat frame_gray;

	//转化为灰度图
	cvtColor(frame, frame_gray, CV_BGR2GRAY);

	//均衡化
	equalizeHist(frame_gray, frame_gray);

	//
	face_cascade.detectMultiScale(frame_gray, facesRect, 1.1, 2, 0 | CV_HAAR_SCALE_IMAGE, Size(60, 60), Size(600, 600));

	//todo:test
	imwrite("graytest.jpg", frame_gray);

	if (facesRect.size() > 0){

		//todo:extract the bigest face
		// get face region
		dectedFace = frame(facesRect[0]);

		//resizeImage
		IplImage faceImage = IplImage(dectedFace);
		IplImage * resizedImage;
		resizedImage = cvCreateImage(size, IPL_DEPTH_8U, 3);
		cvResize(&faceImage, resizedImage, CV_INTER_LINEAR);
		dectedFace = Mat(resizedImage);
		isFaceGet = 1;

		//todo for testing
		imwrite("detected.jpg", dectedFace);
	}

	return isFaceGet;

}

void makecolorwheel(vector<Scalar> &colorwheel)
{
	int RY = 15;
	int YG = 6;
	int GC = 4;
	int CB = 11;
	int BM = 13;
	int MR = 6;

	int i;

	for (i = 0; i < RY; i++) colorwheel.push_back(Scalar(255, 255 * i / RY, 0));
	for (i = 0; i < YG; i++) colorwheel.push_back(Scalar(255 - 255 * i / YG, 255, 0));
	for (i = 0; i < GC; i++) colorwheel.push_back(Scalar(0, 255, 255 * i / GC));
	for (i = 0; i < CB; i++) colorwheel.push_back(Scalar(0, 255 - 255 * i / CB, 255));
	for (i = 0; i < BM; i++) colorwheel.push_back(Scalar(255 * i / BM, 0, 255));
	for (i = 0; i < MR; i++) colorwheel.push_back(Scalar(255, 0, 255 - 255 * i / MR));
}

void motionToColor(Mat flow, Mat &color)
{
	if (color.empty())
		color.create(flow.rows, flow.cols, CV_8UC3);

	static vector<Scalar> colorwheel; //Scalar r,g,b  
	if (colorwheel.empty())
		makecolorwheel(colorwheel);

	// determine motion range:  
	float maxrad = -1;

	// Find max flow to normalize fx and fy  
	for (int i = 0; i < flow.rows; ++i)
	{
		for (int j = 0; j < flow.cols; ++j)
		{
			Vec2f flow_at_point = flow.at<Vec2f>(i, j);
			float fx = flow_at_point[0];
			float fy = flow_at_point[1];
			if ((fabs(fx) >  UNKNOWN_FLOW_THRESH) || (fabs(fy) >  UNKNOWN_FLOW_THRESH))
				continue;
			float rad = sqrt(fx * fx + fy * fy);
			maxrad = maxrad > rad ? maxrad : rad;
		}
	}

	for (int i = 0; i < flow.rows; ++i)
	{
		for (int j = 0; j < flow.cols; ++j)
		{
			uchar *data = color.data + color.step[0] * i + color.step[1] * j;
			Vec2f flow_at_point = flow.at<Vec2f>(i, j);

			//normalize the distance
			float fx = flow_at_point[0] / maxrad;
			float fy = flow_at_point[1] / maxrad;
			if ((fabs(fx) >  UNKNOWN_FLOW_THRESH) || (fabs(fy) >  UNKNOWN_FLOW_THRESH))
			{
				data[0] = data[1] = data[2] = 0;
				continue;
			}
			float rad = sqrt(fx * fx + fy * fy);

			float angle = atan2(-fy, -fx) / CV_PI;
			float fk = (angle + 1.0) / 2.0 * (colorwheel.size() - 1);
			int k0 = (int)fk;
			int k1 = (k0 + 1) % colorwheel.size();
			float f = fk - k0;
			//f = 0; // uncomment to see original color wheel  

			for (int b = 0; b < 3; b++)
			{
				float col0 = colorwheel[k0][b] / 255.0;
				float col1 = colorwheel[k1][b] / 255.0;
				float col = (1 - f) * col0 + f * col1;
				if (rad <= 1)
					col = 1 - rad * (1 - col); // increase saturation with radius  
				else
					col *= .75; // out of range  
				data[2 - b] = (int)(255.0 * col);
			}
		}
	}
}

void getOFBetween2Frames(Mat frame1, Mat frame2, Mat & outColoredMotion){

	Mat grayed_frame1, grayed_frame2, motion;
	try{
		cvtColor(frame1, grayed_frame1, CV_BGR2GRAY);
		if (frame1.empty()) return;

		cvtColor(frame2, grayed_frame2, CV_BGR2GRAY);
		if (frame2.empty()) return;

		calcOpticalFlowFarneback(grayed_frame1, grayed_frame2, motion, 0.5, 3, 15, 3, 5, 1.2, 0);
		motionToColor(motion, outColoredMotion);

	}
	catch (const char * ex){

		cout << ex << endl;
	}

}

void getAverageOFM(vector<Mat>frames, Mat & averageFrame){
	try{
		vector<Mat> opticalFlowFrames;

		for (int i = 0; i < frames.size() - 1; i++){

			//get optical flow map between two face image
			Mat opticalFrame;
			getOFBetween2Frames(frames[i], frames[i + 1], opticalFrame);
			Mat gray_OpticalFrame = cvCreateImage(opticalFrame.size(), IPL_DEPTH_8U, 1);
			IplImage * image_face;
			cvtColor(opticalFrame, gray_OpticalFrame, CV_RGB2GRAY);
			opticalFlowFrames.push_back(gray_OpticalFrame);


			//todo:for testing
			char opticalFileName[40];
			sprintf(opticalFileName, "./opticalFrame/Frame%dface%d.jpg", frameNIndex, i);
			imwrite(opticalFileName, gray_OpticalFrame);
			frameNIndex++;
		}

		averageFrame = cvCreateImage(opticalFlowFrames[0].size(), IPL_DEPTH_8U, 1);
		//averageFrame.create(opticalFlowFrames[0].size(), CV_8UC1,Scalar(0));
		averageFrame = averageFrame * 0;

		for (int i = 0; i < opticalFlowFrames.size(); i++){

			averageFrame = averageFrame + opticalFlowFrames[i] / opticalFlowFrames.size();
		}
	}
	catch (const char * ex){

		cout << ex << endl;
	}
}

HSCNameSpace::Vector<float> * getShearletFeatures(Mat faceImg)
{
	HSCNameSpace::Vector<float> *feat;
	HSC hsc;
	IplImage *img = &IplImage(faceImg);

	/* set parameters */
	hsc.SetParameter("LapFilter", "Burt");			// filter for pyramid decomposition
	hsc.SetParameter("norient", "8");				// number of orientations
	hsc.SetParameter("maskSize", "4");				// convolution mask size
	hsc.SetParameter("nlevels", "4");				// number of decomposition levels
	hsc.SetParameter("useCell", "false");			// split block into four cells 
	hsc.SetParameter("normalization", "global");	// normalization type

	/* initialize feature extraction method */
	hsc.InitializeExtractionMethod();

	/* feature extraction */			// load image
	hsc.AddNewImage(img);	
	// add image to extraction
	feat = hsc.ExtractFeatures(64, 64, 64, 64);	// extract feature descriptors considering block sizes of 256x256 pixels

	return feat;
}

void initImageData(const char * videoPath, vector<Mat> & faceFrames, vector<Mat>&sceneFrames, Mat & shearletFace)
{
	try{
		// read frame from pic file
		//const char* path1 = "Q:/CSharpProject/OpticalFlow/OpticalFlow/car1.jpg";
		//const char* path2 = "Q:/CSharpProject/OpticalFlow/OpticalFlow/car2.jpg";
		//Mat img1 = imread(path1, CV_LOAD_IMAGE_UNCHANGED);
		//if (img1.empty()) return -1;

		//Mat img2 = imread(path2, CV_LOAD_IMAGE_UNCHANGED);
		//if (img2.empty()) return -1;

		int framesToGet = 6;//get 6 frame of faces and scenes
		//load face detector
		face_cascade.load(face_cascade_name);

		// read frame form video
		CvCapture* cap = cvCaptureFromFile(videoPath);//

		//CvCapture* cap = cvCreateCameraCapture(0);
		int numFrames = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FRAME_COUNT);
		int fps = (int)cvGetCaptureProperty(cap, CV_CAP_PROP_FPS);

		//extract 4 frames per second
		int duration_frame = 5;//take 5 as a fixed value
		int duration_millionSeconds = (1000 / (float)fps)*duration_frame;
		//vector<Mat> faceFrames;
		//vector<Mat> sceneFrames;
		//Mat shearletFace;

		int faceCount = 0;
		int frameIndex = 0;
		Mat frame;
		bool isShearletFaceGet = 0;
		while (faceCount < framesToGet && frameIndex < numFrames){

			/*int pos = duration_millionSeconds*frameIndex;
			cvSetCaptureProperty(cap, CV_CAP_PROP_POS_MSEC, pos);*/
			frame = cvQueryFrame(cap);

			//todo for testing
			char extractFrame[40];
			sprintf(extractFrame, "./extractedframe/extracted%d.jpg", frameIndex);
			imwrite(extractFrame, frame);

			Mat detectedFaceFrame;

			if (detectFace(frame, detectedFaceFrame, cvSize(256, 256))){

				faceCount++;

				//get resized face pic
				faceFrames.push_back(detectedFaceFrame);

				// resize sense pic
				resizePic(frame, cvSize(256, 256));
				sceneFrames.push_back(frame);

				//Get face pic for shearlet
				if (!isShearletFaceGet)
				{
					shearletFace = detectedFaceFrame.clone();
					resizePic(shearletFace, cvSize(256, 256));
					isShearletFaceGet = 1;
				}


				//todo:for testing
				char faceFileName[40];
				char sceneFileName[40];
				sprintf(faceFileName, "./face/face%d.jpg", faceCount);
				sprintf(sceneFileName, "./scene/scene%d.jpg", faceCount);
				imwrite(faceFileName, detectedFaceFrame);
				imwrite(sceneFileName, frame);
			}

			//skip frames duration
			for (int i = 0; i < duration_frame; i++){

				cvQueryFrame(cap);
				frameIndex++;
			}

		}

	}
	catch (char * ex){

		cout << "Init error" << ex << endl;

	}
}

int main(int, char**)
{
	try{
		//init to get image data
		vector<Mat> faceFrames;
		vector<Mat> sceneFrames;
		Mat shearletFace;

		initImageData("../Video/2fliped.wmv", faceFrames, sceneFrames, shearletFace);

		// get face average OFM
		Mat averageOFM_Face;
		getAverageOFM(faceFrames, averageOFM_Face);
		resizePic_Gray(averageOFM_Face, cvSize(32, 32));

		// get scene average OFM
		Mat averageOFM_Scene;
		getAverageOFM(sceneFrames, averageOFM_Scene);
		resizePic_Gray(averageOFM_Scene, cvSize(32, 32));

		//get shearlet features
		HSCNameSpace::Vector<float>shearletFeature = getShearletFeatures(shearletFace);


		//todo:for testing
		char faceOpticalAverageFileName[40];
		char sceneOpticalAverageFileName[40];
		sprintf(faceOpticalAverageFileName, "./opticalFrame/averageface%d.jpg", 1);
		sprintf(sceneOpticalAverageFileName, "./opticalFrame/averagescene%d.jpg", 1);
		imwrite(faceOpticalAverageFileName, averageOFM_Face);
		imwrite(sceneOpticalAverageFileName, averageOFM_Scene);

	}
	catch (const char * ex){

		cout << ex << endl;
	}
}