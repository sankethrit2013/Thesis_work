/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/

#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include <iostream>
#include <string>
#include <stdio.h>
// Software Guide : EndCodeSnippet

int main(int argc, char* argv[])
{
	// dicom series folder path
	//argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/dicom_ear_male/1.2.826.0.1.3417726.3.450309712/1.2.826.0.1.3417726.3.101337.20150827141831505";
	//argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/Pump_data_12_09_2015/Skull_base";
	//argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/Pump_data_12_09_2015/Skulll_120min";
	//argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/test_03_7_2016_new_mod";
	//argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/1_7_2016-2016-01-08/1_7_2016/first_base_line/1.2.826.0.1.3417726.3.1160966081/1.2.826.0.1.3417726.3.491044.20160107130027338";
	/*argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/1_7_2016-2016-01-08/1_7_2016/second_base_line/1.2.826.0.1.3417726.3.1161829028/1.2.826.0.1.3417726.3.71413.20160107134216713";*/
	//argv[1] = "c:/users/zx4403/dropbox/research/borkholder - micropump/images/1_7_2016-2016-01-08/1_7_2016/third_line_injection/1.2.826.0.1.3417726.3.1162784237/1.2.826.0.1.3417726.3.73357.20160107144940525";
	//argv[1] = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/1_7_2016-2016-01-08/1_7_2016/fourth_line_injection/1.2.826.0.1.3417726.3.1163653999/1.2.826.0.1.3417726.3.28521.20160107153507838";	
	//argv[1] = "N:/Zhenlin_Xu/Images/Calibrations/sample_1_18.98mg_ml/1.2.826.0.1.3417726.3.3442824753/1.2.826.0.1.3417726.3.342043.20160406162917645";

	//output file path
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20151209/base_ear_matched.nii";
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20151209/120min_ear_matched.nii";
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/first_base_line_ear_matched.nii";
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/second_base_line_ear.nii";
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/third_line_injection_ear.nii";
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160307/RightEar.nii.gz";
	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/Calibration/sample_1_18.98mg_ml.nii.gz";

	///*-----------------------------------------------
	//20160719 seriers
	//-------------------------------------------------*/
	//int numScan = 4; // number of scans

	//std::string *path_imgSeries = new std::string[4]; 
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////scan 1 
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/minutes0/1.2.826.0.1.3417726.3.8852.20160906131731375";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/minutes0_head.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/minutes0.nii";

	////baseline
	//path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/baseline512/512_base";
	//filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/baseline_head.nii";
	//filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/baseline.nii";

	////scan 2
	//path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/minutes30/1.2.826.0.1.3417726.3.516163.20160906134122750";
	//filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/minutes30_head.nii";
	//filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/minutes30.nii";

	////scan 3
	//path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/minutes60";
	//filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/minutes60_head.nii";
	//filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIbaseline/minutes60.nii";

	/*scan 1 
	path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/minutes0/1.2.826.0.1.3417726.3.8852.20160906131731375";
	filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/minutes0_head.nii";
	filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/minutes0.nii";

	baseline
	path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/baseline512/512_base";
	filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/baseline_head.nii";
	filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/baseline.nii";

	scan 2
	path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/minutes30/1.2.826.0.1.3417726.3.516163.20160906134122750";
	filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/minutes30_head.nii";
	filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/minutes30.nii";

	scan 3
	path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160906/minutes60";
	filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/minutes60_head.nii";
	filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROIraw/minutes60.nii";*/

	/*
	//scan 4
	path_imgSeries[4] = "N:/Zhenlin_Xu/Images/scans/20160719/scan_4/1.2.826.0.1.3417726.3.1681083285/1.2.826.0.1.3417726.3.41633.20160629144516250";
	filename_fullHead[4] = "N:/Zhenlin_Xu/Images/ImageVolume/20160719/scan4.nii";
	filename_ear[4] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan4_leftEar_raw.nii";
	*/

	///*-----------------------------------------------
	//20160629 series
	//-------------------------------------------------*/
	//int numScan = 5; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160629/scan1/1.2.826.0.1.3417726.3.821925.20160629120502718";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan1/scan1_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan1/scan1.nii";

	////scan 2
	//path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160629/scan2/1.2.826.0.1.3417726.3.182878.20160629123538140";
	//filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan2/scan2_raw.nii";
	//filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan2/scan2.nii";

	////scan 3
	//path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160629/scan3/1.2.826.0.1.3417726.3.962602.20160629130656281";
	//filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan3/scan3_raw.nii";
	//filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan3/scan3.nii";

	////scan 4
	//path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160629/scan4/1.2.826.0.1.3417726.3.107586.20160629135439468";
	//filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan4/scan4_raw.nii";
	//filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan4/scan4.nii";

	////scan 5
	//path_imgSeries[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20160629/scan5/1.2.826.0.1.3417726.3.41633.20160629144516250";
	//filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan5/scan5_raw.nii";
	//filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan5/scan5.nii";


	///*-------------------------------------------------*/
//
//	/*-----------------------------------------------
//	20170123 series
//	-------------------------------------------------*/
//int numScan = 6; // number of scans
//
//	std::string *path_imgSeries = new std::string[numScan];
//	std::string *filename_fullHead = new std::string[numScan];
//	std::string *filename_ear = new std::string[numScan];
//
//	//baseline / scan1
//	path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170123/scan1/1.2.826.0.1.3417726.3.659867.20170124131234625";
//	/*filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/scan1_raw.nii";
//	filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/scan1.nii";*/
//	filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/scan1_other_raw.nii";
//	filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/scan1_other.nii";
//
//	//scan 2
//	path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170123/scan2/1.2.826.0.1.3417726.3.520836.20170124131234515";
//	/*filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan2/scan2_raw.nii";
//	filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan2/scan2.nii";*/
//	filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan2/scan2_other_raw.nii";
//	filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan2/scan2_other.nii";
//
//	//scan 3
//	path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170123/scan3/1.2.826.0.1.3417726.3.317494.20170124134448078";
//	/*filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan3/scan3_raw.nii";
//	filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan3/scan3.nii";*/
//	filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan3/scan3_other_raw.nii";
//	filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan3/scan3_other.nii";
//
//	//scan 4
//	path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170123/scan4/1.2.826.0.1.3417726.3.292509.20170124141514390";
//	/*filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan4/scan4_raw.nii";
//	filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan4/scan4.nii";*/
//	filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan4/scan4_other_raw.nii";
//	filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan4/scan4_other.nii";
//
//	//scan 5
//	path_imgSeries[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170123/scan5/1.2.826.0.1.3417726.3.632562.20170124144500718";
//	/*filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan5/scan5_raw.nii";
//	filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan5/scan5.nii";*/
//	filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan5/scan5_other_raw.nii";
//	filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan5/scan5_other.nii";
//
//	//scan 6
//	path_imgSeries[5] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170123/scan6/1.2.826.0.1.3417726.3.432160.20170124163840046";
//	/*filename_fullHead[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan6/scan6_raw.nii";
//	filename_ear[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan6/scan6.nii";*/
//	filename_fullHead[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan6/scan6_other_raw.nii";
//	filename_ear[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan6/scan6_other.nii";


	/*-------------------------------------------------*/
	///*-----------------------------------------------
	//20170606 series
	//-------------------------------------------------*/
	//int numScan = 5; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170606/scan1/1.2.826.0.1.3417726.3.506539.20170606143844296";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/scan1_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/scan1.nii";

	////scan 2
	//path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170606/scan2/1.2.826.0.1.3417726.3.913720.20170606150754281";
	//filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan2/scan2_raw.nii";
	//filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan2/scan2.nii";

	////scan 3
	//path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170606/scan3/1.2.826.0.1.3417726.3.590497.20170606154235500";
	//filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan3/scan3_raw.nii";
	//filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan3/scan3.nii";

	////scan 4
	//path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170606/scan4/1.2.826.0.1.3417726.3.797731.20170606160634687";
	//filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan4/scan4_raw.nii";
	//filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan4/scan4.nii";

	////scan 5
	//path_imgSeries[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170606/scan5/1.2.826.0.1.3417726.3.649699.20170606164709625";
	//filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan5/scan5_raw.nii";
	//filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan5/scan5.nii";


	///*-------------------------------------------------*/
	///*-----------------------------------------------
	//20170612 series
	//-------------------------------------------------*/
	//int numScan = 5; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170612/scan1/1.2.826.0.1.3417726.3.338100.20170612112134031";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/scan1_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/scan1.nii";

	////scan 2
	//path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170612/scan2/1.2.826.0.1.3417726.3.734518.20170612114912656";
	//filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan2/scan2_raw.nii";
	//filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan2/scan2.nii";

	////scan 3
	//path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170612/scan3/1.2.826.0.1.3417726.3.841108.20170612122334828";
	//filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan3/scan3_raw.nii";
	//filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan3/scan3.nii";

	////scan 4
	//path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170612/scan4/1.2.826.0.1.3417726.3.452262.20170612125540125";
	//filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan4/scan4_raw.nii";
	//filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan4/scan4.nii";

	////scan 5
	//path_imgSeries[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170612/scan5/1.2.826.0.1.3417726.3.558429.20170612132911468";
	//filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan5/scan5_raw.nii";
	//filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan5/scan5.nii";


	///*-------------------------------------------------*/

    ///*-----------------------------------------------
	//20170619 series
	//-------------------------------------------------*/
	//int numScan = 6; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170619/scan1/1.2.826.0.1.3417726.3.445068.20170619135029265";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan1/scan1_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan1/scan1.nii";

	////scan 2
	//path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170619/scan2/1.2.826.0.1.3417726.3.5847.20170619142010343";
	//filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan2/scan2_raw.nii";
	//filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan2/scan2.nii";

	////scan 3
	//path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170619/scan3/1.2.826.0.1.3417726.3.191608.20170619150109531";
	//filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan3/scan3_raw.nii";
	//filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan3/scan3.nii";

	////scan 4
	//path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170619/scan4/1.2.826.0.1.3417726.3.34819.20170619153322250";
	//filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan4/scan4_raw.nii";
	//filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan4/scan4.nii";

	////scan 5
	//path_imgSeries[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170619/scan5/1.2.826.0.1.3417726.3.706253.20170619160239718";
	//filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan5/scan5_raw.nii";
	//filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan5/scan5.nii";

	////scan 6
	//path_imgSeries[5] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20170619/scan6/1.2.826.0.1.3417726.3.752116.20170619162436296";
	//filename_fullHead[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan6/scan6_raw.nii";
	//filename_ear[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/scan6/scan6.nii";

	///*-------------------------------------------------*/

	///*-----------------------------------------------
	//calibration isovue
	//-------------------------------------------------*/
	//int numScan = 1; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration_isovue/1.2.826.0.1.3417726.3.725346.20170706140851798";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/Isovue/Conc50_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/Isovue/Conc50.nii";

	///*-----------------------------------------------
	//calibration GnP
	//-------------------------------------------------*/
	//int numScan = 1; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration_gold/1.2.826.0.1.3417726.3.346254.20170706140851674";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/Gold/Conc100_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/Gold/Conc100.nii";

	////calibration GnP 20170829
	///*-------------------------------------------------*/
	//int numScan = 1; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration_gold_20170829/08_29_calib";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/20170829/Conc130_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/20170829/Conc130.nii";

	////calibration GnP 20171011
	///*-------------------------------------------------*/
	//int numScan = 1; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration_20171011";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/20171011/Conc100_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/20171011/Conc100.nii";
	///*-------------------------------------------------*/
	
	///*-------------------------------------------------*/
	///*-----------------------------------------------
	//20171017 series
	//-------------------------------------------------*/
	//int numScan = 4; // number of scans

	//std::string *path_imgSeries = new std::string[numScan];
	//std::string *filename_fullHead = new std::string[numScan];
	//std::string *filename_ear = new std::string[numScan];

	////baseline / scan1
	//path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171017/scan1/1.2.826.0.1.3417726.3.570864.20171017150642331";
	//filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan1/scan1_raw.nii";
	//filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan1/scan1.nii";

	////scan 2
	//path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171017/scan2/1.2.826.0.1.3417726.3.790181.20171017150641722";
	//filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan2/scan2_raw.nii";
	//filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan2/scan2.nii";

	////scan 3
	//path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171017/scan3/1.2.826.0.1.3417726.3.13422.20171017154128519";
	//filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan3/scan3_raw.nii";
	//filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan3/scan3.nii";

	////scan 4
	//path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171017/scan4/1.2.826.0.1.3417726.3.412125.20171017160117956";
	//filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan4/scan4_raw.nii";
	//filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan4/scan4.nii";

	///*-------------------------------------------------*/

	/*-------------------------------------------------*/
	/*-----------------------------------------------
	20171017 series
	-------------------------------------------------*/
	int numScan = 6; // number of scans

	std::string *path_imgSeries = new std::string[numScan];
	std::string *filename_fullHead = new std::string[numScan];
	std::string *filename_ear = new std::string[numScan];

	//baseline / scan1
	path_imgSeries[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171026/scan1/1.2.826.0.1.3417726.3.764447.20171026131023324";
	filename_fullHead[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan1/scan1_raw.nii";
	filename_ear[0] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan1/scan1.nii";

	//scan 2
	path_imgSeries[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171026/scan2/1.2.826.0.1.3417726.3.981882.20171026134025355";
	filename_fullHead[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan2/scan2_raw.nii";
	filename_ear[1] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan2/scan2.nii";

	//scan 3
	path_imgSeries[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171026/scan3/1.2.826.0.1.3417726.3.352511.20171026141205589";
	filename_fullHead[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan3/scan3_raw.nii";
	filename_ear[2] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan3/scan3.nii";

	//scan 4
	path_imgSeries[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171026/scan4/1.2.826.0.1.3417726.3.753296.20171026144249152";
	filename_fullHead[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan4/scan4_raw.nii";
	filename_ear[3] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan4/scan4.nii";

	//scan 5
	path_imgSeries[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171026/scan5/1.2.826.0.1.3417726.3.880144.20171026155606183";
	filename_fullHead[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan5/scan5_raw.nii";
	filename_ear[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan5/scan5.nii";

	//scan 6 - faulty scan
	path_imgSeries[5] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/20171026/scan6/1.2.826.0.1.3417726.3.513206.20171026152505527";
	filename_fullHead[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan6/scan6_raw.nii";
	filename_ear[5] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171026/scan6/scan6.nii";

	/*-------------------------------------------------*/

	for (size_t i = 5; i < 6; i++) //Change to 5 depending on number of scans
	{
		std::cout << i << std::endl;
		const char * DicomSeriesDic = path_imgSeries[i].c_str();
		const char * OutputFileName_fullHead = filename_fullHead[i].c_str();
		const char * OutputFileName_ear = filename_ear[i].c_str();

		typedef float    PixelType;
		const unsigned int      Dimension = 3;

		typedef itk::Image< PixelType, Dimension >         ImageType;

		typedef itk::ImageSeriesReader< ImageType >        ReaderType;
		ReaderType::Pointer reader = ReaderType::New();

		typedef itk::GDCMImageIO       ImageIOType;
		ImageIOType::Pointer dicomIO = ImageIOType::New();

		reader->SetImageIO(dicomIO);

		typedef itk::GDCMSeriesFileNames NamesGeneratorType;
		NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

		nameGenerator->SetUseSeriesDetails(true);
		//nameGenerator->AddSeriesRestriction("0008|0021");

		nameGenerator->SetDirectory(DicomSeriesDic);

		try
		{
			std::cout << std::endl << "The directory: " << std::endl;
			std::cout << std::endl << DicomSeriesDic << std::endl << std::endl;
			std::cout << "Contains the following DICOM Series: ";
			std::cout << std::endl << std::endl;


			typedef std::vector< std::string >    SeriesIdContainer;

			const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

			SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
			SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
			while (seriesItr != seriesEnd)
			{
				std::cout << seriesItr->c_str() << std::endl;
				++seriesItr;
			}


			std::string seriesIdentifier;

			if (argc > 3) // If no optional series identifier
			{
				seriesIdentifier = argv[1];
			}
			else
			{
				seriesIdentifier = seriesUID.begin()->c_str();
			}

			std::cout << std::endl << std::endl;
			std::cout << "Now reading series: " << std::endl << std::endl;
			std::cout << seriesIdentifier << std::endl;
			std::cout << std::endl << std::endl;

			typedef std::vector< std::string >   FileNamesContainer;
			FileNamesContainer fileNames;

			fileNames = nameGenerator->GetFileNames(seriesIdentifier);

			reader->SetFileNames(fileNames);

			// Comment bottom two lines if we want whole head
			typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
			FilterType::Pointer regionFilter = FilterType::New();

			
#pragma region "Size for each image seriers"
			////ROI for dicom_ear_male
			//ImageType::IndexType start;
			//start[0] = 700;
			//start[1] = 500;
			//start[2] = 8; 
			//ImageType::SizeType size;
			//size[0] = 250;
			//size[1] = 220;
			//size[2] = 164;

			////ROI for skull base
			//ImageType::IndexType start;
			//start[0] = 130;
			//start[1] = 200;
			//start[2] = 100;
			//ImageType::SizeType size;
			//size[0] = 260;
			//size[1] = 275;
			//size[2] = 300;


			////ROI for skull 120min
			//ImageType::IndexType start;
			//start[0] = 140;
			//start[1] = 250;
			//start[2] = 10;
			//ImageType::SizeType size;
			//size[0] = 240;
			//size[1] = 230;
			//size[2] = 200;

			////ROI for 160108data firstline/second
			//ImageType::IndexType start;
			//start[0] = 310;
			//start[1] = 460;
			//start[2] = 100;
			//ImageType::SizeType size;
			//size[0] = 240;
			//size[1] = 160;
			//size[2] = 220;

			////ROI for 160108data third
			//ImageType::IndexType start;
			//start[0] = 240;
			//start[1] = 150;
			//start[2] = 110;
			//ImageType::SizeType size;
			//size[0] = 270;
			//size[1] = 230;
			//size[2] = 190;

			////ROI for 160307 data 
			//ImageType::IndexType start;
			//start[0] = 490;
			//start[1] = 240;
			//start[2] = 60;
			//ImageType::SizeType size;
			//size[0] = 200;
			//size[1] = 170;
			//size[2] = 210;

			////ROI for 160719 data v1
			//ImageType::IndexType start;
			//start[0] = 130;
			//start[1] = 366;
			//start[2] = 60;
			//ImageType::SizeType size;
			//size[0] = 270;
			//size[1] = 270;
			//size[2] = 230;

			//ROI for 160719 data v2
			//ImageType::IndexType start;
			//start[0] = 183;
			//start[1] = 396;
			//start[2] = 63;
			//ImageType::SizeType size;
			//size[0] = 217;
			//size[1] = 183;
			//size[2] = 260;
			
			////ROI for 160906 data left ear baseline
			//ImageType::IndexType start;
			//start[0] = 315;
			//start[1] = 54;
			//start[2] = 1;
			//ImageType::SizeType size;
			//size[0] = 220;
			//size[1] = 217;
			//size[2] = 252;

			////ROI for 160906 data left ear other scans
			//ImageType::IndexType start;
			//start[0] = 312;
			//start[1] = 171;
			//start[2] = 1;
			//ImageType::SizeType size;
			//size[0] = 216;
			//size[1] = 220;
			//size[2] = 251;

			////ROI for 160629 data left ear other scans
			//ImageType::IndexType start;
			//start[0] = 189;
			//start[1] = 370;
			//start[2] = 84;
			//ImageType::SizeType size;
			//size[0] = 211;
			//size[1] = 216;
			//size[2] = 233;

			////ROI for 170111 data left ear baseline scans
			//ImageType::IndexType start;
			//start[0] = 354;
			//start[1] = 265;
			//start[2] = 222;
			//ImageType::SizeType size;
			//size[0] = 212;
			//size[1] = 212;
			//size[2] = 235;

			////ROI for 170111 data left ear other scans
			//ImageType::IndexType start;
			//start[0] = 343;
			//start[1] = 334;
			//start[2] = 149;
			//ImageType::SizeType size;
			//size[0] = 212;
			//size[1] = 212;
			//size[2] = 235;
			
			////ROI for 170123 data left ear baseline/all scans
			//ImageType::IndexType start;
			//start[0] = 272;
			//start[1] = 201;
			//start[2] = 187;
			//ImageType::SizeType size;
			//size[0] = 229;
			//size[1] = 218;
			//size[2] = 193;


			////ROI for 170123 data contra-lateral ear
			//ImageType::IndexType start;
			//start[0] = 626;
			//start[1] = 204;
			//start[2] = 187;
			//ImageType::SizeType size;
			//size[0] = 229;
			//size[1] = 218;
			//size[2] = 193;

			////ROI for 170606 data left ear baseline
			//ImageType::IndexType start;
			//start[0] = 131;
			//start[1] = 245;
			//start[2] = 60;
			//ImageType::SizeType size;
			//size[0] = 229;
			//size[1] = 218;
			//size[2] = 218;

			////ROI for 170606 data left ear others
			//ImageType::IndexType start;
			//start[0] = 133;
			//start[1] = 251;
			//start[2] = 178;
			//ImageType::SizeType size;
			//size[0] = 229;
			//size[1] = 218;
			//size[2] = 218;

			////ROI for 170612 data left ear others
			//ImageType::IndexType start;
			//start[0] = 230;
			//start[1] = 292;
			//start[2] = 69;
			//ImageType::SizeType size;
			//size[0] = 229;
			//size[1] = 218;
			//size[2] = 214;

			////ROI for 170619 data left ear others
			//ImageType::IndexType start;
			//start[0] = 149;
			//start[1] = 313;
			//start[2] = 63;
			//ImageType::SizeType size;
			//size[0] = 260;
			//size[1] = 260;
			//size[2] = 268;

			////ROI for 171011 data left ear others
			//ImageType::IndexType start;
			//start[0] = 152;
			//start[1] = 304;
			//start[2] = 115;
			//ImageType::SizeType size;
			//size[0] = 238;
			//size[1] = 188;
			//size[2] = 221;

			////ROI for calibration data
			//ImageType::IndexType start;
			//start[0] = 571;
			//start[1] = 674;
			//start[2] = 271;
			//ImageType::SizeType size;
			//size[0] = 37;
			//size[1] = 36;
			//size[2] = 225;

			//ROI for 171026 data left ear others. scan 1,2,3,4; scan 5; scan 6
			ImageType::IndexType start;
			start[0] = 169;//190; //189 //169
			start[1] = 360;//386; //389 //360
			start[2] = 180;//135; //189 //180
			ImageType::SizeType size;
			size[0] = 192; //183;
			size[1] = 182; //154;
			size[2] = 195; //173;

			

#pragma endregion

			ImageType::RegionType ear;
			ear.SetSize(size);
			ear.SetIndex(start);

			regionFilter->SetInput(reader->GetOutput());
			regionFilter->SetRegionOfInterest(ear);

			try
			{
				reader->Update();
			}
			catch (itk::ExceptionObject &ex)
			{
				std::cout << ex << std::endl;
				getchar();
				return EXIT_FAILURE;
			}


			const ImageType::DirectionType& direct = reader->GetOutput()->GetDirection();
			std::cout << "Direction = " << std::endl; std::cout << direct << std::endl;
			

			///*
#pragma region "Ways of contrast adjustment"

			//histogram matching image file path
			const std::string matchImageFile = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/atlasVolume_unlabeled.nii";

			typedef itk::ImageFileReader<ImageType> MatchImageReaderType;
			MatchImageReaderType::Pointer matchImageReader = MatchImageReaderType::New();
			matchImageReader->SetFileName(matchImageFile);
			matchImageReader->Update();

			typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HEFilterType;
			HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();
			IntensityEqualizeFilter->SetReferenceImage(matchImageReader->GetOutput());
			IntensityEqualizeFilter->SetInput(regionFilter->GetOutput());
			//IntensityEqualizeFilter->SetInput(reader->GetOutput());
			IntensityEqualizeFilter->SetNumberOfHistogramLevels(500);
			IntensityEqualizeFilter->SetNumberOfMatchPoints(20);
			IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();
			IntensityEqualizeFilter->Update();

			typedef itk::NormalizeImageFilter< ImageType, ImageType > NormalizeFilterType;
			NormalizeFilterType::Pointer normalizeFilter = NormalizeFilterType::New();
			normalizeFilter->SetInput(regionFilter->GetOutput());
			//normalizeFilter->SetInput(reader->GetOutput());
			normalizeFilter->Update();
			//*/
#pragma endregion

			typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescalerType;
			RescalerType::Pointer intensityRescaler = RescalerType::New();
			intensityRescaler->SetInput(regionFilter->GetOutput());
			//intensityRescaler->SetInput(reader->GetOutput());

			intensityRescaler->SetOutputMinimum(0);
			intensityRescaler->SetOutputMaximum(255);

			//ImageType::Pointer RescaledRegion = intensityRescaler->GetOutput();
			ImageType::Pointer RescaledRegion = regionFilter->GetOutput();//for HU values
			//ImageType::Pointer RescaledRegion = reader->GetOutput();//for HU values

			const double origin[3] = { 0.0, 0.0, 0.0 };
			RescaledRegion->SetOrigin(origin);

			typedef itk::ImageFileWriter< ImageType > WriterType;
			WriterType::Pointer writer = WriterType::New();

#pragma region "write 3D full head image"
			
			writer->SetFileName(OutputFileName_fullHead);
			//writer->SetInput(reader->GetOutput());
			writer->SetInput(RescaledRegion);

			std::cout << "Writing the image as " << std::endl << std::endl;
			std::cout << OutputFileName_fullHead << std::endl << std::endl;
			try
			{
				writer->Update();
			}
			catch (itk::ExceptionObject &ex)
			{
				std::cout << ex << std::endl;
				getchar();
				return EXIT_FAILURE;
			}
			
#pragma endregion

			
			// write 3D ear image
			writer->SetFileName(OutputFileName_ear);
			writer->SetInput(RescaledRegion);
			writer->SetInput(intensityRescaler->GetOutput());

			std::cout << "Writing the image as " << std::endl << std::endl;	
			std::cout << OutputFileName_ear << std::endl << std::endl;
			try
			{
				writer->Update();
			}
			catch (itk::ExceptionObject &ex)
			{
				std::cout << "It does not write!" << std::endl;
				std::cout << ex << std::endl;
				getchar();
				return EXIT_FAILURE;
			}



		}

		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			getchar();
			return EXIT_FAILURE;
		}
	}
	std::cout << "Done!" << std::endl;
	return EXIT_SUCCESS;
}
