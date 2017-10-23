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
#include "itkIntensityWindowingImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include <iostream>
#include <string>
#include <stdio.h>

using std::string;

// Software Guide : EndCodeSnippet

int main(int argc, char* argv[])
{
	argc = 3;
	// dicom series folder path list for first calibration
	//std::string *dicomFolder = new string[9];
	//dicomFolder[0] = "N:/Zhenlin_Xu/Images/Calibrations/sample_1_18.98mg_ml/1.2.826.0.1.3417726.3.3442824753/1.2.826.0.1.3417726.3.342043.20160406162917645";
	//dicomFolder[1] = "N:/Zhenlin_Xu/Images/Calibrations/sample_2_40.68mg_ml/1.2.826.0.1.3417726.3.3442887283/1.2.826.0.1.3417726.3.895282.20160406162917707";
	//dicomFolder[2] = "N:/Zhenlin_Xu/Images/Calibrations/sample_3_75.953mg_ml/1.2.826.0.1.3417726.3.3443713801/1.2.826.0.1.3417726.3.264213.20160406162917769";
	//dicomFolder[3] = "N:/Zhenlin_Xu/Images/Calibrations/sample_4_100.34mg_ml/1.2.826.0.1.3417726.3.3443776142/1.2.826.0.1.3417726.3.650837.20160406162917832";
	//dicomFolder[4] = "N:/Zhenlin_Xu/Images/Calibrations/sample_5_149.16mg_ml/1.2.826.0.1.3417726.3.3444670994/1.2.826.0.1.3417726.3.118749.20160406163944885";
	//dicomFolder[5] = "N:/Zhenlin_Xu/Images/Calibrations/sample_6_200.69mg_ml/1.2.826.0.1.3417726.3.3444733299/1.2.826.0.1.3417726.3.164822.20160406171517663";
	//dicomFolder[6] = "N:/Zhenlin_Xu/Images/Calibrations/sample_7_249.51mg_ml/1.2.826.0.1.3417726.3.3444794802/1.2.826.0.1.3417726.3.174663.20160407105808816";
	//dicomFolder[7] = "N:/Zhenlin_Xu/Images/Calibrations/sample_Optiray_320mg_ml/1.2.826.0.1.3417726.3.3445622316/1.2.826.0.1.3417726.3.689734.20160407110226274";
	//dicomFolder[8] = "N:/Zhenlin_Xu/Images/Calibrations/sample_saline_1g_ml/1.2.826.0.1.3417726.3.32558062/1.2.826.0.1.3417726.3.321899.20160407142617612";
	//int ROI_indexX[] = { 327, 250, 270, 300, 300, 320, 300, 300, 250 };
	//int ROI_indexY[] = { 258, 250, 250, 250, 250, 250, 250, 250, 220 };


	//// dicom series folder path list for second calibration
	//std::string *dicomFolder = new string[7];
	//dicomFolder[0] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample1_18mg_per_ml/1.2.826.0.1.3417726.3.1605811030/1.2.826.0.1.3417726.3.393951.20160608181920937";
	//dicomFolder[1] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample2_40mg_per_ml/1.2.826.0.1.3417726.3.1606617335/1.2.826.0.1.3417726.3.310036.20160609133338137";
	//dicomFolder[2] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample3_75mg_per_ml/1.2.826.0.1.3417726.3.2488552027/1.2.826.0.1.3417726.3.996713.20160609133338199";
	//dicomFolder[3] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample4_100mg_per_ml/1.2.826.0.1.3417726.3.2488634665/1.2.826.0.1.3417726.3.545986.20160609133338277";
	//dicomFolder[4] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample5_149_mg_per_ml/1.2.826.0.1.3417726.3.2488698991/1.2.826.0.1.3417726.3.663595.20160609133338340";
	//dicomFolder[5] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample6_200mg_per_ml/1.2.826.0.1.3417726.3.2489529450/1.2.826.0.1.3417726.3.164421.20160609135326934";
	//dicomFolder[6] = "N:/Zhenlin_Xu/Images/Calibrations/20160714/sample7_240mg_per_ml/1.2.826.0.1.3417726.3.1605660124/1.2.826.0.1.3417726.3.296669.20160608172520718";
	
	// dicom series folder path list for second calibration
	std::string *dicomFolder = new string[7];
	dicomFolder[0] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample1/1.2.826.0.1.3417726.3.393951.20160608181920937";
	dicomFolder[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample2/1.2.826.0.1.3417726.3.310036.20160609133338137";
	dicomFolder[2] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample3/1.2.826.0.1.3417726.3.996713.20160609133338199";
	dicomFolder[3] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample4/1.2.826.0.1.3417726.3.545986.20160609133338277";
	dicomFolder[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample5/1.2.826.0.1.3417726.3.663595.20160609133338340";
	dicomFolder[5] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample6/1.2.826.0.1.3417726.3.164421.20160609135326934";
	dicomFolder[6] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/Calibration/sample7/1.2.826.0.1.3417726.3.296669.20160608172520718";

	
	int ROI_indexX = 335;
	int ROI_indexY = 240;


	////output file name list
	std::string outputFolder1 = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration";
	std::string outputFolder2 = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/Calibration/20160714";

	for (int i = 0; i < 7; i++)
	{
		const char * DicomSeriesDic = dicomFolder[i].c_str();
		const char * OutputFolder1 = outputFolder1.c_str();
		const char * OutputFolder2 = outputFolder2.c_str();

		char originalFileName[300];
		char regionFileName[300];
		char rescaledFileName[300];

		sprintf_s(originalFileName, sizeof(originalFileName), "%s/sample%d/sample%d_orginal.nii.gz", OutputFolder2, i + 1, i + 1);
		sprintf_s(regionFileName, sizeof(regionFileName), "%s/sample%d/sample%d_region.nii", OutputFolder2, i + 1, i + 1);
		sprintf_s(rescaledFileName, sizeof(rescaledFileName), "%s/sample%d/sample%d_rescaledRegion.nii", OutputFolder2, i + 1, i + 1);

		//std::cout << originalFileName << std::endl;
		//std::cout << regionFileName << std::endl;
		//std::cout << rescaledFileName << std::endl;

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

		typedef   itk::CurvatureAnisotropicDiffusionImageFilter <
			ImageType,
			ImageType > SmoothingFilterType;

		SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

		typedef itk::IntensityWindowingImageFilter <ImageType, ImageType> IntensityWindowingImageFilterType;
		IntensityWindowingImageFilterType::Pointer intensityWindowingFilter = IntensityWindowingImageFilterType::New();

		typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescalerType;
		RescalerType::Pointer intensityRescaler = RescalerType::New();


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

			try
			{
				reader->Update();
			}
			catch (itk::ExceptionObject &ex)
			{
				std::cout << ex << std::endl;
				return EXIT_FAILURE;
			}

			typedef itk::ImageFileWriter< ImageType > WriterType;
			WriterType::Pointer writer = WriterType::New();

			/*writer->SetFileName(originalFileName);
			writer->SetInput(reader->GetOutput());

			const ImageType::DirectionType& direct = reader->GetOutput()->GetDirection();
			std::cout << "Direction = " << std::endl; std::cout << direct << std::endl;

			std::cout << "Writing the original image as " << std::endl << std::endl;
			std::cout << originalFileName << std::endl << std::endl;
			try
			{
				writer->Update();
			}
			catch (itk::ExceptionObject &ex)
			{
				std::cout << ex << std::endl;
				return EXIT_FAILURE;
			}*/

				//Set ROI region 

				typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
				FilterType::Pointer regionFilter = FilterType::New();

				ImageType::IndexType start;
				start[0] = 371;
				start[1] = 247;
				start[2] = 0;

				ImageType::SizeType size;
				size[0] = 358;
				size[1] = 377;
				size[2] = 508;

				ImageType::RegionType region;
				region.SetSize(size);
				region.SetIndex(start);

				regionFilter->SetRegionOfInterest(region);
				regionFilter->SetInput(reader->GetOutput());
				
				writer->SetFileName(regionFileName);

				writer->SetInput(regionFilter->GetOutput());
				
				// Software Guide : EndCodeSnippet

				std::cout << "Writing the ROI image as " << std::endl << std::endl;
				std::cout << regionFileName << std::endl << std::endl;



				try
				{
					writer->Update();
				}
				catch (itk::ExceptionObject &ex)
				{
					std::cout << ex << std::endl;
					return EXIT_FAILURE;
				}

				intensityRescaler->SetInput(regionFilter->GetOutput());
				intensityRescaler->SetOutputMinimum(0);
				intensityRescaler->SetOutputMaximum(255);

				//intensityWindowingFilter->SetInput(regionFilter->GetOutput());
				//intensityWindowingFilter->SetWindowMinimum(0);
				//intensityWindowingFilter->SetWindowMaximum(100);
				//intensityWindowingFilter->SetOutputMinimum(0);
				//intensityWindowingFilter->SetOutputMaximum(255);
				//intensityWindowingFilter->Update();
				//smoothing->SetInput(regionFilter->GetOutput());
				//smoothing->SetTimeStep(0.001);//time step should be at or below PixelSpacing/2^(Dimemsion+1)
				//smoothing->SetNumberOfIterations(10);
				//smoothing->SetConductanceParameter(20);

				writer->SetFileName(rescaledFileName);
				writer->SetInput(intensityRescaler->GetOutput());
				std::cout << "Writing the recaled ROI image as " << std::endl << std::endl;
				std::cout << rescaledFileName << std::endl << std::endl;

				try
				{
					writer->Update();
				}
				catch (itk::ExceptionObject &ex)
				{
					std::cout << ex << std::endl;
					return EXIT_FAILURE;
				}

		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			return EXIT_FAILURE;
		}

	}
	return EXIT_SUCCESS;

}