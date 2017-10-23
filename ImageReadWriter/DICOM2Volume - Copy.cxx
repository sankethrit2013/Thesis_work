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
#include "itkIntensityWindowingImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include <iostream>
#include <string>
using std::string;

// Software Guide : EndCodeSnippet

int main(int argc, char* argv[])
{
	argc = 3;
	// dicom series folder path list
	std::string *dicomFolder = new string[9];
	dicomFolder[0] = "sddsdsd";
	dicomFolder[1] = "";
	dicomFolder[2] = "";
	dicomFolder[3] = "";
	dicomFolder[4] = "";
	dicomFolder[5] = "";
	dicomFolder[6] = "";
	dicomFolder[7] = "";
	dicomFolder[8] = "";

	std::cout << dicomFolder[0]<<std::endl;

	//argv[1] = "N:/Zhenlin_Xu/Images/Calibrations/sample_1_18.98mg_ml/1.2.826.0.1.3417726.3.3442824753/1.2.826.0.1.3417726.3.342043.20160406162917645";

	////output file name list

	//argv[2] = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/Calibration/sample_1_18.98mg_ml.nii.gz";



	//for (int i = 1; i < 10; i++)
	//{
	//	const char * DicomSeriesDic = argv[1];
	//	const char * OutputFileName = argv[2];
	//	const char * OutputFileName = argv[2];
	//	typedef float    PixelType;
	//	const unsigned int      Dimension = 3;

	//	typedef itk::Image< PixelType, Dimension >         ImageType;

	//	typedef itk::ImageSeriesReader< ImageType >        ReaderType;
	//	ReaderType::Pointer reader = ReaderType::New();

	//	typedef itk::GDCMImageIO       ImageIOType;
	//	ImageIOType::Pointer dicomIO = ImageIOType::New();

	//	reader->SetImageIO(dicomIO);

	//	typedef itk::GDCMSeriesFileNames NamesGeneratorType;
	//	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

	//	nameGenerator->SetUseSeriesDetails(true);
	//	//nameGenerator->AddSeriesRestriction("0008|0021");

	//	nameGenerator->SetDirectory(DicomSeriesDic);

	//	typedef   itk::CurvatureAnisotropicDiffusionImageFilter <
	//		ImageType,
	//		ImageType > SmoothingFilterType;

	//	SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
	//typedef itk::IntensityWindowingImageFilter <ImageType, ImageType> IntensityWindowingImageFilterType;

	//IntensityWindowingImageFilterType::Pointer filter = IntensityWindowingImageFilterType::New();


	//	try
	//	{
	//		std::cout << std::endl << "The directory: " << std::endl;
	//		std::cout << std::endl << DicomSeriesDic << std::endl << std::endl;
	//		std::cout << "Contains the following DICOM Series: ";
	//		std::cout << std::endl << std::endl;


	//		typedef std::vector< std::string >    SeriesIdContainer;

	//		const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

	//		SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
	//		SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
	//		while (seriesItr != seriesEnd)
	//		{
	//			std::cout << seriesItr->c_str() << std::endl;
	//			++seriesItr;
	//		}


	//		std::string seriesIdentifier;

	//		if (argc > 3) // If no optional series identifier
	//		{
	//			seriesIdentifier = argv[1];
	//		}
	//		else
	//		{
	//			seriesIdentifier = seriesUID.begin()->c_str();
	//		}

	//		std::cout << std::endl << std::endl;
	//		std::cout << "Now reading series: " << std::endl << std::endl;
	//		std::cout << seriesIdentifier << std::endl;
	//		std::cout << std::endl << std::endl;

	//		typedef std::vector< std::string >   FileNamesContainer;
	//		FileNamesContainer fileNames;

	//		fileNames = nameGenerator->GetFileNames(seriesIdentifier);

	//		reader->SetFileNames(fileNames);

	//		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > FilterType;
	//		FilterType::Pointer regionFilter = FilterType::New();


	//		//ROI for skull base
	//		ImageType::IndexType start;
	//		start[0] = 327;
	//		start[1] = 258;
	//		start[2] = 1;

	//		ImageType::SizeType size;
	//		size[0] = 420;
	//		size[1] = 400;
	//		size[2] = 512;

	//		ImageType::RegionType region;
	//		region.SetSize(size);
	//		region.SetIndex(start);

	//		regionFilter->SetRegionOfInterest(region);
	//		regionFilter->SetInput(reader->GetOutput());

	//		try
	//		{
	//			reader->Update();
	//		}
	//		catch (itk::ExceptionObject &ex)
	//		{
	//			std::cout << ex << std::endl;
	//			return EXIT_FAILURE;
	//		}

	//		typedef itk::NormalizeImageFilter< ImageType, ImageType > NormalizeFilterType;
	//		NormalizeFilterType::Pointer normalizeFilter = NormalizeFilterType::New();
	//		normalizeFilter->SetInput(reader->GetOutput());
	//		normalizeFilter->Update();


	//		typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescalerType;
	//		RescalerType::Pointer intensityRescaler = RescalerType::New();
	//		intensityRescaler->SetInput(regionFilter->GetOutput());

	//		intensityRescaler->SetOutputMinimum(0);
	//		intensityRescaler->SetOutputMaximum(255);

	//		ImageType::Pointer RescaledRegion = intensityRescaler->GetOutput();

	//		const double origin[3] = { 0.0, 0.0, 0.0 };
	//		RescaledRegion->SetOrigin(origin);

	//		typedef itk::ImageFileWriter< ImageType > WriterType;
	//		WriterType::Pointer writer = WriterType::New();

	//		writer->SetFileName(OutputFileName);

	//		writer->SetInput(regionFilter->GetOutput());

	//		const ImageType::DirectionType& direct = intensityRescaler->GetOutput()->GetDirection();
	//		std::cout << "Direction = " << std::endl; std::cout << direct << std::endl;

	//		// Software Guide : EndCodeSnippet

	//		std::cout << "Writing the image as " << std::endl << std::endl;
	//		std::cout << OutputFileName << std::endl << std::endl;



	//		try
	//		{
	//			writer->Update();
	//		}
	//		catch (itk::ExceptionObject &ex)
	//		{
	//			std::cout << ex << std::endl;
	//			return EXIT_FAILURE;
	//		}

//filter->SetInput(image);
//filter->SetWindowMinimum(0);
//filter->SetWindowMaximum(100);
//filter->SetOutputMinimum(0);
//filter->SetOutputMaximum(255);
//filter->Update();
	//		smoothing->SetInput(regionFilter->GetOutput());
	//		smoothing->SetTimeStep(0.001);//time step should be at or below PixelSpacing/2^(Dimemsion+1)
	//		smoothing->SetNumberOfIterations(10);
	//		smoothing->SetConductanceParameter(20);
	//		writer->SetFileName(OutputFileName);
	//		writer->SetInput(regionFilter->GetOutput());


	//	}
	//	catch (itk::ExceptionObject &ex)
	//	{
	//		std::cout << ex << std::endl;
	//		return EXIT_FAILURE;
	//	}
	//}
	return EXIT_SUCCESS;

}