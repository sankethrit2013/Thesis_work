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

//  Software Guide : BeginLatex
//
//  This example illustrates how to read a series of 2D slices from independent
//  files in order to compose a volume. The class \doxygen{ImageSeriesReader}
//  is used for this purpose. This class works in combination with a generator
//  of filenames that will provide a list of files to be read. In this
//  particular example we use the \doxygen{NumericSeriesFileNames} class as a
//  filename generator. This generator uses a \code{printf} style of string format
//  with a ``\code{\%d}'' field that will be successively replaced by a number specified
//  by the user. Here we will use a format like ``\code{file\%03d.png}'' for reading
//  PNG files named file001.png, file002.png, file003.png... and so on.
//
//  This requires the following headers as shown.
//
//  \index{itk::ImageSeriesReader!header}
//  \index{itk::NumericSeriesFileNames!header}
//
//  Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkJPEGImageIO.h"
#include "itkPNGImageIO.h"

int main(int argc, char ** argv)
{
	//argc = 5;
	//// input for unlabeled atlas
	//argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/unlabeled/ps%04d.jpg";
	//argv[2] = "2001";
	//argv[3] = "2069";
	////argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/atlasVolume_unlabeled.nii";
	//argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/atlas_trial.nii";

	argc = 5;
	// input for unlabeled atlas
	/*argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/unlabeled/ps%04d.jpg";
	argv[2] = "2001";
	argv[3] = "2069";
	argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/atlas_contralateral.nii";*/
	//argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/all_segments.nii";

	
	//// input for SL Slabel map
	/*argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/SL_labelmap/ps%04d.jpg";
	argv[2] = "2001";
	argv[3] = "2069";
	argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/SL_labelmap_contralateral.nii";*/

	//// input for ST Slabel map
	/*argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/ST_labelmap/ps%04d.jpg";
	argv[2] = "2001";
	argv[3] = "2069";
	argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/ST_labelmap_contralateral.nii";*/

	//// input for SV Slabel map
	/*argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/SV_labelmap/ps%04d.jpg";
	argv[2] = "2001";
	argv[3] = "2069";
	argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/SV_labelmap_contralateral.nii";*/

	////// input for SM Slabel map
	//argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/SM_labelmap/ps%04d.jpg";
	//argv[2] = "2001";
	//argv[3] = "2069";
	//argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/SM_labelmap_contralateral.nii";

	//// input for SM Slabel map
	argv[1] = "E:/Sanketh_Moudgalya_Directory/Nate/Images/3D Atlas/ps2000/CA_labelmap/ps%04d.jpg";
	argv[2] = "2001";
	argv[3] = "2069";
	argv[4] = "E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/CA_labelmap.nii";
	

	// Verify the number of parameters in the command line
	//if (argc < 5)
	//{
	//	std::cerr << "Usage: " << std::endl;
	//	std::cerr << argv[0] << "SeriesFormat firstSliceValue lastSliceValue  outputImageFile " << std::endl;
	//	return EXIT_FAILURE;
	//}

	typedef unsigned int                       PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image< unsigned char, 3 >  JPEGImageType;
	typedef itk::Image< unsigned char, 3 >  PNGImageType;

	typedef itk::Image< unsigned char, 3 >  NiftiImageType;

	typedef itk::ImageSeriesReader< PNGImageType >  ReaderType; // use JPEGImageType or PNGImageType
	typedef itk::ImageFileWriter<   NiftiImageType >  WriterType;

	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();


	const unsigned int first = atoi(argv[2]);
	const unsigned int last = atoi(argv[3]);

	const char * outputFilename = argv[4];

	//const unsigned int first = 2001;
	//const unsigned int last = 2069;
	//const char * outputFilename = "C:/Users/zx4403/Dropbox/Research/Borkholder - Micropump/Images/3D Atlas/ps2000/unlabeled/atlasVolume_unlabeled.nii";
	
	typedef itk::NumericSeriesFileNames NameGeneratorType;

	NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
	
	//nameGenerator->SetSeriesFormat();
	nameGenerator->SetSeriesFormat(argv[1]);

	nameGenerator->SetStartIndex(first);
	nameGenerator->SetEndIndex(last);
	nameGenerator->SetIncrementIndex(1);
	
	reader->SetImageIO(itk::JPEGImageIO::New());
	//reader->SetImageIO(itk::PNGImageIO::New());

	reader->SetFileNames(nameGenerator->GetFileNames());

	NiftiImageType::Pointer imageVolume = reader->GetOutput();

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Reader ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// set the spacing (units: um)
	NiftiImageType::SpacingType spacing;
	spacing[0] = 0.005435;
	spacing[1] = 0.005435;
	spacing[2] = 0.020;

	/*spacing[0] = 0.01638;
	spacing[1] = 0.01638;
	spacing[2] = 0.01638;*/

	imageVolume->SetSpacing(spacing);

	const double origin[3] = { 0.0, 0.0, 0.0 };
	imageVolume->SetOrigin(origin);

	//const NiftiImageType::DirectionType& direction = imageVolume->GetDirection();
	NiftiImageType::DirectionType direction;
	direction.Fill(0);

	/*
	S/I = direction(2, 1) = -1/1
	R/L = direction(0, 2) = 1/-1
	P/A = direction(1, 0) = -1/1

	*/
	// atlas direction for 20150801 FFP
	//direction(0, 2) = -1;
	//direction(1, 0) = 1;
	//direction(2, 1) = 1;

	// atlas direction for 20151209 HFP
	//direction(0, 2) = -1;
	//direction(1, 0) = -1;
	//direction(2, 1) = -1;

	// atlas direction for 20160108 FFP
	//direction(0, 2) = 1;
	//direction(1, 0) = -1;
	//direction(2, 1) = 1;

	//// atlas direction for 20160307 HFP
	//direction(0, 2) =  1;
	//direction(1, 0) = -1;
	//direction(2, 1) = -1;

	// atlas direction for 20160906 HFP
	direction(0, 2) = -1;
	direction(1, 0) = -1;
	direction(2, 1) = -1;

	imageVolume->SetDirection(direction);
	//std::cout << "Direction = " << std::endl; std::cout << direction << std::endl;

	writer->SetFileName(outputFilename);
	
	writer->SetInput(imageVolume);
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Writer ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}
	// Software Guide : EndCodeSnippet
	return EXIT_SUCCESS;
}
