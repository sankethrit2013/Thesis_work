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
//  The following code is an implementation of a small ITK
//  program. It tests including header files and linking with ITK
//  libraries.
//
//  Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkNiftiImageIO.h"

#include <iostream>
using namespace std;
int main(int, char * argv[])
{

	typedef itk::Image< float, 3 > ImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;

	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	ImageType::Pointer atlas = ImageType::New();
	ImageType::Pointer atlas_reoriented = ImageType::New();

	typedef itk::NiftiImageIO       ImageIOType;
	ImageIOType::Pointer niftiIO = ImageIOType::New();


	/*const char * inputFilename = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/atlasVolume_preRegistered.nii";
	const char * outputFilename = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/atlasVolume_preRegistered_reoriented.nii";*/

	// input for unlabel atlas
	/*const char * inputFilename = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/atlasVolume_unlabeled.nii";
	const char * outputFilename = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/atlasVolume_unlabeled_reoriented.nii";*/ 

	// input for label map
	const char * inputFilename = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/atlasVolume_SL_labelmap.nii";
	const char * outputFilename = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/atlasVolume_SL_labelmap_reoriented.nii";

	reader->SetFileName(inputFilename);
	reader->SetImageIO(niftiIO);

	atlas = reader->GetOutput();
	
	typedef itk::PermuteAxesImageFilter <ImageType>
		PermuteAxesImageFilterType;

	itk::FixedArray<unsigned int, 3> order;
	order[0] = 2;
	order[1] = 0;
	order[2] = 1;
	PermuteAxesImageFilterType::Pointer permuteAxesFilter
		= PermuteAxesImageFilterType::New();
	permuteAxesFilter->SetInput(atlas);
	permuteAxesFilter->SetOrder(order);
	permuteAxesFilter->Update();


	
	typedef itk::FlipImageFilter< ImageType >  FlipFilterType;
	FlipFilterType::Pointer flipfilter = FlipFilterType::New();
	typedef FlipFilterType::FlipAxesArrayType     FlipAxesArrayType;
	FlipAxesArrayType flipArray;
	//flipArray[0] = 0;
	//flipArray[1] = 0;
	//flipArray[2] = 1;

	flipArray[0] = 1;
	flipArray[1] = 1;
	flipArray[2] = 1;
	flipfilter->SetFlipAxes(flipArray);
	flipfilter->SetInput(permuteAxesFilter->GetOutput());

	atlas_reoriented = flipfilter->GetOutput();
	const double origin[3] = { 0.0, 0.0, 0.0 };	

	ImageType::DirectionType direction;
	direction.Fill(0);
	direction(0, 0) = -1;
	direction(1, 1) = 1;
	direction(2, 2) = -1;

	atlas_reoriented->SetOrigin(origin);
	atlas_reoriented->SetDirection(direction);
	direction=atlas_reoriented->GetDirection();
	std::cout << "Direction = " << std::endl; std::cout << direction << std::endl;

	writer->SetFileName(outputFilename);
	writer->SetInput(atlas_reoriented);
	writer->SetImageIO(niftiIO);
	
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		system("pause");

		return EXIT_FAILURE;
	}

	system("pause");
	return EXIT_SUCCESS;

}
// Software Guide : EndCodeSnippet

//  Software Guide : BeginLatex
//
//  This code instantiates a $3D$ image\footnote{Also known as a
//  \emph{volume}.} whose pixels are represented with type \code{unsigned
//  short}. The image is then constructed and assigned to a
//  \doxygen{SmartPointer}. Although later in the text we will discuss
//  \code{SmartPointer}s in detail, for now think of it as a handle on an
//  instance of an object (see section \ref{sec:SmartPointers} for more
//  information). The \doxygen{Image} class will be described in
//  Section~\ref{sec:ImageSection}.
//
//  Software Guide : EndLatex
