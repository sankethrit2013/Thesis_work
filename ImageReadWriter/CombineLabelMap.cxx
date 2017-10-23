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
#include "itkNiftiImageIO.h"
#include <iostream>
#include "itkBinaryThresholdImageFilter.h"
#include "itkAddImageFilter.h"
#include <string>
#include <itkNaryAddImageFilter.h>

using namespace std;
int main(int, char * argv[])
{

	//char* labelmapfolder = "C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan1/labelmap_Deformable_Grid10_Iter100/";
	char* labelmapfolder = "N:\Sanketh_Moudgalya_Directory\Nate\ImageVolume\20160719";

	char labelmap_SV[300];
	char labelmap_SM[300];
	char labelmap_SL[300];
	char labelmap_ST[300];
	char labelmap_combined[300];

	//sprintf_s(labelmap_SV, "%sSV_labelmap_registrationOutput.nii", labelmapfolder);
	//sprintf_s(labelmap_SM, "%sSM_labelmap_registrationOutput.nii", labelmapfolder);
	//sprintf_s(labelmap_SL, "%sSL_labelmap_registrationOutput.nii", labelmapfolder);
	//sprintf_s(labelmap_ST, "%sST_labelmap_registrationOutput.nii", labelmapfolder);

	sprintf_s(labelmap_SV, "%satlasVolume_SV_labelmap.nii", labelmapfolder);
	sprintf_s(labelmap_SM, "%satlasVolume_SM_labelmap.nii", labelmapfolder);
	sprintf_s(labelmap_SL, "%satlasVolume_SL_labelmap.nii", labelmapfolder);
	sprintf_s(labelmap_ST, "%satlasVolume_ST_labelmap.nii", labelmapfolder);


	sprintf_s(labelmap_combined, "%slabelmap_combined.nii", labelmapfolder);



	typedef itk::Image< unsigned char, 3 > ImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileWriter< ImageType > WriterType;

	ReaderType::Pointer SLreader = ReaderType::New();
	ReaderType::Pointer SVreader = ReaderType::New();
	ReaderType::Pointer SMreader = ReaderType::New();
	ReaderType::Pointer STreader = ReaderType::New();

	WriterType::Pointer writer = WriterType::New();

	ImageType::Pointer SV = ImageType::New();
	ImageType::Pointer ST = ImageType::New();
	ImageType::Pointer SM = ImageType::New();
	ImageType::Pointer SL = ImageType::New();

	typedef itk::NiftiImageIO       ImageIOType;
	ImageIOType::Pointer niftiIO = ImageIOType::New();

	SVreader->SetImageIO(niftiIO);
	SVreader->SetFileName(labelmap_SV);
	SVreader->Update();
	SV = SVreader->GetOutput();

	SMreader->SetImageIO(niftiIO);
	SMreader->SetFileName(labelmap_SM);
	SMreader->Update();
	SM = SMreader->GetOutput();

	STreader->SetImageIO(niftiIO);
	STreader->SetFileName(labelmap_ST);
	STreader->Update();
	ST = STreader->GetOutput();

	SLreader->SetImageIO(niftiIO);
	SLreader->SetFileName(labelmap_SL);
	SLreader->Update();
	SL = SLreader->GetOutput();

	typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > BinaryFilterType;
	BinaryFilterType::Pointer SM_bin = BinaryFilterType::New();
	BinaryFilterType::Pointer SV_bin = BinaryFilterType::New();
	BinaryFilterType::Pointer ST_bin = BinaryFilterType::New();
	BinaryFilterType::Pointer SL_bin = BinaryFilterType::New();

	SM_bin->SetInput(SM);
	SM_bin->SetOutsideValue(0);
	SM_bin->SetInsideValue(1);
	SM_bin->SetLowerThreshold(128);
	SM_bin->SetUpperThreshold(255);
	
	ST_bin->SetInput(ST);
	ST_bin->SetOutsideValue(0);
	ST_bin->SetInsideValue(2);
	ST_bin->SetLowerThreshold(128);
	ST_bin->SetUpperThreshold(255);

	SV_bin->SetInput(SV);
	SV_bin->SetOutsideValue(0);
	SV_bin->SetInsideValue(3);
	SV_bin->SetLowerThreshold(128);
	SV_bin->SetUpperThreshold(255);

	SL_bin->SetInput(SL);
	SL_bin->SetOutsideValue(0);
	SL_bin->SetInsideValue(4);
	SL_bin->SetLowerThreshold(128);
	SL_bin->SetUpperThreshold(255);
	typedef itk::NaryAddImageFilter<
		ImageType,
		ImageType  >  addFilterType;

	addFilterType::Pointer adder = addFilterType::New();
	adder->SetInput(0, SV_bin->GetOutput());
	adder->SetInput(1, SM_bin->GetOutput());
	adder->SetInput(2, ST_bin->GetOutput());
	adder->SetInput(3, SL_bin->GetOutput());

	adder->Update();

	writer->SetInput(adder->GetOutput());
	writer->SetFileName(labelmap_combined);
	writer->Update();
    //system("pause");
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
