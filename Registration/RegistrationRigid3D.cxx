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



#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkCorrelationImageToImageMetricv4.h"

#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkCenteredVersorTransformInitializer.h"

#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkQuasiNewtonOptimizerv4.h"
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"

class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef  itk::SmartPointer<Self>  Pointer;
	itkNewMacro(Self);

protected:
	CommandIterationUpdate() : m_CumulativeIterationIndex(0) {};

public:
	typedef   itk::GradientDescentOptimizerv4Template<double>  OptimizerType;
	typedef   const OptimizerType *                            OptimizerPointer;

	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *)caller, event);
	}

		void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		OptimizerPointer optimizer = static_cast<OptimizerPointer>(object);
		if (optimizer == ITK_NULLPTR)
		{
			return; // in this unlikely context, just do nothing.
		}
		if (!(itk::IterationEvent().CheckEvent(&event)))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << "  ";
		std::cout << m_CumulativeIterationIndex++ << std::endl;
	}

private:
	unsigned int m_CumulativeIterationIndex;
};


int main(int argc, char *argv[])
{

	const unsigned int                          Dimension = 3;
	typedef double                             PixelType;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;
	typedef itk::Image< PixelType, Dimension >  ResampleInputImageType;


	typedef itk::VersorRigid3DTransform< double > TransformType;
	
	typedef itk::RegularStepGradientDescentOptimizerv4<double>    OptimizerType;

	typedef itk::CorrelationImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType,
		TransformType >           RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);

	// Reader and writer
	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
	typedef itk::ImageFileReader< ResampleInputImageType > ResampleInputImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	ResampleInputImageReaderType::Pointer resampleInputImageReader = ResampleInputImageReaderType::New();
	fixedImageReader->SetFileName(argv[1]);
	movingImageReader->SetFileName(argv[2]);
	resampleInputImageReader->SetFileName(argv[3]);

	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());


	// Set initial transform
	TransformType::Pointer  transform = TransformType::New();

	typedef itk::CenteredVersorTransformInitializer<FixedImageType, MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer = TransformInitializerType::New();

	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());

	/*initializer->MomentsOn();*/
	initializer->InitializeTransform();
	transform->SetIdentity();

	//typedef TransformType::VersorType  VersorType;
	//typedef VersorType::VectorType     VectorType;
	//VersorType     rotation;
	//VectorType     axis;
	//axis[0] = 0.0;
	//axis[1] = 0.0;
	//axis[2] = 1.0;
	//const double angle = 0;
	//rotation.Set(axis, angle);
	//transform->SetRotation(rotation);
	//
	//typedef TransformType::MatrixType MatrixType;
	//typedef TransformType::OffsetType OffsetType;
	//
	//MatrixType initialMatrix;
	////OffsetType initialOffset;
	//
	//initialMatrix(0, 0) = 0.292642789523456;
	//initialMatrix(0, 1) = -0.0720624140477694;
	//initialMatrix(0, 2) = -0.953502598958985;
	//
	//initialMatrix(1, 0) = -0.956177858394543;
	//initialMatrix(1, 1) = -0.0124897827434311;
	//initialMatrix(1, 2) = -0.292519928283611;
	//
	//initialMatrix(2, 0) = -0.00917065188290264;	
	//initialMatrix(2, 1) = -0.997321920850349;
	//initialMatrix(2, 2) = 0.0725595295975165;
	//
	//
	//MatrixType rotationMatrix;
	//rotationMatrix.Fill(0);
	//rotationMatrix(0, 1) = 1;
	//rotationMatrix(1, 0) = 1;
	//rotationMatrix(2, 2) = 1;
	//MatrixType rotationMatrixTran = rotationMatrix.GetTranspose();
	//initialMatrix = rotationMatrixTran*initialMatrix*rotationMatrix;
	//
	//initialOffset[0] = 0.0287;
	//initialOffset[1] = 0.0491;
	//initialOffset[2] = 0.0407;
	//
	//std::cout << "InitialMatrix = " << std::endl << initialMatrix << std::endl;
	//
	//transform->SetMatrix(initialMatrix);
	//
	//transform->SetOffset(initialOffset);

	registration->SetInitialTransform(transform);
	registration->InPlaceOn();

	typedef OptimizerType::ScalesType       OptimizerScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());

	// default optimizer parameters
	double translationScale = 1.0 / 5.0;
	unsigned int numOfIterations = 300;
	double learningRate = 1.0;
	double miminumStepLength = 0.00001;

	// change optimizer parameters according to the arguments
	if (argc > 7){
		translationScale = std::atof(argv[7]);
	}
	if (argc > 8){
		numOfIterations = std::atoi(argv[8]);
	}
	if (argc > 9){
		learningRate = std::atof(argv[9]);
	}
	if (argc > 10){
		miminumStepLength = std::atof(argv[10]);
	}
	std::cout << "translationScale:" << translationScale << std::endl;
	std::cout << "numOfIterations:" << numOfIterations << std::endl;
	std::cout << "learningRate:" << learningRate << std::endl;
	std::cout << "miminumStepLength:" << miminumStepLength << std::endl;

	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizer->SetScales(optimizerScales);

	optimizer->SetNumberOfIterations(numOfIterations);
	optimizer->SetLearningRate(learningRate);
	optimizer->SetMinimumStepLength(miminumStepLength);

	optimizer->SetReturnBestParametersAndValue(true);

	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	const unsigned int numberOfLevels = 1;

	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(1);
	shrinkFactorsPerLevel[0] = 2;

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(1);
	smoothingSigmasPerLevel[0] = 0;

	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);

	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition: "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

#pragma region recording
	ofstream logFile;
	logFile.open("E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170619/log.txt", ios::app);

	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	logFile <<
		"========================================================\n" <<
		asctime(timeinfo) <<
		"Scan-Scan Registration:\n" <<
		"FixedImage: " << argv[1] << "\n" <<
		"MovingImage: " << argv[2] << "\n" <<
		"OutputImage: " << argv[5] << "\n" << "\n" <<


		"Rigid Transform parameters after registration: " << "\n" <<
		transform->GetParameters() << "\n" <<
		"translations: " << transform->GetTranslation() << "\n" <<
		"center: " << transform->GetCenter() << "\n" <<
		"matrix: " << transform->GetMatrix() <<
		"offset: " << transform->GetOffset() << "\n" <<
		"versor: " << transform->GetVersor() << "\n" <<
		"Last Step Length: " << optimizer->GetCurrentStepLength() << "\n" <<
		"Last Metric Vuale: " << metric->GetCurrentValue() << "\n" << endl;

	logFile.close();

	std::cout << "Rigid Transform parameters after registration: " << "\n" <<
		transform->GetParameters() << "\n" <<
		"translations: " << transform->GetTranslation() << "\n" <<
		"center: " << transform->GetCenter() << "\n" <<
		"matrix: " << transform->GetMatrix() <<
		"offset: " << transform->GetOffset() << "\n" <<
		"versor: " << transform->GetVersor() << "\n" <<
		"Last Step Length: " << optimizer->GetCurrentStepLength() << "\n" <<
		"Last Metric Vuale: " << metric->GetCurrentValue() << "\n" << endl;
#pragma endregion


	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();

	resampler->SetTransform(transform);
	resampler->SetInput(movingImageReader->GetOutput());

	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
	FixedImageType::Pointer movingImage = movingImageReader->GetOutput();

	FixedImageType::DirectionType fixedDirection = fixedImage->GetDirection();
	FixedImageType::DirectionType movingDirection = movingImage->GetDirection();

	resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resampler->SetOutputOrigin(fixedImage->GetOrigin());
	resampler->SetOutputSpacing(fixedImage->GetSpacing());
	resampler->SetOutputDirection(fixedDirection);
	resampler->SetDefaultPixelValue(100);

	typedef  double                                       OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension >                OutputImageType;
	typedef itk::CastImageFilter< FixedImageType, OutputImageType > CastFilterType;
	typedef itk::ImageFileWriter< OutputImageType >                 WriterType;

	typedef itk::SubtractImageFilter<
		FixedImageType,
		FixedImageType,
		FixedImageType > DifferenceFilterType;
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	typedef itk::RescaleIntensityImageFilter<
		FixedImageType,
		OutputImageType >   RescalerType;
	RescalerType::Pointer intensityRescaler = RescalerType::New();

	intensityRescaler->SetInput(difference->GetOutput());
	intensityRescaler->SetOutputMinimum(0);
	intensityRescaler->SetOutputMaximum(255);

	difference->SetInput1(fixedImageReader->GetOutput());
	difference->SetInput2(resampler->GetOutput());

	resampler->SetDefaultPixelValue(1);

	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(intensityRescaler->GetOutput());

	// Compute the difference image between the scans
	writer->SetFileName(argv[5]);
	writer->Update();

	//Sanketh code
	std::cout << "difference image written" << std::endl;

	// resample deformed atlas
	typedef itk::ResampleImageFilter<
		ResampleInputImageType,
		FixedImageType >    ResampleFilterType2;
	//Sanketh code
	std::cout << "debug 1" << std::endl;

	ResampleFilterType2::Pointer resampler2 = ResampleFilterType2::New();

	//Sanketh code
	std::cout << "debug 2" << std::endl;

	resampler2->SetTransform(transform);
	resampler2->SetInput(resampleInputImageReader->GetOutput());

	//Sanketh code
	std::cout << "debug 3" << std::endl;

	resampler2->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resampler2->SetOutputOrigin(fixedImage->GetOrigin());
	resampler2->SetOutputSpacing(fixedImage->GetSpacing());
	resampler2->SetOutputDirection(fixedDirection);
	resampler2->SetDefaultPixelValue(100);

	//Sanketh code
	std::cout << "debug 4" << std::endl;

	WriterType::Pointer      writer2 = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();
	
	//Sanketh code
	std::cout << "debug 5" << std::endl;

	//writer2->SetFileName(argv[4]);
	//Sanketh code
	std::cout << "debug 6" << std::endl;

	caster->SetInput(resampler2->GetOutput());

	writer2->SetInput(caster->GetOutput());
	writer2->SetFileName(argv[4]);

	//Sanketh code
	std::cout << "debug 7" << std::endl;

	//sanketh code
	try{
		std::cout << "writing bspline image: " << std::endl;
		writer2->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	//writer2->Update();

	//Sanketh code
	std::cout << "debug 8" << std::endl;

	//Sanketh code
	std::cout << "resample bspline image written" << std::endl;

	intensityRescaler->SetInput(difference->GetOutput());
	intensityRescaler->SetOutputMinimum(0);
	intensityRescaler->SetOutputMaximum(255);

	difference->SetInput1(fixedImageReader->GetOutput());
	difference->SetInput2(resampler2->GetOutput());

	resampler2->SetDefaultPixelValue(1);

	WriterType::Pointer writer3 = WriterType::New();
	writer3->SetInput(intensityRescaler->GetOutput());

	// Compute the difference image between the scans
	if (argc > 6)
	{
		writer3->SetFileName(argv[6]);
		writer3->Update();
		//Sanketh code
		std::cout << "atlas diff image written" << std::endl;
	}

	//resample labelmap
	typedef  unsigned char labelMapPixelType;
	typedef itk::Image< labelMapPixelType, Dimension > LabelMapImageType;

	typedef itk::ImageFileReader< LabelMapImageType  >       LabelMapImageReaderType;

	LabelMapImageReaderType::Pointer labelMapReader = LabelMapImageReaderType::New();

	typedef itk::ResampleImageFilter<
		LabelMapImageType,
		LabelMapImageType >    LabelMapResampleFilterType;

	LabelMapResampleFilterType::Pointer labelMapResampler = LabelMapResampleFilterType::New();

	labelMapResampler->SetTransform(transform);
	labelMapResampler->SetInput(labelMapReader->GetOutput());
	labelMapResampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	labelMapResampler->SetOutputOrigin(fixedImage->GetOrigin());
	labelMapResampler->SetOutputSpacing(fixedImage->GetSpacing());
	labelMapResampler->SetOutputDirection(fixedImage->GetDirection());
	labelMapResampler->SetDefaultPixelValue(0);

	typedef itk::Image< labelMapPixelType, Dimension >                labelMapImageType;
	typedef itk::CastImageFilter< labelMapImageType, labelMapImageType > labelMapCastFilterType;
	typedef itk::ImageFileWriter< labelMapImageType >                 labelMapWriterType;

	labelMapWriterType::Pointer      labelMapWriter = labelMapWriterType::New();
	labelMapCastFilterType::Pointer  labelMapCaster = labelMapCastFilterType::New();


	labelMapCaster->SetInput(labelMapResampler->GetOutput());
	labelMapWriter->SetInput(labelMapCaster->GetOutput());


	for (size_t j = 0; j < 4; j++)
	{
		labelMapReader->SetFileName(argv[11 + 2 * j]);
		labelMapWriter->SetFileName(argv[12 + 2 * j]);
		try{
			cout << "writing label map: " << argv[12 + 2 * j] << endl;
			labelMapWriter->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cout << "ExceptionObject caught !" << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}

	}
	//Sanketh code
	std::cout << "labelmap image written" << std::endl;
	return EXIT_SUCCESS;


}
