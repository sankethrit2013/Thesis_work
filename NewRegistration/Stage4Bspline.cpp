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
// Software Guide : BeginLatex
//
// This example illustrates the use of \code{SpatialObject}s as masks for selecting the
// pixels that should contribute to the computation of Image Metrics. This
// example is almost identical to ImageRegistration6 with the exception that
// the \code{SpatialObject} masks are created and passed to the image metric.
//
//
// Software Guide : EndLatex
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizer.h"

//#include "itkConjugateGradientLineSearchOptimizerv4.h"
#include "itkRegistrationParameterScalesFromIndexShift.h"
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCenteredTransformInitializer.h"
#include <itkMutualInformationImageToImageMetric.h>
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkTranslationTransform.h"
#include "itkImageMomentsCalculator.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCommand.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBSplineTransform.h"
#include "itkBSplineTransformInitializer.h"


#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;
//  Software Guide : BeginLatex
//
//  The most important header in this example is the one corresponding to the
//  \doxygen{ImageMaskSpatialObject} class.
//
//  \index{itk::ImageMaskSpatialObject!header}
//
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
#include "itkImageMaskSpatialObject.h"
// Software Guide : EndCodeSnippet
//
//  The following section of code implements a command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
	typedef  CommandIterationUpdate   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>   Pointer;
	itkNewMacro(Self);
protected:
	CommandIterationUpdate() {};
public:
	typedef itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
	typedef const OptimizerType *                               OptimizerPointer;
	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *)caller, event);
	}
	void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >(object);
		if (!itk::IterationEvent().CheckEvent(&event))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << std::endl;
	}
};
int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " fixedImageFile  movingImageFile outputImagefile ";
		std::cerr << "  [differenceOutputfile] [differenceBeforeRegistration] ";
		std::cerr << " [numberOfGridNodesInOneDimension] [deformationField] ";
		return EXIT_FAILURE;
	}
	
	const    unsigned int    Dimension = 3;
	typedef  double           PixelType;

	const PixelType backgroundGrayLevel = 80;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;

	const unsigned int SpaceDimension = Dimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
		CoordinateRepType,
		SpaceDimension,
		SplineOrder >           BSplineTransformType;
	BSplineTransformType::Pointer BSplineTransform = BSplineTransformType::New();


	typedef itk::CorrelationImageToImageMetricv4<
		FixedImageType,
		MovingImageType >   MetricType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType >     RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();


	registration->SetMetric(metric);

	typedef itk::ImageFileReader< FixedImageType  >       FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType >       MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(argv[1]);
	movingImageReader->SetFileName(argv[2]);

	fixedImageReader->Update();
	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

	unsigned int numberOfGridNodesInOneDimension = 8;
	if (argc >= 7)
	{
		numberOfGridNodesInOneDimension = std::atof(argv[7]);
	}

	typedef itk::BSplineTransformInitializer<
		BSplineTransformType,
		FixedImageType> InitializerType;
	InitializerType::Pointer transformInitializer = InitializerType::New();

	BSplineTransformType::MeshSizeType  meshSize;
	meshSize.Fill(numberOfGridNodesInOneDimension - SplineOrder);
	transformInitializer->SetTransform(BSplineTransform);
	transformInitializer->SetImage(fixedImage);
	transformInitializer->SetTransformDomainMeshSize(meshSize);
	transformInitializer->InitializeTransform();
	BSplineTransform->SetIdentity();

	registration->SetInitialTransform(BSplineTransform);
	registration->InPlaceOn();


	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImageReader->GetOutput());


	unsigned int numberofIterations = 10;
	if (argc >= 8)
	{
		numberofIterations = std::atof(argv[8]);
	}

	typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	registration->SetOptimizer(optimizer);
	optimizer->SetMinimumStepLength(1e-6);
	optimizer->SetRelaxationFactor(0.6);//0.6
	optimizer->SetNumberOfIterations(numberofIterations);
	optimizer->SetLearningRate(0.1); //0.5
	optimizer->SetGradientMagnitudeTolerance(1e-6);
	optimizer->SetReturnBestParametersAndValue(true);

	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	//  Here we instantiate the type of the \doxygen{ImageMaskSpatialObject}
	//  using the same dimension of the images to be registered.
	//
	//  \index{itk::ImageMaskSpatialObject!Instantiation}

	typedef itk::ImageMaskSpatialObject< Dimension >   MaskType;

	//  Then we use the type for creating the spatial object mask that will
	//  restrict the registration to a reduced region of the image.
	//
	//  \index{itk::ImageMaskSpatialObject!New}
	//  \index{itk::ImageMaskSpatialObject!Pointer}


	MaskType::Pointer  spatialObjectMask = MaskType::New();

	//  The mask in this case is read from a binary file using the
	//  \code{ImageFileReader} instantiated for an \code{unsigned char} pixel
	//  type.

	typedef itk::Image< unsigned char, Dimension >   ImageMaskType;
	typedef itk::ImageFileReader< ImageMaskType >    MaskReaderType;

	//  The reader is constructed and a filename is passed to it.

	MaskReaderType::Pointer  maskReader = MaskReaderType::New();
	maskReader->SetFileName(argv[6]);

	//  As usual, the reader is triggered by invoking its \code{Update()} method.
	//  Since this may eventually throw an exception, the call must be placed in
	//  a \code{try/catch} block. Note that a full fledged application will place
	//  this \code{try/catch} block at a much higher level, probably under the
	//  control of the GUI.

	try
	{
		maskReader->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	//  The output of the mask reader is connected as input to the
	//  \code{ImageMaskSpatialObject}.
	//
	//  \index{itk::ImageMaskSpatialObject!SetImage()}

	spatialObjectMask->SetImage(maskReader->GetOutput());

	//  Finally, the spatial object mask is passed to the image metric.
	//
	//  \index{itk::ImageToImageMetricv4!SetFixedImageMask()}

	metric->SetFixedImageMask(spatialObjectMask);


	typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
	ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
	scalesEstimator->SetMetric(metric);
	scalesEstimator->SetTransformForward(true);
	scalesEstimator->SetSmallParameterVariation(1.0);
	optimizer->SetScalesEstimator(scalesEstimator);
	//optimizer->EstimateLearningRate();

	//sanketh code, change number of levels to 4
	//const unsigned int numberOfLevels = 3;
	const unsigned int numberOfLevels = 1;

	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(numberOfLevels);
	/*---Multilevel---*/
	//shrinkFactorsPerLevel[0] = 3;
	//shrinkFactorsPerLevel[1] = 2;
	//shrinkFactorsPerLevel[2] = 1;
	////shrinkFactorsPerLevel[3] = 1;
	/*----------------*/
	shrinkFactorsPerLevel[0] = 2;

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(numberOfLevels);
	smoothingSigmasPerLevel[0] = 0;
	/*---Multilevel---*/
	//smoothingSigmasPerLevel[0] = 2;
	//smoothingSigmasPerLevel[1] = 1;
	//smoothingSigmasPerLevel[2] = 0;
	////smoothingSigmasPerLevel[3] = 0;
	/*----------------*/

	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);

	
	//Add time and memory probes
	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	cout << "This prints" << endl;

	std::cout << std::endl << "Starting Registration" << std::endl;

	try
	{
		memorymeter.Start("Registration");
		chronometer.Start("Registration");

		registration->Update();

		chronometer.Stop("Registration");
		memorymeter.Stop("Registration");

		std::cout << "Optimizer stop condition = "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// Report the time and memory taken by the registration
	chronometer.Report(std::cout);
	memorymeter.Report(std::cout);

	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(BSplineTransform);
	resample->SetInput(movingImageReader->GetOutput());

	resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedImage->GetOrigin());
	resample->SetOutputSpacing(fixedImage->GetSpacing());
	resample->SetOutputDirection(fixedImage->GetDirection());
	resample->SetDefaultPixelValue(128);

	typedef  unsigned char  OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	typedef itk::CastImageFilter<
		FixedImageType,
		OutputImageType > CastFilterType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();

	writer->SetFileName(argv[3]);
	caster->SetInput(resample->GetOutput());
	writer->SetInput(caster->GetOutput());


	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

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
	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput(intensityRescaler->GetOutput());

	//Compute the difference image between the fixed and resampled moving image.
	if (argc >= 5)
	{
		difference->SetInput1(fixedImageReader->GetOutput());
		difference->SetInput2(resample->GetOutput());

		writer2->SetFileName(argv[4]);
		try
		{
			writer2->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}

	
	// //resample the lable map

	typedef  unsigned char labelMapPixelType;
	typedef itk::Image< labelMapPixelType, Dimension > LabelMapImageType;

	typedef itk::ImageFileReader< LabelMapImageType  >       LabelMapImageReaderType;

	LabelMapImageReaderType::Pointer labelMapReader = LabelMapImageReaderType::New();

	typedef itk::ResampleImageFilter<
		LabelMapImageType,
		LabelMapImageType >    LabelMapResampleFilterType;

	LabelMapResampleFilterType::Pointer labelMapResampler = LabelMapResampleFilterType::New();

	labelMapResampler->SetTransform(BSplineTransform);
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


	for (size_t j = 0; j < 5; j++)
	{
		labelMapReader->SetFileName(argv[10 + 2 * j]);
		labelMapWriter->SetFileName(argv[11 + 2 * j]);
		try {
			cout << "writing label map: " << argv[11 + 2 * j] << endl;
			labelMapWriter->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cout << "ExceptionObject caught !" << std::endl;
			std::cout << err << std::endl;
			return EXIT_FAILURE;
		}

	}
	return EXIT_SUCCESS;
}
