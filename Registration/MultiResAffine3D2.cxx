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

// MultiResoluition registration: Affine, MattesMutualInformation, RegularStepGradientDescentOptimizer
#include "itkImageRegistrationMethodv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


#include "itkCommand.h"
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command
{

public:
	typedef  RegistrationInterfaceCommand   Self;
	typedef  itk::Command                   Superclass;
	typedef  itk::SmartPointer<Self>        Pointer;
	itkNewMacro(Self);

protected:
	RegistrationInterfaceCommand() {};

public:
	typedef   TRegistration      RegistrationType;
	typedef   RegistrationType * RegistrationPointer;
	typedef   itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
	typedef   OptimizerType * OptimizerPointer;
	void Execute(itk::Object * object,
		const itk::EventObject & event) ITK_OVERRIDE
	{
		if (!(itk::MultiResolutionIterationEvent().CheckEvent(&event)))
		{
			return;
		}

		RegistrationPointer registration =
			static_cast<RegistrationPointer>(object);
		OptimizerPointer optimizer = static_cast< OptimizerPointer >(
			registration->GetModifiableOptimizer());

		unsigned int currentLevel = registration->GetCurrentLevel();
		typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors =
			registration->GetShrinkFactorsPerDimension(currentLevel);
		typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas =
			registration->GetSmoothingSigmasPerLevel();

		std::cout << "-------------------------------------" << std::endl;
		std::cout << " Current level = " << currentLevel << std::endl;
		std::cout << "    shrink factor = " << shrinkFactors << std::endl;
		std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
		std::cout << std::endl;

		if (registration->GetCurrentLevel() == 0)
		{
			optimizer->SetLearningRate(1.00);
			optimizer->SetMinimumStepLength(0.0001);
		}
		else
		{
			optimizer->SetLearningRate(optimizer->GetCurrentStepLength());
			optimizer->SetMinimumStepLength(
				optimizer->GetMinimumStepLength() * 0.2);
		}
	}

		void Execute(const itk::Object *, const itk::EventObject &) ITK_OVERRIDE
	{
		return;
	}
};

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
	typedef   itk::RegularStepGradientDescentOptimizerv4<double>  OptimizerType;
	typedef   const OptimizerType *                               OptimizerPointer;

	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *)caller, event);
	}

		void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >(object);
		if (!(itk::IterationEvent().CheckEvent(&event)))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetValue() << "   ";
		std::cout << optimizer->GetCurrentPosition() << "   ";
		std::cout << m_CumulativeIterationIndex++ << std::endl;
	}
private:
	unsigned int m_CumulativeIterationIndex;
};



int main(int argc, char *argv[])
{


	if (argc < 4)
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << "   fixedImageFile  movingImageFile " << std::endl;
		std::cerr << "   outputImagefile  [differenceBeforeRegistration] " << std::endl;
		std::cerr << "   [differenceAfterRegistration] " << std::endl;
		std::cerr << "   [maxNumberOfIterations] [translationScale] " << std::endl;
		return EXIT_FAILURE;
	}

	const std::string fixedImageFile = argv[1];
	const std::string movingImageFile = argv[2];
	const std::string outImagefile = argv[3];
	const std::string differenceBeforeRegistration = (argc >4) ? argv[4] : "";
	const std::string differenceAfterRegistration = (argc >5) ? argv[5] : "";
	const int maxNumberOfIterations = (argc >6) ? atoi(argv[6]) : 200;
	const double translationScale = (argc >7) ? atof(argv[7]) : 0.2;
	const int numberOfLevels = (argc >8) ? atof(argv[8]) : 3;
	const int numberOfBins = (argc >9) ? atoi(argv[9]) : 0;
	const double relaxingFactor = 0.8;
	const std::string labelMapInputFile1 = (argc >10) ? argv[10] : "";
	const std::string labelMapOutputFile1 = (argc >11) ? argv[11] : "";
	const std::string labelMapInputFile2 = (argc >12) ? argv[12] : "";
	const std::string labelMapOutputFile2 = (argc >13) ? argv[13] : "";
	const std::string labelMapInputFile3 = (argc >14) ? argv[14] : "";
	const std::string labelMapOutputFile3 = (argc >15) ? argv[15] : "";
	const std::string labelMapInputFile4 = (argc >16) ? argv[16] : "";
	const std::string labelMapOutputFile4 = (argc >17) ? argv[17] : "";


	const    int    Dimension = 3;
	typedef   unsigned char  PixelType;

	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;

	typedef itk::AffineTransform< double, Dimension  > TransformType;


	typedef itk::RegularStepGradientDescentOptimizerv4<double>       OptimizerType;

	typedef itk::MattesMutualInformationImageToImageMetricv4<
		FixedImageType,
		MovingImageType > MetricType;

	//typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<
	//	FixedImageType,	MovingImageType > MetricType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType,
		TransformType >                   RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);

	metric->SetNumberOfHistogramBins(numberOfBins);
	/*MetricType::HistogramType::SizeType histogramSize;
	histogramSize.SetSize(2);
	histogramSize[0] = numberOfBins;
	histogramSize[1] = numberOfBins;
	metric->SetHistogramSize(histogramSize);*/

	TransformType::Pointer  transform = TransformType::New();

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	fixedImageReader->SetFileName(fixedImageFile);
	movingImageReader->SetFileName(movingImageFile);

	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());


	//  Software Guide : BeginLatex
	//
	//  In this example, we again use the
	//  \doxygen{CenteredTransformInitializer} helper class in order to compute
	//  reasonable values for the initial center of rotation and the
	//  translations. The initializer is set to use the center of mass of each
	//  image as the initial correspondence correction.
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	typedef itk::CenteredTransformInitializer<
		TransformType,
		FixedImageType,
		MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer
		= TransformInitializerType::New();
	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	initializer->MomentsOn();
	initializer->InitializeTransform();
	
	//transform->SetIdentity();
	registration->SetInitialTransform(transform);
	registration->InPlaceOn();

	

	typedef OptimizerType::ScalesType       OptimizerScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());

	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = 1.0;
	optimizerScales[4] = 1.0;
	optimizerScales[5] = 1.0;
	optimizerScales[6] = 1.0;
	optimizerScales[7] = 1.0;
	optimizerScales[8] = 1.0;

	optimizerScales[9] = translationScale;
	optimizerScales[10] = translationScale;
	optimizerScales[11] = translationScale;


	optimizer->SetScales(optimizerScales);
	


	//optimizer->SetMinimumStepLength(0.0001);
	optimizer->SetNumberOfIterations(maxNumberOfIterations);
	optimizer->ReturnBestParametersAndValueOn();
	//std::cout << "translationScale:" << translationScale << std::endl;
	//std::cout << "numOfIterations:" << maxNumberOfIterations << std::endl;
	//std::cout << "learningRate:" << steplength << std::endl;
	//std::cout << "miminumStepLength:" << 0.0001 << std::endl;
	optimizer->SetRelaxationFactor(relaxingFactor);
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);


	// One level registration process without shrinking and smoothing.
	//

	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(3);
	shrinkFactorsPerLevel[0] = 2;
	shrinkFactorsPerLevel[1] = 2;
	shrinkFactorsPerLevel[2] = 1;

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(3);
	smoothingSigmasPerLevel[0] = 0;
	smoothingSigmasPerLevel[1] = 0;
	smoothingSigmasPerLevel[2] = 0;

	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);

	typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
	CommandType::Pointer command = CommandType::New();

	registration->AddObserver(itk::MultiResolutionIterationEvent(), command);

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
	// Software Guide : EndCodeSnippet


	//  Software Guide : BeginLatex
	//
	//  Once the optimization converges, we recover the parameters from the
	//  registration method. We can also recover the
	//  final value of the metric with the \code{GetValue()} method and the
	//  final number of iterations with the \code{GetCurrentIteration()}
	//  method.
	//
	//  \index{itk::RegistrationMethodv4!GetValue()}
	//  \index{itk::RegistrationMethodv4!GetCurrentIteration()}
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	const TransformType::ParametersType finalParameters =
		registration->GetOutput()->Get()->GetParameters();

	const double finalRotationCenterX = transform->GetCenter()[0];
	const double finalRotationCenterY = transform->GetCenter()[1];
	const double finalRotationCenterZ = transform->GetCenter()[2];

	const double finalTranslationX = finalParameters[9];
	const double finalTranslationY = finalParameters[10];
	const double finalTranslationZ = finalParameters[11];

	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();
	// Software Guide : EndCodeSnippet


	// Print out results
	//
	std::cout << "Result = " << std::endl;
	std::cout << " Center X      = " << finalRotationCenterX << std::endl;
	std::cout << " Center Y      = " << finalRotationCenterY << std::endl;
	std::cout << " Center Z     = " << finalRotationCenterZ << std::endl;
	std::cout << " Translation X = " << finalTranslationX << std::endl;
	std::cout << " Translation Y = " << finalTranslationY << std::endl;
	std::cout << " Translation Z = " << finalTranslationZ << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue << std::endl;

	//Compute the rotation angle and scaling from SVD of the matrix
	// \todo Find a way to figure out if the scales are along X or along Y.
	// VNL returns the eigenvalues ordered from largest to smallest.

	//vnl_matrix<double> p(3, 3);
	//p[0][0] = (double)finalParameters[0];
	//p[0][1] = (double)finalParameters[1];
	//p[0][2] = (double)finalParameters[2];
	//p[1][0] = (double)finalParameters[3];
	//p[1][1] = (double)finalParameters[4];
	//p[1][2] = (double)finalParameters[5];
	//p[2][0] = (double)finalParameters[6];
	//p[2][1] = (double)finalParameters[7];
	//p[2][2] = (double)finalParameters[8];

	//vnl_svd<double> svd(p);
	//vnl_matrix<double> r(3, 3);
	//r = svd.U() * vnl_transpose(svd.V());
	//double angle = std::asin(r[1][0]);

	//const double angleInDegrees = angle * 180.0 / vnl_math::pi;

	//std::cout << " Scale 1         = " << svd.W(0) << std::endl;
	//std::cout << " Scale 2         = " << svd.W(1) << std::endl;
	//std::cout << " Angle (degrees) = " << angleInDegrees << std::endl;



	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	//TransformType::ParametersType parameters(Dimension * Dimension + Dimension);
	//parameters[1] = 0.7883421180709784;
	//parameters[2] = -0.11161831859622057;
	//parameters[3] = -0.4044797227228365;
	//parameters[4] = 0.184849923280962;
	//parameters[5] = 0.9028515840651747;
	//parameters[6] = 0.17379856465734245;
	//parameters[7] = 0.3691846493706708;
	//parameters[8] = -0.23667646154851427;
	//parameters[9] = 0.8139244516196104;
	//parameters[10] = -5.278204448242631;
	//parameters[11] = 2.9914186603661292;
	//parameters[12] = 1.3206958229901822;
	//transform->SetParameters(parameters);
	resampler->SetTransform(transform);
	resampler->SetInput(movingImageReader->GetOutput());

	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

	resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resampler->SetOutputOrigin(fixedImage->GetOrigin());
	resampler->SetOutputSpacing(fixedImage->GetSpacing());
	resampler->SetOutputDirection(fixedImage->GetDirection());
	resampler->SetDefaultPixelValue(80);

	typedef  unsigned char  OutputPixelType;

	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

	typedef itk::CastImageFilter<
		FixedImageType,
		OutputImageType > CastFilterType;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;


	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();


	writer->SetFileName(outImagefile);


	caster->SetInput(resampler->GetOutput());
	writer->SetInput(caster->GetOutput());
	writer->Update();

	//resample the lable map
	MovingImageReaderType::Pointer labelMapReader = MovingImageReaderType::New();

	ResampleFilterType::Pointer labelMapResampler = ResampleFilterType::New();

	labelMapResampler->SetTransform(transform);
	labelMapResampler->SetInput(labelMapReader->GetOutput());
	labelMapResampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	labelMapResampler->SetOutputOrigin(fixedImage->GetOrigin());
	labelMapResampler->SetOutputSpacing(fixedImage->GetSpacing());
	labelMapResampler->SetOutputDirection(fixedImage->GetDirection());
	labelMapResampler->SetDefaultPixelValue(0);

	typedef  float	labelMapPixelType;
	typedef itk::Image< labelMapPixelType, Dimension >                labelMapImageType;
	typedef itk::CastImageFilter< FixedImageType, labelMapImageType > labelMapCastFilterType;
	typedef itk::ImageFileWriter< labelMapImageType >                 labelMapWriterType;

	labelMapWriterType::Pointer      labelMapWriter = labelMapWriterType::New();
	labelMapCastFilterType::Pointer  labelMapCaster = labelMapCastFilterType::New();

	labelMapReader->SetFileName(labelMapInputFile1);
	labelMapWriter->SetFileName(labelMapOutputFile1);
	labelMapCaster->SetInput(labelMapResampler->GetOutput());
	labelMapWriter->SetInput(labelMapCaster->GetOutput());
	if (labelMapInputFile1 != "" & labelMapOutputFile1 != "")
	{
		labelMapWriter->Update();
	}
	labelMapReader->SetFileName(labelMapInputFile2);
	labelMapWriter->SetFileName(labelMapOutputFile2);
	if (labelMapInputFile2 != "" & labelMapOutputFile2 != "")
	{
		labelMapWriter->Update();
	}
	labelMapReader->SetFileName(labelMapInputFile3);
	labelMapWriter->SetFileName(labelMapOutputFile3);
	if (labelMapInputFile3 != "" & labelMapOutputFile3 != "")
	{
		labelMapWriter->Update();
	}
	labelMapReader->SetFileName(labelMapInputFile4);
	labelMapWriter->SetFileName(labelMapOutputFile4);
	if (labelMapInputFile4 != "" & labelMapOutputFile4 != "")
	{
		labelMapWriter->Update();
	}
	typedef itk::SubtractImageFilter<
		FixedImageType,
		FixedImageType,
		FixedImageType > DifferenceFilterType;

	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	difference->SetInput1(fixedImageReader->GetOutput());
	difference->SetInput2(resampler->GetOutput());

	WriterType::Pointer writer2 = WriterType::New();

	typedef itk::RescaleIntensityImageFilter<
		FixedImageType,
		OutputImageType >   RescalerType;

	RescalerType::Pointer intensityRescaler = RescalerType::New();

	intensityRescaler->SetInput(difference->GetOutput());
	intensityRescaler->SetOutputMinimum(0);
	intensityRescaler->SetOutputMaximum(255);

	writer2->SetInput(intensityRescaler->GetOutput());
	resampler->SetDefaultPixelValue(1);

	// Compute the difference image between the
	// fixed and resampled moving image.
	if (argc > 5)
	{
		writer2->SetFileName(differenceAfterRegistration);
		writer2->Update();
	}


	typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
	IdentityTransformType::Pointer identity = IdentityTransformType::New();

	// Compute the difference image between the
	// fixed and moving image before registration.
	if (argc > 4)
	{
		resampler->SetTransform(identity);
		writer2->SetFileName(differenceBeforeRegistration);
		writer2->Update();
	}

	return EXIT_SUCCESS;
}
