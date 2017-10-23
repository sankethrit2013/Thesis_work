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

#include "itkImageMomentsCalculator.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkCommand.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

//  The following section of code implements a Command observer
//  that will monitor the configurations of the registration process
//  at every change of stage and resolution level.
//
template <typename SRegistration>
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
	typedef   SRegistration                          RegistrationType;

	// The Execute function simply calls another version of the \code{Execute()}
	// method accepting a \code{const} input object
	void Execute(itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *) object, event);
	}

		void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		if (!(itk::MultiResolutionIterationEvent().CheckEvent(&event)))
		{
			return;
		}

		std::cout << "\nObserving from class " << object->GetNameOfClass();
		if (!object->GetObjectName().empty())
		{
			std::cout << " \"" << object->GetObjectName() << "\"" << std::endl;
		}

		const RegistrationType * registration = static_cast<const RegistrationType *>(object);
		if (registration == 0)
		{
			itkExceptionMacro(<< "Dynamic cast failed, object of type " << object->GetNameOfClass());
		}

		unsigned int currentLevel = registration->GetCurrentLevel();
		typename RegistrationType::ShrinkFactorsPerDimensionContainerType shrinkFactors =
			registration->GetShrinkFactorsPerDimension(currentLevel);
		typename RegistrationType::SmoothingSigmasArrayType smoothingSigmas =
			registration->GetSmoothingSigmasPerLevel();

		std::cout << "-------------------------------------" << std::endl;
		std::cout << " Current multi-resolution level = " << currentLevel << std::endl;
		std::cout << "    shrink factor = " << shrinkFactors << std::endl;
		std::cout << "    smoothing sigma = " << smoothingSigmas[currentLevel] << std::endl;
		std::cout << std::endl;
	}
};

//  The following section of code implements an observer
//  that will monitor the evolution of the registration process.
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
		std::cout << optimizer->GetCurrentPosition() << "  " <<
			m_CumulativeIterationIndex++ << std::endl;
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
		std::cerr << " fixedImageFile  movingImageFile ";
		std::cerr << " outputImagefile [backgroundGrayLevel]";
		std::cerr << " [checkerboardbefore] [CheckerBoardAfter]";
		std::cerr << " [numberOfBins] " << std::endl;
		return EXIT_FAILURE;
	}

	const std::string fixedImageFile = argv[1];
	const std::string movingImageFile = argv[2];
	const std::string outImagefile = argv[3];
	const std::string differenceBeforeRegistration = (argc >4) ? argv[4] : "";
	const std::string differenceAfterRegistration = (argc > 5) ? argv[5] : "";
	const int NumberOfIterations1 = (argc > 6) ? atoi(argv[6]) : 200;
	const double translationScale = (argc > 7) ? atof(argv[7]) : 0.2;
	const double rescalefactor = (argc > 8) ? atof(argv[8]) : 2;
	const int NumberOfIterations2 = (argc > 9) ? atoi(argv[9]) : 0;
	//const int numberOfLevels = (argc >8) ? atof(argv[8]) : 3;
	const int numberOfHistogramBins = 64;
	const int NumberOfIterations3 = (argc > 10) ? atoi(argv[10]) : 0;
	const double translationScale2 = (argc > 11) ? atof(argv[11]) : 0.2;
	const double rescalefactor2 = (argc > 12) ? atof(argv[12]) : 2;
	const std::string scan2ImageFile = argv[13];

	
	//const int nn = 9;//number of argument except labelmap names
	//const std::string labelMapInputFile1 = (argc > nn + 1) ? argv[nn + 1] : "";
	//const std::string labelMapOutputFile1 = (argc > nn + 2) ? argv[nn + 2] : "";
	//const std::string labelMapInputFile2 = (argc > nn + 3) ? argv[nn + 3] : "";
	//const std::string labelMapOutputFile2 = (argc > nn + 4) ? argv[nn + 4] : "";
	//const std::string labelMapInputFile3 = (argc > nn + 5) ? argv[nn + 5] : "";
	//const std::string labelMapOutputFile3 = (argc > nn + 6) ? argv[nn + 6] : "";
	//const std::string labelMapInputFile4 = (argc > nn + 7) ? argv[nn + 7] : "";
	//const std::string labelMapOutputFile4 = (argc > nn + 8) ? argv[nn + 8] : "";

	std::cout << "Started..." << std::endl;
	const    unsigned int    Dimension = 3;
	typedef  float           PixelType;

	const PixelType backgroundGrayLevel = 80;


	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;

	typedef itk::Similarity3DTransform<double>     STransformType;
	typedef itk::RegularStepGradientDescentOptimizerv4<double> SOptimizerType;
	//typedef itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType >  MetricType;
	typedef itk::CorrelationImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
	//typedef itk::MutualInformationImageToImageMetric<FixedImageType, MovingImageType> MetricType;
	//typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<FixedImageType,MovingImageType> MetricType;


	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(fixedImageFile);
	movingImageReader->SetFileName(movingImageFile);

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType >  SRegistrationType;
	
	/*----------------------------------------------------
	registration between scans
	----------------------------------------------------*/
	
	SOptimizerType::Pointer      s2sOptimizer = SOptimizerType::New();
	MetricType::Pointer         s2sMetric = MetricType::New();
	SRegistrationType::Pointer   s2sRegistration = SRegistrationType::New();
	STransformType::Pointer s2sTx = STransformType::New();

	s2sRegistration->SetOptimizer(s2sOptimizer);
	s2sRegistration->SetMetric(s2sMetric);


	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	FixedImageReaderType::Pointer  scan2ImageReader = FixedImageReaderType::New();
	scan2ImageReader->SetFileName(scan2ImageFile);
	scan2ImageReader->Update();
	FixedImageType::Pointer scan2Image = scan2ImageReader->GetOutput();

	typedef itk::CenteredTransformInitializer<STransformType, FixedImageType, FixedImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer2
		= TransformInitializerType::New();

	initializer2->SetTransform(s2sTx);
	initializer2->SetFixedImage(scan2ImageReader->GetOutput());
	initializer2->SetMovingImage(fixedImageReader->GetOutput());
	initializer2->MomentsOn();
	initializer2->InitializeTransform();
	

	s2sRegistration->SetInitialTransform(s2sTx);
	s2sRegistration->InPlaceOn();

	s2sRegistration->SetFixedImage(scan2ImageReader->GetOutput());
	s2sRegistration->SetMovingImage(fixedImageReader->GetOutput());
	s2sRegistration->SetObjectName("subject2subjectRegistration");

	typedef itk::RegistrationParameterScalesFromPhysicalShift<
		MetricType> ScalesEstimatorType;

	const unsigned int numberOfLevels3 = 1;

	SRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel3;
	shrinkFactorsPerLevel3.SetSize(numberOfLevels3);
	shrinkFactorsPerLevel3[0] = 2;
	//shrinkFactorsPerLevel1[1] = 1;

	SRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel3;
	smoothingSigmasPerLevel3.SetSize(numberOfLevels3);
	smoothingSigmasPerLevel3[0] = 0;
	//smoothingSigmasPerLevel1[1] = 0;

	s2sRegistration->SetNumberOfLevels(numberOfLevels3);
	s2sRegistration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel3);
	s2sRegistration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel3);

	s2sOptimizer->SetNumberOfIterations(NumberOfIterations3);
	//s2sOptimizer->SetRelaxationFactor(0.5);
	s2sOptimizer->SetLearningRate(0.8);
	//similarityOptimizer->SetDoEstimateLearningRateOnce(true);
	s2sOptimizer->SetDoEstimateScales(true);
	s2sOptimizer->SetMinimumStepLength(0.00001);
	s2sOptimizer->SetGradientMagnitudeTolerance(0.00001);
	//similarityOptimizer->SetReturnBestParametersAndValue(true);
	//Set optimizer scales(movestep = 1/scaling factor * OptimizerStepLenght)
	typedef SOptimizerType::ScalesType       SOptimizerScalesType;
	SOptimizerScalesType optimizerScales2(s2sTx->GetNumberOfParameters());
	optimizerScales2[0] = 1.0;
	optimizerScales2[1] = 1.0;
	optimizerScales2[2] = 1.0;
	optimizerScales2[3] = translationScale2;
	optimizerScales2[4] = translationScale2;
	optimizerScales2[5] = translationScale2;
	optimizerScales2[6] = rescalefactor2;
	s2sOptimizer->SetScales(optimizerScales2);

	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer3 = CommandIterationUpdate::New();
	s2sOptimizer->AddObserver(itk::IterationEvent(), observer3);

	// Create the Command interface observer and register it with the optimizer.
	//
	typedef RegistrationInterfaceCommand<SRegistrationType> SimilarityCommandType;
	SimilarityCommandType::Pointer command3 = SimilarityCommandType::New();
	s2sRegistration->AddObserver(itk::MultiResolutionIterationEvent(), command3);

	try
	{
		s2sRegistration->Update();
		std::cout << "Optimizer stop condition: "
			<< s2sRegistration->
			GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	/*--------------------------------
		Subject-atlas Registration: Similarity Transform Stage
	---------------------------------*/

	SOptimizerType::Pointer      similarityOptimizer = SOptimizerType::New();
	MetricType::Pointer         similarityMetric = MetricType::New();
	SRegistrationType::Pointer   similarityRegistration = SRegistrationType::New();

	similarityRegistration->SetOptimizer(similarityOptimizer);
	similarityRegistration->SetMetric(similarityMetric);

	STransformType::Pointer similarityTx = STransformType::New();



	typedef itk::CenteredTransformInitializer<STransformType, FixedImageType, MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer
		= TransformInitializerType::New();

	initializer->SetTransform(similarityTx);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	initializer->MomentsOn();
	initializer->InitializeTransform();
	
	//STransformType::ParametersType sparameters = similarityTx->GetParameters();
	//std::cout << "\nInitial transform" << sparameters << std::endl;
	//sparameters[3] = sparameters[3] + 0.9;
	//sparameters[4] = sparameters[4] + 0.7;

	//sparameters[5] = sparameters[5] + 0.9;
	//sparameters[6] = 0.9;
	////sparameters[5] = sparameters[0] + 1;
	//similarityTx->SetParameters(sparameters);
	//std::cout << "\nModified Initial transform" << similarityTx->GetParameters() << std::endl;
	
	similarityRegistration->SetInitialTransform(similarityTx);
	similarityRegistration->InPlaceOn();

	similarityRegistration->SetFixedImage(fixedImageReader->GetOutput());
	similarityRegistration->SetMovingImage(movingImageReader->GetOutput());
	similarityRegistration->SetObjectName("SimilarityRegistration");

	// Set ScaleEstimator
	typedef itk::RegistrationParameterScalesFromPhysicalShift<
		MetricType> ScalesEstimatorType;
	//ScalesEstimatorType::Pointer scalesEstimator1 = ScalesEstimatorType::New();

	//scalesEstimator1->SetMetric(similarityMetric);
	//scalesEstimator1->SetTransformForward(true);
	//scalesEstimator1->SetScalesSamplingStrategy();
	////scalesEstimator1->SetSmallParameterVariation(0.1);
	//similarityOptimizer->SetScalesEstimator(scalesEstimator1);
	//similarityOptimizer->SetDoEstimateLearningRateOnce(true);
	//similarityOptimizer->SetDoEstimateLearningRateAtEachIteration(false);
	//std::cout << "scales:"<<similarityOptimizer->GetScales() << std::endl;

	const unsigned int numberOfLevels1 = 1;

	SRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel1;
	shrinkFactorsPerLevel1.SetSize(numberOfLevels1);
	shrinkFactorsPerLevel1[0] = 2;
	//shrinkFactorsPerLevel1[1] = 2;

	SRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel1;
	smoothingSigmasPerLevel1.SetSize(numberOfLevels1);
	smoothingSigmasPerLevel1[0] = 1; //changed from 0 to 1
	//smoothingSigmasPerLevel1[1] = 0;

	similarityRegistration->SetNumberOfLevels(numberOfLevels1);
	similarityRegistration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel1);
	similarityRegistration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel1);

	/*similarityMetric->SetNumberOfHistogramBins(numberOfHistogramBins);
	MetricType::HistogramType::SizeType histogramSize;
	histogramSize.SetSize(2);
	histogramSize[0] = numberOfHistogramBins;
	histogramSize[1] = numberOfHistogramBins;
	similarityMetric->SetHistogramSize(histogramSize);*/

	similarityOptimizer->SetNumberOfIterations(NumberOfIterations1);
	similarityOptimizer->SetRelaxationFactor(0.5);
	similarityOptimizer->SetLearningRate(0.8);
	//similarityOptimizer->SetDoEstimateLearningRateOnce(true);
	similarityOptimizer->SetDoEstimateScales(true);
	similarityOptimizer->SetMinimumStepLength(0.00001);
	//similarityOptimizer->SetReturnBestParametersAndValue(true);
	//Set optimizer scales(movestep = 1/scaling factor * OptimizerStepLenght)
	typedef SOptimizerType::ScalesType       SOptimizerScalesType;
	SOptimizerScalesType optimizerScales(similarityTx->GetNumberOfParameters());
	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizerScales[5] = translationScale;
	optimizerScales[6] = rescalefactor;
	similarityOptimizer->SetScales(optimizerScales);

	//std::cout << "OptimizerScale: " << similarityOptimizer->GetScales()<< std::endl;
	//std::cout << "LearningRate: " << similarityOptimizer->GetCurrentStepLength() << std::endl;


	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer1 = CommandIterationUpdate::New();
	similarityOptimizer->AddObserver(itk::IterationEvent(), observer1);

	// Create the Command interface observer and register it with the optimizer.
	//
	typedef RegistrationInterfaceCommand<SRegistrationType> SimilarityCommandType;
	SimilarityCommandType::Pointer command1 = SimilarityCommandType::New();
	similarityRegistration->AddObserver(itk::MultiResolutionIterationEvent(), command1);

	/*------------------
	Affine stage
	------------------*/
	typedef itk::AffineTransform< double, Dimension > ATransformType;

	//typedef itk::ConjugateGradientLineSearchOptimizerv4Template<double>  AOptimizerType;
	typedef itk::RegularStepGradientDescentOptimizerv4<double>  AOptimizerType;

	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType > ARegistrationType;

	AOptimizerType::Pointer      affineOptimizer = AOptimizerType::New();
	MetricType::Pointer          affineMetric = MetricType::New();
	ARegistrationType::Pointer   affineRegistration = ARegistrationType::New();

	affineRegistration->SetOptimizer(affineOptimizer);
	affineRegistration->SetMetric(affineMetric);

	//affineMetric->SetNumberOfHistogramBins(50);

	fixedImageReader->Update();
	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

	ATransformType::Pointer   affineTx = ATransformType::New();

	affineRegistration->SetInitialTransform(affineTx);
	affineRegistration->InPlaceOn();

	affineRegistration->SetFixedImage(fixedImageReader->GetOutput());
	affineRegistration->SetMovingImage(movingImageReader->GetOutput());
	affineRegistration->SetObjectName("AffineRegistration");

	affineRegistration->SetMovingInitialTransformInput(
		similarityRegistration->GetTransformOutput());


	ScalesEstimatorType::Pointer scalesEstimator2 =
		ScalesEstimatorType::New();
	scalesEstimator2->SetMetric(affineMetric);
	scalesEstimator2->SetTransformForward(true);

	affineOptimizer->SetScalesEstimator(scalesEstimator2);
	affineOptimizer->SetDoEstimateLearningRateOnce(true);
	affineOptimizer->SetDoEstimateLearningRateAtEachIteration(false);
	//affineOptimizer->SetLowerLimit(0);
	//affineOptimizer->SetUpperLimit(2);
	//affineOptimizer->SetEpsilon(0.2);
	affineOptimizer->SetNumberOfIterations(NumberOfIterations2);
	affineOptimizer->SetMinimumConvergenceValue(1e-6);
	//affineOptimizer->SetConvergenceWindowSize(10);

	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdate::Pointer observer2 = CommandIterationUpdate::New();
	affineOptimizer->AddObserver(itk::IterationEvent(), observer2);


	const unsigned int numberOfLevels2 = 1;

	ARegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel2;
	shrinkFactorsPerLevel2.SetSize(numberOfLevels2);
	shrinkFactorsPerLevel2[0] = 2;
	//shrinkFactorsPerLevel2[1] = 1;

	ARegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel2;
	smoothingSigmasPerLevel2.SetSize(numberOfLevels2);
	smoothingSigmasPerLevel2[0] = 0;
   // smoothingSigmasPerLevel2[1] = 0;

	affineRegistration->SetNumberOfLevels(numberOfLevels2);
	affineRegistration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel2);
	affineRegistration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel2);

	// Create the Command interface observer and register it with the optimizer.
	//
	typedef RegistrationInterfaceCommand<ARegistrationType> AffineCommandType;
	AffineCommandType::Pointer command2 = AffineCommandType::New();
	affineRegistration->AddObserver(itk::MultiResolutionIterationEvent(), command2);

	//  Software Guide : BeginLatex
	//
	//  Once all the registration components are in place,
	//  finally we trigger the whole registration process, including two cascaded
	//  registration stages, by calling \code{Update()} on the registration
	//  filter of the last stage, which causes both stages be updated.
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	try
	{
		affineRegistration->Update();
		std::cout << "Optimizer stop condition: "
			<< affineRegistration->
			GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}



	/*---------------------------------
	combine the transforms
	---------------------------------*/
	std::cout << "Start combine transforms" << std::endl;
	typedef itk::CompositeTransform< double,
		Dimension >  CompositeTransformType;
	CompositeTransformType::Pointer   compositeTransform =
		CompositeTransformType::New();
	CompositeTransformType::Pointer   compositeTransform2 =
		CompositeTransformType::New();

	compositeTransform->AddTransform(similarityTx);
	compositeTransform->AddTransform(affineTx);
	//compositeTransform->AddTransform(s2sTx);

	//compositeTransform->AddTransform(s2sTx);
	//compositeTransform->AddTransform(similarityTx);
	//compositeTransform->AddTransform(affineTx);

	// Software Guide : EndCodeSnippet

	std::cout << " Similaritytransform parameters after registration: " << std::endl
		<< similarityOptimizer->GetCurrentPosition() << std::endl
		<< " Last LearningRate: " << similarityOptimizer->GetCurrentStepLength() << std::endl;
	std::cout << " Affine transform parameters after registration: " << std::endl
		<< affineOptimizer->GetCurrentPosition() << std::endl
		<< " Last LearningRate: " << affineOptimizer->GetLearningRate() << std::endl;

	std::cout << " subject2subject transform parameters after registration: " << std::endl
		<< s2sOptimizer->GetCurrentPosition() << std::endl
		<< " Last LearningRate: " << s2sOptimizer->GetCurrentStepLength() << std::endl;


	// Resmaple 

	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;
	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(compositeTransform);
	resample->SetInput(movingImageReader->GetOutput());

	resample->SetSize(scan2Image->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(scan2Image->GetOrigin());
	resample->SetOutputSpacing(scan2Image->GetSpacing());
	resample->SetOutputDirection(scan2Image->GetDirection());
	resample->SetDefaultPixelValue(backgroundGrayLevel);

	typedef  unsigned char                           OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	typedef itk::CastImageFilter<
		FixedImageType,
		OutputImageType >          CastFilterType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();

	writer->SetFileName(outImagefile);

	caster->SetInput(resample->GetOutput());
	writer->SetInput(caster->GetOutput());
	writer->Update();

	typedef itk::SubtractImageFilter<
		FixedImageType,
		FixedImageType,
		FixedImageType > DifferenceFilterType;
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	difference->SetInput1(scan2Image);
	difference->SetInput2(resample->GetOutput());

	WriterType::Pointer writer2 = WriterType::New();

	typedef itk::RescaleIntensityImageFilter<
		FixedImageType,
		OutputImageType >   RescalerType;

	RescalerType::Pointer intensityRescaler = RescalerType::New();

	intensityRescaler->SetInput(difference->GetOutput());
	intensityRescaler->SetOutputMinimum(0);
	intensityRescaler->SetOutputMaximum(255);

	writer2->SetInput(intensityRescaler->GetOutput());
	resample->SetDefaultPixelValue(1);

	// Compute the difference image between the
	// fixed and resampled moving image.
	if (differenceAfterRegistration != "")
	{
		writer2->SetFileName(differenceAfterRegistration);
		writer2->Update();
	}


	typedef itk::IdentityTransform< double, Dimension > IdentityTransformType;
	IdentityTransformType::Pointer identity = IdentityTransformType::New();

	// Compute the difference image between the
	// fixed and moving image before registration.
	if (differenceBeforeRegistration != "")
	{
		resample->SetTransform(identity);
		writer2->SetFileName(differenceBeforeRegistration);
		writer2->Update();
	}

	////resample the lable map
	//MovingImageReaderType::Pointer labelMapReader = MovingImageReaderType::New();

	//ResampleFilterType::Pointer labelMapResampler = ResampleFilterType::New();

	//labelMapResampler->SetTransform(compositeTransform);
	//labelMapResampler->SetInput(labelMapReader->GetOutput());
	//labelMapResampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	//labelMapResampler->SetOutputOrigin(fixedImage->GetOrigin());
	//labelMapResampler->SetOutputSpacing(fixedImage->GetSpacing());
	//labelMapResampler->SetOutputDirection(fixedImage->GetDirection());
	//labelMapResampler->SetDefaultPixelValue(0);

	//typedef  float	labelMapPixelType;
	//typedef itk::Image< labelMapPixelType, Dimension >                labelMapImageType;
	//typedef itk::CastImageFilter< FixedImageType, labelMapImageType > labelMapCastFilterType;
	//typedef itk::ImageFileWriter< labelMapImageType >                 labelMapWriterType;

	//labelMapWriterType::Pointer      labelMapWriter = labelMapWriterType::New();
	//labelMapCastFilterType::Pointer  labelMapCaster = labelMapCastFilterType::New();

	//labelMapCaster->SetInput(labelMapResampler->GetOutput());
	//labelMapWriter->SetInput(labelMapCaster->GetOutput());
	//if (labelMapInputFile1 != "" && labelMapOutputFile1 != "")
	//{
	//	labelMapReader->SetFileName(labelMapInputFile1);
	//	labelMapWriter->SetFileName(labelMapOutputFile1);
	//	labelMapWriter->Update();
	//}
	//if (labelMapInputFile2 != "" && labelMapOutputFile2 != "")
	//{
	//	labelMapReader->SetFileName(labelMapInputFile2);
	//	labelMapWriter->SetFileName(labelMapOutputFile2);
	//	labelMapWriter->Update();
	//}
	//if (labelMapInputFile3 != "" && labelMapOutputFile3 != "")
	//{
	//	labelMapReader->SetFileName(labelMapInputFile3);
	//	labelMapWriter->SetFileName(labelMapOutputFile3);
	//	labelMapWriter->Update();
	//}
	//if (labelMapInputFile4 != "" && labelMapOutputFile4 != "")
	//{
	//	labelMapReader->SetFileName(labelMapInputFile4);
	//	labelMapWriter->SetFileName(labelMapOutputFile4);
	//	labelMapWriter->Update();
	//}
	//
	
	return EXIT_SUCCESS;
}
