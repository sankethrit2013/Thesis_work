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
#include "itkSimilarity2DTransform.h"
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
#include <itkTIFFImageIO.h>
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

//  The following section of code implements a Command observer
//  that will monitor the configurations of the registration process
//  at every change of stage and resolution level.

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
		std::cout << optimizer->GetCurrentPosition() << "  ";
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
	/*const double translationScale = (argc > 7) ? atof(argv[7]) : 0.2;
	const double rescalefactor = (argc > 8) ? atof(argv[8]) : 2;*/
	/*const double similarityLearnRate = (argc > 9) ? atof(argv[9]) : 0.01;
	const double affineLearnRate = (argc > 10) ? atof(argv[10]) : 0.01;*/
	const int NumberOfIterations2 = (argc > 9) ? atoi(argv[9]) : 0;


	// Change dimendion to 2 for Masoumeh's code
	const    unsigned int    Dimension = 2;
	typedef  double           PixelType;

	const PixelType backgroundGrayLevel = 80;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;

	typedef itk::Similarity2DTransform<double>     STransformType;
	typedef itk::RegularStepGradientDescentOptimizerv4<double> SOptimizerType;
	//typedef itk::MattesMutualInformationImageToImageMetricv4<FixedImageType, MovingImageType >  MetricType;
	typedef itk::CorrelationImageToImageMetricv4< FixedImageType, MovingImageType >   MetricType;
	//typedef itk::MutualInformationImageToImageMetric<FixedImageType, MovingImageType> MetricType;
	//typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<FixedImageType,MovingImageType> MetricType;
	typedef itk::ImageRegistrationMethodv4<
		FixedImageType,
		MovingImageType >  SRegistrationType;

	SOptimizerType::Pointer      similarityOptimizer = SOptimizerType::New();
	MetricType::Pointer         similarityMetric = MetricType::New();
	SRegistrationType::Pointer   similarityRegistration = SRegistrationType::New();

	similarityRegistration->SetOptimizer(similarityOptimizer);
	similarityRegistration->SetMetric(similarityMetric);

	STransformType::Pointer similarityTx = STransformType::New();

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	try{
		//fixedImageReader->SetImageIO(itk::TIFFImageIO::New());
		fixedImageReader->SetFileName(fixedImageFile);
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "Fixed Image ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	try{
		//movingImageReader->SetImageIO(itk::TIFFImageIO::New());
		movingImageReader->SetFileName(movingImageFile);
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "Moving Image ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

	typedef itk::CenteredTransformInitializer<STransformType, FixedImageType, MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer
		= TransformInitializerType::New();

	initializer->SetTransform(similarityTx);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	initializer->MomentsOn();
	try{
		initializer->InitializeTransform();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "Moving Image ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	cout << "Similarity starting" << endl;



	//similarityTx->SetScale(0.7);

	//similarityRegistration->SetMovingInitialTransform(similarityTx);
	similarityRegistration->SetInitialTransform(similarityTx);
	std::cout << "\nInitial transform" << similarityRegistration->GetInitialTransform() << std::endl;

	//similarityRegistration->InPlaceOn();
	similarityRegistration->SetFixedImage(fixedImageReader->GetOutput());
	similarityRegistration->SetMovingImage(movingImageReader->GetOutput());
	similarityRegistration->SetObjectName("SimilarityRegistration");

	// Set ScaleEstimator
	typedef itk::RegistrationParameterScalesFromPhysicalShift<
		MetricType> ScalesEstimatorType;
	ScalesEstimatorType::Pointer scalesEstimator1 = ScalesEstimatorType::New();

	scalesEstimator1->SetMetric(similarityMetric);
	scalesEstimator1->SetTransformForward(true);
	////scalesEstimator1->SetSmallParameterVariation(0.1);
	//similarityOptimizer->SetScalesEstimator(scalesEstimator1);
	//similarityOptimizer->SetDoEstimateScales(true);
	//similarityOptimizer->SetDoEstimateLearningRateOnce(true);
	//similarityOptimizer->SetDoEstimateLearningRateAtEachIteration(true);

	const unsigned int numberOfLevels1 = 1;

	SRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel1;
	shrinkFactorsPerLevel1.SetSize(numberOfLevels1);
	shrinkFactorsPerLevel1[0] = 3;
	//shrinkFactorsPerLevel1[1] = 2;

	SRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel1;
	smoothingSigmasPerLevel1.SetSize(numberOfLevels1);
	smoothingSigmasPerLevel1[0] = 0;
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
	similarityOptimizer->SetRelaxationFactor(0.7); //0.7
	similarityOptimizer->SetLearningRate(0.1);
	//similarityOptimizer->SetLearningRate(similarityLearnRate);
	similarityOptimizer->SetMinimumStepLength(1e-6);
	similarityOptimizer->SetGradientMagnitudeTolerance(1e-6);
	similarityOptimizer->SetReturnBestParametersAndValue(true);



	//Set optimizer scales(movestep = 1/scaling factor * OptimizerStepLength)
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

	try
	{
		similarityRegistration->Update();
		std::cout << "Optimizer stop condition: "
			<< similarityRegistration->
			GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}

#pragma region AffineStage
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
	//// Sanketh line of code
	//affineOptimizer->SetLearningRate(0.1);
	////affineOptimizer->SetLearningRate(affineLearnRate);


	affineOptimizer->SetMinimumConvergenceValue(1e-6);
	affineOptimizer->SetGradientMagnitudeTolerance(1e-6);
	//affineOptimizer->SetConvergenceWindowSize(10);
	affineOptimizer->SetReturnBestParametersAndValue(true);
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

#pragma endregion

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


	typedef itk::CompositeTransform< double,
		Dimension >  CompositeTransformType;
	CompositeTransformType::Pointer   compositeTransform =
		CompositeTransformType::New();
	//compositeTransform->AddTransform(TranslationTransform);
	compositeTransform->AddTransform(similarityTx);
	compositeTransform->AddTransform(affineTx);


	// Software Guide : EndCodeSnippet


	ofstream logFile;
	logFile.open("E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/log.txt", ios::app);
	//logFile.open("E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/log.txt", ios::app);

	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	logFile <<
		"========================================================\n" <<
		asctime(timeinfo) <<
		"MultiStage Registration:\n" <<
		"FixedImage: " << fixedImageFile << "\n" <<
		"MovingImage: " << movingImageFile << "\n" <<
		"OutputImage: " << differenceAfterRegistration << "\n" << "\n" <<

		/*"Stage 1: \n" <<
		"Scales: " << similarityOptimizer->GetScales() << "\n" <<
		"Learning Rate: " << similarityOptimizer->GetLearningRate() << "\n" <<
		"Iterations" <<*/
		"Similaritytransform parameters after registration: " << "\n" <<
		similarityTx->GetParameters() << "\n" <<
		"translations: " << similarityTx->GetTranslation() << "\n" <<
		"center: " << similarityTx->GetCenter() << "\n" <<
		"matrix: " << similarityTx->GetMatrix() <<
		"scale: " << similarityTx->GetScale() << "\n" <<
		"offset: " << similarityTx->GetOffset() << "\n" <<
		//"versor: " << similarityTx->GetVersor() << "\n" <<
		"Last Step Length: " << similarityOptimizer->GetCurrentStepLength() << "\n" <<
		"Last Metric Vuale: " << similarityMetric->GetCurrentValue() << "\n" << "\n" <<
		"Similarity Learning Rate: " << similarityOptimizer->GetLearningRate() << "\n" <<

		"Affine transform parameters after registration: " << "\n" <<
		affineTx->GetParameters() << "\n" <<
		"translations: " << affineTx->GetTranslation() << "\n" <<
		"center: " << affineTx->GetCenter() << "\n" <<
		"matrix: " << affineTx->GetMatrix() <<
		"offset: " << affineTx->GetOffset() << "\n" <<
		"Last Step Length: " << affineOptimizer->GetCurrentStepLength() << "\n" <<
		"Last Metric Vuale: " << affineMetric->GetCurrentValue() << "\n" <<
		"Affine Learning Rate: " << affineOptimizer->GetLearningRate() << "\n" <<
		"The composite transform: " << compositeTransform->GetParameters() << "\n" <<
		"The composite transform: " << compositeTransform->GetFixedParameters() << std::endl;

	logFile.close();

	std::cout << " Similaritytransform parameters after registration: " << std::endl
		<< similarityTx->GetParameters() << "\n" <<
		"translations: " << similarityTx->GetTranslation() << "\n" <<
		"center: " << similarityTx->GetCenter() << "\n" <<
		"matrix: " << similarityTx->GetMatrix() << "\n" <<
		"scale: " << similarityTx->GetScale() << "\n" <<
		"offset: " << similarityTx->GetOffset() << "\n" <<
		//"versor: " << similarityTx->GetVersor() << "\n" <<
		" Last Step Length: " << similarityOptimizer->GetCurrentStepLength() << "\n" <<
		"Similarity Learning Rate: " << similarityOptimizer->GetLearningRate() << std::endl;

	std::cout << " Affine transform parameters after registration: " << std::endl
		<< affineTx->GetParameters() << "\n" <<
		"translations: " << affineTx->GetTranslation() << "\n" <<
		"center: " << affineTx->GetCenter() << "\n" <<
		"matrix: " << affineTx->GetMatrix() << "\n" <<
		"offset: " << affineTx->GetOffset() << "\n" <<
		" Last Step Length: " << affineOptimizer->GetCurrentStepLength() << "\n" <<
		"Affine Learning Rate: " << affineOptimizer->GetLearningRate() << "\n" << std::endl;

	std::cout << "The composite transform: " << compositeTransform->GetParameters() << std::endl;
	// Resmaple 

	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;
	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(compositeTransform);
	resample->SetInput(movingImageReader->GetOutput());

	resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedImage->GetOrigin());
	resample->SetOutputSpacing(fixedImage->GetSpacing());
	resample->SetOutputDirection(fixedImage->GetDirection());
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

	difference->SetInput1(fixedImageReader->GetOutput());
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
	
	return EXIT_SUCCESS;
}
