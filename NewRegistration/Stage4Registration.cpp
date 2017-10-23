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
		std::cerr << " fixedImageFile  movingImageFile fixedImageMaskFile";
		std::cerr << " outputImagefile  [differenceOutputfile] ";
		std::cerr << " [differenceBeforeRegistration] " << std::endl;
		return EXIT_FAILURE;
	}

	const std::string fixedImageFile = argv[1];
	const std::string movingImageFile = argv[2];
	const std::string maskImageFile = argv[3];
	const std::string outImagefile = argv[4];
	const std::string differenceAfterRegistration = (argc > 5) ? argv[5] : "";
	const int NumberOfIterations1 = (argc > 6) ? atoi(argv[6]) : 100;
	const int nn = 6;//number of argument except labelmap names
					  //const int nn = 11;//number of argument except labelmap names
	const std::string labelMapInputFile1 = (argc > nn + 1) ? argv[nn + 1] : "";
	const std::string labelMapOutputFile1 = (argc > nn + 2) ? argv[nn + 2] : "";
	const std::string labelMapInputFile2 = (argc > nn + 3) ? argv[nn + 3] : "";
	const std::string labelMapOutputFile2 = (argc > nn + 4) ? argv[nn + 4] : "";
	const std::string labelMapInputFile3 = (argc > nn + 5) ? argv[nn + 5] : "";
	const std::string labelMapOutputFile3 = (argc > nn + 6) ? argv[nn + 6] : "";
	const std::string labelMapInputFile4 = (argc > nn + 7) ? argv[nn + 7] : "";
	const std::string labelMapOutputFile4 = (argc > nn + 8) ? argv[nn + 8] : "";
	const std::string labelMapInputFile5 = (argc > nn + 9) ? argv[nn + 9] : "";
	const std::string labelMapOutputFile5 = (argc > nn + 10) ? argv[nn + 10] : "";

	const    unsigned int    Dimension = 3;
	typedef  double           PixelType;

	const PixelType backgroundGrayLevel = 80;
	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;

	typedef itk::AffineTransform< double, Dimension > ATransformType;
	typedef itk::RegularStepGradientDescentOptimizerv4<double> AOptimizerType;
	typedef itk::CorrelationImageToImageMetricv4<FixedImageType,MovingImageType > MetricType;
	typedef itk::ImageRegistrationMethodv4<FixedImageType,MovingImageType,ATransformType > RegistrationType;
	
	MetricType::Pointer metric = MetricType::New();

	AOptimizerType::Pointer optimizer = AOptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	ATransformType::Pointer  transform = ATransformType::New();
	registration->SetInitialTransform(transform);
	registration->InPlaceOn();

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	fixedImageReader->SetFileName(fixedImageFile);
	movingImageReader->SetFileName(movingImageFile);
	registration->SetFixedImage(fixedImageReader->GetOutput());
	registration->SetMovingImage(movingImageReader->GetOutput());
	fixedImageReader->Update();

	typedef itk::CenteredTransformInitializer<ATransformType,FixedImageType,MovingImageType >  TransformInitializerType;
	TransformInitializerType::Pointer initializer = TransformInitializerType::New();
	initializer->SetTransform(transform);
	initializer->SetFixedImage(fixedImageReader->GetOutput());
	initializer->SetMovingImage(movingImageReader->GetOutput());
	initializer->MomentsOn();
	initializer->InitializeTransform();
	//transform->SetAngle(0.0);

	typedef AOptimizerType::ScalesType OptimizerScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
	const double translationScale = 1.0 / 1000.0;
	optimizerScales[0] = 1.0;
	optimizerScales[1] = translationScale;
	optimizerScales[2] = translationScale;
	optimizerScales[3] = translationScale;
	optimizerScales[4] = translationScale;
	optimizer->SetScales(optimizerScales);
	optimizer->SetLearningRate(0.001);
	optimizer->SetMinimumStepLength(0.0001);
	optimizer->SetNumberOfIterations(NumberOfIterations1); // argument 4
	
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
	maskReader->SetFileName(maskImageFile);
	
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
	
	// One level registration process without shrinking and smoothing.
	
	const unsigned int numberOfLevels = 1;
	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(1);
	shrinkFactorsPerLevel[0] = 1;
	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(1);
	smoothingSigmasPerLevel[0] = 0;
	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
	try
	{
		registration->Update();
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
	AOptimizerType::ParametersType finalParameters = transform->GetParameters();
	const double finalAngle = finalParameters[0];
	const double finalRotationCenterX = finalParameters[1];
	const double finalRotationCenterY = finalParameters[2];
	const double finalTranslationX = finalParameters[3];
	const double finalTranslationY = finalParameters[4];
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();
	// Print out results
	//
	const double finalAngleInDegrees = finalAngle * 45.0 / std::atan(1.0);
	std::cout << "Result = " << std::endl;
	std::cout << " Angle (radians) " << finalAngle << std::endl;
	std::cout << " Angle (degrees) " << finalAngleInDegrees << std::endl;
	std::cout << " Center X      = " << finalRotationCenterX << std::endl;
	std::cout << " Center Y      = " << finalRotationCenterY << std::endl;
	std::cout << " Translation X = " << finalTranslationX << std::endl;
	std::cout << " Translation Y = " << finalTranslationY << std::endl;
	std::cout << " Iterations    = " << numberOfIterations << std::endl;
	std::cout << " Metric value  = " << bestValue << std::endl;

	//  Software Guide : BeginLatex
	//
	//  Let's execute this example over some of the images provided in
	//  \code{Examples/Data}, for example:
	//
	//  \begin{itemize}
	//  \item \code{BrainProtonDensitySliceBorder20.png}
	//  \item \code{BrainProtonDensitySliceR10X13Y17.png}
	//  \end{itemize}
	//
	//  The second image is the result of intentionally rotating the first
	//  image by $10$ degrees and shifting it $13mm$ in $X$ and $17mm$ in
	//  $Y$. Both images have unit-spacing and are shown in Figure
	//  \ref{fig:FixedMovingImageRegistration5}.
	//
	//  The registration converges after $23$ iterations and produces the
	//  following results:
	//
	//  \begin{verbatim}
	//
	//  Angle (radians) 0.174407
	//  Angle (degrees) 9.99281
	//  Center X      = 111.172
	//  Center Y      = 131.563
	//  Translation X = 12.4584
	//  Translation Y = 16.0726
	//
	//  \end{verbatim}
	//
	//  These values are a very close match to the true misalignments
	//  introduced in the moving image.
	//
	//  Now we resample the moving image using the transform resulting from the
	//  registration process.
	//
	//  Software Guide : EndLatex
	// Software Guide : BeginCodeSnippet
	ATransformType::MatrixType matrix = transform->GetMatrix();
	ATransformType::OffsetType offset = transform->GetOffset();
	std::cout << "Matrix = " << std::endl << matrix << std::endl;
	std::cout << "Offset = " << std::endl << offset << std::endl;
	// Software Guide : EndCodeSnippet
	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;
	ATransformType::Pointer finalTransform = ATransformType::New();
	finalTransform->SetParameters(finalParameters);
	finalTransform->SetFixedParameters(transform->GetFixedParameters());
	ResampleFilterType::Pointer resample = ResampleFilterType::New();
	resample->SetTransform(transform); //finalTransform
	resample->SetInput(movingImageReader->GetOutput());
	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
	resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedImage->GetOrigin());
	resample->SetOutputSpacing(fixedImage->GetSpacing());
	resample->SetOutputDirection(fixedImage->GetDirection());
	resample->SetDefaultPixelValue(100);
	typedef  unsigned char  OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	typedef itk::CastImageFilter<
		FixedImageType,
		OutputImageType > CastFilterType;
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
		OutputImageType > DifferenceFilterType;
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();
	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput(difference->GetOutput());
	// Compute the difference image between the
	// fixed and resampled moving image.
	if (argc >= 6)
	{
		difference->SetInput1(fixedImageReader->GetOutput());
		difference->SetInput2(resample->GetOutput());
		writer2->SetFileName(differenceAfterRegistration);
		writer2->Update();
	}

	ofstream logFile;
	logFile.open("E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/log.txt", ios::app);
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	logFile <<
		"========================================================\n" <<
		asctime(timeinfo) <<
		"Stage 4 Registration:\n" <<
		"FixedImage: " << fixedImageFile << "\n" <<
		"MovingImage: " << movingImageFile << "\n" <<
		"OutputImage: " << differenceAfterRegistration << "\n" << "\n" <<

		"Affine transform parameters after registration: " << "\n" <<
		finalTransform->GetParameters() << "\n" <<
		"translations: " << finalTransform->GetTranslation() << "\n" <<
		"center: " << finalTransform->GetCenter() << "\n" <<
		"matrix: " << finalTransform->GetMatrix() <<
		"offset: " << finalTransform->GetOffset() << "\n" <<
		"Last Step Length: " << optimizer->GetCurrentStepLength() << "\n" <<
		"Last Metric Vuale: " << metric->GetCurrentValue() << "\n" <<
		"Affine Learning Rate: " << optimizer->GetLearningRate() << "\n" << endl;

	logFile.close();

	std::cout << " Affine transform parameters after registration: " << std::endl
		<< finalTransform->GetParameters() << "\n" <<
		"translations: " << finalTransform->GetTranslation() << "\n" <<
		"center: " << finalTransform->GetCenter() << "\n" <<
		"matrix: " << finalTransform->GetMatrix() << "\n" <<
		"offset: " << finalTransform->GetOffset() << "\n" <<
		" Last Step Length: " << optimizer->GetCurrentStepLength() << "\n" <<
		"Affine Learning Rate: " << optimizer->GetLearningRate() << "\n" << std::endl;

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

	labelMapCaster->SetInput(labelMapResampler->GetOutput());
	labelMapWriter->SetInput(labelMapCaster->GetOutput());

	// Labelmap writer

	if (labelMapInputFile1 != "" && labelMapOutputFile1 != "")
	{
		labelMapReader->SetFileName(labelMapInputFile1);
		labelMapWriter->SetFileName(labelMapOutputFile1);
		labelMapWriter->Update();
	}
	if (labelMapInputFile2 != "" && labelMapOutputFile2 != "")
	{
		labelMapReader->SetFileName(labelMapInputFile2);
		labelMapWriter->SetFileName(labelMapOutputFile2);
		labelMapWriter->Update();
	}
	if (labelMapInputFile3 != "" && labelMapOutputFile3 != "")
	{
		labelMapReader->SetFileName(labelMapInputFile3);
		labelMapWriter->SetFileName(labelMapOutputFile3);
		labelMapWriter->Update();
	}
	if (labelMapInputFile4 != "" && labelMapOutputFile4 != "")
	{
		labelMapReader->SetFileName(labelMapInputFile4);
		labelMapWriter->SetFileName(labelMapOutputFile4);
		labelMapWriter->Update();
	}
	if (labelMapInputFile5 != "" && labelMapOutputFile5 != "")
	{
		labelMapReader->SetFileName(labelMapInputFile5);
		labelMapWriter->SetFileName(labelMapOutputFile5);
		labelMapWriter->Update();
	} //uncomment till here
	//// Compute the difference image between the
	//// fixed and moving image before registration.
	//if (argc >= 7)
	//{
	//	writer2->SetFileName(argv[7]);
	//	difference->SetInput1(fixedImageReader->GetOutput());
	//	difference->SetInput2(movingImageReader->GetOutput());
	//	writer2->Update();
	//}
	return EXIT_SUCCESS;
}
