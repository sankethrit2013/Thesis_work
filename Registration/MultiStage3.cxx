
#include "itkImageRegistrationMethodv4.h"
#include "itkCorrelationImageToImageMetricv4.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSOptimizerv4.h"
#include "itkLBFGSBOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"

#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkSubtractImageFilter.h"
using namespace std;

// NOTE: the LBFGSOptimizerv4 does not invoke events
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
	typedef itk::LBFGSBOptimizerv4     OptimizerType;
	typedef   const OptimizerType *    OptimizerPointer;
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
		std::cout << optimizer->GetCurrentMetricValue() << "   ";
		std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
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
		getchar();
		return EXIT_FAILURE;
	}

	const    unsigned int    ImageDimension = 3;
	typedef  double           PixelType;

	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
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
	try{ 
		fixedImageReader->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		getchar();
		return EXIT_FAILURE;
	}

	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

	unsigned int numberOfGridNodesInOneDimension = 8;
	if (argc >= 7)
	{
		numberOfGridNodesInOneDimension = std::atof(argv[6]);
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

	/*
	typedef BSplineTransformType::ParametersType BSplineParametersType;
	BSplineParametersType BSplineParameters(BSplineTransform->GetNumberOfParameters());
	BSplineParameters.Fill(0.0);
	BSplineTransform->SetParameters(BSplineParameters);*/


	registration->SetInitialTransform(BSplineTransform);
	registration->InPlaceOn();


	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImageReader->GetOutput());




	unsigned int numberofIterations = 10;
	if (argc >= 8)
	{
		numberofIterations = std::atof(argv[7]);
	}





	//setting for LBFGS optimizer
	typedef itk::LBFGSOptimizerv4                         OptimizerType;
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	registration->SetOptimizer(optimizer);
	typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
	ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
	scalesEstimator->SetMetric(metric);
	scalesEstimator->SetTransformForward(true);
	scalesEstimator->SetSmallParameterVariation(1.0);
	optimizer->SetGradientConvergenceTolerance( 1e-10 );
	optimizer->SetLineSearchAccuracy(0.9);
	optimizer->SetDefaultStepLength( 1 );
	optimizer->TraceOn();
	optimizer->SetNumberOfIterations(numberofIterations);
	//optimizer->SetMaximumNumberOfFunctionEvaluations(MaximumNumberOfFunctionEvaluations);
	optimizer->SetScalesEstimator( scalesEstimator );


	////setting for LBFGSB optimizer
	//typedef itk::LBFGSBOptimizerv4       OptimizerType;
	//OptimizerType::Pointer      optimizer = OptimizerType::New();
	//registration->SetOptimizer(optimizer);
	//const unsigned int numParameters =
	// BSplineTransform->GetNumberOfParameters();
	//OptimizerType::BoundSelectionType boundSelect(numParameters);
	//OptimizerType::BoundValueType upperBound(numParameters);
	//OptimizerType::BoundValueType lowerBound(numParameters);
	//boundSelect.Fill(OptimizerType::UNBOUNDED);
	//upperBound.Fill(0.0);
	//lowerBound.Fill(0.0);
	//optimizer->SetBoundSelection(boundSelect);
	//optimizer->SetUpperBound(upperBound);
	//optimizer->SetLowerBound(lowerBound);
	//optimizer->SetCostFunctionConvergenceFactor(1e+12);
	//optimizer->SetGradientConvergenceTolerance(1.0e-8);
	//optimizer->SetNumberOfIterations(numberofIterations);
	//optimizer->SetMaximumNumberOfFunctionEvaluations(500);
	//optimizer->SetMaximumNumberOfCorrections(100);

	//typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
	//OptimizerType::Pointer      optimizer = OptimizerType::New();
	//registration->SetOptimizer(optimizer);
	//optimizer->SetMinimumStepLength(1e-6);
	//optimizer->SetRelaxationFactor(0.7);
	//optimizer->SetNumberOfIterations(10);


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


	//Add time and memory probes
	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;

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
		getchar();
		return EXIT_FAILURE;
	}

	// Report the time and memory taken by the registration
	chronometer.Report(std::cout);
	memorymeter.Report(std::cout);

	/*
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);
	std::cout << "Starting Registration "
	<< std::endl;
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
	*/
	std::cout << "Last BSplineTransform Parameters" << std::endl;
	std::cout << BSplineTransform->GetParameters() << std::endl;

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
	resample->SetDefaultPixelValue(80);

	typedef  unsigned char  OutputPixelType;

	typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;

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
		getchar();
		return EXIT_FAILURE;
	}

	typedef itk::SubtractImageFilter<
		FixedImageType,
		FixedImageType,
		OutputImageType > DifferenceFilterType;

	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput(difference->GetOutput());


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
			getchar();
			return EXIT_FAILURE;
		}
	}


	// Compute the difference image between the
	// fixed and moving image before registration.
	//if( argc >= 6 )
	//  {
	//  writer2->SetFileName( argv[5] );
	//  difference->SetInput1( fixedImageReader->GetOutput() );
	//  difference->SetInput2( movingImageReader->GetOutput() );
	//  try
	//    {
	//    writer2->Update();
	//    }
	//  catch( itk::ExceptionObject & err )
	//    {
	//    std::cerr << "ExceptionObject caught !" << std::endl;
	//    std::cerr << err << std::endl;
	//    return EXIT_FAILURE;
	//    }
	//  }

	// //resample the lable map
	MovingImageReaderType::Pointer labelMapReader = MovingImageReaderType::New();

	ResampleFilterType::Pointer labelMapResampler = ResampleFilterType::New();

	labelMapResampler->SetTransform(BSplineTransform);
	labelMapResampler->SetInput(labelMapReader->GetOutput());
	labelMapResampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	labelMapResampler->SetOutputOrigin(fixedImage->GetOrigin());
	labelMapResampler->SetOutputSpacing(fixedImage->GetSpacing());
	labelMapResampler->SetOutputDirection(fixedImage->GetDirection());
	labelMapResampler->SetDefaultPixelValue(0);

	typedef  float	labelMapPixelType;
	typedef itk::Image< labelMapPixelType, ImageDimension >                labelMapImageType;
	typedef itk::CastImageFilter< FixedImageType, labelMapImageType > labelMapCastFilterType;
	typedef itk::ImageFileWriter< labelMapImageType >                 labelMapWriterType;

	labelMapWriterType::Pointer      labelMapWriter = labelMapWriterType::New();
	labelMapCastFilterType::Pointer  labelMapCaster = labelMapCastFilterType::New();

	//labelMapReader->SetFileName(argv[9]);
	//labelMapWriter->SetFileName(argv[10]);
	//labelMapCaster->SetInput(labelMapResampler->GetOutput());
	//labelMapWriter->SetInput(labelMapCaster->GetOutput());
	//labelMapWriter->Update();

	//labelMapReader->SetFileName(argv[11]);
	//labelMapWriter->SetFileName(argv[12]);
	//labelMapWriter->Update();

	//labelMapReader->SetFileName(argv[13]);
	//labelMapWriter->SetFileName(argv[14]);
	//labelMapWriter->Update();

	//labelMapReader->SetFileName(argv[15]);
	//labelMapWriter->SetFileName(argv[16]);
	//labelMapWriter->Update();

	//// Generate the explicit deformation field resulting from
	//// the registration.
	//typedef itk::Vector< float, ImageDimension >          VectorPixelType;
	//typedef itk::Image< VectorPixelType, ImageDimension > DisplacementFieldImageType;

	//typedef itk::TransformToDisplacementFieldFilter<
	//                      DisplacementFieldImageType,
	//                      CoordinateRepType >             DisplacementFieldGeneratorType;

	///** Create an setup displacement field generator. */
	//DisplacementFieldGeneratorType::Pointer dispfieldGenerator =
	//                                                DisplacementFieldGeneratorType::New();
	//dispfieldGenerator->UseReferenceImageOn();
	//dispfieldGenerator->SetReferenceImage( fixedImage );
	//dispfieldGenerator->SetTransform( BSplineTransform );

	/*try
	{
	dispfieldGenerator->Update();
	}
	catch ( itk::ExceptionObject & err )
	{
	std::cerr << "Exception detected while generating deformation field";
	std::cerr << " : "  << err << std::endl;
	return EXIT_FAILURE;
	}*/

	//typedef itk::ImageFileWriter< DisplacementFieldImageType >  FieldWriterType;
	//FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

	//fieldWriter->SetInput( dispfieldGenerator->GetOutput() );

	/*if( argc >= 9 )
	{
	fieldWriter->SetFileName( argv[8] );
	try
	{
	fieldWriter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
	std::cerr << "Exception thrown " << std::endl;
	std::cerr << excp << std::endl;
	return EXIT_FAILURE;
	}
	}*/

	getchar();
	return EXIT_SUCCESS;
}
