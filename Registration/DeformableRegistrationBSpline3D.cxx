
#include "itkImageRegistrationMethodv4.h"
#include "itkCorrelationImageToImageMetricv4.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSOptimizerv4.h"
#include "itkLBFGSBOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"

#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
using namespace std;

// NOTE: the LBFGSOptimizerv4 does not invoke events
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
		std::cout << optimizer->GetLearningRate()<< "  ";
		std::cout << m_CumulativeIterationIndex++ << std::endl;
	}

private:
	unsigned int m_CumulativeIterationIndex;
};

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile outputImagefile ";
    std::cerr << "  [differenceOutputfile] [differenceBeforeRegistration] ";
    std::cerr << " [numberOfGridNodesInOneDimension] [deformationField] ";
    return EXIT_FAILURE;
    }

  const    unsigned int    ImageDimension = 3;
  typedef  double           PixelType;

  typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
  typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

  const unsigned int SpaceDimension = ImageDimension;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;
#pragma region addPreviousStage

  typedef itk::Similarity3DTransform<double> Similarity3DTransformType;
  Similarity3DTransformType::Pointer SimilarityTransform = Similarity3DTransformType::New();
  Similarity3DTransformType::ParametersType SimilarityParameters(SimilarityTransform->GetNumberOfParameters());
  Similarity3DTransformType::CenterType SCenter = SimilarityTransform->GetCenter();

  typedef itk::AffineTransform<double, ImageDimension> AffineTransformType;
  AffineTransformType::Pointer AffineTransform = AffineTransformType::New();
  AffineTransformType::ParametersType AffineParameters(AffineTransform->GetNumberOfParameters());
  AffineTransformType::CenterType ACenter = AffineTransform->GetCenter();

#pragma region parameterScan1
  // change from log file - sanketh
  SimilarityParameters[0] = -0.018570519712191588;
  SimilarityParameters[1] = -0.0038260363690700897;
  SimilarityParameters[2] = 0.07280084421522154;
  SimilarityParameters[3] = -0.2664431362988693;
  SimilarityParameters[4] = -0.36015435795723477;
  SimilarityParameters[5] = 0.22954843271829656;
  SimilarityParameters[6] = 0.965282633860967;
  SCenter[0] = -2.00214;
  SCenter[1] = -1.55205;
  SCenter[2] = 1.67857;

  AffineParameters[0] = 0.9485629222905053;
  AffineParameters[1] = 0.06871110174135162;
  AffineParameters[2] = -0.00640578592279396;
  AffineParameters[3] = 0.0026359840492128353;
  AffineParameters[4] = 1.0341588561548753;
  AffineParameters[5] = 0.010372887712632008;
  AffineParameters[6] = 0.008983963690218226;
  AffineParameters[7] = -0.02126857634803973;
  AffineParameters[8] = 1.001568959145986;
  AffineParameters[9] = -0.012037008765033977;
  AffineParameters[10] = 0.004438024854140598;
  AffineParameters[11] = -0.003009373867658991;
  ACenter[0] = -2.00214;
  ACenter[1] = -1.55205;
  ACenter[2] = 1.67857;

  //SimilarityParameters[0] = -0.1436948652872146;
  //SimilarityParameters[1] = 0.03157415693267113;
  //SimilarityParameters[2] = 0.1820820683634619;
  //SimilarityParameters[3] = -0.11189565299615623;
  //SimilarityParameters[4] = -0.0853683469364656;
  //SimilarityParameters[5] = 0.015195276647691595;
  //SimilarityParameters[6] = 0.9975690042262584;
  //SCenter[0] = -1.85716;
  //SCenter[1] = -1.89245;
  //SCenter[2] = 1.58383;

  //AffineParameters[0] = 1.0041000386596557;
  //AffineParameters[1] = -0.0155977747088994;
  //AffineParameters[2] = -0.00011447819396342225;
  //AffineParameters[3] = -0.01417572977582712;
  //AffineParameters[4] = 0.991160028135731;
  //AffineParameters[5] = -0.0056367411982091815;
  //AffineParameters[6] = -0.015897718722077483;
  //AffineParameters[7] = -0.028900609861898536;
  //AffineParameters[8] = 0.9891848115955435;
  //AffineParameters[9] = -0.049312842664449213;
  //AffineParameters[10] = -0.051683798376623125;
  //AffineParameters[11] = -0.04527019426983195;
  //ACenter[0] = -1.85716;
  //ACenter[1] = -1.89245;
  //ACenter[2] = 1.58383;
#pragma endregion


  SimilarityTransform->SetParameters(SimilarityParameters);
  SimilarityTransform->SetCenter(SCenter);
  AffineTransform->SetParameters(AffineParameters);
  AffineTransform->SetCenter(ACenter);

#pragma endregion

  //cout << AffineParameters << endl;


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

  MetricType::Pointer         metric        = MetricType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();


  registration->SetMetric(        metric        );
  
  typedef itk::ImageFileReader< FixedImageType  >       FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType >       MovingImageReaderType;

  FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );

  fixedImageReader->Update();
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
  //cout<<BSplineTransform->GetNumberOfParameters()<<endl;
  /*
  typedef BSplineTransformType::ParametersType BSplineParametersType;
  BSplineParametersType BSplineParameters(BSplineTransform->GetNumberOfParameters());
  BSplineParameters.Fill(0.0);
  BSplineTransform->SetParameters(BSplineParameters);*/

  typedef itk::CompositeTransform< double, ImageDimension> CompositeType;
  CompositeType::Pointer compositeTransform = CompositeType::New();
  compositeTransform->AddTransform(SimilarityTransform);
  compositeTransform->AddTransform(AffineTransform);
/*
  std::cout << " Similaritytransform parameters after registration: " << std::endl
	  << SimilarityTransform->GetParameters() << "\n" <<
	  "translations: " << SimilarityTransform->GetTranslation() << "\n" <<
	  "center: " << SimilarityTransform->GetCenter() << "\n" <<
	  "matrix: " << SimilarityTransform->GetMatrix() << "\n" <<
	  "scale: " << SimilarityTransform->GetScale() << "\n" <<
	  "offset: " << SimilarityTransform->GetOffset() << "\n" <<
	  "versor: " << SimilarityTransform->GetVersor() << endl;

  std::cout << " Affine transform parameters after registration: " << std::endl
	  << AffineTransform->GetParameters() << "\n" <<
	  "translations: " << AffineTransform->GetTranslation() << "\n" <<
	  "center: " << AffineTransform->GetCenter() << "\n" <<
	  "matrix: " << AffineTransform->GetMatrix() << "\n" <<
	  "offset: " << AffineTransform->GetOffset() << endl;*/

  std::cout << "The initial composite transform: " << compositeTransform->GetParameters() << std::endl;
  registration->SetMovingInitialTransform(compositeTransform);
  registration->SetInitialTransform(BSplineTransform);
  registration->InPlaceOn();


  registration->SetFixedImage( fixedImage );
  registration->SetMovingImage( movingImageReader->GetOutput() );


  

  unsigned int numberofIterations = 10;
  if (argc >= 8)
  {
	  numberofIterations = std::atof(argv[7]);
  }





  ////setting for LBFGS optimizer
  //typedef itk::LBFGSOptimizerv4                         OptimizerType;
  //OptimizerType::Pointer      optimizer = OptimizerType::New();
  //registration->SetOptimizer(optimizer);
  //typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
  //ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
  //scalesEstimator->SetMetric(metric);
  //scalesEstimator->SetTransformForward(true);
  //scalesEstimator->SetSmallParameterVariation(1.0);
  //optimizer->SetGradientConvergenceTolerance( 1e-6 );
  //optimizer->SetLineSearchAccuracy( 2);
  //optimizer->SetDefaultStepLength( 0.5);
  //optimizer->TraceOn();
  //optimizer->SetNumberOfIterations(numberofIterations);
  ////optimizer->SetMaximumNumberOfFunctionEvaluations(MaximumNumberOfFunctionEvaluations);
  //optimizer->SetScalesEstimator( scalesEstimator );


  //setting for LBFGSB optimizer
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
  ////optimizer->SetCostFunctionConvergenceFactor(1e+12);
  //optimizer->SetGradientConvergenceTolerance(1.0e-8);
  //optimizer->SetNumberOfIterations(numberofIterations);
  //optimizer->TraceOn();
  //optimizer->SetMaximumNumberOfFunctionEvaluations(500);
  //optimizer->SetMaximumNumberOfCorrections(5);

  typedef itk::RegularStepGradientDescentOptimizerv4<double> OptimizerType;
  OptimizerType::Pointer      optimizer = OptimizerType::New();
  registration->SetOptimizer(optimizer);
  optimizer->SetMinimumStepLength(1e-6);
  optimizer->SetRelaxationFactor(0.6);//0.6
  optimizer->SetNumberOfIterations(numberofIterations);
  optimizer->SetLearningRate(0.1); //0.5
  optimizer->SetGradientMagnitudeTolerance(1e-6);
  optimizer->SetReturnBestParametersAndValue(true);

  typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
  ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
  scalesEstimator->SetMetric(metric);
  scalesEstimator->SetTransformForward(true);
  scalesEstimator->SetSmallParameterVariation(1.0);
  optimizer->SetScalesEstimator( scalesEstimator );
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
  smoothingSigmasPerLevel.SetSize( numberOfLevels );
  smoothingSigmasPerLevel[0] = 0;
  /*---Multilevel---*/
  //smoothingSigmasPerLevel[0] = 2;
  //smoothingSigmasPerLevel[1] = 1;
  //smoothingSigmasPerLevel[2] = 0;
  ////smoothingSigmasPerLevel[3] = 0;
  /*----------------*/

  registration->SetNumberOfLevels( numberOfLevels );
  registration->SetSmoothingSigmasPerLevel( smoothingSigmasPerLevel );
  registration->SetShrinkFactorsPerLevel( shrinkFactorsPerLevel );

 // run optimization
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver(itk::IterationEvent(), observer);
   //Add time and memory probes
  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;

  std::cout << std::endl << "Starting Registration" << std::endl;

  try
    {
    memorymeter.Start( "Registration" );
    chronometer.Start( "Registration" );

    registration->Update();

    chronometer.Stop( "Registration" );
    memorymeter.Stop( "Registration" );

    std::cout << "Optimizer stop condition = "
              << registration->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Report the time and memory taken by the registration
  chronometer.Report( std::cout );
  memorymeter.Report( std::cout );


  /*std::cout << "Last BSplineTransform Parameters" << std::endl;
  std::cout << BSplineTransform->GetParameters() << std::endl;
*/
  compositeTransform->AddTransform(BSplineTransform);

  typedef itk::ResampleImageFilter<
                            MovingImageType,
                            FixedImageType >    ResampleFilterType;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( compositeTransform );
  resample->SetInput( movingImageReader->GetOutput() );

  resample->SetSize(   fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin(  fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );
  resample->SetDefaultPixelValue( 128 );

  typedef  unsigned char  OutputPixelType;
  typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;
  typedef itk::CastImageFilter<
                        FixedImageType,
                        OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer      writer =  WriterType::New();
  CastFilterType::Pointer  caster =  CastFilterType::New();

  writer->SetFileName( argv[3] );
  caster->SetInput( resample->GetOutput() );
  writer->SetInput( caster->GetOutput()   );


  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
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
  if( argc >= 5 )
    {
    difference->SetInput1( fixedImageReader->GetOutput() );
    difference->SetInput2( resample->GetOutput() );

    writer2->SetFileName( argv[4] );
    try
      {
      writer2->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
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

  typedef  unsigned char labelMapPixelType;
  typedef itk::Image< labelMapPixelType, ImageDimension > LabelMapImageType;

  typedef itk::ImageFileReader< LabelMapImageType  >       LabelMapImageReaderType;

  LabelMapImageReaderType::Pointer labelMapReader = LabelMapImageReaderType::New();

  typedef itk::ResampleImageFilter<
	  LabelMapImageType,
	  LabelMapImageType >    LabelMapResampleFilterType;

  LabelMapResampleFilterType::Pointer labelMapResampler = LabelMapResampleFilterType::New();

  labelMapResampler->SetTransform(compositeTransform);
  labelMapResampler->SetInput(labelMapReader->GetOutput());
  labelMapResampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
  labelMapResampler->SetOutputOrigin(fixedImage->GetOrigin());
  labelMapResampler->SetOutputSpacing(fixedImage->GetSpacing());
  labelMapResampler->SetOutputDirection(fixedImage->GetDirection());
  labelMapResampler->SetDefaultPixelValue(0);

  typedef itk::Image< labelMapPixelType, ImageDimension >                labelMapImageType;
  typedef itk::CastImageFilter< labelMapImageType, labelMapImageType > labelMapCastFilterType;
  typedef itk::ImageFileWriter< labelMapImageType >                 labelMapWriterType;

  labelMapWriterType::Pointer      labelMapWriter = labelMapWriterType::New();
  labelMapCastFilterType::Pointer  labelMapCaster = labelMapCastFilterType::New();


  labelMapCaster->SetInput(labelMapResampler->GetOutput());
  labelMapWriter->SetInput(labelMapCaster->GetOutput());


  for (size_t j = 0; j < 5; j++)
  {
	  labelMapReader->SetFileName(argv[9 + 2 * j]);
	  labelMapWriter->SetFileName(argv[10 + 2 * j]);
	  try{
		  cout << "writing label map: " << argv[10 + 2 * j] << endl;
		  labelMapWriter->Update(); 
	  }
	  catch (itk::ExceptionObject & err)
	  {
		  std::cout << "ExceptionObject caught !" << std::endl;
		  std::cout << err << std::endl;
		  return EXIT_FAILURE;
	  }

  }

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

  return EXIT_SUCCESS;
}
