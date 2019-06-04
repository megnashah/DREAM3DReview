/*
 * Your License or Copyright can go here
 */

#include "ITKPairwiseImageRegistration.h"
#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/IntFilterParameter.h"
#include "SIMPLib/FilterParameters/FloatFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/ITK/itkInPlaceImageToDream3DDataFilter.h"
#include "itkBSplineTransform.h"
#include "itkBSplineTransformInitializer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkEuler2DTransform.h"
#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkLBFGSOptimizerv4.h"
#include "itkAmoebaOptimizerv4.h"
#include "itkGradientDescentLineSearchOptimizerv4.h"
#include "itkMattesMutualInformationImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkCorrelationImageToImageMetricv4.h"
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"
#include "SIMPLib/Common/TemplateHelpers.h"

#include "SIMPLib/Geometry/ImageGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

 // -----------------------------------------------------------------------------
 //
 // -----------------------------------------------------------------------------
ITKPairwiseImageRegistration::ITKPairwiseImageRegistration()
	: m_FixedImageArrayPath("", "", "")
	, m_MovingImageArrayPath("", "", "")
	, m_TransformType(0)
	, m_BSplineOrder(3)
	, m_MeshSize(3)
	, m_MaximumNumberOfFunctionEvaluations(1000)
	, m_NumberOfIterations(100)
	, m_GradientConvergenceTolerance(1e-5)
	, m_MovingImage(nullptr)
	, m_FixedImage(nullptr)
{
	initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKPairwiseImageRegistration::~ITKPairwiseImageRegistration() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::initialize()
{
	clearErrorCode();
	clearWarningCode();
	setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::setupFilterParameters()
{
	FilterParameterVectorType parameters;
	{
		LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
		parameter->setHumanLabel("Transformation Type");
		parameter->setPropertyName("TransformType");
		parameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, TransformType));
		parameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, TransformType));

		QVector<QString> choices;
		choices.push_back("BSpline");
		choices.push_back("Rigid");
		choices.push_back("Affine");
		parameter->setChoices(choices);

		QStringList linkedProps;
		linkedProps << "MeshSize"
			<< "BSplineOrder";
		parameter->setLinkedProperties(linkedProps);
		parameter->setEditable(false);
		parameter->setCategory(FilterParameter::Parameter);
		parameters.push_back(parameter);
	}

	ChoiceFilterParameter::Pointer metricparameter = ChoiceFilterParameter::New();
	metricparameter->setHumanLabel("Evaluation Metric");
	metricparameter->setPropertyName("Metric");
	metricparameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, Metric));
	metricparameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, Metric));

	QVector<QString> metricchoices;
	metricchoices.push_back("Mattes Mutual Information");
	metricchoices.push_back("Mean Squares");
	metricchoices.push_back("Normalized Cross Correlation");
	metricparameter->setChoices(metricchoices);
	metricparameter->setEditable(false);
	metricparameter->setCategory(FilterParameter::Parameter);
	parameters.push_back(metricparameter);

	LinkedChoicesFilterParameter::Pointer optimizerparameter = LinkedChoicesFilterParameter::New();
	optimizerparameter->setHumanLabel("Optimizer");
	optimizerparameter->setPropertyName("Optimizer");
	optimizerparameter->setSetterCallback(SIMPL_BIND_SETTER(ITKPairwiseImageRegistration, this, Optimizer));
	optimizerparameter->setGetterCallback(SIMPL_BIND_GETTER(ITKPairwiseImageRegistration, this, Optimizer));

	QVector<QString> optimizerchoices;
	optimizerchoices.push_back("LBFGS");
	optimizerchoices.push_back("Amoeba");
	optimizerchoices.push_back("Gradient Descent");
	optimizerparameter->setChoices(optimizerchoices);

	QStringList optimizerlinkedProps;
	optimizerlinkedProps << "GradientConvergenceTolerance"
		<< "MaximumNumberOfFunctionEvaluations"
		<<"LearningRate";
	optimizerparameter->setLinkedProperties(optimizerlinkedProps);
	optimizerparameter->setEditable(false);
	optimizerparameter->setCategory(FilterParameter::Parameter);
	parameters.push_back(optimizerparameter);



	DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
	parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Fixed Image", FixedImageArrayPath, FilterParameter::RequiredArray, ITKPairwiseImageRegistration, dasReq));
	parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Moving Image", MovingImageArrayPath, FilterParameter::RequiredArray, ITKPairwiseImageRegistration, dasReq));
	parameters.push_back(SIMPL_NEW_INTEGER_FP("B-Spline Order", BSplineOrder, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
	parameters.push_back(SIMPL_NEW_INTEGER_FP("Mesh Size", MeshSize, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
	parameters.push_back(SIMPL_NEW_FLOAT_FP("Gradient Convergence Tolerance", GradientConvergenceTolerance, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
	parameters.push_back(SIMPL_NEW_INTEGER_FP("Num Max Function Evaluations", MaximumNumberOfFunctionEvaluations, FilterParameter::Parameter, ITKPairwiseImageRegistration, 0));
	parameters.push_back(SIMPL_NEW_INTEGER_FP("Num Iterations", NumberOfIterations, FilterParameter::Parameter, ITKPairwiseImageRegistration));
	parameters.push_back(SIMPL_NEW_FLOAT_FP("Learning Rate", LearningRate, FilterParameter::Parameter, ITKPairwiseImageRegistration, 2));





	setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::dataCheck()
{
	clearErrorCode();
	clearWarningCode();

	QVector<size_t> cDims(1, 1);
	ImageGeom::Pointer fixedImage = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getFixedImageArrayPath().getDataContainerName());
	ImageGeom::Pointer movingImage = getDataContainerArray()->getPrereqGeometryFromDataContainer<ImageGeom, AbstractFilter>(this, getMovingImageArrayPath().getDataContainerName());

	m_FixedImagePtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getFixedImageArrayPath());
	m_MovingImagePtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */


	if (m_BSplineOrder > 3)
	{
		QString ss = QObject::tr("The largest B-spline order is 3rd Order");
		setErrorCondition(-5555, ss);
	}

	if (m_MeshSize < m_BSplineOrder)
	{
		QString ss = QObject::tr("The mesh size cannot be smaller than the B-spline order");
		setErrorCondition(-5555, ss);
	}
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::preflight()
{
	// These are the REQUIRED lines of CODE to make sure the filter behaves correctly
	setInPreflight(true);              // Set the fact that we are preflighting.
	emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
	emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
	dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
	emit preflightExecuted();          // We are done preflighting this filter
	setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void ITKPairwiseImageRegistration::convertToDouble(std::weak_ptr<IDataArray> imagePtr, DataArrayPath dataArrayPath, QString dataArrayName)
{
	typename DataArray<T>::Pointer FixedImagePtr = std::dynamic_pointer_cast<DataArray<T>>(imagePtr.lock());
	//typename DataArray<T>::Pointer MovingImagePtr = std::dynamic_pointer_cast<DataArray<T>>(m_MovingImagePtr.lock());

	DataArray<double>::Pointer FixedImageDoublePtr = DataArray<double>::CreateArray(FixedImagePtr->getNumberOfTuples(), dataArrayName);
	//DataArray<double>::Pointer MovingImageDoublePtr = DataArray<double>::CreateArray(MovingImagePtr->getNumberOfTuples(), "_INTERNAL_USE_ONLY_MovingImageDouble");

	T* FixedImageDataPtr = FixedImagePtr->getPointer(0);
	//T* MovingImageDataPtr = MovingImagePtr->getPointer(0);

	double* FixedImageDoubleDataPtr = FixedImageDoublePtr->getPointer(0);
	//double* MovingImageDoubleDataPtr = MovingImageDoublePtr->getPointer(0);

	for (size_t i = 0; i < FixedImagePtr->getNumberOfTuples(); i++)
	{
		FixedImageDoubleDataPtr[i] = static_cast<double>(FixedImageDataPtr[i]);
	}

	//for (size_t i = 0; i < MovingImagePtr->getNumberOfTuples(); i++)
	//{
	//	MovingImageDoubleDataPtr[i] = static_cast<double>(MovingImageDataPtr[i]);
	//}

	getDataContainerArray()->getDataContainer(dataArrayPath.getDataContainerName())->getAttributeMatrix(dataArrayPath.getAttributeMatrixName())->addOrReplaceAttributeArray(FixedImageDoublePtr);
		
	//getDataContainerArray()->getDataContainer(m_MovingImageArrayPath.getDataContainerName())->getAttributeMatrix(m_MovingImageArrayPath.getAttributeMatrixName())->addAttributeArray("_INTERNAL_USE_ONLY_MovingImageDouble", MovingImageDoublePtr);

}

#define CHOOSE_LBFGS_OPTIMIZER()											\
typedef itk::LBFGSOptimizerv4 OptimizerType;							\
OptimizerType::Pointer optimizer = OptimizerType::New();				\
registration->SetOptimizer(optimizer);									\
typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;\
ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New(); \
scalesEstimator->SetMetric(metric); \
scalesEstimator->SetTransformForward(true);\
scalesEstimator->SetSmallParameterVariation(1.0);\
optimizer->SetGradientConvergenceTolerance(m_GradientConvergenceTolerance);\
optimizer->SetMaximumNumberOfFunctionEvaluations(m_MaximumNumberOfFunctionEvaluations);\
optimizer->SetScalesEstimator(scalesEstimator);\
optimizer->SetNumberOfIterations(m_NumberOfIterations);

#define CHOOSE_AMOEBA_OPTIMIZER() \
typedef itk::AmoebaOptimizerv4 OptimizerType;			\
OptimizerType::Pointer optimizer = OptimizerType::New();\
registration->SetOptimizer(optimizer);\
typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;\
ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();\
scalesEstimator->SetMetric(metric);\
scalesEstimator->SetTransformForward(true);\
scalesEstimator->SetSmallParameterVariation(1.0);\
optimizer->SetScalesEstimator(scalesEstimator);\
optimizer->SetNumberOfIterations(m_NumberOfIterations);

#define CHOOSE_GRADIENT_LINE_SEARCH_OPTIMIZER()\
typedef itk::GradientDescentLineSearchOptimizerv4 OptimizerType;\
OptimizerType::Pointer optimizer = OptimizerType::New();\
registration->SetOptimizer(optimizer);\
typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;\
ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();\
scalesEstimator->SetMetric(metric);\
scalesEstimator->SetTransformForward(true);\
scalesEstimator->SetSmallParameterVariation(1.0);\
optimizer->SetScalesEstimator(scalesEstimator);\
optimizer->SetNumberOfIterations(m_NumberOfIterations);

#define CHOOSE_METRIC(metricv4)\
typedef itk::metricv4<ImageType, ImageType> MetricType;\
MetricType::Pointer metric = MetricType::New();\
registration->SetMetric(metric);\
if (m_Optimizer == 0)\
{\
	CHOOSE_LBFGS_OPTIMIZER();\
}\
if (m_Optimizer == 1)\
{\
	CHOOSE_AMOEBA_OPTIMIZER();\
}\
if (m_Optimizer == 2)\
{\
	CHOOSE_GRADIENT_LINE_SEARCH_OPTIMIZER();\
}










// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::registerImagePair2D()
{


	const unsigned int ImageDimension = 2;
	const unsigned int SpaceDimension = 2;


	typedef itk::Dream3DImage<double, ImageDimension> ImageType;

	typedef itk::InPlaceDream3DDataToImageFilter<double, ImageDimension> ToITKType;
	ToITKType::Pointer fixedtoITK = ToITKType::New();
	ToITKType::Pointer movingtoITK = ToITKType::New();


	fixedtoITK->SetDataArrayName("_INTERNAL_USE_ONLY_FixedImageDouble");
	fixedtoITK->SetAttributeMatrixArrayName(m_FixedImageArrayPath.getAttributeMatrixName().toStdString());
	fixedtoITK->SetInput(getDataContainerArray()->getDataContainer(m_FixedImageArrayPath.getDataContainerName()));
	fixedtoITK->InPlaceOn();
	fixedtoITK->Update();

	ImageType::Pointer itkFixedImage = fixedtoITK->GetOutput();

	movingtoITK->SetDataArrayName("_INTERNAL_USE_ONLY_MovingImageDouble");
	movingtoITK->SetAttributeMatrixArrayName(m_MovingImageArrayPath.getAttributeMatrixName().toStdString());
	movingtoITK->SetInput(getDataContainerArray()->getDataContainer(m_MovingImageArrayPath.getDataContainerName()));
	movingtoITK->InPlaceOn();
	movingtoITK->Update();

	ImageType::Pointer itkMovingImage = movingtoITK->GetOutput();

	//////////////////////////////// SETUP ORIGIN _ RESOLUTION //////////////////////////////////////////////////

	ImageType::PointType origin = (0, 0);
	itkFixedImage->SetSpacing(1);
	itkFixedImage->SetOrigin(origin);

	itkMovingImage->SetSpacing(1);
	itkMovingImage->SetOrigin(origin);

	//////////////////////////////// REGISTRATION //////////////////////////////////////////////////

	typedef itk::ImageRegistrationMethodv4<ImageType, ImageType> ImageRegistrationType;
	ImageRegistrationType::Pointer registration = ImageRegistrationType::New();

	registration->InPlaceOn();
	registration->SetFixedImage(itkFixedImage);
	registration->SetMovingImage(itkMovingImage);

	//////////////////////////////// TRANSFORM //////////////////////////////////////////////////

	typedef double CoordinateRepType;
	typedef itk::Dream3DImage<double, ImageDimension> ImageType;

	if (m_TransformType == 0)
	{
		const unsigned int SplineOrder = 3;
		typedef itk::BSplineTransform<CoordinateRepType, ImageDimension, SplineOrder> TransformType;
		TransformType::Pointer transform = TransformType::New();
		typedef itk::BSplineTransformInitializer<TransformType, ImageType> TransformInitializerType;
		TransformInitializerType::Pointer transformInitializer = TransformInitializerType::New();

		// unsigned int numberOfGridNodesInOneDimension = 6;
		TransformType::MeshSizeType meshSize;
		meshSize.Fill(m_MeshSize);

		transformInitializer->SetTransform(transform);
		transformInitializer->SetImage(itkFixedImage);
		transformInitializer->SetTransformDomainMeshSize(meshSize);
		transformInitializer->InitializeTransform();

		transform->SetIdentity();

		registration->SetInitialTransform(transform);
	}

	if (m_TransformType == 1)
	{
		typedef itk::Euler2DTransform<double> TransformType; 
		TransformType::Pointer transform = TransformType::New(); 
		typedef itk::CenteredTransformInitializer<TransformType, ImageType, ImageType> TranformInitializerType; 
		TranformInitializerType::Pointer transformInitializer = TranformInitializerType::New();
		transformInitializer->SetTransform(transform);
		transformInitializer->SetFixedImage(itkFixedImage);
		transformInitializer->SetMovingImage(itkMovingImage);
		transformInitializer->GeometryOn(); 
		transformInitializer->InitializeTransform();

		registration->SetInitialTransform(transform);

	}

	if (m_TransformType == 2)
	{
		typedef itk::TranslationTransform<double, ImageDimension> TransformType; 
		TransformType::Pointer transform = TransformType::New();

		TransformType::ParametersType initialParameters(transform->GetNumberOfParameters()); 
		initialParameters[0] = 0; 
		initialParameters[1] = 0;
		
	}







	//////////////////////////////// METRIC //////////////////////////////////////////////////

	//For each metric, an optimizer needs to be set. The optimizer needs to know the type of the metric, so there is no easy way to abstract this. 
	//Instead each case must be explicitly written out. Here we use macros to help. 

	typedef itk::Dream3DImage<double, 2> ImageType;


	if (m_Metric == 0)
	{
		CHOOSE_METRIC(MattesMutualInformationImageToImageMetricv4);
	}


	if (m_Metric == 1)
	{
		CHOOSE_METRIC(MeanSquaresImageToImageMetricv4);

	}

	if (m_Metric == 2)
	{
		CHOOSE_METRIC(CorrelationImageToImageMetricv4);
	}


	//////////////////////////////// MUTLI-RESOLUTION //////////////////////////////////////////////////

	int numberOfLevels = 3;
	ImageRegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(3);
	shrinkFactorsPerLevel[0] = 6;
	shrinkFactorsPerLevel[1] = 2;
	shrinkFactorsPerLevel[2] = 1;
	ImageRegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(3);
	smoothingSigmasPerLevel[0] = 6;
	smoothingSigmasPerLevel[1] = 2;
	smoothingSigmasPerLevel[2] = 1;
	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);


	//////////////////////////////// RUN REGISTRATION //////////////////////////////////////////////////

	registration->Update();

	TransformType::ParametersType finalParameters = registration->GetTransform()->GetParameters();





	//typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
	//ResampleFilterType::Pointer resample = ResampleFilterType::New();
	//resample->SetTransform(transform);
	//resample->SetInput(itkMovingImage);
	//resample->SetSize(itkFixedImage->GetLargestPossibleRegion().GetSize());
	//resample->SetOutputOrigin(itkFixedImage->GetOrigin());
	//resample->SetOutputSpacing(itkFixedImage->GetSpacing());
	//resample->SetOutputDirection(itkFixedImage->GetDirection());
	//resample->SetDefaultPixelValue(0);

	//typedef unsigned char OutputPixelType;

	//typedef itk::Image<OutputPixelType, ImageDimension> OutputImageType;

	//typedef itk::CastImageFilter<ImageType, OutputImageType> CastFilterType;

	//typedef itk::ImageFileWriter<OutputImageType> OutputWriterType;

	//OutputWriterType::Pointer writer = OutputWriterType::New();
	//CastFilterType::Pointer caster = CastFilterType::New();

	//writer->SetFileName("C:/Users/shahmn/Documents/test.png");

	//caster->SetInput(resample->GetOutput());
	//writer->SetInput(caster->GetOutput());

	//writer->Update();



}





// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKPairwiseImageRegistration::execute()
{
	initialize();
	dataCheck();
	if (getErrorCode() < 0)
	{
		return;
	}

	if (getCancel())
	{
		return;
	}
	EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, convertToDouble, m_FixedImagePtr.lock(), m_FixedImagePtr, m_FixedImageArrayPath, "_INTERNAL_USE_ONLY_FixedImageDouble");
	EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, convertToDouble, m_MovingImagePtr.lock(), m_MovingImagePtr, m_MovingImageArrayPath, "_INTERNAL_USE_ONLY_MovingImageDouble");
	registerImagePair2D();





}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ITKPairwiseImageRegistration::newFilterInstance(bool copyFilterParameters) const
{
	ITKPairwiseImageRegistration::Pointer filter = ITKPairwiseImageRegistration::New();
	if (copyFilterParameters)
	{
		copyFilterParameterInstanceVariables(filter.get());
	}
	return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getCompiledLibraryName() const
{
	return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getBrandingString() const
{
	return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getFilterVersion() const
{
	QString version;
	QTextStream vStream(&version);
	vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
	return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getGroupName() const
{
	return SIMPL::FilterGroups::Unsupported;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getSubGroupName() const
{
	return DREAM3DReviewConstants::FilterSubGroups::RegistrationFilters;

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKPairwiseImageRegistration::getHumanLabel() const
{
	return "ITKPairwiseImageRegistration";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ITKPairwiseImageRegistration::getUuid()
{
	return QUuid("{0a283663-e1c8-596d-bcd0-fd27b747551d}");
}
