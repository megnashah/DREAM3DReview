/*
 * Your License or Copyright can go here
 */

#include "ITKResampleImage.h"

#include "SIMPLib/Common/Constants.h"

#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/ChoiceFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedPathCreationFilterParameter.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"
#include "SIMPLib/ITK/itkInPlaceImageToDream3DDataFilter.h"
#include "itkResampleImageFilter.h"
#include "SIMPLib/Geometry/TransformContainer.h"
#include "H5Support/H5ScopedSentinel.h"
#include "H5Support/QH5Utilities.h"
#include "itkBSplineTransform.h"
#include "itkBSplineTransformInitializer.h"
#include "itkCenteredTransformInitializer.h"
#include "itkAffineTransform.h"
#include "itkEuler2DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"
#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"
#include "DREAM3DReview/DREAM3DReviewFilters/util/ITKTransformHelpers.h"



enum createdPathID : RenameDataPath::DataID_t
{
	AttributeMatrixID21 = 21,
	
	DataArrayID31 = 31,
	
	DataContainerID = 1
};



// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKResampleImage::ITKResampleImage() 
:  m_TransformFileName("")
,  m_InterpolationType(0)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ITKResampleImage::~ITKResampleImage() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::initialize()
{
  clearErrorCode();
  clearWarningCode();
  setCancel(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  parameters.push_back(SIMPL_NEW_INPUT_FILE_FP("Transform File Name", TransformFileName, FilterParameter::Parameter, ITKResampleImage, "*.hdf5"));
  {
    ChoiceFilterParameter::Pointer parameter = ChoiceFilterParameter::New();
    parameter->setHumanLabel("Interpolation Type");
    parameter->setPropertyName("InterpolationType");
    
    QVector<QString> choices;          // Please add choices to the choices QVector to finish this widget
	choices.push_back("Linear");
	choices.push_back("Nearest Neighbor");
	parameter->setChoices(choices);
    parameter->setChoices(choices);
    parameter->setCategory(FilterParameter::Parameter);
    parameter->setSetterCallback(SIMPL_BIND_SETTER(ITKResampleImage, this, InterpolationType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(ITKResampleImage, this, InterpolationType));
    parameters.push_back(parameter);
  }

  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, SIMPL::Defaults::AnyComponentSize, AttributeMatrix::Type::Cell, IGeometry::Type::Image);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Moving Image", MovingImageArrayPath, FilterParameter::RequiredArray, ITKResampleImage, dasReq, 0));

  parameters.push_back(SIMPL_NEW_DC_CREATION_FP("Data Container", DataContainerName, FilterParameter::CreatedArray, ITKResampleImage));
  parameters.push_back(SeparatorFilterParameter::New("Cell Data", FilterParameter::CreatedArray));
  parameters.push_back(SIMPL_NEW_AM_WITH_LINKED_DC_FP("Cell Attribute Matrix", CellAttributeMatrixName, DataContainerName, FilterParameter::CreatedArray, ITKResampleImage));
  parameters.push_back(SIMPL_NEW_DA_WITH_LINKED_AM_FP("Image Data", ImageDataArrayName, DataContainerName, CellAttributeMatrixName, FilterParameter::CreatedArray, ITKResampleImage));


  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::dataCheck()
{
  clearErrorCode();
  clearWarningCode();

  m_MovingImagePtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath());
  if (getErrorCode() < 0)
  {
	  return;
  }

  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getDataContainerName(), DataContainerID);
  if (getErrorCode() < 0)
  {
	  return;
  }

  // Create the Image Geometry
  ImageGeom::Pointer image = ImageGeom::CreateGeometry(SIMPL::Geometry::ImageGeometry);
  m->setGeometry(image);


  QVector<size_t> tDims(3, 0);
  AttributeMatrix::Pointer cellAttrMat = m->createNonPrereqAttributeMatrix(this, getCellAttributeMatrixName(), tDims, AttributeMatrix::Type::Cell, AttributeMatrixID21);
  if (getErrorCode() < 0)
  {
	  return;
  }

  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, createCompatibleArrays, m_MovingImagePtr.lock())


  
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true); // Set the fact that we are preflighting.
  emit preflightAboutToExecute(); // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck(); // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted(); // We are done preflighting this filter
  setInPreflight(false); // Inform the system this filter is NOT in preflight mode anymore.
}

template <typename T> void ITKResampleImage::createCompatibleArrays()
{
	

	DataArrayPath path(m_DataContainerName.getDataContainerName(), m_CellAttributeMatrixName, m_ImageDataArrayName);


	m_FinalImagePtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<T>, AbstractFilter, T>(this, path, 0, m_MovingImagePtr.lock()->getComponentDimensions());
	
}


#define BSPLINESETUP()\
typedef itk::BSplineTransform<double, ImageDimension, SplineOrder> TransformType; \
TransformType::Pointer transform = TransformType::New();\
TransformType::MeshSizeType meshSize;\
meshSize.Fill(3);\
transform->SetTransformDomainMeshSize(meshSize);\
int numFixedElements = d3dtransform->getFixedParameters().size();\
TransformType::FixedParametersType fixedParameters(numFixedElements);\
for (int i = 0; i < numFixedElements; i++)\
{\
	fixedParameters.SetElement(i, d3dtransform->getFixedParameters()[i]);\
}\
transform->SetFixedParameters(fixedParameters);\
int numLearnedElements = d3dtransform->getParameters().size();\
TransformType::ParametersType learnedParameters(numLearnedElements);\
for (int i = 0; i < numLearnedElements; i++)\
{\
	learnedParameters.SetElement(i, d3dtransform->getParameters()[i]);\
}\
transform->SetParameters(learnedParameters);\
resample->SetTransform(transform);





template <typename T> void ITKResampleImage::Resample2D()
{
	const unsigned int ImageDimension = 2;
	const unsigned int SpaceDimension = 2;

	using PixelType = itk::Vector<T, 3>;

	typedef itk::InPlaceDream3DDataToImageFilter<PixelType, ImageDimension> ToITKType;
	ToITKType::Pointer movingtoITK = ToITKType::New();


	movingtoITK->SetDataArrayName(m_MovingImageArrayPath.getDataArrayName().toStdString());
	movingtoITK->SetAttributeMatrixArrayName(m_MovingImageArrayPath.getAttributeMatrixName().toStdString());
	movingtoITK->SetInput(getDataContainerArray()->getDataContainer(m_MovingImageArrayPath.getDataContainerName()));
	movingtoITK->InPlaceOn();
	movingtoITK->Update();
	
	typedef itk::Dream3DImage<PixelType, ImageDimension> ImageType;

	ImageType::Pointer itkMovingImage = movingtoITK->GetOutput();

	typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	
	ITKTransformHelpers transformhelper = getTransformAndFixedParams("0"); 
	TransformContainer::Pointer d3dtransform = transformhelper.transform;

	d3dtransform->getFixedParameters();
	QString stringTransformType = QString::fromStdString(d3dtransform->getTransformTypeAsString()).split("_")[0];

	if (stringTransformType == "BSplineTransform")
	{
		QString transformOrderString = QString::fromStdString(d3dtransform->getTransformTypeAsString()); 

		if (transformOrderString == "BSplineTransform_double_2_2")
		{
			const unsigned int SplineOrder = 3;
			BSPLINESETUP()
		}

		if (transformOrderString == "BSplineTransform_double_2_2_2")
		{
			const unsigned int SplineOrder = 2;
			BSPLINESETUP()
		}

		if (transformOrderString == "BSplineTransform_double_2_2_1")
		{
			const unsigned int SplineOrder = 1;
			BSPLINESETUP()

		}

		resample->SetInput(itkMovingImage);

		resample->SetSize(transformhelper.FixedSize);
		resample->SetOutputOrigin(transformhelper.FixedOrigin);
		resample->SetOutputSpacing(transformhelper.FixedSpacing);
		resample->SetOutputDirection(transformhelper.FixedDirection);
		//resample->SetDefaultPixelValue(0);

		if (m_InterpolationType == 1)
		{
			typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
			InterpolatorType::Pointer interpolator = InterpolatorType::New();
			resample->SetInterpolator(interpolator);
		}


		resample->Update();
		ImageType::Pointer output = resample->GetOutput();
		ImageType::SizeType size = output->GetLargestPossibleRegion().GetSize();

		typedef itk::InPlaceImageToDream3DDataFilter<PixelType, ImageDimension> ToD3DType;
		ToD3DType::Pointer movedD3D = ToD3DType::New();

		DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(getDataContainerName());


		AttributeMatrix::Pointer am = dc->getAttributeMatrix(m_CellAttributeMatrixName); 
		QVector<size_t> tdims(3, 0); 
		tdims[0] = size[0]; 
		tdims[1] = size[1]; 
		tdims[2] = 1; 
		am->setTupleDimensions(tdims);
		
		movedD3D->SetDataContainer(dc); 
		movedD3D->SetDataArrayName(m_ImageDataArrayName.toStdString());
		movedD3D->SetAttributeMatrixArrayName(m_CellAttributeMatrixName.toStdString()); 
		movedD3D->InPlaceOn();
		movedD3D->SetInput(output);
		movedD3D->Update();
	}

	if (stringTransformType == "Affine")
	{

	}

	if (stringTransformType == "Euler2DTransform_double_2_2")
	{

	}






}



ITKTransformHelpers ITKResampleImage::getTransformAndFixedParams(QString sliceNo)
{
	hid_t fileId = QH5Utilities::openFile(m_TransformFileName, true);
	ITKTransformHelpers transformhelper(fileId, sliceNo.toStdString());
	transformhelper.readFixedInformation();
	transformhelper.readTransform();

	return transformhelper;
}





// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ITKResampleImage::execute()
{
  initialize();
  dataCheck();
  if(getErrorCode() < 0) { return; }

  if (getCancel()) { return; }

  if (getWarningCode() < 0)
  {
    QString ss = QObject::tr("Some warning message");
    setWarningCondition(-88888888, ss);
  }

  if (getErrorCode() < 0)
  {
    QString ss = QObject::tr("Some error message");
    setErrorCondition(-99999999, ss);
    return;
  }

  m_MovingImagePtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getMovingImageArrayPath()); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */



  EXECUTE_FUNCTION_TEMPLATE_NO_BOOL(this, Resample2D, m_MovingImagePtr.lock())
  




  
  	 




  //typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  //ResampleFilterType::Pointer resample = ResampleFilterType::New();
  //resample->SetTransform(transform);
  //resample->SetInput(m_MovingITKImage);
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
AbstractFilter::Pointer ITKResampleImage::newFilterInstance(bool copyFilterParameters) const
{
  ITKResampleImage::Pointer filter = ITKResampleImage::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getCompiledLibraryName() const
{ 
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream <<  DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getGroupName() const
{ 
  return SIMPL::FilterGroups::Unsupported; 
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getSubGroupName() const
{ 
  return "DREAM3DReview"; 
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ITKResampleImage::getHumanLabel() const
{ 
  return "ITKResampleImage"; 
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ITKResampleImage::getUuid()
{
  return QUuid("{7d478bf6-1acc-5e49-86b6-45c94776bc48}");
}

