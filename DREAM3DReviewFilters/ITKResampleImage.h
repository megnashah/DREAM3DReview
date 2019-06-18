/*
 * Your License or Copyright can go here
 */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"
#include "SIMPLib/Geometry/TransformContainer.h"

#include "DREAM3DReview/DREAM3DReviewFilters/util/ITKTransformHelpers.h"
#include "DREAM3DReview/DREAM3DReviewPlugin.h"

/**
 * @brief The ITKResampleImage class. See [Filter documentation](@ref itkresampleimage) for details.
 */
class DREAM3DReview_EXPORT ITKResampleImage : public AbstractFilter
{
  Q_OBJECT
    // clang-format off
    PYB11_CREATE_BINDINGS(ITKResampleImage SUPERCLASS AbstractFilter)
    PYB11_PROPERTY(QString TransformFileName READ getTransformFileName WRITE setTransformFileName)
    PYB11_PROPERTY(int InterpolationType READ getInterpolationType WRITE setInterpolationType)
	PYB11_PROPERTY(DataArrayPath MovingImageArrayPath READ getMovingImageArrayPath WRITE setMovingImageArrayPath)
	PYB11_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)
	PYB11_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)
	PYB11_PROPERTY(QString ImageDataArrayName READ getImageDataArrayName WRITE setImageDataArrayName)


    // clang-format on

  public:
    SIMPL_SHARED_POINTERS(ITKResampleImage)
    SIMPL_FILTER_NEW_MACRO(ITKResampleImage)
    SIMPL_TYPE_MACRO_SUPER(ITKResampleImage, AbstractFilter)

    ~ITKResampleImage() override;

    SIMPL_FILTER_PARAMETER(QString, TransformFileName)
    Q_PROPERTY(QString TransformFileName READ getTransformFileName WRITE setTransformFileName)

    SIMPL_FILTER_PARAMETER(int, InterpolationType)
    Q_PROPERTY(int InterpolationType READ getInterpolationType WRITE setInterpolationType)

	SIMPL_FILTER_PARAMETER(DataArrayPath, MovingImageArrayPath)
	Q_PROPERTY(DataArrayPath MovingImageArrayPath READ getMovingImageArrayPath WRITE setMovingImageArrayPath)

	SIMPL_FILTER_PARAMETER(DataArrayPath, DataContainerName)
	Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

	SIMPL_FILTER_PARAMETER(QString, CellAttributeMatrixName)
	Q_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)

	SIMPL_FILTER_PARAMETER(QString, ImageDataArrayName)
	Q_PROPERTY(QString ImageDataArrayName READ getImageDataArrayName WRITE setImageDataArrayName)


    /**
     * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
     */
    const QString getCompiledLibraryName() const override;

    /**
     * @brief getBrandingString Returns the branding string for the filter, which is a tag
     * used to denote the filter's association with specific plugins
     * @return Branding string
    */
    const QString getBrandingString() const override;

    /**
     * @brief getFilterVersion Returns a version string for this filter. Default
     * value is an empty string.
     * @return
     */
    const QString getFilterVersion() const override;

    /**
     * @brief newFilterInstance Reimplemented from @see AbstractFilter class
     */
    AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

    /**
     * @brief getGroupName Reimplemented from @see AbstractFilter class
     */
    const QString getGroupName() const override;

    /**
     * @brief getSubGroupName Reimplemented from @see AbstractFilter class
     */
    const QString getSubGroupName() const override;

    /**
     * @brief getUuid Return the unique identifier for this filter.
     * @return A QUuid object.
     */
    const QUuid getUuid() override;
  
    /**
     * @brief getHumanLabel Reimplemented from @see AbstractFilter class
     */
    const QString getHumanLabel() const override;

    /**
     * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
     */
    void setupFilterParameters() override;

    /**
     * @brief execute Reimplemented from @see AbstractFilter class
     */
    void execute() override;

    /**
    * @brief preflight Reimplemented from @see AbstractFilter class
    */
    void preflight() override;

  signals:
    /**
     * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
     * be pushed from a user-facing control (such as a widget)
     * @param filter Filter instance pointer 
     */
    void updateFilterParameters(AbstractFilter* filter);

    /**
     * @brief parametersChanged Emitted when any Filter parameter is changed internally
     */
    void parametersChanged();

    /**
     * @brief preflightAboutToExecute Emitted just before calling dataCheck()
     */
    void preflightAboutToExecute();

    /**
     * @brief preflightExecuted Emitted just after calling dataCheck()
     */
    void preflightExecuted();

  protected:
    ITKResampleImage();

    /**
    * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
    */
    void dataCheck();

    /**
    * @brief Initializes all the private instance variables.
    */
    void initialize();
  private:
	  DEFINE_IDATAARRAY_VARIABLE(MovingImage)
	  DEFINE_IDATAARRAY_VARIABLE(FinalImage)

	  template <typename T> void Resample2D();
	  ITKTransformHelpers getTransformAndFixedParams(QString sliceNo);
	  template <typename T> void createCompatibleArrays();
	  
  public:
    /* Rule of 5: All special member functions should be defined if any are defined.
    * https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#c21-if-you-define-or-delete-any-default-operation-define-or-delete-them-all
    */
    ITKResampleImage(const ITKResampleImage&) = delete;             // Copy Constructor Not Implemented
    ITKResampleImage& operator=(const ITKResampleImage&) = delete;  // Copy Assignment Not Implemented
    ITKResampleImage(ITKResampleImage&&) = delete;                  // Move Constructor Not Implemented
    ITKResampleImage& operator=(ITKResampleImage&&) = delete;       // Move Assignment Not Implemented

};

