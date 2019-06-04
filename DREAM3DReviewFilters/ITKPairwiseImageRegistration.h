/*
 * Your License or Copyright can go here
 */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "itkImageRegistrationMethodv4.h"
#include "SIMPLib/ITK/itkInPlaceDream3DDataToImageFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"




 /**
  * @brief The ITKPairwiseImageRegistration class. See [Filter documentation](@ref itkpairwiseimageregistration) for details.
  */



class DREAM3DReview_EXPORT ITKPairwiseImageRegistration : public AbstractFilter
{
	Q_OBJECT
		// clang-format off
		PYB11_CREATE_BINDINGS(ITKPairwiseImageRegistration SUPERCLASS AbstractFilter)
		PYB11_PROPERTY(DataArrayPath FixedImageArrayPath READ getFixedImageArrayPath WRITE setFixedImageArrayPath)
		PYB11_PROPERTY(DataArrayPath MovingImageArrayPath READ getMovingImageArrayPath WRITE setMovingImageArrayPath)

		// clang-format on

public:
	SIMPL_SHARED_POINTERS(ITKPairwiseImageRegistration)
		SIMPL_FILTER_NEW_MACRO(ITKPairwiseImageRegistration)
		SIMPL_TYPE_MACRO_SUPER(ITKPairwiseImageRegistration, AbstractFilter)

		~ITKPairwiseImageRegistration() override;

	SIMPL_FILTER_PARAMETER(DataArrayPath, FixedImageArrayPath)
		Q_PROPERTY(DataArrayPath FixedImageArrayPath READ getFixedImageArrayPath WRITE setFixedImageArrayPath)

		SIMPL_FILTER_PARAMETER(DataArrayPath, MovingImageArrayPath)
		Q_PROPERTY(DataArrayPath MovingImageArrayPath READ getMovingImageArrayPath WRITE setMovingImageArrayPath)

		SIMPL_FILTER_PARAMETER(int, TransformType)
		Q_PROPERTY(int TransformType READ getTransformType WRITE setTransformType)

		SIMPL_FILTER_PARAMETER(int, Metric)
		Q_PROPERTY(int Metric READ getMetric WRITE setMetric)

		SIMPL_FILTER_PARAMETER(int, Optimizer)
		Q_PROPERTY(int Optimizer READ getOptimizer WRITE setOptimizer)

		SIMPL_FILTER_PARAMETER(int, BSplineOrder)
		Q_PROPERTY(int BSplineOrder READ getBSplineOrder WRITE setBSplineOrder)

		SIMPL_FILTER_PARAMETER(int, MeshSize)
		Q_PROPERTY(int MeshSize READ getMeshSize WRITE setMeshSize)

		SIMPL_FILTER_PARAMETER(int, MaximumNumberOfFunctionEvaluations)
		Q_PROPERTY(int MaximumNumberOfFunctionEvaluations READ getMaximumNumberOfFunctionEvaluations WRITE setMaximumNumberOfFunctionEvaluations)

		SIMPL_FILTER_PARAMETER(int, NumberOfIterations)
		Q_PROPERTY(int NumberOfIterations READ getNumberOfIterations WRITE setNumberOfIterations)

		SIMPL_FILTER_PARAMETER(float, GradientConvergenceTolerance)
		Q_PROPERTY(float GradientConvergenceTolerance READ getGradientConvergenceTolerance WRITE setGradientConvergenceTolerance)

		SIMPL_FILTER_PARAMETER(float, LearningRate)
		Q_PROPERTY(float LearningRate READ getLearningRate WRITE setLearningRate)




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
	ITKPairwiseImageRegistration();

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
	DEFINE_IDATAARRAY_VARIABLE(FixedImage)

	void registerImagePair2D();
	template <typename T> void convertToDouble(std::weak_ptr<IDataArray> imagePtr, DataArrayPath dataArrayPath, QString dataArrayName);

public:
	/* Rule of 5: All special member functions should be defined if any are defined.
	* https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#c21-if-you-define-or-delete-any-default-operation-define-or-delete-them-all
	*/
	ITKPairwiseImageRegistration(const ITKPairwiseImageRegistration&) = delete;             // Copy Constructor Not Implemented
	ITKPairwiseImageRegistration& operator=(const ITKPairwiseImageRegistration&) = delete;  // Copy Assignment Not Implemented
	ITKPairwiseImageRegistration(ITKPairwiseImageRegistration&&) = delete;                  // Move Constructor Not Implemented
	ITKPairwiseImageRegistration& operator=(ITKPairwiseImageRegistration&&) = delete;       // Move Assignment Not Implemented

};

