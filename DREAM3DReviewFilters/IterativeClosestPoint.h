/* ============================================================================
 * Software developed by US federal government employees (including military personnel)
 * as part of their official duties is not subject to copyright protection and is
 * considered "public domain" (see 17 USC Section 105). Public domain software can be used
 * by anyone for any purpose, and cannot be released under a copyright license
 * (including typical open source software licenses).
 *
 * This source code file was originally written by United States DoD employees. The
 * original source code files are released into the Public Domain.
 *
 * Subsequent changes to the codes by others may elect to add a copyright and license
 * for those changes.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The IterativeClosestPoint class. See [Filter documentation](@ref iterativeclosestpoint) for details.
 */
class DREAM3DReview_EXPORT IterativeClosestPoint : public AbstractFilter
{
  Q_OBJECT

  PYB11_BEGIN_BINDINGS(IterativeClosestPoint SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(IterativeClosestPoint)
  PYB11_FILTER_NEW_MACRO(IterativeClosestPoint)
  PYB11_PROPERTY(DataArrayPath MovingVertexGeometry READ getMovingVertexGeometry WRITE setMovingVertexGeometry)
  PYB11_PROPERTY(DataArrayPath TargetVertexGeometry READ getTargetVertexGeometry WRITE setTargetVertexGeometry)
  PYB11_PROPERTY(int Iterations READ getIterations WRITE setIterations)
  PYB11_PROPERTY(bool ApplyTransform READ getApplyTransform WRITE setApplyTransform)
  PYB11_PROPERTY(QString TransformAttributeMatrixName READ getTransformAttributeMatrixName WRITE setTransformAttributeMatrixName)
  PYB11_PROPERTY(QString TransformArrayName READ getTransformArrayName WRITE setTransformArrayName)
  PYB11_PROPERTY(bool SaveSummedDistance READ getSaveSummedDistance WRITE setSaveSummedDistance)
  PYB11_PROPERTY(QString SummedDistanceAttributeMatrixName READ getSummedDistanceAttributeMatrixName WRITE setSummedDistanceAttributeMatrixName)
  PYB11_PROPERTY(QString SummedDistanceArrayName READ getSummedDistanceArrayName WRITE setSummedDistanceArrayName)

  PYB11_END_BINDINGS()

public:
  using Self = IterativeClosestPoint;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static Pointer New();

  /**
   * @brief Returns the name of the class for IterativeClosestPoint
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for IterativeClosestPoint
   */
  static QString ClassName();

  ~IterativeClosestPoint() override;

  /**
   * @brief Setter property for MovingVertexGeometry
   */
  void setMovingVertexGeometry(const DataArrayPath& value);

  /**
   * @brief Getter property for MovingVertexGeometry
   * @return Value of MovingVertexGeometry
   */
  DataArrayPath getMovingVertexGeometry() const;
  Q_PROPERTY(DataArrayPath MovingVertexGeometry READ getMovingVertexGeometry WRITE setMovingVertexGeometry)

  /**
   * @brief Setter property for TargetVertexGeometry
   */
  void setTargetVertexGeometry(const DataArrayPath& value);

  /**
   * @brief Getter property for TargetVertexGeometry
   * @return Value of TargetVertexGeometry
   */
  DataArrayPath getTargetVertexGeometry() const;
  Q_PROPERTY(DataArrayPath TargetVertexGeometry READ getTargetVertexGeometry WRITE setTargetVertexGeometry)

  /**
   * @brief Setter property for Iterations
   */
  void setIterations(const int& value);

  /**
   * @brief Getter property for Iterations
   * @return Value of Iterations
   */
  int getIterations() const;
  Q_PROPERTY(int Iterations READ getIterations WRITE setIterations)

  /**
   * @brief Setter property for ApplyTransform
   */
  void setApplyTransform(const bool& value);

  /**
   * @brief Getter property for ApplyTransform
   * @return Value of ApplyTransform
   */
  bool getApplyTransform() const;
  Q_PROPERTY(bool ApplyTransform READ getApplyTransform WRITE setApplyTransform)

  /**
   * @brief Setter property for ApplyTransform
   */
  void setTransformAttributeMatrixName(const QString& value);

  /**
   * @brief Getter property for ApplyTransform
   * @return Value of ApplyTransform
   */
  QString getTransformAttributeMatrixName() const;
  Q_PROPERTY(QString TransformAttributeMatrixName READ getTransformAttributeMatrixName WRITE setTransformAttributeMatrixName)

  /**
   * @brief Setter property for ApplyTransform
   */
  void setTransformArrayName(const QString& value);

  /**
   * @brief Getter property for ApplyTransform
   * @return Value of ApplyTransform
   */
  QString getTransformArrayName() const;
  Q_PROPERTY(QString TransformArrayName READ getTransformArrayName WRITE setTransformArrayName)

  /**
   * @brief Setter property for SaveSummedDistance
   */
    void setSaveSummedDistance(const bool& value);

  /**
   * @brief Getter property for SaveSummedDistance
   * @return Value of SaveSummedDistance
   */
  bool getSaveSummedDistance() const;
  Q_PROPERTY(bool SaveSummedDistance READ getSaveSummedDistance WRITE setSaveSummedDistance)

    /**
     * @brief Setter property for SummedDistanceAttributeMatrixName
     */
    void setSummedDistanceAttributeMatrixName(const QString& value);

  /**
   * @brief Getter property for SummedDistanceAttributeMatrixName
   * @return Value of SummedDistanceAttributeMatrixName
   */
  QString getSummedDistanceAttributeMatrixName() const;
  Q_PROPERTY(QString SummedDistanceAttributeMatrixName READ getSummedDistanceAttributeMatrixName WRITE setSummedDistanceAttributeMatrixName)

    /**
     * @brief Setter property for SummedDistanceArrayName
     */
    void setSummedDistanceArrayName(const QString& value);

  /**
   * @brief Getter property for SummedDistanceArrayName
   * @return Value of SummedDistanceArrayName
   */
  QString getSummedDistanceArrayName() const;
  Q_PROPERTY(QString SummedDistanceArrayName READ getSummedDistanceArrayName WRITE setSummedDistanceArrayName)




  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

protected:
  IterativeClosestPoint();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DataArrayPath m_MovingVertexGeometry = {"", "", ""};
  DataArrayPath m_TargetVertexGeometry = {"", "", ""};
  int m_Iterations = {100};
  bool m_ApplyTransform = {false};
  QString m_TransformAttributeMatrixName = {"TransformAttributeMatrix"};
  QString m_TransformArrayName = {"Transform"};
  bool m_SaveSummedDistance = { false };
  QString m_SummedDistanceAttributeMatrixName = { "SummedDistanceAttributeMatrix" };
  QString m_SummedDistanceArrayName = { "SummedDistance" };

  std::weak_ptr<DataArray<float>> m_SummedDistPtr;
  float* m_SummedDist = nullptr;

public:
  IterativeClosestPoint(const IterativeClosestPoint&) = delete;            // Copy Constructor Not Implemented
  IterativeClosestPoint& operator=(const IterativeClosestPoint&) = delete; // Copy Assignment Not Implemented
  IterativeClosestPoint(IterativeClosestPoint&&) = delete;                 // Move Constructor Not Implemented
  IterativeClosestPoint& operator=(IterativeClosestPoint&&) = delete;      // Move Assignment Not Implemented
};
