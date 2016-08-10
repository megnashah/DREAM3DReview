/* ============================================================================
* Copyright (c) 2012 Michael A. Jackson (BlueQuartz Software)
* Copyright (c) 2012 Dr. Michael A. Groeber (US Air Force Research Laboratories)
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
*
* Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright notice, this
* list of conditions and the following disclaimer in the documentation and/or
* other materials provided with the distribution.
*
* Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
* BlueQuartz Software nor the names of its contributors may be used to endorse
* or promote products derived from this software without specific prior written
* permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
* USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*  This code was written under United States Air Force Contract number
*                           FA8650-07-D-5800
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _readmicdata_h_
#define _readmicdata_h_

#include <QtCore/QString>
#include <QtCore/QScopedPointer>
#include <QtCore/QDateTime>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataArrays/StringDataArray.hpp"
#include "SIMPLib/Common/AbstractFilter.h"
#include "SIMPLib/DataContainers/DataContainer.h"

#include "HEDMAnalysisFilters/HEDM/MicPhase.h"

#include "HEDMAnalysis/HEDMAnalysisConstants.h"

class MicReader;

struct Mic_Private_Data
{
    QVector<size_t> dims;
    QVector<float> resolution;
    QVector<float> origin;
    QVector<MicPhase::Pointer> phases;
};

enum MIC_READ_FLAG
{
  MIC_FULL_FILE,
  MIC_HEADER_ONLY
};

// our PIMPL private class
class ReadMicDataPrivate;

/**
* @class ReadMicData ReadMicData.h /FilterCategoryFilters/ReadMicData.h
* @brief
* @author
* @date
* @version 1.0
*/
class ReadMicData : public AbstractFilter
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(ReadMicData)

  public:
    SIMPL_SHARED_POINTERS(ReadMicData)
    SIMPL_STATIC_NEW_MACRO(ReadMicData)
    SIMPL_TYPE_MACRO_SUPER(ReadMicData, AbstractFilter)

    virtual ~ReadMicData();

    SIMPL_FILTER_PARAMETER(QString, DataContainerName)
    Q_PROPERTY(QString DataContainerName READ getDataContainerName WRITE setDataContainerName)

    SIMPL_FILTER_PARAMETER(QString, CellEnsembleAttributeMatrixName)
    Q_PROPERTY(QString CellEnsembleAttributeMatrixName READ getCellEnsembleAttributeMatrixName WRITE setCellEnsembleAttributeMatrixName)

    SIMPL_FILTER_PARAMETER(QString, CellAttributeMatrixName)
    Q_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)

    SIMPL_FILTER_PARAMETER(bool, FileWasRead)
    Q_PROPERTY(bool FileWasRead READ getFileWasRead)

    SIMPL_INSTANCE_STRING_PROPERTY(PhaseNameArrayName)
    SIMPL_INSTANCE_STRING_PROPERTY(MaterialNameArrayName)

    SIMPL_FILTER_PARAMETER(QString, InputFile)
    Q_PROPERTY(QString InputFile READ getInputFile WRITE setInputFile)

    SIMPL_INSTANCE_PROPERTY(QString, CellEulerAnglesArrayName)
    SIMPL_INSTANCE_PROPERTY(QString, CellPhasesArrayName)
    SIMPL_INSTANCE_PROPERTY(QString, CrystalStructuresArrayName)
    SIMPL_INSTANCE_PROPERTY(QString, LatticeConstantsArrayName)

    /**
  * @brief This returns the group that the filter belonds to. You can select
  * a different group if you want. The string returned here will be displayed
  * in the GUI for the filter
  */
    virtual const QString getCompiledLibraryName();
    virtual AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters);
    virtual const QString getGroupName();

    /**
  * @brief This returns the group that the filter belonds to. You can select
  * a different group if you want. The string returned here will be displayed
  * in the GUI for the filter
  */
    virtual const QString getSubGroupName();

    /**
  * @brief This returns a string that is displayed in the GUI. It should be readable
  * and understandable by humans.
  */
    virtual const QString getHumanLabel();

    /**
  * @brief This method will instantiate all the end user settable options/parameters
  * for this filter
  */
    virtual void setupFilterParameters();

    /**
  * @brief This method will read the options from a file
  * @param reader The reader that is used to read the options from a file
  */
    virtual void readFilterParameters(AbstractFilterParametersReader* reader, int index);

    /**
  * @brief Reimplemented from @see AbstractFilter class
  */
    virtual void execute();

    /**
  * @brief This function runs some sanity checks on the DataContainer and inputs
  * in an attempt to ensure the filter can process the inputs.
  */
    virtual void preflight();

    /* These are non-exposed to the user through the GUI. Manual Pipelines are OK to set them */
    SIMPL_INSTANCE_PROPERTY(uint32_t, RefFrameZDir)
    SIMPL_INSTANCE_PROPERTY(int, Manufacturer)

    SIMPL_PIMPL_PROPERTY_DECL(QString, InputFile_Cache)
    SIMPL_PIMPL_PROPERTY_DECL(QDateTime, TimeStamp_Cache)
    SIMPL_PIMPL_PROPERTY_DECL(Mic_Private_Data, Data)
    Q_PROPERTY(Mic_Private_Data Data READ getData WRITE setData)


    public slots:
      void flushCache();


  signals:
    void updateFilterParameters(AbstractFilter* filter);
    void parametersChanged();
    void preflightAboutToExecute();
    void preflightExecuted();

  protected:
    ReadMicData();

    /**
     * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
     */
    void dataCheck();

    /**
     * @brief Initializes all the private instance variables.
     */
    void initialize();

    /**
    * @brief readMicFile This reads the HEDM .mic file and puts the data into the Voxel Data container
    */
    void readMicFile();

    /**
    * @brief This method reads the values for the phase type, crystal structure
    * and precipitate fractions from the EBSD file.
    * @param reader The EbsdReader instance
    * @return Zero/Positive on Success - Negative on error.
    */
    int loadMaterialInfo(MicReader* reader);

    /**
     * @brief populateMicData
     * @param reader
     * @param m
     * @param dims
     */
    void populateMicData(MicReader* reader, DataContainer::Pointer m, QVector<size_t> dims, MIC_READ_FLAG = MIC_FULL_FILE);

  private:
    QScopedPointer<ReadMicDataPrivate> const d_ptr;

    DEFINE_DATAARRAY_VARIABLE(int32_t, CellPhases)
    DEFINE_DATAARRAY_VARIABLE(float, CellEulerAngles)
    DEFINE_DATAARRAY_VARIABLE(uint32_t, CrystalStructures)
    DEFINE_DATAARRAY_VARIABLE(float, LatticeConstants)

    ReadMicData(const ReadMicData&); // Copy Constructor Not Implemented
    void operator=(const ReadMicData&); // Operator '=' Not Implemented
};

Q_DECLARE_METATYPE(Mic_Private_Data)


#endif /* _ReadMicData_H_ */





