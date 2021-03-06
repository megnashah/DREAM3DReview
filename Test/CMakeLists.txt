set(PLUGIN_NAME "DREAM3DReview")

set(${PLUGIN_NAME}Test_SOURCE_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Test)
set(${PLUGIN_NAME}Test_BINARY_DIR ${${PLUGIN_NAME}_BINARY_DIR}/Test)

#------------------------------------------------------------------------------
# List all the source files here. They will NOT be compiled but instead
# be directly included in the main test source file. We list them here so that
# they will show up in IDEs
set(TEST_NAMES
  AnisotropyFilterTest
  ImportMASSIFDataTest
  FFTHDFWriterFilterTest
)

#------------------------------------------------------------------------------
# Include this file from the CMP Project
include(${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake)
include(${SIMPLProj_SOURCE_DIR}/Source/SIMPLib/SIMPLibMacros.cmake)


get_filename_component(${PLUGIN_NAME}_PARENT_SOURCE_DIR "${${PLUGIN_NAME}_SOURCE_DIR}" DIRECTORY)
get_filename_component(${PLUGIN_NAME}_PARENT_BINARY_DIR "${${PLUGIN_NAME}_BINARY_DIR}" DIRECTORY)

set(Anisotropy_Data_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Data/Anisotropy)
list(APPEND ${PLUGIN_NAME}_LINK_LIBS Qt5::Core H5Support SIMPLib ${ITK_LIBRARIES})

SIMPL_GenerateUnitTestFile(PLUGIN_NAME ${PLUGIN_NAME}
                           TEST_DATA_DIR ${${PLUGIN_NAME}_SOURCE_DIR}/Test/Data
                           SOURCES ${TEST_NAMES}
                           LINK_LIBRARIES ${${PLUGIN_NAME}_LINK_LIBS} ${plug_target_name}
                           INCLUDE_DIRS ${${PLUGIN_NAME}_PARENT_SOURCE_DIR}
                                        ${${PLUGIN_NAME}Test_SOURCE_DIR}
                                        ${${PLUGIN_NAME}Test_BINARY_DIR}
                                        ${${PLUGIN_NAME}_PARENT_BINARY_DIR}
)


#------------------------------------------------------------------------------
# If Python is enabled, then enable the Python unit tests for this plugin
if(SIMPL_ENABLE_PYTHON)
    get_property(SIMPL_ANACONDA_OUTPUT_DIR GLOBAL PROPERTY SIMPL_ANACONDA_OUTPUT_DIR)
    get_property(SIMPL_PYTHON_TEST_SCRIPT GLOBAL PROPERTY SIMPL_PYTHON_TEST_SCRIPT)
    get_property(PYTHON_SITE_PACKAGES_NAME GLOBAL PROPERTY PYTHON_SITE_PACKAGES_NAME)

    #------------------------------
    # Initialize the PYTHON_TEST_INPUT_DIR variable to point to the "Testing Directory"
    file(TO_NATIVE_PATH "${${PLUGIN_NAME}Test_SOURCE_DIR}/Python" PYTHON_TEST_INPUT_DIR)
    #------------------------------
    # These names should match the names "EXACTLY" (including capitalization).
    # NO Spaces in the names (which means no spaces in the variable names)
    set(PLUGIN_PYTHON_TESTS
        Apply_Transformation_To_Geometry
        Average_EdgeFaceCell_Array_To_Vertex_Array
        Average_Vertex_Array_To_EdgeFaceCell_Array
        Compute_Umeyama_Transform
        DBSCAN
        Extract_Internal_Surfaces_From_Triangle_Geometry
        Find_Attribute_Array_Statistics
        Find_Element_Centroids
        Find_Norm
        K_Medoids
        K_Distance_Graph
        K_Means_Test
        Normalize_Attribute_Arrays
        Point_Sample_Triangle_Geometry
        Potts_Model
        Principal_Component_Analysis
        Remove_Flagged_Vertices
        Robust_Automatic_Threshold
        # 01_Adaptive_Alignment_Import_Data # The conversion is already complete so this does not need to run
        02_Adaptive_Alignment_Misorientation_Zero_Shifts
        03_Adaptive_Alignment_Mutual_Information_SEM_Images
        04_Adaptive_Alignment_Steiner_Compact
        InsertTransformationPhase
    )

    set(TEST_SCRIPT_FILE_EXT "sh")
    if(WIN32)
      set(TEST_SCRIPT_FILE_EXT "bat")
    endif()
    foreach(test ${PLUGIN_PYTHON_TESTS})
        set(SIMPL_PYTHON_TEST_DRIVER "${SIMPL_ANACONDA_OUTPUT_DIR}/UnitTest/${PLUGIN_NAME}/${test}.${TEST_SCRIPT_FILE_EXT}")

        configure_file(${SIMPL_PYTHON_TEST_SCRIPT}
                        "${SIMPL_PYTHON_TEST_DRIVER}" @ONLY)

        add_test(NAME PY_${PLUGIN_NAME}_${test} COMMAND "${SIMPL_PYTHON_TEST_DRIVER}" )
    endforeach(test ${PLUGIN_PYTHON_TESTS})
    #------------------------------
    # Also setup a unit test for the base python unit test file that is generated as part
    # of the Pybind11 generated codes
    set(PYTHON_TEST_INPUT_DIR "${SIMPL_ANACONDA_OUTPUT_DIR}/UnitTest/${PLUGIN_NAME}")
    set(SIMPL_PYTHON_TEST_DRIVER "${SIMPL_ANACONDA_OUTPUT_DIR}/UnitTest/${PLUGIN_NAME}/${PLUGIN_NAME}.${TEST_SCRIPT_FILE_EXT}")
    set(test "${PLUGIN_NAME}_UnitTest")
    configure_file(${SIMPL_PYTHON_TEST_SCRIPT}  "${SIMPL_PYTHON_TEST_DRIVER}" @ONLY)
    add_test(NAME PY_${PLUGIN_NAME}_UnitTest COMMAND "${SIMPL_PYTHON_TEST_DRIVER}" )
endif()
