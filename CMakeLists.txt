CMAKE_MINIMUM_REQUIRED (VERSION 3.11)
PROJECT(VKLab VERSION 0.25.0 LANGUAGES CXX)
#include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
#conan_basic_setup(NO_OUTPUT_DIRS)
SET (CMAKE_CXX_STANDARD 11)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
INCLUDE(GNUInstallDirs)

ADD_SUBDIRECTORY(src)

set(EXAMPLE_DIR examples)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/${EXAMPLE_DIR}
        DESTINATION Examples
		COMPONENT Examples)
		
set(SCRIPTS_DIR scripts)
file(GLOB SCRIPT_FILES "${CMAKE_SOURCE_DIR}/${SCRIPTS_DIR}/*")
install(FILES ${SCRIPT_FILES}
        DESTINATION ${CMAKE_INSTALL_BINDIR}
		COMPONENT applications)
file(WRITE setpath.bat 	"@echo off\necho Setting up environment for OpenMKM\n%comspec% /k \"cd /d \"%HOMEPATH%\" && setx PATH \"%~dp0\"\"")
install(FILES setpath.bat
        DESTINATION ${CMAKE_INSTALL_BINDIR}
		COMPONENT applications)

#
# CPack related stuff to install OpenMKM project defined in src folder
#
set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP TRUE)
include(InstallRequiredSystemLibraries)
install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
        DESTINATION ${CMAKE_INSTALL_BINDIR}
		COMPONENT applications)

set(CPACK_PACKAGE_NAME "OpenMKMInstaller")
set(CPACK_PACKAGE_VENDOR "University of Delaware")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "C++ software for microkinetic modeling")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "3")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenMKM")
set(CPACK_NSIS_MODIFY_PATH "ON")

set(CPACK_COMPONENTS_ALL applications Examples)# CTData instructionData)
set(CPACK_COMPONENT_DATA_GROUP "Examples")
#set(CPACK_COMPONENT_INSTRUCTIONDATA_GROUP "data")

set(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "applications (omkm)")
set(CPACK_COMPONENT_GROUP_DATA_DESCRIPTION "Examples (OpenMKM Examples)")


set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}\\openmkm.png")
set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/openmkm.ico")
set(CPACK_NSIS_MUI_UNIICON "${CMAKE_SOURCE_DIR}/openmkm.ico")

set(CPACK_NSIS_MENU_LINKS "bin/omkm" "OpenMKM"
						  "bin/setpath.bat" "OpenMKM CmdLine" #
                          "Examples" "Examples"
                          #${CT_DATA_FILE_DIR} "CT Data"
                          #"data\\example.pdf" "documentation"
						  )
#list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS " ExecWait ./tmp/vcredist_x64.exe")
#list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS " ExecWait ./tmp/some-other-installer.exe")
#list(APPEND CPACK_NSIS_EXTRA_INSTALL_COMMANDS " RMDir /r ./tmp/")
include(CPack)