#
# Rather than bundling the code, we fetch it from github.
#

INCLUDE(FetchContent)

FETCHCONTENT_DECLARE(
	yaml_cpp
	GIT_REPOSITORY	https://github.com/jbeder/yaml-cpp.git
)

FETCHCONTENT_GETPROPERTIES(yaml_cpp)
IF (NOT yaml_cpp_POPULATED)
	FETCHCONTENT_POPULATE(yaml_cpp)

	# Set yaml-cpp flags:
	OPTION(YAML_CPP_BUILD_TESTS FALSE)
	OPTION(YAML_CPP_BUILD_TOOLS FALSE)
	OPTION(YAML_CPP_INSTALL FALSE)

	ADD_SUBDIRECTORY(${yaml_cpp_SOURCE_DIR} ${yaml_cpp_BINARY_DIR})

	# Let other parts of this project know about the header dirs:
	SET(yaml_cpp_INCLUDE_DIRS "${yaml_cpp_SOURCE_DIR}/include" PARENT_SCOPE)
    SET(yaml_cpp_LIBRARIES "${yaml_cpp_BINARY_DIR}/libyaml-cpp.a" PARENT_SCOPE)
ENDIF()
