/**
 * @file hctexceptions.h
 *   Definitions for the classes that are
 *   thrown when %HeteroCt experiences an error condition
 *   (also contains errorhandling module text - see \ref errorhandling).
 */

// This file is part of HeteroCt. See License.txt in the top-level directory.

#ifndef HCT_CTEXCEPTIONS_H
#define HCT_CTEXCEPTIONS_H

#include "cantera/base/fmt.h"
#include "cantera/base/ctexceptions.h"
#include <exception>

namespace OpenMKM
{

/*!
 * @defgroup errorhandling Error Handling
 *
 * \brief These classes and related functions are used to handle errors and
 *        unknown events within Heteroct.
 *
 * The general idea is that exceptions arising within Cantera are thrown 
 * using the common base class called CanteraError. However errors arising
 * out of tube driver and tube yaml input parsing are thrown with
 * derivative classes of CanteraError.  A list of all of the thrown errors 
 * is kept in the Application class.
 *
 * Any exceptions which are not caught cause a fatal error exit from the
 * program.
 */


//! YAML parser error.
/*!
 * This error is thrown if a supplied tube parameter is invalid or if a
 * required parameter is not supplied.
 *
 * @ingroup errorhandling
 */
class YAMLParserError : public Cantera::CanteraError
{
public:
    //! Constructor
    /*!
     * The length needed is supplied by the argument, reqd, and the
     * length supplied is given by the argument sz.
     *
     * @param procedure String name for the function within which the error was
     *             generated.
     * @param node_name    This is the node where parser failed.
     * @param failure_mode This is the failure mode (Node not found, value null etc.)
     */
    YAMLParserError(const std::string& procedure, const std::string& node_name, 
                    const std::string& failure_mode) :
        Cantera::CanteraError(procedure), node_(node_name), failure_mode_(failure_mode) {}

    virtual std::string getMessage() const;
    virtual std::string getClass() const {
        return "YAMLParserError";
    }

private:
    std::string node_, failure_mode_;
};

}

#endif
