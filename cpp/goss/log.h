// Copyright (C) 2003-2013 Anders Logg and Jim Tilander
//
// This file is part of GOSS.
//
// GOSS is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GOSS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GOSS. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Ola Skavhaug 2007, 2009
//
// First added:  2003-03-13
// Last changed: 2013-01-07

#ifndef __LOG_H
#define __LOG_H

#include <string>

#include "LogLevel.h"

namespace goss
{

  class Parameters;

  /// The GOSS log system provides the following set of functions for
  /// uniform handling of log messages, warnings and errors. In addition,
  /// macros are provided for debug messages and goss_assertions.
  ///
  /// Only messages with a debug level higher than or equal to the current
  /// log level are printed (the default being zero). Logging may also be
  /// turned off by calling set_log_active(false).

  /// Print message
  void info(std::string msg, ...);

  /// Print message to stream
  void info_stream(std::ostream& out, std::string msg);

  /// Print underlined message
  void info_underline(std::string msg, ...);

  /// Print warning
  void warning(std::string msg, ...);

  /// Print error message and throw an exception.
  /// Note to developers: this function should not be used internally
  /// in GOSS. Use the more informative goss_error instead.
  void error(std::string msg, ...);

  /// Print error message. Prefer this to the above generic error message.
  ///
  /// *Arguments*
  ///     location (std::string)
  ///         Name of the file from which the error message was generated.
  ///     task (std::string)
  ///         Name of the task that failed.
  ///         Note that this string should begin with lowercase.
  ///         Note that this string should not be punctuated.
  ///     reason (std::string)
  ///         A format string explaining the reason for the failure.
  ///         Note that this string should begin with uppercase.
  ///         Note that this string should not be punctuated.
  ///         Note that this string may contain printf style formatting.
  ///     ... (primitive types like int, uint, double, bool)
  ///         Optional arguments for the format string.
  ///
  /// Developers should read the file goss/log/README in the GOSS
  /// source tree for further notes about the use of this function.
  void goss_error(std::string location,
		  std::string task,
		  std::string reason, ...);

  /// Issue deprecation warning for removed feature
  ///
  /// *Arguments*
  ///     feature (std::string)
  ///        Name of the feature that has been removed.
  ///     version (std::string)
  ///        Version number of the release in which the feature was removed.
  ///     message (std::string)
  ///        A format string explaining the deprecation.
  void deprecation(std::string feature,
                   std::string version,
                   std::string message, ...);

  /// Print message at given debug level
  void log(int debug_level, std::string msg, ...);

  /// Begin task (increase indentation level)
  void begin(std::string msg, ...);

  /// Begin task (increase indentation level)
  void begin(int debug_level, std::string msg, ...);

  /// End task (decrease indentation level)
  void end();

  /// Turn logging on or off
  void set_log_active(bool active=true);

  /// Set log level
  void set_log_level(int level);

  /// Set output stream
  void set_output_stream(std::ostream& out);

  /// Get log level
  int get_log_level();

  /// Monitor memory usage. Call this function at the start of a
  /// program to continuously monitor the memory usage of the process.
  void monitor_memory_usage();

  // Helper function for goss_debug macro
  void __debug(std::string file,
               unsigned long line,
               std::string function,
               std::string format, ...);

  // Helper function for goss_goss_assert macro
  void __goss_assert(std::string file,
                unsigned long line,
                std::string function,
                std::string check);

  // Helper function to indent stings
  std::string indent(std::string block);

}

// The following three macros are the only "functions" in GOSS
// named goss_foo. Other functions can be placed inside the
// GOSS namespace and therefore don't require a prefix.

// Debug macros (with varying number of arguments)
#define goss_debug(msg)                  do { goss::__debug(__FILE__, __LINE__, __FUNCTION__, msg); } while (false)
#define goss_debug1(msg, a0)             do { goss::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0); } while (false)
#define goss_debug2(msg, a0, a1)         do { goss::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1); } while (false)
#define goss_debug3(msg, a0, a1, a2)     do { goss::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1, a2); } while (false)
#define goss_debug4(msg, a0, a1, a2, a3) do { goss::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1, a2, a3); } while (false)

// Not implemented error, reporting function name and line number
#define goss_not_implemented() \
  do { \
    goss::goss_error("log.h", \
                 "perform call to GOSS function", \
                 "The function %s has not been implemented (in %s line %d)", \
                 __FUNCTION__, __FILE__, __LINE__); \
  } while (false)

// Assertion, only active if DEBUG is defined
#ifdef DEBUG
#define goss_assert(check) \
  do { \
    if (!(check)) \
    { \
      goss::__goss_assert(__FILE__, __LINE__, __FUNCTION__, #check);    \
    } \
  } while (false)
#else
#define goss_assert(check)
#endif

#endif
