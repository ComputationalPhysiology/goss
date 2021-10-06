// Copyright (C) 2003-2013 Anders Logg
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
// Thanks to Jim Tilander for many helpful hints.
//
// Modified by Ola Skavhaug 2007, 2009
//
// First added:  2003-03-13
// Last changed: 2013-01-07

#ifndef __LOGGER_H
#define __LOGGER_H

#include <map>
#include <ostream>
#include <string>
#include <boost/scoped_ptr.hpp>
#include "Table.h"
#include "LogLevel.h"

// Forward declarations
namespace boost { class thread; }

namespace goss
{

  class Logger
  {
  public:

    /// Constructor
    Logger();

    /// Destructor
    ~Logger();

    /// Print message
    void log(std::string msg, int log_level=INFO) const;

    /// Print underlined message
    void log_underline(std::string msg, int log_level=INFO) const;

    /// Print warning
    void warning(std::string msg) const;

    /// Print error message and throw exception
    void error(std::string msg) const;

    /// Print error message, prefer this to the above generic error message
    void goss_error(std::string location,
                      std::string task,
                      std::string reason) const;

    /// Issue deprecation warning for removed feature
    void deprecation(std::string feature,
                     std::string version,
                     std::string message) const;

    /// Begin task (increase indentation level)
    void begin(std::string msg, int log_level=INFO);

    /// End task (decrease indentation level)
    void end();

    /// Draw progress bar
    void progress (std::string title, double p) const;

    /// Set output stream
    void set_output_stream(std::ostream& stream);

    /// Get output stream
    std::ostream& get_output_stream() { return *logstream; }

    /// Turn logging on or off
    void set_log_active(bool active);

    /// Return true iff logging is active
    inline bool is_active() { return _active; }

    /// Set log level
    void set_log_level(int log_level);

    /// Get log level
    inline int get_log_level() const { return _log_level; }

    /// Register timing (for later summary)
    void register_timing(std::string task, double elapsed_time);

    /// Return a summary of timings and tasks as a Table, optionally clearing
    /// stored timings
    Table timings(bool reset=false);

    /// Print summary of timings and tasks, optionally clearing stored timings
    void list_timings(bool reset=false);

    /// Return timing (average) for given task, optionally clearing timing for task
    double timing(std::string task, bool reset=false);

    /// Monitor memory usage. Call this function at the start of a
    /// program to continuously monitor the memory usage of the process.
    void monitor_memory_usage();

    /// Helper function for reporting memory usage
    void _report_memory_usage(size_t num_mb);

    /// Helper function for goss_debug macro
    void __debug(std::string msg) const;

    /// Helper function for goss_goss_assert macro
    void __goss_assert(std::string file, unsigned long line,
                  std::string function, std::string check) const;

  private:

    // Write message
    void write(int log_level, std::string msg) const;

    // True iff logging is active
    bool _active;

    // Current log level
    int _log_level;

    // Current indentation level
    int indentation_level;

    // Optional stream for logging
    std::ostream* logstream;

    // List of timings for tasks, map from string to (num_timings, total_time)
    std::map<std::string, std::pair<std::size_t, double> > _timings;

    // MPI data (initialized to 0)
    mutable std::size_t num_processes;
    mutable std::size_t process_number;

    // Thread used for monitoring memory usage
    boost::scoped_ptr<boost::thread> _thread_monitor_memory_usage;

    // Maximum memory usage so far
    long int _maximum_memory_usage;

  };

}

#endif
