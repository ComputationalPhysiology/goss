// Copyright (C) 2003-2012 Anders Logg
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
// Modified by Ola Skavhaug, 2007, 2009.
// Modified by Garth N. Wells, 2011.
//
// First added:  2003-03-13
// Last changed: 2011-11-15

#include <unistd.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <string>

#include <boost/thread.hpp>
#include <boost/bind/bind.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#ifdef __linux__
#include <sys/types.h>
#include <unistd.h>
#endif

#include "constants.h"
#include "LogLevel.h"
#include "Logger.h"
#include "log.h"

using namespace goss;

typedef std::map<std::string, std::pair<std::size_t, double> >::iterator map_iterator;
typedef std::map<std::string, std::pair<std::size_t, double> >::const_iterator const_map_iterator;

// Function for monitoring memory usage, called by thread
#ifdef __linux__
void _monitor_memory_usage(goss::Logger* logger)
{
  assert(logger);

  // Open statm
  //std::fstream

  // Get process ID and page size
  const std::size_t pid = getpid();
  const size_t page_size = getpagesize();

  // Print some info
  std::stringstream s;
  s << "Initializing memory monitor for process " << pid << ".";
  logger->log(s.str());

  // Prepare statm file
  std::stringstream filename;
  filename << "/proc/" << pid << "/statm";
  std::ifstream statm;

  // Enter loop
  while (true)
  {
    // Sleep for a while
    boost::this_thread::sleep(boost::posix_time::seconds(1));

    // Read number of pages from statm
    statm.open(filename.str().c_str());
    if (!statm)
      logger->error("Unable to open statm file for process.");
    size_t num_pages;
    statm >> num_pages;
    statm.close();

    // Convert to MB and report memory usage
    const size_t num_mb = num_pages*page_size / (1024*1024);
    logger->_report_memory_usage(num_mb);
  }
}
#endif

//-----------------------------------------------------------------------------
Logger::Logger() : _active(true), _log_level(INFO), indentation_level(0),
  logstream(&std::cout), num_processes(0), process_number(0),
  _maximum_memory_usage(-1)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Logger::~Logger()
{
  // Join memory monitor thread if it exists
  if (_thread_monitor_memory_usage)
    _thread_monitor_memory_usage->join();
}
//-----------------------------------------------------------------------------
void Logger::log(std::string msg, int log_level) const
{
  write(log_level, msg);
}
//-----------------------------------------------------------------------------
void Logger::log_underline(std::string msg, int log_level) const
{
  if (msg.empty())
    log(msg, log_level);

  std::stringstream s;
  s << msg;
  s << "\n";
  for (int i = 0; i < indentation_level; i++)
    s << "  ";
  for (std::size_t i = 0; i < msg.size(); i++)
    s << "-";

  log(s.str(), log_level);
}
//-----------------------------------------------------------------------------
void Logger::warning(std::string msg) const
{
  std::string s = std::string("*** Warning: ") + msg;
  write(WARNING, s);
}
//-----------------------------------------------------------------------------
void Logger::error(std::string msg) const
{
  std::string s = std::string("*** Error: ") + msg;
  throw std::runtime_error(s);
}
//-----------------------------------------------------------------------------
void Logger::goss_error(std::string location,
                          std::string task,
                          std::string reason) const
{
  std::stringstream s;
  s << std::endl << std::endl
    << "*** "
    << "-------------------------------------------------------------------------"
    << std::endl
    << "*** GOSS encountered an error. If you are not able to resolve this issue"
    << std::endl
    << "*** using the information listed below, you can ask for help at"
    << std::endl
    << "***" << std::endl
    << "***     hake@gmail.org"
    << std::endl
    << "***" << std::endl
    << "*** Remember to include the error message listed below and, if possible,"
    << std::endl
    << "*** include a *minimal* running example to reproduce the error."
    << std::endl
    << "***" << std::endl
    << "*** "
    << "-------------------------------------------------------------------------"
    << std::endl
    << "*** " << "Error:   Unable to " << task << "." << std::endl
    << "*** " << "Reason:  " << reason << "." << std::endl
    << "*** " << "Where:   This error was encountered inside " << location << "."
    << std::endl
    << "*** "
    << "-------------------------------------------------------------------------"
    << std::endl;

  throw std::runtime_error(s.str());
}
//-----------------------------------------------------------------------------
void Logger::deprecation(std::string feature,
                         std::string version,
                         std::string message) const
{
  std::stringstream s;
  s << "*** "
    << "-------------------------------------------------------------------------"
    << std::endl
    << "*** Warning: " << feature << " has been deprecated in GOSS version "
    << version << "." << std::endl
    << "*** " << message << std::endl
    << "*** "
    << "-------------------------------------------------------------------------"
    << std::endl;

  write(WARNING, s.str());
}
//-----------------------------------------------------------------------------
void Logger::begin(std::string msg, int log_level)
{
  // Write a message
  log(msg, log_level);
  indentation_level++;
}
//-----------------------------------------------------------------------------
void Logger::end()
{
  indentation_level--;
}
//-----------------------------------------------------------------------------
void Logger::progress(std::string title, double p) const
{
  std::stringstream line;
  line << title << " [";

  const int N = GOSS_TERM_WIDTH - title.size() - 12 - 2*indentation_level;
  const int n = static_cast<int>(p*static_cast<double>(N));

  for (int i = 0; i < n; i++)
    line << '=';
  if (n < N)
    line << '>';
  for (int i = n+1; i < N; i++)
    line << ' ';

  line << std::setiosflags(std::ios::fixed);
  line << std::setprecision(1);
  line << "] " << 100.0*p << '%';

  write(PROGRESS, line.str());
}
//-----------------------------------------------------------------------------
void Logger::set_output_stream(std::ostream& ostream)
{
  logstream = &ostream;
}
//-----------------------------------------------------------------------------
void Logger::set_log_active(bool active)
{
  _active = active;
}
//-----------------------------------------------------------------------------
void Logger::set_log_level(int log_level)
{
  _log_level = log_level;
}
//-----------------------------------------------------------------------------
void Logger::register_timing(std::string task, double elapsed_time)
{
  // Remove small or negative numbers
  if (elapsed_time < GOSS_EPS)
    elapsed_time = 0.0;

  // Print a message
  std::stringstream line;

  // Store values for summary
  map_iterator it = _timings.find(task);
  if (it == _timings.end())
  {
    std::pair<std::size_t, double> timing_(1, elapsed_time);
    _timings[task] = timing_;
  }
  else
  {
    it->second.first += 1;
    it->second.second += elapsed_time;
  }
}
//-----------------------------------------------------------------------------
void Logger::list_timings(bool reset)
{
  // Check if timings are empty
  if (_timings.empty())
  {
    log("Timings: no timings to report.");
    return;
  }
  else
  {
    log("");
    log(timings(reset).str(true));
  }

  // Print maximum memory usage if available
  if (_maximum_memory_usage >= 0)
  {
    std::stringstream s;
    s << "\nMaximum memory usage: " << _maximum_memory_usage << " MB";
    log(s.str());
  }

}
//-----------------------------------------------------------------------------
Table Logger::timings(bool reset)
{
  // Generate timing table
  Table table("Summary of timings");
  for (const_map_iterator it = _timings.begin(); it != _timings.end(); ++it)
  {
    const std::string task    = it->first;
    const std::size_t num_timings    = it->second.first;
    const double total_time   = it->second.second;
    const double average_time = total_time / static_cast<double>(num_timings);

    table(task, "Average time") = average_time;
    table(task, "Total time")   = total_time;
    table(task, "Reps")         = num_timings;

  }

  // Clear timings
  if (reset)
    _timings.clear();

  return table;
}
//-----------------------------------------------------------------------------
double Logger::timing(std::string task, bool reset)
{
  // Find timing
  map_iterator it = _timings.find(task);
  if (it == _timings.end())
  {
    std::stringstream line;
    line << "No timings registered for task \"" << task << "\".";
    goss_error("Logger.cpp",
                 "extract timing for task",
                 line.str());
  }

  // Compute average
  const std::size_t num_timings  = it->second.first;
  const double total_time   = it->second.second;
  const double average_time = total_time / static_cast<double>(num_timings);

  // Clear timing
  if (reset)
    _timings.erase(it);

  return average_time;
}
//-----------------------------------------------------------------------------
void Logger::monitor_memory_usage()
{
  #ifndef __linux__
  warning("Unable to initialize memory monitor; only available on GNU/Linux.");
  return;

  #else

  // Check that thread has not alrady been started
  if (_thread_monitor_memory_usage)
  {
    log("Memory monitor already initialize.");
    return;
  }

  // Create thread
  _thread_monitor_memory_usage.reset(new boost::thread(boost::bind(&_monitor_memory_usage, this)));

  #endif
}
//-----------------------------------------------------------------------------
void Logger::_report_memory_usage(size_t num_mb)
{
  std::stringstream s;
  s << "Memory usage: " << num_mb << " MB";
  log(s.str());
  _maximum_memory_usage = std::max(_maximum_memory_usage,
                                   static_cast<long int>(num_mb));
}
//-----------------------------------------------------------------------------
void Logger::__debug(std::string msg) const
{
  std::string s = std::string("Debug: ") + msg;
  write(DBG, s);
}
//-----------------------------------------------------------------------------
void Logger::__goss_assert(std::string file, unsigned long line,
                      std::string function, std::string check) const
{
  std::stringstream location;
  location << file << " (line " << line << ")";
  std::stringstream task;
  task << "complete call to function " << function << "()";
  std::stringstream reason;
  reason << "Assertion " << check << " failed";
  goss_error(location.str(), task.str(), reason.str());
}
//-----------------------------------------------------------------------------
void Logger::write(int log_level, std::string msg) const
{
  // Check log level
  if (!_active || log_level < _log_level)
    return;

  // Add indentation
  for (int i = 0; i < indentation_level; i++)
    msg = "  " + msg;

  // Write to stream
  *logstream << msg << std::endl;
}
//----------------------------------------------------------------------------
