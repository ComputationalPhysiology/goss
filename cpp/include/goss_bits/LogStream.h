// Copyright (C) 2003-2009 Anders Logg
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
// Modified by Garth N. Wells 2005.
//
// First added:  2003-03-13
// Last changed: 2009-08-11

#ifndef __LOG_STREAM_H
#define __LOG_STREAM_H

#include <string>
#include <sstream>

namespace goss
{

  /// This class provides functionality similar to standard C++
  /// streams (std::cout, std::endl) for output but working through
  /// the GOSS log system.

  class LogStream
  {
  public:

    /// Stream types
    enum Type {COUT, ENDL};

    /// Create log stream of given type
    LogStream(Type type);

    /// Destructor
    ~LogStream();

    /// Output for log stream
    LogStream& operator<< (const LogStream& stream);

    /// Output for string
    LogStream& operator<< (const std::string& s);

    /// Output for int
    LogStream& operator<< (int a);

    /// Output for unsigned int
    LogStream& operator<< (unsigned int a);

    /// Output for long int
    LogStream& operator<< (long int a);

    /// Output for long int
    LogStream& operator<< (long unsigned int a);

    /// Output for double
    LogStream& operator<< (double a);

    void setprecision(std::streamsize n);

  private:

    // Type of stream
    Type _type;

    // Buffer
    std::stringstream buffer;

  };

  /// goss::cout
  extern LogStream cout;

  /// goss::endl;
  extern LogStream endl;

}

#endif
