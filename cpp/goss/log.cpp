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
// Modified by Ola Skavhaug 2007,
// Modified by Garth N. Wells 2009
//
// First added:  2003-03-13
// Last changed: 2013-01-07

#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <vector>

#include "constants.h"
#include "LogManager.h"
#include "log.h"


using namespace goss;

static std::vector<char> goss_buffer(GOSS_LINELENGTH);
static unsigned int goss_buffer_size = GOSS_LINELENGTH;

// Goss_Buffer allocation
void allocate_goss_buffer(std::string msg)
{
  // va_list, start, end require a char pointer of fixed size so we
  // need to allocate the buffer here. We allocate twice the size of
  // the format string and at least GOSS_LINELENGTH. This should be
  // ok in most cases.
  unsigned int new_size = std::max(static_cast<unsigned int>(2*msg.size()),
                                   static_cast<unsigned int>(GOSS_LINELENGTH));
  //static_cast<unsigned int>(GOSS_LINELENGTH));
  if (new_size > goss_buffer.size())
  {
    goss_buffer.resize(new_size);
    goss_buffer_size = new_size;
  }
}

// Macro for parsing arguments
#define read(local_buffer, msg) \
  allocate_goss_buffer(msg); \
  va_list aptr; \
  va_start(aptr, msg); \
  vsnprintf(local_buffer, goss_buffer_size, msg.c_str(), aptr); \
  va_end(aptr);

//-----------------------------------------------------------------------------
void goss::info(std::string msg, ...)
{
  if (!LogManager::logger.is_active() || INFO < LogManager::logger.get_log_level())
    return; // optimization
  read(goss_buffer.data(), msg);
  LogManager::logger.log(goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::info_stream(std::ostream& out, std::string msg)
{
  if (!LogManager::logger.is_active())
    return; // optimization
  std::ostream& old_out = LogManager::logger.get_output_stream();
  LogManager::logger.set_output_stream(out);
  LogManager::logger.log(msg);
  LogManager::logger.set_output_stream(old_out);
}
//-----------------------------------------------------------------------------
void goss::info_underline(std::string msg, ...)
{
  if (!LogManager::logger.is_active() || INFO < LogManager::logger.get_log_level())
    return; // optimization
  read(goss_buffer.data(), msg);
  LogManager::logger.log_underline(goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::warning(std::string msg, ...)
{
  if (!LogManager::logger.is_active() || WARNING < LogManager::logger.get_log_level())
    return; // optimization
  read(goss_buffer.data(), msg);
  LogManager::logger.warning(goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::error(std::string msg, ...)
{
  read(goss_buffer.data(), msg);
  LogManager::logger.error(goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::goss_error(std::string location,
		      std::string task,
		      std::string reason, ...)
{
  read(goss_buffer.data(), reason);
  LogManager::logger.goss_error(location, task, goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::deprecation(std::string feature,
                         std::string version,
                         std::string message, ...)
{
  read(goss_buffer.data(), message);
  LogManager::logger.deprecation(feature, version, goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::log(int log_level, std::string msg, ...)
{
  if (!LogManager::logger.is_active() || log_level < LogManager::logger.get_log_level())
    return; // optimization
  read(goss_buffer.data(), msg);
  LogManager::logger.log(goss_buffer.data(), log_level);
}
//-----------------------------------------------------------------------------
void goss::begin(std::string msg, ...)
{
  if (!LogManager::logger.is_active())
    return; // optimization
  read(goss_buffer.data(), msg);
  LogManager::logger.begin(goss_buffer.data());
}
//-----------------------------------------------------------------------------
void goss::begin(int log_level, std::string msg, ...)
{
  if (!LogManager::logger.is_active())
    return; // optimization
  read(goss_buffer.data(), msg);
  LogManager::logger.begin(goss_buffer.data(), log_level);
}
//-----------------------------------------------------------------------------
void goss::end()
{
  if (!LogManager::logger.is_active())
    return; // optimization
  LogManager::logger.end();
}
//-----------------------------------------------------------------------------
void goss::set_log_active(bool active)
{
  LogManager::logger.set_log_active(active);
}
//-----------------------------------------------------------------------------
void goss::set_log_level(int level)
{
  LogManager::logger.set_log_level(level);
}
//-----------------------------------------------------------------------------
void goss::set_output_stream(std::ostream& out)
{
  LogManager::logger.set_output_stream(out);
}
//-----------------------------------------------------------------------------
int goss::get_log_level()
{
  return LogManager::logger.get_log_level();
}
//-----------------------------------------------------------------------------
void goss::monitor_memory_usage()
{
  LogManager::logger.monitor_memory_usage();
}
//-----------------------------------------------------------------------------
void goss::__debug(std::string file, unsigned long line,
		   std::string function, std::string format, ...)
{
  read(goss_buffer.data(), format);
  std::ostringstream ost;
  ost << file << ":" << line << " in " << function << "()";
  std::string msg = std::string(goss_buffer.data()) + " [at " + ost.str() + "]";
  LogManager::logger.__debug(msg);
}
//-----------------------------------------------------------------------------
void goss::__goss_assert(std::string file, unsigned long line,
                      std::string function, std::string check)
{
  LogManager::logger.__goss_assert(file, line, function, check);
}
//-----------------------------------------------------------------------------
std::string goss::indent(std::string block)
{
  std::string indentation("  ");
  std::stringstream s;

  s << indentation;
  for (std::size_t i = 0; i < block.size(); ++i)
  {
    s << block[i];
    if (block[i] == '\n' && i < block.size() - 1)
      s << indentation;
  }

  return s.str();
}
//-----------------------------------------------------------------------------
