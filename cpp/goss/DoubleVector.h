// Copyright (C) 2006-2012 Ola Skavhaug
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
// Modified by Johan Hake 2012

#ifndef DOUBLEVECTOR_H_IS_INCLUDED
#define DOUBLEVECTOR_H_IS_INCLUDED

#include <boost/scoped_array.hpp>

#include "types.h"

namespace goss {

  // A small wrapper around a double pointer
  struct DoubleVector
  {

    boost::scoped_array<double> data;
    uint n;

  };

  // A small wrapper around a double pointer
  struct DoubleVector2D
  {

    boost::scoped_array<double> data;
    uint m;
    uint n;

  };

}

#endif
