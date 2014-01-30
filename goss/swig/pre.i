/* -*- C -*- */
// Copyright (C) 2012 Johan Hake
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

%ignore goss::DoubleVector::data;
%ignore goss::DoubleVector2D::data;

%ignore goss::ODESystemSolver::states;

//-----------------------------------------------------------------------------
// Need to ignore these dues to SWIG confusion of overloaded functions
//-----------------------------------------------------------------------------
%ignore goss::Table::set(std::string,std::string,std::size_t);

//-----------------------------------------------------------------------------
// Ignore operators so SWIG stop complaining
//-----------------------------------------------------------------------------
%ignore goss::TableEntry::operator std::string;
%ignore goss::Progress::operator++;
%ignore goss::Progress::operator=;
%ignore goss::Table::operator=;
%ignore goss::TableEntry::operator=;

//-----------------------------------------------------------------------------
// Ignore GOSS C++ stream handling
//-----------------------------------------------------------------------------
%ignore goss::LogStream;
%ignore goss::cout;
%ignore goss::endl;
