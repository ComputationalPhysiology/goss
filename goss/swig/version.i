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

//-----------------------------------------------------------------------------
// Include code to generate a __swigversion__ and a __gossversion__ 
// attributes, from defines during compile time, to the cpp module
//-----------------------------------------------------------------------------

%inline %{
int goss_swigversion() { return SWIGVERSION; }
std::string goss_version() { return GOSS_VERSION; }
%}

%pythoncode %{
tmp = hex(goss_swigversion())
__swigversion__ = "%d.%d.%d"%(tuple(map(int, [tmp[-5], tmp[-3], tmp[-2:]])))
__gossversion__ = goss_version()
del tmp, goss_swigversion, goss_version
%}
