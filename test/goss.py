__author__ = "Johan Hake (hake.dev@gmail.com)"
__date__ = "2012-09-20 -- 2012-10-10"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU LGPL Version 3.0 or later"

from gotran2 import *
from gotran2.codegeneration.gosscodegenerator import GossCodeGenerator
from gotran2.codegeneration.codegenerator import ODERepresentation

ode = load_ode("winslow")
oderepr = ODERepresentation(ode, keep_intermediates=True, use_cse=False)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

oderepr = ODERepresentation(ode, ode.name+"NoIntermediates", keep_intermediates=False)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

oderepr = ODERepresentation(ode, ode.name+"CSE", keep_intermediates=False, use_cse=True)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

oderepr = ODERepresentation(ode, ode.name+"CSEArray", keep_intermediates=False,
                            use_cse=True, use_state_names=False)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

print

ode = load_ode("panfilov")
oderepr = ODERepresentation(ode, keep_intermediates=True, use_cse=False)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

oderepr = ODERepresentation(ode, ode.name+"NoIntermediates", keep_intermediates=False)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

oderepr = ODERepresentation(ode, ode.name+"CSE", keep_intermediates=False, use_cse=True)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

oderepr = ODERepresentation(ode, ode.name+"CSEArray", keep_intermediates=False,
                            use_cse=True, use_state_names=False)
gossgen = GossCodeGenerator(oderepr)
open(gossgen.name+".h", "w").write(gossgen.generate())

