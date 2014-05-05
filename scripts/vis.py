from pymol.cgo import *
from pymol import cmd

cmd.do("load traj.xyz, traj")
cmd.do("hide all")
cmd.do("show spheres")
cmd.do("alter elem h, vdw=0.1")
cmd.do("rebuild")
#cmd.mset("1 -100")
#cmd.do("set cache_frames, 1 ")
#cmd.do("viewport 400, 280")
#cmd.do("set ray_trace_frames=0")
#cmd.do("pngseq traj")
#cmd.do("mpng traj")

