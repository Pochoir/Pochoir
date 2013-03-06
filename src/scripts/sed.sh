#!/bin/sed

s/Pochoir_Shape<2>/Pochoir_Shape_2D/g
s/Pochoir_Shape<3>/Pochoir_Shape_3D/g
s/Pochoir_Shape <2>/Pochoir_Shape_2D/g
s/Pochoir_Shape <3>/Pochoir_Shape_3D/g
s/Pochoir<2>/Pochoir_2D/g
s/Pochoir<3>/Pochoir_3D/g
s/Pochoir <2>/Pochoir_2D/g
s/Pochoir <3>/Pochoir_3D/g
s/registerBV/Register_Boundary/g
s/registerArray/Register_Array/g
s/Pochoir_kernel_end/Pochoir_Kernel_End/g
s/Pochoir_kernel/Pochoir_Kernel/g
s/Pochoir_kernel/Pochoir_Kernel/g
s/Pochoir_kernel/Pochoir_Kernel/g
s/Pochoir_Boundary_end/Pochoir_Boundary_End/g
s/run_obase/Run_Obase/g
s/run/Run/g
s/registerShape/Register_Shape/g
s/registerDomain/Register_Domain/g

