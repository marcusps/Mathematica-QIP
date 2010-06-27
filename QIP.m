BeginPackage["QIP`"]

ComplexMatrixPlot::usage =
  "ComplexMatrixPlot[m] plots the entries in a complex-valued matrix as an array
   of notched disks.The radius of the disk is proportional to the absolute value
   of the matrix element, and the radial notch corresponds to the 
   phase/angle/argument of the complex number. By default, the disks are scaled 
   so that the disks never overlap or touch. In other words, the entry with the
   largest absolute values sets the scaling of all entries in the plot."
  
RealMatrixPlot::usage =
  "RealMatrixPlot[m] plots the entries in a real-valued matrix as an array of 
   filles squares. The perimeter of the square is proportional to the absolute 
   value of the entry, and color indicates whether it is negative (Red) or 
   positive (Black). All squares are scaled so that no entries overlap or 
   touch.  In other words, the entry with the largest absolute value sets the
   scaling"
  
Begin["`Private`"]

End[]

EndPackage[]
