BeginPackage["QIP`"]

    (*
    QIP is a Mathematic package of useful functions for quantum information processing
    Copyright (C) 2010  Marcus P. da Silva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    *)

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
  
RealRepHermitianMatrix::usage = 
  "RealRepHermitianMatrix[m] computes a real-valued matrix representation of the
   Hermitian matrix m. This representation consists the real part of the upper
   triangular section, the imaginary part of the lower triangular section, and the
   diagonal entries unchanged. This is possible due symmetry of the entries of a
   Hermitian matrix."

Begin["`Private`"]

Options[ComplexToDisk] = {DiskColor -> Black, LineColor -> White, 
   MaxScale -> .95};

ComplexToDisk[scale_, c_, index_, OptionsPattern[]] :=
 Module[{r, i, idx},
  r = Re[c];
  i = Im[c];
  idx = RotateLeft[index];
  {
   Opacity[1],
   OptionValue[DiskColor],
   Disk[idx, .5 OptionValue[MaxScale] scale Abs[c]],
   Thickness[.003], OptionValue[LineColor],
   Line[{idx, idx + .5 OptionValue[MaxScale] scale {r, i}}],
   Opacity[0],
   Disk[idx, .5 OptionValue[MaxScale]]
   }
  ]

 (*
  The conventions for matrix indices and graphical coordinates are
  not consistent, however. For example, index [0,0] in a matrix
  is actually in the top left corner of a matrix, while coordinates
  (0,0) is at the bottom right of a plot. To compensate for
  that in the plots, we must reverse the order of the row.
  *)
 
Flip[m_] := Reverse[m]
 
Options[ComplexMatrixPlot] = {DiskColor -> Black, LineColor -> White, 
   MaxScale -> .95};

ComplexMatrixPlot[rho_, opts : OptionsPattern[]] :=
 Module[{fm, scale},
  fm = Flip[rho];
  scale = 1/Max[Abs[rho]];
  Graphics[
   MapIndexed[
    ComplexToDisk[scale, #1, #2, 
      FilterRules[{opts}, Options[ComplexMatrixPlot]]] &, fm, {2}], 
   FilterRules[{opts}, Options[Graphics]
    ]
   ]
  ]

 (*
  We can plot real valued matrices by representing positive valued by \
  black squares, and negative values by red squares. The sizes of the \
  squares are scaled such that the largest absolute value has size \
  one, and all other sizes are scaled accordingly.
  *)

CenteredRectangle[s_, index_] :=
 Rectangle[index - {s/2, s/2}, index + {s/2, s/2}]

Options[RealToSquare] = {PosColor -> Black, NegColor -> Red, MaxScale -> .95};

RealToSquare[scale_, r_, index_, OptionsPattern[]] :=
 Module[{idx},
  idx = RotateLeft[index];
  {
   If[r < 0, OptionValue[NegColor], OptionValue[PosColor]],
   CenteredRectangle[OptionValue[MaxScale] scale r, idx ],
   Opacity[0],
   CenteredRectangle[OptionValue[MaxScale], idx ]
   }
  ]
 
Options[RealMatrixPlot] = {PosColor -> Black, NegColor -> Red, 
   MaxScale -> .95};

RealMatrixPlot[rho_, opts : OptionsPattern[]] :=
 Module[{fm, scale, xsize, ysize},
  fm = Flip[rho];
  scale = 1/Max[Abs[rho]];
  Graphics[
   GraphicsGroup[
    {
     MapIndexed[
      RealToSquare[scale, #1, #2, 
        FilterRules[{opts}, Options[RealToSquare]]] &, fm, {2}]
     }
    ]
   , FilterRules[{opts}, Options[Graphics]]]
  ]

UpperTriangularPart[rho_] := 
 MapIndexed[If[#2[[1]] < #2[[2]], #1, 0] &, rho, {2}]

LowerTriangularPart[rho_] := 
 MapIndexed[If[#2[[1]] > #2[[2]], #1, 0] &, rho, {2}]

RealRepHermitianMatrix[rho_] :=
  Re[UpperTriangularPart[rho]] + Im[LowerTriangularPart[rho]] + 
  Re[DiagonalMatrix[Diagonal[rho]]]

End[]

EndPackage[]
