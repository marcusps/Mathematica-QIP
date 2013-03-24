(* ::Package:: *)

BeginPackage["QIP`"];


(*
    QIP is a Mathematica package of useful functions for quantum information processing
    Copyright (C) 2011  Marcus P. da Silva

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


Id::usage = "Id[n] return an identity matrix of dimension n.";


Ket::usage = "Ket[n,d] returns the unit column vector of index 0 <= n < d in a d dimensional Hilbert space.";


Bra::usage = "Bra[n,d] returns the unit row vector of index 0 <= n < d in a d dimensional Hilbert space.";


KetBra::usage="KetBra[i,j,d] returns a matrix which maps the Ket[j,d] to Ket[i,d].";


Projector::usage="Returns the projector correcponding to a given ket.";


TrA::usage = "TrA[rho,da] returns the partial trace of rho, by tracing over the first da dimensions.";


TrB::usage = "TrB[rho,db] returns the partial trace of rho, by tracing over the last db dimensions.";


CircleTimes::usage = "CircleTime is the binary KroneckerProduct between matrices or vectors.";


Pauli::usage =
  "Pauli[k] returns the Pauli matrix corresponding to index k.";

Tensor::usage =
  "Tensor[l] computes the tensor product of a list of matrices.";


\[Sigma]x::usage="Pauli matrix x, acting on 2 levels."
\[Sigma]y::usage="Pauli matrix y, acting on 2 levels."
\[Sigma]z::usage="Pauli matrix z, acting on 2 levels."


ComplexMatrixPlot::usage =
  "ComplexMatrixPlot[m] plots the entries in a complex-valued matrix as an array
   of notched disks. The radius of the disk is proportional to the absolute value
   of the matrix element, and the radial notch corresponds to the 
   phase/angle/argument of the complex number. By default, the disks are scaled 
   so that the disks never overlap or touch. In other words, the entry with the
   largest absolute values sets the scaling of all entries in the plot.";
  
RealMatrixPlot::usage =
  "RealMatrixPlot[m] plots the entries in a real-valued matrix as an array of 
   filles squares. The perimeter of the square is proportional to the absolute 
   value of the entry, and color indicates whether it is negative (Red) or 
   positive (Black). All squares are scaled so that no entries overlap or 
   touch.  In other words, the entry with the largest absolute value sets the
   scaling";
  
RealRepHermitianMatrix::usage = 
  "RealRepHermitianMatrix[m] computes a real-valued matrix representation of the
   Hermitian matrix m. This representation consists the real part of the upper
   triangular section, the imaginary part of the lower triangular section, and the
   diagonal entries unchanged. This is possible due symmetry of the entries of a
   Hermitian matrix.";


Vec::usage="Vec[m] takes a matrix m and returns a column vector corresponding to the stacked columns of the matrix.";
VecInv::usage="VecInv[m,r,c] takes a vectorized matrix and returns the corresponding matrix with r rows and c columns.";
HermitianPart::usage="HermitianPart[m] returns the Hermitian part of the matrix m.";
SkewHermitianPart::usage="SkewHermitianPart[m] returns the skew-Hermitian part of the matrix m.";
ChoiLiouvilleInvolution::usage="ChoiLiouvilleInvolution[m] rearranges the matrix elements of m according to the involution mapping between the Choi matrix of a process and the Liouville (or natural) representation of the same process assuming column major ordering.";
MaximallyEntangledDensityMatrix::usage="MaximallyEntangledDensityMatrix[d] returns the maximally entangled state of two d-dimensional subsystems, with perfect correlations between the basis states.";
HaarRandomUnitary::usage="HarrRandomUnitary[d] returns a Harr-distributed unitary acting on a Hilber space of dimension d.";
SquareMatrixQ::usage="SquareMatrixQ[m] returns True if m is a square matrix and False otherwise.";
Reshape::usage="Reshape[m,r,c] reshapes an arbitrarily nested list of elements into a matrix with r rows and c columns. It is not compatible with Vec and VecInv.";
Liou::usage=""
LiouvilleRepresentation::usage=""
Dissipator::usage=""
Hamiltonian::usage=""
Purity::usage="Computes the purity of a density matrix."
MultiPauli::usage="Return the matrix corresponding to the tensor product of the single qubit Pauli operators with indices given by a list of integers between 0 and 3"


Begin["`Private`"];


Purity[m_]:=Tr[MatrixPower[m,2]]


InsertAt[a_,j_,L_]:=MapAt[a&,Table[0,{i,L}],j]


CirclePlus[x__]:=ArrayFlatten[MapIndexed[InsertAt[#1,#2[[1]],Length[List[x]]]&,List[x]]]


BlockMatrix[m_]:=ArrayFlatten[m]


LiouvilleRepresentation[left_?MatrixQ,right_?MatrixQ]:=Transpose[right]\[CircleTimes]left;


Liou=LiouvilleRepresentation;


Hamiltonian[H_?SquareMatrixQ]:=-I Liou[H,Id[Dimensions[H][[1]]]]+I Liou[Id[Dimensions[H][[1]]],H]


Dissipator[m_?SquareMatrixQ]:=Liou[m,m\[ConjugateTranspose]]-1/2 Liou[m\[ConjugateTranspose].m,Id[Dimensions[m][[1]]]]-1/2 Liou[Id[Dimensions[m][[1]]],m\[ConjugateTranspose].m]


CircleTimes[x__]:=KroneckerProduct[x]


Id[d_]:=IdentityMatrix[d]


\[Sigma]x={
 {0, 1},
 {1, 0}
};
\[Sigma]y={
 {0, -I},
 {I, 0}
};
\[Sigma]z={
 {1, 0},
 {0, -1}
};


Pauli[x_] := Switch[x,0,Id[2],1,\[Sigma]x,2,\[Sigma]y,3,\[Sigma]z];


MultiPauli[l_]:=Fold[KroneckerProduct,{{1}},Map[Pauli,l]];


\[Sigma]=Pauli;


Tensor[l_] := Apply[KroneckerProduct,l];


Ket[n_,d_]:={UnitVector[d,n+1]}\[Transpose];


Bra[n_,d_]:=Ket[n,d]\[ConjugateTranspose];


KetBra[i_,j_,d_]:=Ket[i,d].Bra[j,d]


Projector[\[Psi]_]:=\[Psi].\[Psi]\[ConjugateTranspose];


Reshape[m_,_,l_]:=Partition[Flatten[m],l]


Vec[mtx_?MatrixQ]:=Partition[Flatten[mtx\[Transpose]],1]


VecInv[vec_,m_/;m>0&&m\[Element]Integers,n_/;n>0&&n\[Element]Integers]:=Partition[Flatten[vec],n]\[Transpose]


HermitianPart[m_?SquareMatrixQ]:=1/2 (m+m\[ConjugateTranspose])


SkewHermitianPart[m_?SquareMatrixQ]:=1/2 (m-m\[ConjugateTranspose])


TrA[m_?SquareMatrixQ,da_/;da>0&&da\[Element]Integers]:=Sum[(Bra[i,da]\[CircleTimes]Id[Dimensions[m][[2]]/da]).m.(Ket[i,da]\[CircleTimes]Id[Dimensions[m][[2]]/da]),{i,0,da-1}]


TrB[m_?SquareMatrixQ,db_/;db>0&&db\[Element]Integers]:=Sum[(Id[Dimensions[m][[2]]/db]\[CircleTimes]Bra[i,db]).m.(Id[Dimensions[m][[2]]/db]\[CircleTimes]Ket[i,db]),{i,0,db-1}]


ChoiLiouvilleInvolution[m_?SquareMatrixQ]:=
Block[{d=Sqrt[Dimensions[m][[1]]]},
  Sum[
    (KetBra[n,i,d]\[CircleTimes]KetBra[l,l,d])
    .m.
    (KetBra[j,j,d]\[CircleTimes]KetBra[n,i,d]),
    {i,0,d-1},
    {j,0,d-1},
    {l,0,d-1},
    {n,0,d-1}
  ]
]


MaximallyEntangledDensityMatrix[d_/;d>0&&d\[Element]Integers]:=
  1/d Projector[Sum[Ket[i,d]\[CircleTimes]Ket[i,d],{i,0,d-1}]]


HaarRandomUnitary[d_/;d>0&&d\[Element]Integers]:=
  Orthogonalize[
    RandomReal[NormalDistribution[0,1],{d,d}]+
    I RandomReal[NormalDistribution[0,1],{d,d}]
  ]


SquareMatrixQ = MatrixQ @ # && Equal @@ Dimensions @ # &;



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



End[];

EndPackage[];
