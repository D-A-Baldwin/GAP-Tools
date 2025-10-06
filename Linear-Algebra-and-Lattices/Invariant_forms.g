##########
# Functions to compute the space of symmetric forms invariant under the action
# of a matrix or matrix group.
##########

Read("Linear_tools.g");

InvariantFormsOfMatrix := function(M)
	local n, i, bas_vec, SymMatBasis, off_diag_inds, ind, mat, M_SymRep, row, ind2, j, transf_SymMats,
			nspace;
	
	n := Length(M);
	SymMatBasis := []; #generate standard basis of symmetric nxn matrices
	
	for i in [1..n] do
		bas_vec := NullMat(n,n);
		bas_vec[i][i] := 1;
		Add(SymMatBasis, bas_vec);
	od;
	
	off_diag_inds := Combinations([1..n], 2);
	
	for ind in off_diag_inds do
		bas_vec := NullMat(n,n);
		bas_vec[ind[1]][ind[2]] := 1;
		bas_vec[ind[2]][ind[1]] := 1;
		Add(SymMatBasis, bas_vec);
	od;
		
	transf_SymMats := [];
	for mat in SymMatBasis do
		Add(transf_SymMats, TransposedMat(M) * mat * M);
	od;
	
	M_SymRep := [];
	for mat in transf_SymMats do
		row := [];
		for j in [1..n] do
			Add(row, mat[j][j]);
		od;
		for ind2 in off_diag_inds do
			Add(row, mat[ind2[1]][ind2[2]]);
		od;
		Add(M_SymRep, row);
	od;
	
	nspace := NullspaceMat(M_SymRep - IdentityMat(Length(SymMatBasis)));
	
	return nspace;
end;


InvariantFormsOfMatrixGroup := function(matgrp)
	local matgens, InvFormsList, int, i;
	
	matgens := GeneratorsOfGroup(matgrp);
	
	InvFormsList := List(matgens, gen -> InvariantFormsOfMatrix(gen));
	
	int := InvFormsList[1];
	for i in [1..Length(matgens)] do
		int := IntersectionOfSubspaces(int, InvFormsList[i]);
	od;
	
	return int;
end;


VecToSymMat := function(vec, n)
	local SymMatBasis, i, j, bas_vec, off_diag_inds, ind;
	
	SymMatBasis := []; #generate standard basis of symmetric nxn matrices
	
	for i in [1..n] do
		bas_vec := NullMat(n,n);
		bas_vec[i][i] := 1;
		Add(SymMatBasis, bas_vec);
	od;
	
	off_diag_inds := Combinations([1..n], 2);
	
	for ind in off_diag_inds do
		bas_vec := NullMat(n,n);
		bas_vec[ind[1]][ind[2]] := 1;
		bas_vec[ind[2]][ind[1]] := 1;
		Add(SymMatBasis, bas_vec);
	od;
	
	for j in [1..Length(SymMatBasis)] do
		SymMatBasis[j] := vec[j]*SymMatBasis[j];
	od;
	
	return Sum(SymMatBasis);
end;
