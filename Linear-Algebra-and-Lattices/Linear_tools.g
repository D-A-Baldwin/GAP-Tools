###########
# A variety a small-ish functions related to matrices and linear algebra. Including several functions
# relating to lattices.
###########

ReflectionMatrix := function(vec, gram_mat)
  # This returns the matrix corresponding to the reflection through vec
  # with respect to the inner product given by gram_mat.
  # GAP actually already has a function that does this.

	local n;
	
	n := Length(gram_mat);
	
	return IdentityMat(n) -  List(vec, i-> i*(2/(vec*gram_mat*vec))*vec*gram_mat)  ;
end;


TransformationMatrixByFunction := function(func, n, m)
	# Given an R-linear function, func, that takes in a length-n vector
	# as an input and outputs a length-m vector, this function
	# outputs the corresponding transformation mxn matrix that maps R^n -> R^m
	# using the standard orthonormal basis for R^n.
	
	local M, bas_vec;
	
	M := [];
	
	for bas_vec in IdentityMat(n) do
		Add(M, func(bas_vec));
	od;
	
	return TransposedMat(M);
end;

IntersectionOfSubspaces := function(M1, M2)
	#Algorithm 2.3.9 in Henri Cohen's Computational Number theory book 
	#M1 and M2 are matrices whose rows each span a subspace 
	#THEY ARE ASSUMED TO EACH HAVE LINEARLY INDEPENDENT rows

	local n, M1copy, nullmat, nullmatreduced;
	if M1 = [] or M2 = [] then
		return [];
	else 
		n :=  DimensionsMat(M1)[1];
	
		M1copy := ShallowCopy(M1);
		Append(M1copy, M2);
		#M1copy is now the transpose of what is called M1 in Cohen, since GAP likes right multiplication
		
		nullmat := NullspaceMat(M1copy);
		
		if nullmat = [] then
			return [];
		else
		
			nullmatreduced := nullmat{[1..Length(nullmat)]}{[1..n]};#find a basis of the kernel, take first n columns
			#note that GAP computes the kernel using right multiplication
			return nullmatreduced*M1;
		fi;
	fi;
end;

CommonFixedSpaceOfMatrixGroup := function(matgrp)
	local matgens, FixedSpaces, int, i;
	
	matgens := GeneratorsOfGroup(matgrp);
	
	FixedSpaces := List(matgens, gen -> NullspaceMat(IdentityMat(Length(gen)) - TransposedMat(gen)));
	
	int := FixedSpaces[1];
	for i in [1..Length(matgens)] do
		int := IntersectionOfSubspaces(int, FixedSpaces[i]);
	od;
	
	return int;
end;

NFormRep := function(mat, N)
	# Given a square matrix acting on R^n, this function
	# returns a matrix which gives the action on the N-th 
	# exterior product of R^n, i.e. the action on N-forms
	
	local mat_nform, indexlist, mat_nform_row, minor, inds1, inds2;
	mat_nform := [];
	indexlist := Combinations([1..Length(mat)], N);
	
	for inds1 in indexlist do
		mat_nform_row := [];
		for inds2 in indexlist do
			minor := List(inds1, i -> mat[i]);
			minor := TransposedMat(List(inds2, i -> TransposedMat(minor)[i]));
			Add(mat_nform_row, Determinant(minor));
		od;
		Add(mat_nform, mat_nform_row);
	od;
	return TransposedMat(mat_nform);
end;


Gram_Schmidt := function(basis, orthonormal_ind)
	# The rows of basis span some vector space, this function returns an 
	# orthogonal set of vectors spanning the same space. If orthonormal_ind =
	# true then the function will normalise the output vectors so that they have
	# norm 1. Otherwise it won't.
	
	local n, basisnew, i, bold, mu_sum, j, bnew;
	
	n := Length(basis);
	
	if orthonormal_ind = true then
		basisnew := [basis[1] * (1/Sqrt(basis[1]*basis[1]))];
	else
		basisnew := [basis[1]];
	fi;
	
	for i in [2..n] do
		bold := basis[i];
		
		mu_sum := [];
		for j in [1..i-1] do
			Add(mu_sum, (bold*basisnew[j]/(basisnew[j]*basisnew[j])) * basisnew[j]);
		od;
		
		bnew := bold - Sum(mu_sum);
		
		if orthonormal_ind = true then
			bnew := bnew * (1/Sqrt(bnew*bnew));
		
			Add(basisnew, bnew);
		else	
			Add(basisnew, bnew);
		fi;
	od;
	
	return basisnew;
end;


MatrixDirectSum := function(mat1, mat2)
	# The matrices must be square. This gives the block diagonal sum
	# of mat1 and mat 2.
	# GAP's own BlockMatrix function requires the blocks to have 
	# the same shape.
	
	local n1,n2,row,rnew,newmat;
	
	if ( not DimensionsMat(mat1)[1] = DimensionsMat(mat1)[2]) or 
			( not DimensionsMat(mat1)[1] = DimensionsMat(mat1)[2]) then
		return "One or both of the matrices are not square";
	else
		n1 := Length(mat1);
		n2 := Length(mat2);
	
		newmat := [];
		for row in mat1 do
			rnew := ShallowCopy(row);
			Append(rnew, NullMat(1,n2)[1]);
			Add(newmat, rnew);
		od;
		
		for row in mat2 do
			rnew := NullMat(1,n1)[1];
			Append(rnew, row);
			Add(newmat, rnew);
		od;
	
		return newmat;
	fi;
end;


LatticeDirectSum := function(genmat_1, autgens_1, genmat_2, autgens_2)
	# genmat_1 and genmat_2 are the generator matrices of the lattices 
	# (i.e. the matrix whose columns are the basis of the lattice)
	# and autgens_1 and autgens_2 are lists of the generators of the automorphism
	# groups of each lattices in the orthogonal basis.
	# The output is the generator matrix of the direct sum of the lattices and
	# the list of generators of its automorphism group.
	
	local n1,n2,bas_1,bas_2,vec,newbas,t,newautgens, gen;
	
	n1 := DimensionsMat(genmat_1)[1];
	n2 := DimensionsMat(genmat_2)[1];
	
	bas_1 := TransposedMat(genmat_1);
	newbas := [];
	for vec in bas_1 do
		t := ShallowCopy(vec);
		Append(t, NullMat(1,n2)[1]);
		Add(newbas, t);
	od;
	
	bas_2 := TransposedMat(genmat_2);
	for vec in bas_2 do
		t := NullMat(1,n1)[1];
		Append(t, vec);
		Add(newbas,t);
	od;
	
	newautgens := [];
	for gen in autgens_1 do
		Add(newautgens, MatrixDirectSum(gen, IdentityMat(n2)));
	od;
	
	for gen in autgens_2 do
		Add(newautgens, MatrixDirectSum(IdentityMat(n1), gen));
	od;
	
	return [TransposedMat(newbas), newautgens];
end;
