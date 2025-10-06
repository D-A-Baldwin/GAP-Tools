#########
# Why do these functions work:
# TO ADD
#########


GlueGroupOfLattice := function(B_lat, IntMat)
  # This returns the "glue group" of a lattice, that is, the quotient 
  # of the dual of an integral lattice by the lattice.
  # B_lat is a matrix whose columns are the basis of the lattice and
  # IntMat is the matrix defining the inner product on the space in 
  # which the lattice lives. In the orthogonal basis this is just the 
  # Identity matrix. The output is a list of numbers [d1, d2,...,dn]
  # which corresponds to Z_d1 x Z_d2 x ... x Z_dn. A 1 corresponds 
  # to the trivial group and a 0 corresponds to a factor of Z.
  
	local GramMat, S;
	GramMat := TransposedMat(B_lat) *IntMat* B_lat;
	S := SmithNormalFormIntegerMatTransforms(GramMat).normal;;
	return DiagonalOfMatrix(S);
end;


QuotientGroupOfLattices := function(B_lat_super, B_lat_sub, IntMat)
  # This is a more general function than the above one. This returns the quotient
  # of a lattice by some sublattice.
  
	local GramMat_super, S, T, n, B_lat_sub_new, i;
	
	n := DimensionsMat(B_lat_super)[1];;
	
	if Length(TransposedMat(B_lat_sub)) < Length(TransposedMat(B_lat_super)) then
		B_lat_sub_new := TransposedMatMutable(B_lat_sub);;
		for i in [1..Length(TransposedMat(B_lat_super))-
						Length(TransposedMat(B_lat_sub))] do
			Add(B_lat_sub_new, NullMat(1,n)[1]);
		od;
		B_lat_sub := TransposedMat(B_lat_sub_new);
	fi;

	GramMat_super := TransposedMat(B_lat_super) *IntMat* B_lat_super;
	T := GramMat_super^-1 * TransposedMat(B_lat_super) * IntMat * B_lat_sub;
	S := SmithNormalFormIntegerMatTransforms(T).normal;;
	return DiagonalOfMatrix(S);
end;
