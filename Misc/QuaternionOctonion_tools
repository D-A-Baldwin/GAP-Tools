# quaternion and octonion multiplication and conjugation


qmult := function(x,y)
  # This takes two quaternions x and y which are lists of four numbers and outputs
  # quaternionic product.
	local a,b,c,d;
	
	a := x[1]*y[1] - x[2]*y[2] - x[3]*y[3] - x[4]*y[4];
	b := x[1]*y[2] + x[2]*y[1] + x[3]*y[4] - x[4]*y[3];
	c := x[1]*y[3] - x[2]*y[4] + x[3]*y[1] + x[4]*y[2];
	d := x[1]*y[4] + x[2]*y[3] - x[3]*y[2] + x[4]*y[1];
	
	return [a,b,c,d];
end;

qconj := function(v)
  # The conjugate of a quaternion v, written as a list of four numbers.
	local conj;
	conj := [1,-1,-1,-1];
	return List([1..4], i-> conj[i]*v[i]);
end;

oconj := function(v)
  # The conjugate of an octonion v, written as a list of eight numbers.
	local conj;
	conj := [1,-1,-1,-1,-1,-1,-1,-1];
	return List([1..8], i-> conj[i]*v[i]);
end;

omult := function(v,w)
  # This takes two quaternions v and w which are lists of eight numbers and outputs
  # octonionic product, defined via Dixon doubling by the quaternionic product.
	local a1,a2,b1,b2,vw1,vw2;
	
	a1 := v{[1..4]};
	a2 := v{[5..8]};
	b1 := w{[1..4]};
	b2 := w{[5..8]};
	
	vw1 := qmult(a1,b1) - qmult(qconj(b2), a2);
	vw2 := qmult(b2,a1) + qmult(a2,qconj(b1));
	Append(vw1,vw2);
	
	return vw1;
end;
