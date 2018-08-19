//------------------------------------------------------------------------------
// Options
//------------------------------------------------------------------------------

// Set path
path:="suhome/Magma/db/mixed_degree_one/";

//------------------------------------------------------------------------------
// Useful functions
//------------------------------------------------------------------------------

// Returns an integer sequence as a string
function IntegerSequenceToString(s)
	string:="[";
	for i in [1..#s] do
		string:=string cat IntegerToString(s[i]);
		if i lt #s then
			string:=string cat ",";
		else
			string:=string cat "]";
		end if;
	end for;
	return string;
end function;

// Returns a sequence of vertices as a string
function VerticesToString(s)
	string:="[";
	for i in [1..#s] do
		string:=string cat IntegerSequenceToString(s[i]);
		if i lt #s then
			string:=string cat ",";
		else
			string:=string cat "]";
		end if;
	end for;
	return string;
end function;

// Embed P in a dim(P)+1 dimensional space, at height 1
function homogenize(P)
	return Polytope([[1] cat Eltseq(v) : v in Vertices(P)]);
end function;

// Coordinate projection of P on the last dim(P)-1 entries
function dehomogenize(P)
	s:=[Eltseq(v) : v in Vertices(P)];
	return Polytope([Remove(v,1) : v in s]);
end function;

// Check if P and Q are the same polytope
function are_same_polytope(P,Q)
	return {Eltseq(v) : v in Vertices(P)} eq {Eltseq(v) : v in Vertices(Q)};
end function;

// Check if P and Q differs only by a translation
function are_translations(P,Q)
	v:=Representative(Vertices(P));
	for u in Vertices(Q) do
			if are_same_polytope(P-v,Q-u) then
				return true;
			end if;
	end for;
	return false;
end function;

// Check if two vectors have the same direction
function are_same_direction(u,v)
	p_u:=u/GCD(Eltseq(u));
	p_v:=v/GCD(Eltseq(v));
	return p_u eq p_v or p_u eq -p_v;
end function;

//------------------------------------------------------------------------------
// Important functions
//------------------------------------------------------------------------------

// Returns all the lattice preserving isomorphisms mapping 
// A x {1} to C x {1} (one needs to dehomogenize later).
function all_maps(A,C)
	Ap:=homogenize(A);
	Cp:=homogenize(C);
	bool,phi:=IsEquivalent(Ap,Cp);
	if not bool then
		return false,_,_;
	end if;
	L:=Ambient(Cp);
	return true,phi,[LatticeMap(L,[L ! p : p in Rows(m)]) : m in AutomorphismGroup(Cp)];
end function;

// Given a vector p in R^3, returns a map R^3 -> R^2 projecting along p
function project_along_direction(p)
	M:=Ambient(Polytope([p]));
	l:=LinearSpanEquations([p]);
	assert #l eq #Eltseq(p)-1;
	m:=Matrix(l);
	L,pi:=ToricLattice([Eltseq(c) : c in Rows(Transpose(m))]);
	Z2:=ToricLattice(2);
	phi:=LatticeMap(M,Z2,Matrix(([Z2 ! Eltseq(p) : p in Rows(Transpose(m))])));
	B:=Basis(M);	
	phi_basis:=Image(phi,B);
	phi_basis:=[Codomain(pi) ! a : a in phi_basis];
	ipi_phi_basis:=[Eltseq(p) : p in Preimage(pi,phi_basis)];
	projection:=LatticeMap(M,L,Matrix(ipi_phi_basis));
	return projection;
end function;

// Returns the Cayley sum of a family of polytopes
function cayley(list)
	n:=#list;
	assert n gt 1;
	if n eq 2 then
		return Polytope([Eltseq(v) cat [0] : v in Vertices(list[1])] cat [Eltseq(v) cat [1] : v in Vertices(list[2])] );
	else
		return cayley([cayley(list[1..n-1]),Polytope([Eltseq(v) cat [0 : i in [1..n-2]] : v in Vertices(list[n])])]);
	end if;
end function;

//------------------------------------------------------------------------------
// Extract all the hollow three dimensional polytopes with width >= 1 contained
// in one of the 12 hollow maximal polytopes
//------------------------------------------------------------------------------
D:=3;

A:=[
	[[0,0,0],[1,0,0],[1,2,0],[2,2,0],[1,0,2],[2,0,2],[2,2,2],[3,2,2]],
	[[0,0,0],[2,0,0],[2,4,0],[2,0,4]],
	[[0,0,0],[1,0,0],[2,4,0],[3,0,4]],
	[[0,0,0],[1,0,0],[2,3,0],[1,0,3],[-1,-6,3]],
	[[0,0,0],[3,0,0],[0,3,0],[0,0,3]],
	[[0,0,0],[1,0,0],[1,4,0],[3,0,4],[-1,4,-4]],
	[[0,0,0],[1,0,0],[4,6,0],[4,0,6]],
	[[0,0,0],[1,0,0],[2,3,0],[5,3,9]],
	[[0,0,0],[1,0,0],[3,4,0],[7,4,8]],
	[[0,0,0],[1,0,0],[2,3,0],[1,0,3],[2,0,3],[3,3,3]],
	[[0,0,0],[1,0,0],[2,5,0],[3,0,5]],
	[[0,0,0],[1,0,0],[1,2,0],[3,2,4],[2,2,0],[4,2,4]]
];

N:=Max([NumberOfPoints(Polytope(s)) : s in A]);

list:=[[] : i in [1..N]];

for s in A do
	P:=Polytope(s);
	n:=NumberOfPoints(P);
	Append(~list[n],Polytope([Eltseq(p) : p in AffineNormalForm(P)]));
end for;

for i in [N..1 by -1] do
	if #list[i] eq 0 then break; end if;
	for P in list[i] do
		for v in Vertices(P) do
			Q:=Polytope(Exclude(Points(P),v));
			if Dimension(Q) eq D and Width(Q) gt 1 then
				QQ:=Polytope([Eltseq(p) : p in AffineNormalForm(Q)]);
				Append(~list[NumberOfPoints(QQ)],QQ);
			end if;
		end for;
	end for;
	set_p:=SequenceToSet(list[i-1]);
	list[i-1]:=SetToSequence(set_p);
	printf "%o - %o\n",i,#list[i];
end for;

set:=SequenceToSet(&cat(list));
#set;

//------------------------------------------------------------------------------
// ... or load them from the file
//------------------------------------------------------------------------------
set:={};
file:=path cat "hollow_width_gt1.txt";
fh:=Open(file,"r");
while true do
	sl:=Gets(fh);
	if IsEof(sl) or #sl lt 2 then
		break;
	end if;
	while sl[#sl] eq "\\" do
		sl:=Prune(sl) cat Gets(fh);
	end while;
	s:=eval(sl);
	Include(~set,Polytope(s));
end while;

//------------------------------------------------------------------------------
// Pick the ones having a nontrivial Minkowski decomposition
//------------------------------------------------------------------------------

list_nd:=[];
i:=0;
for P in set do
	i+:=1;
	i;
	MDs:=[MD : MD in MinkowskiDecomposition(P) | #MD gt 1];
	if #MDs gt 0 then
		"*";
		Append(~list_nd,P);
	end if;
end for;

//------------------------------------------------------------------------------
// Check the different types (= dimensions of Minkowski summands). After
// that, we select the ones that can lead to a couple of three dimensional
// polytopes
//------------------------------------------------------------------------------

possible_types:=&cat[[[Dimension(S) : S in MD] : MD in MinkowskiDecomposition(P)] : P in list_nd];
possible_types:={Sort(a) : a in possible_types};

// possible_types:
//    [ 3, 3, 3 ],
//    [ 2, 2, 3 ],
//    [ 2, 2, 2 ],
//    [ 2, 3 ],
//    [ 2, 3, 3 ],
//    [ 1, 1, 1 ],
//    [ 2, 2 ],
//    [ 1, 2, 3 ],
//    [ 1, 2, 2 ],
//    [ 1, 3 ],
//    [ 1, 2 ],
//    [ 3, 3 ],
//    [ 1, 3, 3 ]

good_types:=[
	[ 3, 3, 3 ],
	[ 2, 2, 3 ],
	[ 2, 3, 3 ],
	[ 1, 2, 3 ],
	[ 3, 3 ],
	[ 1, 3, 3 ]
];

//------------------------------------------------------------------------------
// Keep the ones with a "good" Minkowski Decomposition
//------------------------------------------------------------------------------

good_MDs:=[];
for P in list_nd do
	MDs:=MinkowskiDecomposition(P);
	for MD in MDs do
		type:=Sort([Dimension(S) : S in MD]);
		if type in good_types then
			Append(~good_MDs,MD);
		end if;
	end for;
end for;

//------------------------------------------------------------------------------
// Create all the possible candidates for couples [P1,P2]
//------------------------------------------------------------------------------

couples:={};
for P in list_nd do
	MDs:=MinkowskiDecomposition(P);
	for MD in MDs do
		type:=Sort([Dimension(S) : S in MD]);
		if type eq [3,3,3] then // there is just one! ;)
			Include(~couples,[MD[1]+MD[2],MD[3]]);
		end if;
		if type eq [2,2,3] then
			P1:=&+[S : S in MD | Dimension(S) eq 2];
			P2:=[S : S in MD | Dimension(S) eq 3][1];
			if Dimension(P1) ge 3 then
				Include(~couples,[P1,P2]);
			end if;
		end if;
		if type eq [2,3,3] then
			P1:=[S : S in MD | Dimension(S) eq 2][1];
			P2:=[S : S in MD | Dimension(S) eq 3][1];
			P3:=[S : S in MD | Dimension(S) eq 3][2];
			Include(~couples,[P1+P2,P3]);
			Include(~couples,[P1+P3,P2]);
		end if;
		if type eq [1,2,3] then
			P1:=[S : S in MD | Dimension(S) eq 1][1];
			P2:=[S : S in MD | Dimension(S) eq 2][1];
			P3:=[S : S in MD | Dimension(S) eq 3][1];
			if Dimension(P1+P2) eq 3 then
				Include(~couples,[P1+P2,P3]);
			end if;
		end if;
		if type eq [3,3] then
			Include(~couples,[MD[1],MD[2]]);
		end if;
		if type eq [1,3,3] then
			P1:=[S : S in MD | Dimension(S) eq 1][1];
			P2:=[S : S in MD | Dimension(S) eq 3][1];
			P3:=[S : S in MD | Dimension(S) eq 3][2];
			Include(~couples,[P1+P2,P3]);
			Include(~couples,[P1+P3,P2]);
		end if;
	end for;
end for;

//------------------------------------------------------------------------------
// Filter the couples
//------------------------------------------------------------------------------

exceptional_couples:=couples;
for c in couples do
	P1:=c[1];
	P2:=c[2];
	M:=Ambient(P1);
	for e in Edges(P1) do
		endpoints:=[v : v in Vertices(e)];
		p:=endpoints[1]-endpoints[2];
		proj:=project_along_direction(p);
		p1:=Image(proj,P1);
		if Volume(p1) eq 1 then
			p2:=Image(proj,P2);
			if are_translations(p1,p2) then
				Exclude(~exceptional_couples,c);
			end if;
		end if;
	end for;
end for;

cayleys:={};
exceptional_couple_bis:={};
i:=0;
for c in exceptional_couples do
	i+:=1;
	i;
	C:=cayley(c);
	anfC:=AffineNormalForm(C);
	if not anfC in cayleys then
		Include(~cayleys,AffineNormalForm(C));
		Include(~exceptional_couple_bis,c);
		"*";
	end if;
end for;

couples:=exceptional_couple_bis;

//------------------------------------------------------------------------------
// ... or load them from a file
//------------------------------------------------------------------------------
couples:={};
file:=path cat "couples.txt";
fh:=Open(file,"r");
while true do
	sl:=Gets(fh);
	if IsEof(sl) or #sl lt 2 then
		break;
	end if;
	while sl[#sl] eq "\\" do
		sl:=Prune(sl) cat Gets(fh);
	end while;
	s1:=eval(sl);
	sl:=Gets(fh);
	while sl[#sl] eq "\\" do
		sl:=Prune(sl) cat Gets(fh);
	end while;
	s2:=eval(sl);
	Include(~couples,[Polytope(s1),Polytope(s2)]);
	sl:=Gets(fh);
end while;

//------------------------------------------------------------------------------
// Create all the triples
//------------------------------------------------------------------------------

couples_l:=SetToSequence(couples);
set_couples_of_couples:={};
for i in [1..#couples_l] do
	c:=couples_l[i];
	A:=c[1];
	B:=c[2];
	for j in [i..#couples_l] do
		d:=couples_l[j];
		C:=d[1];
		D:=d[2];
		if AffineNormalForm(A) eq AffineNormalForm(C) then
			Include(~set_couples_of_couples,[[A,B],[C,D]]);
		end if;
		if AffineNormalForm(B) eq AffineNormalForm(C) then
			Include(~set_couples_of_couples,[[B,A],[C,D]]);
		end if;
		if AffineNormalForm(A) eq AffineNormalForm(D) then
			Include(~set_couples_of_couples,[[A,B],[D,C]]);
		end if;
		if AffineNormalForm(B) eq AffineNormalForm(D) then
			Include(~set_couples_of_couples,[[B,A],[D,C]]);
		end if;
	end for;
end for;

triples:={};
i:=0;
for cc in set_couples_of_couples do
	i+:=1;
	c:=cc[1];
	d:=cc[2];
	A:=c[1];
	B:=c[2];
	C:=d[1];
	D:=d[2];
	bool,phi,all_psi:=all_maps(A,C);
	assert bool;
	Bp:=homogenize(B);
	for psi in all_psi do 
		imageB:=dehomogenize(Image(psi,Image(phi,Bp)));
		assert Dimension(imageB) eq 3;
		if NumberOfInteriorPoints(imageB + D) eq 0 then
			Include(~triples,[C,imageB,D]);
		end if;
	end for;
	printf "%o / %o\n",#triples,i;
end for;

// Check for equivalences in the good triples using Cayley
cayleys:={};
perfect_triples:={};
i:=0;
for t in triples do
	i+:=1;
	i;
	C:=cayley(t);
	anfC:=AffineNormalForm(C);
	if not anfC in cayleys then
		Include(~cayleys,AffineNormalForm(C));
		Include(~perfect_triples,t);
		"*";
	end if;
end for;

//------------------------------------------------------------------------------
// By counting the common projection we complete cases (a) and (b)
//------------------------------------------------------------------------------

list12:={};
list23:={};
list13:={};
for t in perfect_triples do
	P1:=t[1];
	P2:=t[2];
	P3:=t[3];
	M:=Ambient(P1);
	
	for e in Edges(P1) do
		endpoints:=[v : v in Vertices(e)];
		p:=endpoints[1]-endpoints[2];
		proj:=project_along_direction(p);
		p1:=Image(proj,P1);
		if Volume(p1) eq 1 then
			p2:=Image(proj,P2);
			if are_translations(p1,p2) then
				Include(~list12,t);
			end if;
		end if;
	end for;
	
	for e in Edges(P1) do
		endpoints:=[v : v in Vertices(e)];
		p:=endpoints[1]-endpoints[2];
		proj:=project_along_direction(p);
		p1:=Image(proj,P1);
		if Volume(p1) eq 1 then
			p3:=Image(proj,P3);
			if are_translations(p1,p3) then
				Include(~list13,t);
			end if;
		end if;
	end for;
	
	for e in Edges(P2) do
		endpoints:=[v : v in Vertices(e)];
		p:=endpoints[1]-endpoints[2];
		proj:=project_along_direction(p);
		p2:=Image(proj,P2);
		if Volume(p2) eq 1 then
			p3:=Image(proj,P3);
			if are_translations(p2,p3) then
				Include(~list23,t);
			end if;
		end if;
	end for;

end for;

// IMPORTANT: at this point only cases (a) and (b) are complete

case_a := perfect_triples diff (list12 join list13 join list23);
case_b := (list12 join list13 join list23) diff ((list12 meet list13) join (list12 meet list23) join (list13 meet list23));

#case_a;// 29
#case_b;// 141

// saved on case_a.txt and case_b.txt

//------------------------------------------------------------------------------
// case (c)
//------------------------------------------------------------------------------
lawrence_prisms:={};
i:=0;
for c in couples do
	i+:=1;
	i;
	P1:=c[1];
	P2:=c[2];
	for e1 in Edges(P1) do
		endpoints:=[v : v in Vertices(e1)];
		p1:=endpoints[1]-endpoints[2];
		proj1:=project_along_direction(p1);
		if Volume(Image(proj1,P1)) eq 1 then
			for e2 in Edges(P2) do
				endpoints:=[v : v in Vertices(e2)];
				p2:=endpoints[1]-endpoints[2];
				proj2:=project_along_direction(p2);
				if Volume(Image(proj2,P2)) eq 1 then
					Include(~lawrence_prisms,c);
					break;
				end if;
			end for;
		end if;
	end for;
end for;

case_c:={};
for c in lawrence_prisms do
	P1:=c[1];
	P2:=c[2];
	for e1 in Edges(P1) do
		endpoints:=[v : v in Vertices(e1)];
		p1:=endpoints[1]-endpoints[2];
		proj1:=project_along_direction(p1);
		if Volume(Image(proj1,P1)) eq 1 then
			Q1:=P1 + Cone([p1,-p1]);//line
			for e2 in Edges(P2) do
				endpoints:=[v : v in Vertices(e2)];
				p2:=endpoints[1]-endpoints[2];
				proj2:=project_along_direction(p2);
				if Volume(Image(proj2,P2)) eq 1 and not are_same_direction(p1,p2) then
					Q2:=P2 + Cone([p2,-p2]);
					for v1 in Vertices(P1) do
						new_Q1:=Q1-v1;
						for v2 in Vertices(P2) do
							new_Q2:=Q2-v2;
							I:=new_Q1 meet new_Q2;
							if IsPolytope(I) then
								I:=Polytope(Points(I));
								if Dimension(I) eq 3 then
									Include(~case_c,[P1,P2,I]);
									for u in Vertices(I) do
										Ip:=Polytope(Points(I) diff {u});
										if Dimension(Ip) eq 3 then
											Include(~case_c,[P1,P2,Ip]);
										end if;
									end for;
								end if;
							end if;
						end for;
					end for;
				end if;
			end for;
		end if;
	end for;
end for;
	
cayleys_case_c:={};
perfect_triples_case_c:={};
i:=0;
for t in case_c do
	i+:=1;
	i;
	C:=cayley(t);
	anfC:=AffineNormalForm(C);
	if not anfC in cayleys_case_c then
		Include(~cayleys_case_c,AffineNormalForm(C));
		Include(~perfect_triples_case_c,t);
		"*";
	end if;
end for;

#perfect_triples_case_c;
case_c := perfect_triples_case_c;
#case_c; \\82

// saved on case_c.txt

