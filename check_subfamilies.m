//------------------------------------------------------------------------------
// This script checks that all the families of cases (a), (b) and (c) are 
// are subfamilies of the following six families. This script needs to be run
// after main_script.m
//------------------------------------------------------------------------------
D:=3;

A1:=Polytope([[0,0,0],[2,0,0],[0,2,0],[0,0,2]]);
A2:=Polytope([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]);
A3:=Polytope([[0,0,0],[1,0,0],[0,1,0],[0,0,1]]);

B1:=Polytope([[0,0,0],[2,0,0],[0,2,0],[0,0,1]]);
B2:=Polytope([[0,0,0],[2,0,0],[0,2,0],[0,0,1]]);
B3:=Polytope([[0,0,0],[2,0,0],[0,2,0],[0,0,1]]);

C1:=Polytope([[0,0,0],[2,0,0],[0,1,0],[0,0,1]]);
C2:=Polytope([[0,0,0],[1,0,0],[0,2,0],[0,0,1]]);
C3:=Polytope([[0,0,0],[1,0,0],[0,1,0],[0,0,2]]);

D1:=Polytope([[1,0,0],[0,1,0],[0,-1,0],[0,0,1],[1,0,1]]);
D2:=Polytope([[0,0,0],[1,0,0],[0,1,0],[0,1,1]]);
D3:=Polytope([[0,0,0],[1,0,0],[0,-1,0],[0,-1,1]]);

E1:=Polytope([[0,0,0],[-2,2,0],[0,0,1],[1,0,1]]);
E2:=Polytope([[0,0,0],[-1,1,0],[-1,0,0],[1,0,1]]);
E3:=Polytope([[0,0,0],[0,-1,0],[-1,0,0],[1,-2,1]]);

F1:=Polytope([[0,0,0],[0,2,0],[1,0,1],[0,0,1]]);
F2:=Polytope([[0,0,0],[0,-1,0],[-1,0,0],[1,-2,1]]);
F3:=Polytope([[0,0,0],[-1,0,0],[-1,-1,0],[-1,-2,1]]);

subfamilies:=[[[A1,A2,A3],[B1,B2,B3],[C1,C2,C3],[D1,D2,D3],[E1,E2,E3],[F1,F2,F3]]];
cayleys:={};
i:=1;
while true do
    Append(~subfamilies,[]); // here we'll put the 'good' subfamilies of the ones
                             // created in the previous iteration
    for t in subfamilies[i] do
        for P in t do
            for v in Vertices(P) do
                Q:=Polytope(Exclude(Points(P),v));
                if Dimension(Q) eq D then
                    nt:= Exclude(t,P) cat [Q]; // new triple
                    assert #nt eq 3;
                    C:=cayley(nt);
                    anfC:=AffineNormalForm(C);
                    if not anfC in cayleys then
                        Include(~cayleys,anfC);
                        Include(~subfamilies[i+1],nt);
                    end if;
                end if;
            end for;
        end for;
    end for;
    if #subfamilies[i+1] eq 0 then // if there are no new families, we stop
        break;
    end if;
    i:=i+1;
    [#step : step in subfamilies];
end while;

all:=SequenceToSet(&cat(subfamilies));
list12:={};
list23:={};
list13:={};
for t in all do
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

case_a2 := all diff (list12 join list13 join list23);
case_b2 := (list12 join list13 join list23) diff ((list12 meet list13) join (list12 meet list23) join (list13 meet list23));
case_c2 := ((list12 meet list13) join (list12 meet list23) join (list13 meet list23)) diff (list12 meet list13 meet list23);
[#case_a2,#case_b2,#case_c2];

// Now we check that we got the same results

cayleys_a:={AffineNormalForm(cayley(t)) : t in case_a};
cayleys_a2:={AffineNormalForm(cayley(t)) : t in case_a2};

cayleys_b:={AffineNormalForm(cayley(t)) : t in case_b};
cayleys_b2:={AffineNormalForm(cayley(t)) : t in case_b2};

cayleys_c:={AffineNormalForm(cayley(t)) : t in case_c};
cayleys_c2:={AffineNormalForm(cayley(t)) : t in case_c2};

assert cayleys_a2 eq cayleys_a;
assert cayleys_b2 eq cayleys_b;
assert cayleys_c2 eq cayleys_c;
