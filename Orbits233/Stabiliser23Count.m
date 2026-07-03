//============================================
//   Stabiliser23Count.m
//============================================
load "MAGMAFUNCTIONS.m";

function Space(q,str)
    /*
        Given q a prime power and str a string among:
          "0","1","2","3","4","4T","5","6","7","7T","8",
          "9","10","11","11T","12","13","14","15","16","17"
        return a representative of the orbit numbered by str in GL(3,q) x GL(3,q) 
        acting on the space of 2x3 matrices over GF(q) by (P,Q).A = P*A*Q^(-1).
    */
    Amb := KMatrixSpace(GF(q),2,3);
    if str eq "0" then
        return sub<Amb |{} >;


    elif str eq "1" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0]} >;
    elif str eq "4" then //====================================================
        return sub<Amb | {[1,0,0, 0,1,0]} >;


    elif str eq "2" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,1,0, 0,0,0]} >;
    elif str eq "4T" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,0,0, 0,1,0]} >;
    elif str eq "5" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,0,0, 0,1,0]} >;
    elif str eq "6" then //====================================================
        return sub<Amb | {[1,0,0, 0,1,0], [0,0,0, 1,0,0]} >;
    elif str eq "7" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,1], [0,1,0, 0,0,0]} >;
    elif str eq "10" then //====================================================
        for uu in GF(q), vv in GF(q) do
            if vv ne GF(q)!0 and forall{ x : x in GF(q) | vv*x^2 + uu*vv*x -1 ne 0 } then
                u := uu;
                v := vv;
                break;
            end if;
        end for;
        return sub<Amb | {[1,u,0, 0,1,0], [0,1,0, v,0,0]} >;
    elif str eq "11" then //====================================================
        return sub<Amb | {[1,0,0, 0,1,0], [0,1,0, 0,0,1]} >;
    

    elif str eq "3" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,1,0, 0,0,0], [0,0,1, 0,0,0]} >;    
    elif str eq "7T" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,0,0, 1,0,0], [0,0,0, 0,1,0]} >;
    elif str eq "8" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,0,0, 0,1,0], [0,0,0, 0,0,1]} >; 
    elif str eq "9" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,1], [0,0,0, 1,0,0], [0,0,0, 0,1,0]} >;
    elif str eq "11T" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,1,0, 1,0,0], [0,0,0, 0,1,0]} >;
    elif str eq "12" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,1], [0,1,0, 0,0,0], [0,0,0, 0,1,0]} >;
    elif str eq "13" then //====================================================
        return sub<Amb | {[1,0,0, 0,1,0], [0,1,0, 0,0,0], [0,0,0, 0,0,1]} >;
    elif str eq "14" then //====================================================
        return sub<Amb | {[1,0,0, 0,0,0], [0,1,0, 0,1,0], [0,0,0, 0,0,1]} >;
    elif str eq "15" then //====================================================
        for uu in GF(q), vv in GF(q) do
            if vv ne GF(q)!0 and forall{ x : x in GF(q) | vv*x^2 + uu*vv*x -1 ne 0 } then
                u := uu;
                v := vv;
                break;
            end if;
        end for;
        return sub<Amb | {[1,u,0, 0,1,0], [0,1,0, v,0,0], [0,0,1, 0,0,0]} >;
    elif str eq "16" then //====================================================
        return sub<Amb | {[1,0,0, 0,1,0], [0,1,0, 0,0,1], [0,0,1, 0,0,0]} >;
    elif str eq "17" then //====================================================
        for aa in GF(q), bb in GF(q) ,cc in GF(q) do
            if forall{ x : x in GF(q) | aa - bb*x + cc*x^2 + x^3 ne 0 } then
                a := aa;
                b := bb;
                c := cc;
                break;
            end if;
        end for;
        return sub<Amb | {[1,0,0, 0,1,0], [0,1,0, 0,0,1], [0,0,1, a,b,c]} >;
    else
        error "Unknown orbit indicator, chose among : 
          \"0\",\"1\",\"2\",\"3\",\"4\",\"4T\",\"5\",\"6\",\"7\",\"7T\",\"8\",
          \"9\",\"10\",\"11\",\"11T\",\"12\",\"13\",\"14\",\"15\",\"16\",\"17\"";
    end if;
end function;

function CountStab(V)
    /*
        Given V a subspace of mxn matrices over GF(q), return the size of the stabilizer of V
        in GL(m,q) x GL(n,q) acting on the space of mxn matrices over GF(q) by (P,Q).A = P*A*Q^(-1).
    */
    m := Nrows(V!0);
    n := Ncols(V!0);
    q := #BaseRing(V);

    EligibleP := {P : P in GL(m,q) | P[1][Min({j : j in [1..m] | P[1][j] ne 0})] eq GF(q)!1};
    EligibleQ := {Q : Q in GL(n,q) | Q[1][Min({j : j in [1..n] | Q[1][j] ne 0})] eq GF(q)!1};

    BasisOfTheSpace := Basis(V);

    /*
        Add other techniques here to reduce the complexity by reducing the number of
        ennumerated matrices : use intersection and sums of right (resp left) kernels of V, use flags, use permutation group for low dimension....
    */

    count := 0;
    for P in EligibleP do
        for Q in EligibleQ do
            if forall{ M : M in BasisOfTheSpace | P*M*(Q^(-1)) in V } then 
                //if (P,Q) works then (aP,bQ) works for all a,b in GF(q)^*, i.e. (q-1)^2 ok.
                count := count + (q-1)^2;
            end if;
        end for;
    end for;
    return count;
end function;

procedure PrintAllStabSizes(q)
    /*
        Given q a prime power, print the sizes of the stabilizers of all orbits 
        of GL(3,q) x GL(3,q) acting on the space of 2x3 matrices over GF(q) 
        by (P,Q).A = P*A*Q^(-1).
    */
    for str in ["0","1","4","2","4T","5","6","7","10","11","3","7T","8","9","11T","12","13","14","15","16","17"] do
        V := Space(q,str);
        printf ("Orbit " cat str cat " : \t");
        printf " %o \n", CountStab(V);
    end for;
end procedure;

//========= This is an example ====================
/*

> V := Space(3,"7T");
> Basis(V);
[
    [1 0 0]
    [0 0 0],

    [0 0 0]
    [1 0 0],

    [0 0 0]
    [0 1 0]
]
> CountStab(V);
2592
> PrintAllStabSizes(3);
Orbit 0 :        539136 
Orbit 1 :        10368 
Orbit 4 :        1728 
Orbit 2 :        10368 
Orbit 4T :       576 
Orbit 5 :        576 
Orbit 6 :        1296 
Orbit 7 :        144 
Orbit 10 :       2304 
Orbit 11 :       96 
Orbit 3 :        134784 
Orbit 7T :       2592 
Orbit 8 :        384 
Orbit 9 :        1296 
Orbit 11T :      1728 
Orbit 12 :       1728 
Orbit 13 :       48 
Orbit 14 :       96 
Orbit 15 :       64 
Orbit 16 :       216 
Orbit 17 :       156 

*/
