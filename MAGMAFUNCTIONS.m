/*
    ===/ Version 02 07 2026 /===
    Author: Lucien François, University College Dublin


    Start your file with 
        load "MAGMAFUNCTIONS.m";
    to use the functions defined in this file.

    Notes: 
       Use the functions MAGMAFUNCTIONS_help() and MAGMAFUNCTIONS_list() to get a list of
       all functions and their descriptions.

*/

MAGMAFUNCTIONS_registery := [[], [], [], []];
// MAGMAFUNCTIONS_registery[1] : list of functions
// MAGMAFUNCTIONS_registery[2] : list of descriptions
// MAGMAFUNCTIONS_registery[3] : list of types (ARITH, CMBTRCS, BIJCTNS)
// MAGMAFUNCTIONS_registery[4] : list of parameters description
function MAGMAFUNCTIONS_add(registery, LISTOFSTUFF)
    THESTUFF := LISTOFSTUFF;
    if #THESTUFF lt 4 then
        THESTUFF := THESTUFF cat <"No info." : i in [#THESTUFF + 1..4]>;
    end if;
    registery[1] := Append(registery[1], THESTUFF[1]);
    registery[2] := Append(registery[2], THESTUFF[2]);
    registery[3] := Append(registery[3], THESTUFF[3]);
    registery[4] := Append(registery[4], THESTUFF[4]);
    return registery;
end function;

function MAGMAFUNCTIONS_FormatError(functionName, parameters, message)
    errorMessage := functionName cat "(\n";
    for i in [1..#parameters] do
        errorMessage cat:= "    " cat parameters[i][1] cat " : " cat Sprint(parameters[i][2]) cat "\n";
    end for;
    errorMessage cat:= ")\n" cat message;
    return errorMessage;
end function;
// Types : ARITH, CMBTRCS, BIJCTNS, TEST, LINALG, DISP





//===========================================================
//          ARITHMETICS
//===========================================================

function MF_PrimePowersUpTo(m)
    /*
        Returns the list of prime powers p^eta that are less or equal than m.
    */
    if m lt 1 then
        return [];
    end if;
    ListPrimes := PrimesUpTo(m);
    Result := [];
    for p in ListPrimes do
        Result := Result cat [p^eta : eta in [1..Floor(Log(p, m))]];
    end for;
    return Sort(Result);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_PrimePowersUpTo(m)", "Returns the list of prime powers p^eta that are less or equal than m.", "ARITH", "\tm: positive integer.">
);

function MF_LongestConsecutiveSequence(SeqInt)
    /*
        Given a sequence of integers, find the longest consecutive sequence.
        Returns <BestFound,BestLength> where BestFound = <start,end> are two elements of SeqInt. For instance, if SeqInt = [1,2,3,5,6,10], the function returns < <1,3>, 3 > since the longest consecutive sequence is 1,2,3 and has length 3.
    */
    if SeqInt ne Sort(SeqInt) then
        return Sort(SeqInt);
    end if;


    if #SeqInt eq 0 then
        return <<0,-1>,0>;
    end if;
    BestFound := <SeqInt[1],SeqInt[1]>;
    BestLength := 1;
    Candidate := BestFound;
    CurrentLength := BestLength;
    for t in [2..#SeqInt] do
        //printf "%o,%o,%o,%o,%o\n",t,BestFound,BestLength,Candidate,CurrentLength;
        //Continue candidate?
        IsConsecutive := (SeqInt[t] eq Candidate[2]+1);
        if IsConsecutive then //A consecutive integer: continue the current sequence
            Candidate := <Candidate[1],SeqInt[t]>;
            CurrentLength := CurrentLength +1;
        end if;
        if (t eq #SeqInt or not IsConsecutive) then //if we have to stop the current sequence, check if best
            if CurrentLength gt BestLength then 
                BestFound := Candidate;
                BestLength := CurrentLength;
            end if;
        end if;
        if not IsConsecutive then //if not connsecutive, reset current sequence.
            Candidate := <SeqInt[t],SeqInt[t]>;
            CurrentLength := 1;
        end if;    
    end for;
    return <BestFound,BestLength>;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_LongestConsecutiveSequence(SeqInt)", "Given an ordered sequence of integers, find the longest consecutive sequence. Returns <BestFound,BestLength> where BestFound = <start,end> are two elements of SeqInt.", "ARITH", "\tSeqInt: sequence of integers.">
);




//===========================================================
//          COMBINATORICS
//===========================================================

function MF_GBC(q, n, r)
    /*
        Returns the number of r-dimensional subspaces in Fq^n.
    */
    if #Factorization(q) ne 1 then
        error MAGMAFUNCTIONS_FormatError("MF_GBC", < <"q", q>, <"n", n>, <"r", r> >, "q must be a prime power.");
    end if;
    if r lt 0 or n lt 0 then
        error MAGMAFUNCTIONS_FormatError("MF_GBC", < <"q", q>, <"n", n>, <"r", r> >, "n and r must be nonnegative integers.");
    end if;
    if r gt n then
        return 0;
    end if;
    return Floor(&*([1 - q^(n - t) : t in [0..r - 1]] cat [1]) / &*([1 - q^t : t in [1..r]] cat [1]));
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_GBC(q,n,r)", "Returns the number of r-dimensional subspaces in Fq^n.", "CMBTRCS", "\tq: prime power,\n\tn,r: non-negative integers.">
);

function MF_NbMatrixRank(n, m, q, r)
    /*
        Counts the number of matrices of rank r of size nxm over Fq.
    */
    if #Factorization(q) ne 1 then
        error MAGMAFUNCTIONS_FormatError("MF_NbMatrixRank", < <"n", n>, <"m", m>, <"q", q>, <"r", r> >, "q must be a prime power.");
    end if;
    if r lt 0 or n lt 0 or m lt 0 then
        error MAGMAFUNCTIONS_FormatError("MF_NbMatrixRank", < <"n", n>, <"m", m>, <"q", q>, <"r", r> >, "n, m, q, and r must be nonnegative integers.");
    end if;
    res := 0;
    for k in [0..r] do
        res := res + (-1)^(r - k) * MF_GBC(q, r, k) * q^(n * k + Binomial(r - k, 2));
    end for;
    return MF_GBC(q, m, r) * res;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_NbMatrixRank(n,m,q,r)", "Counts the number of matrices of rank r of size nxm over Fq.", "CMBTRCS", "\tn,m: non-negative integers,\n\tq: prime power,\n\tr: non-negative integer.">
);

procedure MF_NbMatrixRankTest()
    for test in [1..10] do
        m := Random([1..5]);
        n := Random([1..5]);
        q := Random([2,3,4,5,7]);
        while q^(m*n) gt 10^9 do
            q := Random([2,3,4,5,7]);
            m := Random([1..5]);
            n := Random([1..5]);
        end while;
        r := Random([0..Min(m, n)]);
        printf "Test %o: m = %o, n = %o, q = %o, r = %o, ", test, m, n, q, r;
        count := 0;
        for M in KMatrixSpace(GF(q), m, n) do
            if Rank(M) eq r then
                count +:= 1;
            end if;
        end for;
        if count ne MF_NbMatrixRank(n, m, q, r) then
            error MAGMAFUNCTIONS_FormatError("MF_NbMatrixRankTest", < <"test", test>, <"m", m>, <"n", n>, <"q", q>, <"r", r>, <"count", count>, <"returned", MF_NbMatrixRank(n, m, q, r)> >, "Test failed: counted a different number of matrices of rank r than MF_NbMatrixRank returned.");
        else 
            printf "OK.\n";
        end if;
    end for;
    printf "All tests passed for MF_NbMatrixRank.\n";
end procedure;

function MF_MRDrankdistribution(q,m,n,d,r)    /*
        Given q a prime power, given m,n non-negative integers, given d a positive integer, and given r a non-negative integer,
        returns the number of matrices of rank r in an MRD rank-metric code of minimum distance d in the space of mxn matrices over Fq.
    */

    if #Factorization(q) ne 1 then
        error MAGMAFUNCTIONS_FormatError("MF_MRDrankdistribution", < <"q", q>, <"m", m>, <"n", n>, <"d", d>, <"r", r> >, "q must be a prime power.");
    end if;
    if m lt 0 or n lt 0 or d lt 1 or r lt 0 then
        error MAGMAFUNCTIONS_FormatError("MF_MRDrankdistribution", < <"q", q>, <"m", m>, <"n", n>, <"d", d>, <"r", r> >, "m and n must be nonnegative integers, d must be a positive integer, and r must be a nonnegative integer.");
    end if;

    M := Min(m, n);
    N := Max(m, n);

    if r gt M or (r gt 0 and r lt d) then 
        return 0;
    end if;
    if r eq 0 then
        return 1;
    end if;
    
    res := 0;
    k := N*(M-d+1);
    for j in [0..r-d] do
        res := res + (-1)^j * q^(Binomial(j, 2)) * MF_GBC(q,r,j) *(q^(k - N*(M + j -r)) -1);
    end for;
    return (MF_GBC(q, M, r) * res);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_MRDrankdistribution(q,m,n,d,r)", "Given q a prime power, given m,n non-negative integers, given d a positive integer, and given r a non-negative integer, returns the number of matrices of rank r in an MRD rank-metric code of minimum distance d in the space of mxn matrices over Fq.", "CMBTRCS", "\tq: prime power,\n\tm,n: non-negative integers,\n\td: positive integer,\n\tr: non-negative integer.">
);


function MF_MRDrankenumerator(q,m,n,d,x,y)
    /*
        Given q a prime power, given m,n non-negative integers, given d a positive integer, and given x,y variables,
        returns the rank enumerator of an MRD rank-metric code of minimum distance d in the space of mxn matrices over Fq.
    */
    if #Factorization(q) ne 1 then
        error MAGMAFUNCTIONS_FormatError("MF_MRDrankenumerator", < <"q", q>, <"m", m>, <"n", n>, <"d", d>, <"x", x>, <"y", y> >, "q must be a prime power.");
    end if;
    if m lt 0 or n lt 0 or d lt 1 then
        error MAGMAFUNCTIONS_FormatError("MF_MRDrankenumerator", < <"q", q>, <"m", m>, <"n", n>, <"d", d>, <"x", x>, <"y", y> >, "m and n must be nonnegative integers, and d must be a positive integer.");
    end if;
    return &+ [MF_MRDrankdistribution(q, m, n, d, r) * x^(Min(m, n) - r) * y^r : r in [0..Min(m, n)]];
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_MRDrankenumerator(q,m,n,d,x,y)", "Given q a prime power, given m,n non-negative integers, given d a positive integer, and given x,y variables, returns the rank enumerator of an MRD rank-metric code of minimum distance d in the space of mxn matrices over Fq.", "CMBTRCS", "\tq: prime power,\n\tm,n: non-negative integers,\n\td: positive integer,\n\tx,y: variables.">
);

function MF_MRDeigenvaluetrilinearforms(q,n1,n2,n3,d,mode)
    /*
        Given q a prime power, given n1,n2,n3 non-negative integers, and given d a positive integer, returns 
        the eigenvalue in the cayley graph of tensors Fq^(n1,n2,n3) generated by rank one tensors, associated to 
        a tensor that generates along the mode-th mode an MRD code of parameters [MxN,...,d], where {M,N,mode} = {n1,n2,n3}.
    */
    if #Factorization(q) ne 1 then
        error MAGMAFUNCTIONS_FormatError("MF_MRDeigenvaluetrilinearforms", < <"q", q>, <"n1", n1>, <"n2", n2>, <"n3", n3>, <"d", d>, <"mode", mode> >, "q must be a prime power.");
    end if;
    if n1 lt 0 or n2 lt 0 or n3 lt 0 or d lt 1 then
        error MAGMAFUNCTIONS_FormatError("MF_MRDeigenvaluetrilinearforms", < <"q", q>, <"n1", n1>, <"n2", n2>, <"n3", n3>, <"d", d>, <"mode", mode> >, "n1, n2, and n3 must be nonnegative integers, and d must be a positive integer.");
    end if;
    if Type(n1) ne RngIntElt or Type(n2) ne RngIntElt or Type(n3) ne RngIntElt or Type(d) ne RngIntElt then
        error MAGMAFUNCTIONS_FormatError("MF_MRDeigenvaluetrilinearforms", < <"q", q>, <"n1", n1>, <"n2", n2>, <"n3", n3>, <"d", d>, <"mode", mode> >, "n1, n2, n3, and d must be integers.");
    end if;

    if mode eq 1 then
        M := Min({n2,n3});
        N := Max({n2,n3});
    elif mode eq 2 then
        M := Min({n1,n3});
        N := Max({n1,n3});
    elif mode eq 3 then
        M := Min({n1,n2});
        N := Max({n1,n2});
    else
        error MAGMAFUNCTIONS_FormatError("MF_MRDeigenvaluetrilinearforms", < <"q", q>, <"n1", n1>, <"n2", n2>, <"n3", n3>, <"d", d>, <"mode", mode> >, "mode must be in {1,2,3}.");
    end if;

    k := N*(M-d+1);
    return (q-1)^(-2) * ( (q^(n1)-1)*(q^(n2)-1)*(q^(n3)-1) - q^(n1+n2+n3) + q^(n1+n2+n3-k)*MF_MRDrankenumerator(q, M, N, d, 1, 1/q ) );
end function;

function MF_NbTrankOne(q,N)
    /*
        Counts the number of tensors of trank 1 in (Fq^N1 otimes ... otimes Fq^Nk)
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_NbTrankOne", < <"q", q>, <"N", N> >, "N must be a non-empty sequence of nonnegative integers.");
    end if;
    if exists{j : j in [1..#N] | N[j] lt 0} then
        error MAGMAFUNCTIONS_FormatError("MF_NbTrankOne", < <"q", q>, <"N", N> >, "All entries of N must be nonnegative integers.");
    end if;
    return Floor(&*([q^n - 1 : n in N] cat [1]) / (q-1)^#N);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_NbTrankOne(q,N)", "Counts the number of tensors of trank 1 in (Fq^N1 otimes ... otimes Fq^Nk).", "CMBTRCS", "\tq: prime power,\n\tN = [N1,...,Nk]: non-empty sequence of nonnegative integers.">
);

function MF_NbTrankTwo(q,N)
    /*
        Counts the number of tensors of trank 2 in (Fq^N1 otimes ... otimes Fq^Nk)
        KNOWN ONLY FOR t = 1,2,3.
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_NbTrankTwo", < <"q", q>, <"N", N> >, "N must be a non-empty sequence of nonnegative integers.");
    end if;
    if #N eq 1 then
        return 0;
    end if;
    if #N eq 2 then
        return MF_NbMatrixRank(N[1], N[2], q, 2);
    end if;
    if #N gt 3 then
        error MAGMAFUNCTIONS_FormatError("MF_NbTrankTwo", < <"q", q>, <"N", N> >, "The formula for counting tensors of trank 2 is only known for k = 1, 2, 3.");
    end if;
    //==== Compute q(q-1)(q2-1) (A + BC) ====
    A := (q^N[1]-1)*(q^N[2]-1)*(q^(N[2]-1)-1)*(q^N[3]-1)*(q^(N[3]-1)-1)/((q-1)^3 * (q^2-1)^2);
    B := (q^N[1]-1)*(q^(N[1]-1)-1)*(q^N[2]-1)*(q^N[3]-1)/((q-1)^3 * (q^2-1));
    C := ( (q^(N[2]-1) -1)/(q-1) + (q^(N[3]-1) -1)/(q-1) )/(q+1) 
                    + 0.5*q^2*((q^(N[2]-1)-1)*(q^(N[3]-1)-1))/((q-1)^2);
    return Floor(q*(q-1)*(q^2-1)*(A + B*C));
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_NbTrankTwo(q,N)", "Counts the number of tensors of trank 2 in (Fq^N1 otimes ... otimes Fq^Nk) for k=1,2,3.", "CMBTRCS", "\tq: prime power,\n\tN = [N1,...,Nk]: non-empty sequence of nonnegative integers with k=1,2,3.">
);

function MF_SmallestRectangle(SET,n)
    /*
    Given a SET of plane integer points in [[0,n-1]]², stored as length-2 sequences, find the 
    smallest rectangle modulo n containing all the points.
    */
    SETx := { s[1] : s in SET} join { s[1] - n : s in SET};
    SETy := { s[2] : s in SET} join { s[2] - n : s in SET};
    CPLTx := SetToSequence({-n..n-1} diff SETx);
    Sort(~CPLTx);
    CPLTy := SetToSequence({-n..n-1} diff SETy);
    Sort(~CPLTy);
    RESx := MF_LongestConsecutiveSequence(CPLTx);
    RESy := MF_LongestConsecutiveSequence(CPLTy);
    return <<RESx[1][2]+1,n+RESx[1][1]-1>,<RESy[1][2]+1,n+RESy[1][1]-1>,n-1-RESx[2]+1,n-1-RESy[2]+1>;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_SmallestRectangle(SET,n)", "Given a SET of plane integer points in [[0,n-1]]², stored as length-2 sequences, find the smallest rectangle modulo n containing all the points. Returns <<x_min,x_max>,<y_min,y_max>,width,height>.", "CMBTRCS", "\tSET: set of length-2 sequences representing points in [[0,n-1]]².\n\tn: integer.">
);




//===========================================================
//          BIJECTIONS OF COUNTABLE SETS
//===========================================================

function MF_IntToqary(n, q, size)
    /*
        Converts an integer n in {0..q^size-1} to a sequence in {0,..,q-1} of length size
        such that n = sum_{i=0}^{size-1} a_i * (q^(i)) .
    */
    if n lt 0 or n ge q^size then
        error MAGMAFUNCTIONS_FormatError("MF_IntToqary", < <"n", n>, <"q", q>, <"size", size> >, "n must be in the interval [0..q^size - 1].");
    end if;
    nn := n;
    seq := [0 : i in [1..size]];
    for i in [1..size] do
        seq[i] := nn mod q;
        nn := nn - seq[i];
        nn := Floor(nn / q);
    end for;
    return seq;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_IntToqary(n,q,size)", "Converts an integer n in {0..q^size-1} to a sequence in {0,..,q-1} of length size such that n = sum_{i=0}^{size-1} a_i * (q^(i)) .", "BIJCTNS", "\tn: integer,\n\tq: integer,\n\tsize: integer.">
);

function MF_qarytoInt(seq, q)
    /*
        Converts a sequence in {0,..,q-1} of length size to an integer 0 <= n < q^(size)
        such that n = sum_{i=0}^{size-1} a_i * (q^(i)) .
    */
    if exists{j : j in [1..#seq] | seq[j] lt 0 or seq[j] ge q} then
        error MAGMAFUNCTIONS_FormatError("MF_qarytoInt", < <"seq", seq>, <"q", q> >, "All entries of seq must be in {0,.., q - 1}.");
    end if;
    size := #seq;
    n := 0;
    for i in [1..size] do
        if seq[i] lt 0 or seq[i] ge q then
            error MAGMAFUNCTIONS_FormatError("MF_qarytoInt", < <"seq", seq>, <"q", q> >, "All entries of seq must be in {0,.., q - 1}.");
        end if;
        n +:= seq[i] * (q^(i - 1));
    end for;
    return n;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_qarytoInt(seq,q)", "Converts a sequence in {0,..,q-1} of length size to an integer 0 <= n < q^(size) such that n = sum_{i=0}^{size-1} a_i * (q^(i)) .", "BIJCTNS", "\tseq: sequence of integers,\n\tq: integer.">
);

procedure MF_BijectionqaryTest(q, size)
    /*
        Tests the bijection between integers and q-ary sequences of given size.
    */
    for n in [0..q^size - 1] do
        printf "n = %o: ", n;
        seq := MF_IntToqary(n, q, size);
        printf "seq = %o; ", seq;
        m := MF_qarytoInt(seq, q);
        printf "m = %o\n", m;
        if n ne m then
            error MAGMAFUNCTIONS_FormatError("MF_BijectionqaryTest", < <"q", q>, <"size", size>, <"n", n>, <"m", m> >, "Bijection test failed: the round-trip n -> seq -> m was not exact.");
        end if;
    end for;
    printf "Bijection test passed for all n in [0..%o] over GF(%o) with size = %o.\n", q^size - 1, q, size;
end procedure;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_BijectionqaryTest(q,size)", "Tests the bijections MF_IntToqary(n,q,size) and MF_qarytoInt(seq,q) between integers and q-ary sequences of given size.", "TEST", "\tq: integer,\n\tsize: integer.">
);

function MF_IntToFq(n, q)
    /*
        Converts an integer n to a vector in GF(q).
        Details:
              [0..q-1] -> GF(q)
                     0 -> 0
            non-zero n -> primitive element ^ (n-1)
    */
    if n eq 0 then
        return GF(q)!0;
    end if;
    return PrimitiveElement(GF(q))^(n - 1);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_IntToFq(n,q)", "Converts an integer n to a vector in GF(q) via the function\n\t[0..q-1] -> GF(q),\n\t0 -> 0,\n\tnon-zero n -> primitive element ^ (n-1).", "BIJCTNS", "\tn: integer,\n\tq: prime power.">
);

function MF_FqtoInt(x)
    /*
        Converts a vector x in GF(q) to an integer.
        Details:
              GF(q) -> [0..q-1]
                     0 -> 0
            non-zero x -> discrete log base primitive element + 1
    */
    q := #(Parent(x));
    if x eq GF(q)!0 then
        return 0;
    end if;
    return Log(x) + 1;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_FqtoInt(x)", "Converts a vector x in GF(q) to an integer via the function\n\tGF(q) -> [0..q-1],\n\t0 -> 0,\n\tnon-zero x -> discr log base primitive element + 1.", "BIJCTNS", "\tx: element of GF(q).">
);

procedure MF_BijectionFieldTest(q)
    /*
        Tests the bijection between integers and vectors in GF(q).
    */
    for n in [0..q - 1] do
        printf "n = %o: ", n;
        x := MF_IntToFq(n, q);
        printf "x = %o; ", x;
        m := MF_FqtoInt(x);
        printf "m = %o\n", m;
        if n ne m then
            error MAGMAFUNCTIONS_FormatError("MF_BijectionFieldTest", < <"q", q>, <"n", n>, <"m", m> >, "Bijection test failed: the round-trip n -> x -> m was not exact.");
        end if;
    end for;
    printf "Bijection test passed for all n in [0..%o] over GF(%o).\n", q - 1, q;
end procedure;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_BijectionFieldTest(q)", "Tests the bijections MF_IntToFq(n,q) and MF_FqtoInt(x) between integers and vectors in GF(q) for n in [0..q-1].", "TEST", "\tq: prime power.">
);


function MF_IntToVSpace(n, q, Vdim)
    /*
        Converts an integer 0 <= n < q^(Vdim) to a vector in VectorSpace(GF(q),Vdim).
        Details:
                [0..Vdim-1] -> VectorSpace(GF(q),Vdim)
                n = sum_{i=0}^{Vdim-1} a_i * (q^(i)) -> (...,MF_IntToFq(a_i,q),...)
    */
    seq := MF_IntToqary(n, q, Vdim);
    x := VectorSpace(GF(q), Vdim)!0;
    for i in [1..Vdim] do
        x[i] := MF_IntToFq(seq[i], q);
    end for;
    return x;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_IntToVSpace(n,q,Vdim)", "Converts an integer 0 <= n < q^(Vdim) to a vector in VectorSpace(GF(q),Vdim) using the function\n\t[0..q^Vdim-1] -> VectorSpace(GF(q),Vdim),\n\tn = sum_{i=0}^{Vdim-1} a_i * (q^(i)) -> (...,IntToFq(a_i,q),...).", "BIJCTNS", "\tn: integer,\n\tq: prime power,\n\tVdim: positive integer.">
);

function MF_VSpaceToInt(x)
    /*
        Converts a vector x in Parent(x)=VectorSpace(GF(q),Vdim) to an integer.
        Details:
                VectorSpace(GF(q),Vdim) -> [0..Vdim-1]
                (x_1,...,x_k) -> n = sum_{i=0}^{Vdim-1} FqtoInt(x_i) * (q^(i))
    */
    q := #(Parent(x));
    return &+ [MF_FqtoInt(x[i]) * (q^(i - 1)) : i in [1..Dimension(Parent(x))]];
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_VSpaceToInt(x)", "Converts a vector x in Parent(x)=VectorSpace(GF(q),Vdim) to an integer using the function\n\tVectorSpace(GF(q),Vdim) -> [0..Vdim-1],\n\t(x_1,...,x_k) -> n = sum_{i=0}^{Vdim-1} MF_FqtoInt(x_i) * (q^(i))", "BIJCTNS", "\tx: element of VectorSpace(GF(q),Vdim).">
);

procedure MF_BijectionTest()
    /*
        Tests the bijection between integers and vectors in X.
    */
    theq := Random([2, 3, 4, 5]);
    Vdim := Random([2..5]);
    for n in [0..theq^Vdim - 1] do
        printf "n = %o: ", n;
        x := MF_IntToVSpace(n, theq, Vdim);
        printf "x = %o; ", x;
        m := MF_VSpaceToInt(x);
        printf "m = %o\n", m;
        if n ne m then
            error MAGMAFUNCTIONS_FormatError("MF_BijectionTest", < <"q", theq>, <"Vdim", Vdim>, <"n", n>, <"m", m> >, "Bijection test failed: the round-trip n -> x -> m was not exact.");
        end if;
    end for;
    printf "Bijection test passed for all n in [0..%o].\n", theq^Vdim - 1;
end procedure;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_BijectionTest()", "Tests the bijections MF_IntToVSpace(n,q,Vdim) and MF_VSpaceToInt(x) between integers and vectors in VectorSpace(GF(q),Vdim) for n in [0..q^Vdim-1] with random q in [2..5] and random Vdim in [2..5].", "TEST", "\t*no parameter*">
);

function MF_IntToIntTuple(N,i)
    /*
        Given N = [n1,...,nk] a sequence of positive integers and i an integer in [0..n1*...*nk-1], returns the unique tuple (i1,...,ik) such that i = sum_{j=1}^k (ij-1) * (n1*...*n(j-1)).
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_IntToIntTuple", < <"N", N>, <"i", i> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if exists(j){ j : j in [1..#N] | N[j] le 0 } then
        error MAGMAFUNCTIONS_FormatError("MF_IntToIntTuple", < <"N", N>, <"i", i> >, "All entries of N must be positive integers.");
    end if;
    if i lt 0 or i ge &*N then
        error MAGMAFUNCTIONS_FormatError("MF_IntToIntTuple", < <"N", N>, <"i", i> >, "i must be in the interval [0..n1*...*nk-1] where N = [n1,...,nk].");
    end if;
    res := [];
    for j in [#N..1 by -1] do
        res[j] := (i mod N[j]) + 1;
        i := Floor(i / N[j]);
    end for;
    return res;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_IntToIntTuple(N,i)", "Given N = [n1,...,nk] a sequence of positive integers and i an integer in [0..n1*...*nk-1], returns the unique tuple (i1,...,ik) such that i = sum_{j=1}^k (ij-1) * (n1*...*n(j-1)).", "BIJCTNS", "\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\ti: integer.">
);
function MF_IntTupleToInt(N,i)
    /*
        Given N = [n1,...,nk] a sequence of positive integers and i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the unique integer in [0..n1*...*nk-1] corresponding to the tuple (i1,...,ik) via the function i = sum_{j=1}^k (ij-1) * (n1*...*n(j-1)).
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_IntTupleToInt", < <"N", N>, <"i", i> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if exists(j){ j : j in [1..#N] | N[j] le 0 } then
        error MAGMAFUNCTIONS_FormatError("MF_IntTupleToInt", < <"N", N>, <"i", i> >, "All entries of N must be positive integers.");
    end if;
    if exists(j){ j : j in [1..#i] | i[j] lt 1 or i[j] gt N[j] } then
        error MAGMAFUNCTIONS_FormatError("MF_IntTupleToInt", < <"N", N>, <"i", i> >, "All entries of i must be in the range [1..nj] where N = [n1,...,nk].");
    end if;
    if #N ne #i then
        error MAGMAFUNCTIONS_FormatError("MF_IntTupleToInt", < <"N", N>, <"i", i> >, "N and i must have the same length.");
    end if;
    res := 0;
    i_from_zero := [i[j] - 1 : j in [1..#i]];
    for j in [1..#N] do
        res +:= i_from_zero[j] * &*([1] cat [N[t] : t in [j+1..#N]]); 
    end for;
    return res;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_IntTupleToInt(N,i)", "Given N = [n1,...,nk] a sequence of positive integers and i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the unique integer in [0..n1*...*nk-1] corresponding to the tuple (i1,...,ik) via the function i = sum_{j=1}^k (ij-1) * (n1*...*n(j-1)).", "BIJCTNS", "\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\ti = [i1,...,ik]: sequence of integers.">
);
procedure MF_BijectionIntTupleTest(N)
    /*
        Tests the bijection between integers and tuples of integers defined by MF_IntToIntTuple and MF_IntTupleToInt for a given N = [n1,...,nk].
    */
    total := &*N;
    for i in [0..total - 1] do
        printf "i = %o: ", i;
        tuple := MF_IntToIntTuple(N, i);
        printf "tuple = %o; ", tuple;
        m := MF_IntTupleToInt(N, tuple);
        printf "m = %o\n", m;
        if i ne m then
            error MAGMAFUNCTIONS_FormatError("MF_BijectionIntTupleTest", < <"N", N>, <"i", i>, <"m", m> >, "Bijection test failed: the round-trip i -> tuple -> m was not exact.");
        end if;
    end for;
    printf "Bijection test passed for all i in [0..%o] where N = %o.\n", total - 1, N;
end procedure;

//===========================================================
//          LINEAR ALGEBRA
//===========================================================

function MF_RandomBasisOfTheField(Fn, F, n)
    /*
        Returns a random basis of the extension Fn of degree n over F.
    */
    Fnvs, phi := VectorSpace(Fn, F);
    B := [Random(Fnvs) : i in [1..n]];
    while Dimension(sub<Fnvs | B>) ne n do
        B := [Random(Fnvs) : i in [1..n]];
    end while;
    return Inverse(phi)(B);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_RandomBasisOfTheField(Fn,F,n)", "Returns a random basis of the extension Fn of degree n over F.", "LINALG", "\tFn: field,\n\tF: subfield of Fn,\n\tn: positive integer.">
);

function MF_FirstHyperdet(TensorAsVector,n)
    /*
    Computes the first Cayley hyperdeterminant of an n x n x n tensor given as a vector of size n^3.
    */
    if n lt 1 then
        error MAGMAFUNCTIONS_FormatError("MF_FirstHyperdet", < <"TensorAsVector", TensorAsVector>, <"n", n> >, "n must be a positive integer.");
    end if;
    theRing := Parent(TensorAsVector);
    Indexes := {Integers() | x : x in [1..n]};
    Sn := Sym(Indexes);

    CutSizesOfT := [];
    for i in [1..n] do
        CutSizesOfT[i] := n*n;
    end for;

    Slices := Partition(TensorAsVector,CutSizesOfT);

    res := theRing ! 0;
    for s in Sn do
        for t in Sn do
            Prod := Sign(s*t);
            for i in [1..n] do
                Prod := Prod* (theRing ! Slices[i][n*(i^s-1)+i^t]);
            end for;
            res := res + Prod;
        end for;
    end for;

    return(res);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_FirstHyperdet(TensorAsVector,n)", "Computes the first Cayley hyperdeterminant of an n x n x n tensor given as a vector of size n^3.", "LINALG", "\tTensorAsVector: vector of size n^3 representing the entries of the tensor,\n\tn: positive integer.">
);

function MF_TensorProduct(T)
    /*
        Given a tuple T = <T1,...,Tk> of vectors over the same field, 
        returns the tensor product T1 otimes ... otimes Tk as a vector.
        If Ti is in VectorSpace(F,ni) for i in [1..k], then the output is in VectorSpace(F, n1*...*nk).
        This function is compatible with the other tensor functions in this package.
    */
    if #T eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_TensorProduct", < <"T", T> >, "T must be a non-empty sequence of vectors.");
    end if;
    if #T eq 1 then
        return T[1];
    end if;
    return TensorProduct(T[1], MF_TensorProduct(<T[i] : i in [2..#T]>));
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorProduct(T)", "Given a tuple T = <T1,...,Tk> of vectors over the same field, returns the tensor product T1 otimes ... otimes Tk as a vector. If Ti is in VectorSpace(F,ni) for i in [1..k], then the output is in VectorSpace(F, n1*...*nk). This function is compatible with the other tensor functions in this package.", "LINALG", "\tT: tuple of vectors over the same field.">
);

function MF_Segre(N,q)
    /*
        Given N = [n1,...,nk] a sequence of positive integers and q a prime power, returns the Segre variety as a set of elements in VectorSpace(GF(q), n1*...*nk).
    */
    U := <{1..q^N[1] -1}>;
    for j in [2..#N] do
        Ui_int := {q^(N[j]-1)} join &join[{q^r*(1+q*i) : i in [0..q^(N[j]-r-1)-1]} : r in [0..N[j]-2]];
        Append(~U, Ui_int);
    end for;
    
    S := {};
    for u_int in CartesianProduct(U) do
        x := MF_TensorProduct(<MF_IntToVSpace(u_int[j], q, N[j]) : j in [1..#N]>);
        Include(~S, x);
    end for;

    return S;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_Segre(N,q)", "Given N = [n1,...,nk] a sequence of positive integers and q a prime power, returns the Segre variety as a set of elements in VectorSpace(GF(q), n1*...*nk).", "LINALG", "\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\tq: prime power.">
);

function MF_ElementaryTensor(N,F,i)
    /*
        Given N = [n1,...,nk] a sequence of positive integers,
        given F a field and given i = [i1,...,ik] a sequence of integers 
        such that ij in [1..nj] for j in [1..k], returns the elementary tensor e_{i1} otimes ... otimes e_{ik} where {e_1,...,e_nj} is the canonical basis of VectorSpace(F,nj).
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_ElementaryTensor", < <"N", N>, <"F", F>, <"i", i> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if exists(j){ j : j in [1..#N] | N[j] le 0 } then
        error MAGMAFUNCTIONS_FormatError("MF_ElementaryTensor", < <"N", N>, <"F", F>, <"i", i> >, "All entries of N must be positive integers.");
    end if;
    if exists(j){ j : j in [1..#i] | i[j] lt 1 or i[j] gt N[j] } then
        error MAGMAFUNCTIONS_FormatError("MF_ElementaryTensor", < <"N", N>, <"F", F>, <"i", i> >, "All entries of i must be in the range [1..nj] where N = [n1,...,nk].");
    end if;
    if #N ne #i then
        error MAGMAFUNCTIONS_FormatError("MF_ElementaryTensor", < <"N", N>, <"F", F>, <"i", i> >, "N and i must have the same length.");
    end if;
    return MF_TensorProduct(<VectorSpace(F, N[j]).(i[j] ) : j in [1..#N]>);
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_ElementaryTensor(N,F,i)", "Given N = [n1,...,nk] a sequence of positive integers, given F a field and given i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the elementary tensor e_{i1} otimes ... otimes e_{ik} where {e_1,...,e_nj} is the canonical basis of VectorSpace(F,nj).", "LINALG", "\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\tF: field,\n\ti = [i1,...,ik]: sequence of integers.">
);
procedure MF_PrintTheCoordinates(N)
    /*
        Given N = [n1,...,nk] a sequence of positive integers, returns a sequence of n1*...*nk sequences of the form [i1,...,ik] with ij in [1..nj] for j in [1..k], such that if [i1,...,ik] is the j-th sequence in the output, then MF_ElementaryTensor(N,F,[i1,...,ik]) is the j-th element of the canonical basis of VectorSpace(F,n1*...*nk) for any field F.
    */
    F := GF(2);
    Res := [[-1 : i in [1..#N]] : j in [1..&*N]];
    for i_int in [1..&*N] do
        i := MF_IntToIntTuple(N, i_int - 1);
        elem := MF_ElementaryTensor(N, F, i);
        j := 1;
        while elem[j] eq 0 do
            j +:= 1;
        end while;
        Res[j] := i;
    end for;
    
    for j in [1..#Res] do
        printf "%o (inttuple+1=%o) : %o\n", j, MF_IntTupleToInt(N, Res[j])+1, Res[j];
    end for;
end procedure;

function MF_TensorCoordToInt(N,i)
    /*
        Given N = [n1,...,nk] a sequence of positive integers, given i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the integer in [1..n1*...*nk] such that for each tensor Ten in VectorSpace(F,n1*...*nk) for some field F, its coefficient at the inputed coordinate i is exacly Ten[j].
    */
    return MF_IntTupleToInt(N, i) + 1;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorCoordToInt(N,i)", "Given N = [n1,...,nk] a sequence of positive integers, given i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the integer in [1..n1*...*nk] such that for each tensor Ten in VectorSpace(F,n1*...*nk) for some field F, its coefficient at the inputed coordinate i is exacly Ten[j].", "LINALG", "\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\ti = [i1,...,ik]: sequence of integers.">
);

function MF_TensorCoefficient(Ten,N,i)
    /*
        Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, given i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the coefficient of the elementary tensor e_{i1} otimes ... otimes e_{ik} in the expression of Ten in the basis of such tensors.
    */
    i_int := MF_IntTupleToInt(N, i) +1;
    return Ten[i_int];
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorCoefficient(Ten,N,i)", "Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, given i = [i1,...,ik] a sequence of integers such that ij in [1..nj] for j in [1..k], returns the coefficient of the elementary tensor e_{i1} otimes ... otimes e_{ik} in the expression of Ten in the basis of such tensors.", "LINALG", "\tTen: tensor in VectorSpace(F,n1*...*nk),\n\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\ti = [i1,...,ik]: sequence of integers.">
);
procedure MF_TensorCoefficientTest()
    for test in [1..10] do
        k := Random([1..5]);
        N := [Random([2..5]) : i in [1..k]];
        F := GF(2);
        i := [Random([1..N[j]]) : j in [1..k]];
        Ten := MF_ElementaryTensor(N, F, i);
        if MF_TensorCoefficient(Ten, N, i) ne 1 then
            error MAGMAFUNCTIONS_FormatError("MF_TensorCoefficientTest", < <"N", N>, <"i", i>, <"Ten", Ten> >, "Tensor coefficient test failed: the coefficient of the elementary tensor in its own expression was not 1.");
        end if;
    end for;
    printf "Tensor coefficient test passed for 10 random tests.\n";
end procedure;

function MF_TensorContraction(Ten,N,contr_modes,contr_indices)
    /*
        Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, and given contr_modes and contr_indices two sequences of same length L, returns the contraction space of Ten obtained by contracting the modes contr_modes with the indices contr_indices corresponding.
        The output is an element of VectorSpace(F, nn1*...*nn(k-L)) where nn obtained by removing from N the entries in positions contr_modes. 
        The output Res is such that the coefficient of Res is the coefficient of Ten where we append contr_indices to the tuple of indices corresponding to the coefficient of Res.
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if #contr_modes ne #contr_indices then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "contr_modes and contr_indices must have the same length.");
    end if;
    if exists(j){ j : j in [1..#contr_modes] | contr_modes[j] lt 1 or contr_modes[j] gt #N } then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "All entries of contr_modes must be in the range [1..k] where N = [n1,...,nk].");
    end if;
    if exists(j){ j : j  in [1..#contr_indices] | contr_indices[j] lt 1 or contr_indices[j] gt N[contr_modes[j]] } then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "All entries of contr_indices must be in the range [1..nj] where N = [n1,...,nk] and contr_modes[j] = n. ");
    end if;
    if #contr_modes ne #Seqset(contr_modes) then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "contr_modes must not contain duplicates.");
    end if;
    if Degree(Parent(Ten)) ne &*N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "The degree of the ambient space of Ten must be n1*...*nk where N = [n1,...,nk].");
    end if;
    if #contr_modes gt #N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "The number of modes to contract must be less or equal to the length of N.");
    end if;
    if contr_modes ne Sort(contr_modes) then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContraction", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices> >, "contr_modes must be sorted in increasing order.");
    end if;

    if #contr_modes eq 0 then
        return Ten;
    end if;
    if #contr_modes eq #N then
        return MF_TensorCoefficient(Ten, N, contr_indices);
    end if;

    nn := [N[i] : i in [1..#N] | not i in contr_modes];
    Res := VectorSpace(BaseRing(Parent(Ten)), &*nn)!0;
    for i_int in [1..&*nn] do
        i := MF_IntToIntTuple(nn, i_int - 1);
        i_for_Ten := [];
        jnn := 1;
        jcontr := 1;
        for j in [1..#N] do
            if j in contr_modes then
                i_for_Ten[j] := contr_indices[jcontr];
                jcontr +:= 1;
            else
                i_for_Ten[j] := i[jnn];
                jnn +:= 1;
            end if;
        end for;
        Res[i_int] := MF_TensorCoefficient(Ten, N, i_for_Ten);
    end for;
    return Res;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorContraction(Ten,N,contr_modes,contr_indices)", "Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, and given contr_modes and contr_indices two sequences of same length L, returns the contraction space of Ten obtained by contracting the modes contr_modes with the indices contr_indices corresponding. The output is an element of VectorSpace(F, nn1*...*nn(k-L)) where nn obtained by removing from N the entries in positions contr_modes. The output Res is such that the coefficient of Res is the coefficient of Ten where we append contr_indices to the tuple of indices corresponding to the coefficient of Res.", "LINALG", "\tTen: tensor in VectorSpace(F,n1*...*nk),\n\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\tcontr_modes: sequence of integers,\n\tcontr_indices: sequence of integers.">
);
procedure MF_TestTensorContraction()
    for test in [1..100] do
        k := Random([3..5]);
        N := [Random([2..5]) : i in [1..k]];
        F := GF(2);
        //Take a random elementary tensor, contract it and check if we get the right other elementary tensor.
        i := [Random([1..N[j]]) : j in [1..k]];
        elem := MF_ElementaryTensor(N, F, i);
        nb_contr := Random([1..k-2]);
        contr_modes := Sort(Setseq(Random(Subsets({1..k}, nb_contr))));
        contr_indices := [Random([1..N[contr_modes[j]]]) : j in [1..nb_contr]];
        contr := MF_TensorContraction(elem, N, contr_modes, contr_indices);
        //So, contr should either be zero or an elementary tensor. It is zero if and only if the contraction has been done along a mode j where ij is not equal to the corresponding contr_index. Otherwise, contr should be the elementary tensor corresponding to the indices obtained from i by removing the entries in positions contr_modes and appending contr_indices at the end.
        if exists(j){ j : j in [1..#contr_modes] | i[contr_modes[j]] ne contr_indices[j] } then
            //supposed to get zero
            if contr ne Parent(contr)!0 then
                error MAGMAFUNCTIONS_FormatError("MF_TestTensorContraction", < <"N", N>, <"i", i>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices>, <"contr", contr> >, "Tensor contraction test failed: the contraction was supposed to be zero but it was not.");
            end if;
        else
            //supposed to be the elementary tensor corresponding to the indices obtained from i by removing the entries in positions contr_modes and appending contr_indices at the end.
            nn := [N[t] : t in [1..#N] | not t in contr_modes];
            i_expected := [];
            for j in [1..#N] do
                if j in contr_modes then
                    continue;
                else
                    Append(~i_expected, i[j]);
                end if;
            end for;
            elem_expected := MF_ElementaryTensor(nn, F, i_expected);
            if contr ne elem_expected then
                error MAGMAFUNCTIONS_FormatError("MF_TestTensorContraction", < <"N", N>, <"i", i>, <"contr_modes", contr_modes>, <"contr_indices", contr_indices>, <"contr", contr>, <"elem_expected", elem_expected> >, "Tensor contraction test failed: the contraction was supposed to be the elementary tensor corresponding to the indices obtained from i by removing the entries in positions contr_modes and appending contr_indices at the end, but it was not.");
            end if;
        end if;
    end for;
    printf "Tensor contraction test passed for 100 random tests.\n";
end procedure;

function MF_TensorContractionSpace(Ten,N,contr_modes)
    /*
        Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, and given contr_modes a sequence of integers, returns the contraction space of Ten obtained by contracting the modes contr_modes with all possible indices.
        The output is a vector subspace of VectorSpace(F, nn1*...*nn(k-L)) where nn obtained by removing from N the entries in positions contr_modes.
        The output is the span of all possible contractions along the modes contr_modes with all possible indices.
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContractionSpace", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if exists(j){ j : j in [1..#contr_modes] | contr_modes[j] lt 1 or contr_modes[j] gt #N } then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContractionSpace", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes> >, "All entries of contr_modes must be in the range [1..k] where N = [n1,...,nk].");
    end if;
    if #contr_modes ne #Seqset(contr_modes) then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContractionSpace", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes> >, "contr_modes must not contain duplicates.");
    end if;
    if Degree(Parent(Ten)) ne &*N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContractionSpace", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes> >, "The degree of the ambient space of Ten must be n1*...*nk where N = [n1,...,nk].");
    end if;
    if #contr_modes gt #N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContractionSpace", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes> >, "The number of modes to contract must be less than or equal to the length of N.");
    end if;
    if contr_modes ne Sort(contr_modes) then
        error MAGMAFUNCTIONS_FormatError("MF_TensorContractionSpace", < <"N", N>, <"Ten", Ten>, <"contr_modes", contr_modes> >, "contr_modes must be sorted in increasing order.");
    end if;

    if #contr_modes eq 0 then
        return sub<Parent(Ten) | Ten>;
    end if;
    if #contr_modes eq #N then
        F := Parent(Ten[1]);
        if Ten eq Parent(Ten)!0 then
            return sub<VectorSpace(F, 1) | {} >;
        else
            return sub<VectorSpace(F, 1) | VectorSpace(F, 1).1 >;
        end if;
    end if;

    nn := [N[i] : i in [1..#N] | not i in contr_modes];
    ContrSpace := VectorSpace(BaseRing(Parent(Ten)), &*nn);
    ContrSpaceGenFam := [];
    nn_contr := [N[i] : i in contr_modes];
    for i_int in [0..&*nn_contr -1] do
        contr_indices := MF_IntToIntTuple(nn_contr, i_int);
        Append(~ContrSpaceGenFam, MF_TensorContraction(Ten, N, contr_modes, contr_indices));
    end for;
    return sub<ContrSpace | ContrSpaceGenFam>;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorContractionSpace(Ten,N,contr_modes)", "Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, and given contr_modes a sequence of integers, returns the contraction space of Ten obtained by contracting the modes contr_modes with all possible indices. The output is a vector subspace of VectorSpace(F, nn1*...*nn(k-L)) where nn obtained by removing from N the entries in positions contr_modes. The output is the span of all possible contractions along the modes contr_modes with all possible indices.", "LINALG", "\tTen: tensor in VectorSpace(F,n1*...*nk),\n\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\tcontr_modes: sequence of integers.">
);

function MF_TensorSliceSpaceDimension(Ten,N,j)
    /*
        Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, and given j an integer in [1..k], returns the dimension of the slice space of Ten along the mode j, i.e. the vector subspace of VectorSpace(F, n1*...*n(j-1)*n(j+1)*...*nk) spanned by all possible slices of Ten along the mode j.
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_TensorSliceSpaceDimension", < <"N", N>, <"Ten", Ten>, <"j", j> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if Degree(Parent(Ten)) ne &*N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorSliceSpaceDimension", < <"N", N>, <"Ten", Ten>, <"j", j> >, "The degree of the ambient space of Ten must be n1*...*nk where N = [n1,...,nk].");
    end if;
    if j lt 1 or j gt #N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorSliceSpaceDimension", < <"N", N>, <"Ten", Ten>, <"j", j> >, "j must be in the range [1..k] where N = [n1,...,nk].");
    end if;
    return Dimension(MF_TensorContractionSpace(Ten, N, [j]));
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorSliceSpaceDimension(Ten,N,j)", "Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, and given j an integer in [1..k], returns the dimension of the slice space of Ten along the mode j, i.e. the vector subspace of VectorSpace(F, n1*...*n(j-1)*n(j+1)*...*nk) spanned by all possible slices of Ten along the mode j.", "LINALG", "\tTen: tensor in VectorSpace(F,n1*...*nk),\n\tN = [n1,...,nk]: non-empty sequence of positive integers,\n\tj: integer.">
);
procedure MF_TestTensorSliceSpaceDimension()
    for test in [1..100] do
        k := Random([1..5]);
        N := [Random([2..5]) : i in [1..k]];
        q := Random([2,3,4]);
        F := GF(q);
        T := <Random(VectorSpace(F, N[j])) : j in [1..k]>;
        Ten := MF_TensorProduct(T);
        if exists(j){ j : j in [1..k] | T[j] eq Parent(T[j])!0 } then
            if MF_TensorSliceSpaceDimension(Ten, N, j) ne 0 then
                error MAGMAFUNCTIONS_FormatError("MF_TestTensorSliceSpaceDimension", < <"N", N>, <"Ten", Ten>, <"j", j> >, "Tensor slice space dimension test failed: the dimension of the slice space along some mode was not 0 for a tensor with a zero slice along this mode.");
            end if;
        elif exists(j){ j : j in [1..k] | MF_TensorSliceSpaceDimension(Ten, N, j) ne 1} then
            error MAGMAFUNCTIONS_FormatError("MF_TestTensorSliceSpaceDimension", < <"N", N>, <"Ten", Ten> >, "Tensor slice space dimension test failed: the dimension of the slice space along some mode was not 1 for a rank-1 tensor.");
        end if;
    end for;
    printf "Tensor slice space dimension test passed for 100 random tests.\n";
end procedure;

function MF_TensorIsTensorRankOne(Ten,N)
    /*
        Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, returns true if Ten is of tensor rank 1, i.e. if it can be written as the tensor product of k vectors, and false otherwise.
        The function checks if Ten has every slice space of dimension exactly one.
    */
    if #N eq 0 then
        error MAGMAFUNCTIONS_FormatError("MF_TensorIsTensorRankOne", < <"N", N>, <"Ten", Ten> >, "N must be a non-empty sequence of positive integers.");
    end if;
    if Degree(Parent(Ten)) ne &*N then
        error MAGMAFUNCTIONS_FormatError("MF_TensorIsTensorRankOne", < <"N", N>, <"Ten", Ten> >, "The degree of the ambient space of Ten must be n1*...*nk where N = [n1,...,nk].");
    end if;
    for j in [1..#N] do
        if MF_TensorSliceSpaceDimension(Ten, N, j) ne 1 then
            return false;
        end if;
    end for;
    return true;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_TensorIsTensorRankOne(Ten,N)", "Given Ten a tensor in VectorSpace(F,n1*...*nk) for some field F, given N = [n1,...,nk] a sequence of positive integers, returns true if Ten is of tensor rank 1, i.e. if it can be written as the tensor product of k vectors, and false otherwise. The function checks if Ten has every slice space of dimension exactly one.", "LINALG", "\tTen: tensor in VectorSpace(F,n1*...*nk),\n\tN = [n1,...,nk]: non-empty sequence of positive integers.">
);






//=============================================================
//          DISPLAY, WRITE, STRING HANDLING
//=============================================================
function MF_RealToString(x)
    /*
        Converts a real number to a string with 30 digits after the comma.
    */
    Ent := IntegerToString(Floor(x));
    NbDigInt := #Ent;
    Str := IntegerToString(Floor(x * 10^30));
    if Floor(x) eq 0 then
        mm := 1;
    else
        mm := NbDigInt + 1;
    end if;
    Str := Ent cat "." cat Substring(Str, mm, 30 - NbDigInt);
    return Str;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_RealToString(x)", "Converts a real number to a string with 30 digits after the comma.", "DISP", "\tx: real number.">
);

function MF_IntMatrixToCSV(M)
    L := Nrows(M);
    C := Ncols(M);
    if L * C eq 0 then
        return "";
    end if;
    Str := "";
    for i in [1..(L - 1)] do
        for j in [1..(C - 1)] do
            Str := Str cat IntegerToString(M[i, j]) cat ",";
        end for;
        Str := Str cat IntegerToString(M[i, C]) cat "\n";
    end for;
    for j in [1..(C - 1)] do
        Str := Str cat IntegerToString(M[L, j]) cat ",";
    end for;
    Str := Str cat IntegerToString(M[L, C]);
    return Str;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_IntMatrixToCSV(M)", "Converts an integer matrix to a CSV string.", "DISP", "\tM: integer matrix.">
);

function MF_FieldToString(x, Q)
    Field := Parent(x);
    if x in PrimeField(Field) then
        for k in [0..(#PrimeField(Field) - 1)] do
            if x eq Field!k then
                return "GF(" cat IntegerToString(Q) cat ")!" cat IntegerToString(k);
            end if;
        end for;
    else
        for m in [0..(Q - 1)] do
            if x eq Field.1^m then
                return "GF(" cat IntegerToString(Q) cat ").1^" cat IntegerToString(m);
            end if;
        end for;
    end if;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_FieldToString(x,Q)", "Converts an element x of GF(Q) to a string of the form GF(Q)!k if x is in the prime field, or GF(Q).1^m if x is a power m of the generator of the field.", "DISP", "\tx: element of GF(Q),\n\tQ: prime power.">
);

function MF_AUX_FieldToStringSeq(x, Q)
    Field := Parent(x);
    if x in PrimeField(Field) then
        for k in [0..(#PrimeField(Field) - 1)] do
            if x eq Field!k then
                return IntegerToString(k);
            end if;
        end for;
    else
        for m in [0..(Q - 1)] do
            if x eq Field.1^m then
                return "GF(" cat IntegerToString(Q) cat ").1^" cat IntegerToString(m);
            end if;
        end for;
    end if;
end function;
function MF_SequenceFieldToString(Seq, Q)
    ElemAsString := &cat [MF_AUX_FieldToStringSeq(s, Q) cat "," : s in Seq];
    return "[GF(" cat IntegerToString(Q) cat ")! " cat Substring(ElemAsString, 1, #ElemAsString - 1) cat "]";
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_SequenceFieldToString(Seq,Q)", "Converts a sequence Seq of elements of GF(Q) to a string of the form [GF(Q)! k1, GF(Q)! k2, ...] where each ki is either an integer if the element is in the prime field, or a string of the form GF(Q).1^m if the element is a power m of the generator of the field.", "DISP", "\tSeq: sequence of elements of GF(Q),\n\tQ: prime power.">
);

function MF_MatrixFieldToString(M, Q)
    nr := Nrows(M);
    nc := Ncols(M);
    SeqM := &cat [[M[i, j] : j in [1..nc]] : i in [1..nr]];
    return "Matrix(" cat IntegerToString(nr) cat "," cat IntegerToString(nc) cat "," cat MF_SequenceFieldToString(SeqM, Q) cat ")";
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_MatrixFieldToString(M,Q)", "Converts a matrix M with entries in GF(Q) to a string of the form Matrix(nr,nc,[GF(Q)! k1, GF(Q)! k2, ...]) where nr and nc are the number of rows and columns of M, and each ki is either an integer if the entry is in the prime field, or a string of the form GF(Q).1^m if the entry is a power m of the generator of the field.", "DISP", "\tM: matrix with entries in GF(Q),\n\tQ: prime power.">
);


function MF_RealMatrixToCSV(M)
    L := Nrows(M);
    C := Ncols(M);
    if L * C eq 0 then
        return "";
    end if;
    Str := "";
    for i in [1..(L - 1)] do
        for j in [1..(C - 1)] do
            Str := Str cat RealToString(M[i, j]) cat ",";
        end for;
        Str := Str cat RealToString(M[i, C]) cat "\n";
    end for;
    for j in [1..(C - 1)] do
        Str := Str cat RealToString(M[L, j]) cat ",";
    end for;
    Str := Str cat RealToString(M[L, C]);
    return Str;
end function;
MAGMAFUNCTIONS_registery := MAGMAFUNCTIONS_add(
    MAGMAFUNCTIONS_registery,
    <"MF_RealMatrixToCSV(M)", "Converts a real matrix to a CSV string with 30 digits after the comma for each entry.", "DISP", "\tM: real matrix.">
);






//===========================================================
//          MATHFUNCTIONS DESCRIPTION
//===========================================================
procedure MAGMAFUNCTIONS_help()
    printf "\n==========================\n";
    nb := 1;
    for type in Set(MAGMAFUNCTIONS_registery[3]) do
        for i in [1..#MAGMAFUNCTIONS_registery[1]] do
            if MAGMAFUNCTIONS_registery[3][i] eq type then
                printf " %2o %o\t%o\n", nb, MAGMAFUNCTIONS_registery[3][i], MAGMAFUNCTIONS_registery[1][i];
                nb +:= 1;
            end if;
        end for;
    end for;
    printf "\nUse MAGMAFUNCTIONS_info(x) for more details on a specific function, where x is either the function number or the function name. Use MAGMAFUNCTIONS_help_type(type) for more details on all functions of a specific type, where type is either ARITH, CMBTRCS, BIJCTNS, TEST.\n\n";
end procedure;
procedure MAGMAFUNCTIONS_info(x)
    if Type(x) eq RngIntElt then
        // If the input is an integer, its x-th function when ordered by type.
        nb := 1;
        for type in Set(MAGMAFUNCTIONS_registery[3]) do
            for i in [1..#MAGMAFUNCTIONS_registery[1]] do
                if MAGMAFUNCTIONS_registery[3][i] eq type then
                    if nb eq x then
                        printf "\n";
                        printf "[%o] %o\n%o\n%o\n\n",
                            MAGMAFUNCTIONS_registery[3][i],
                            MAGMAFUNCTIONS_registery[1][i],
                            MAGMAFUNCTIONS_registery[4][i],
                            MAGMAFUNCTIONS_registery[2][i];
                        return;
                    end if;
                    nb +:= 1;
                end if;
            end for;
        end for;
        error MAGMAFUNCTIONS_FormatError("MAGMAFUNCTIONS_info", < <"x", x> >, "No function number was found in MAGMAFUNCTIONS.");
    elif Type(x) eq MonStgElt then
        // If the input is a string, the function with that name.
        for i in [1..#MAGMAFUNCTIONS_registery[1]] do
            if MAGMAFUNCTIONS_registery[1][i] eq x then
                printf "\n";
                printf "[%o] %o\n%o\n%o\n\n",
                    MAGMAFUNCTIONS_registery[3][i],
                    MAGMAFUNCTIONS_registery[1][i],
                    MAGMAFUNCTIONS_registery[4][i],
                    MAGMAFUNCTIONS_registery[2][i];
                return;
            end if;
        end for;
        error MAGMAFUNCTIONS_FormatError("MAGMAFUNCTIONS_info", < <"x", x> >, "No function named was found in MAGMAFUNCTIONS.");
    else
        error MAGMAFUNCTIONS_FormatError("MAGMAFUNCTIONS_info", < <"x", x> >, "Input must be either an integer or a string.");
    end if;
end procedure;
procedure MAGMAFUNCTIONS_help_type(thetype)
    printf "\n==========================\n";
    nb := 1;
    for type in Set(MAGMAFUNCTIONS_registery[3]) do
        for i in [1..#MAGMAFUNCTIONS_registery[1]] do
            if MAGMAFUNCTIONS_registery[3][i] eq type then
                if type eq thetype then
                    printf " %2o %o\t%o\n", nb, MAGMAFUNCTIONS_registery[3][i], MAGMAFUNCTIONS_registery[1][i];
                end if;
                nb +:= 1;
            end if;
        end for;
    end for;
    printf "\nUse MAGMAFUNCTIONS_info(x) for more details on a specific function, where x is either the function number or the function name. Use MAGMAFUNCTIONS_help_type(type) for more details on all functions of a specific type, where type is either ARITH, CMBTRCS, BIJCTNS, TEST.\n\n";
end procedure;
procedure MAGMAFUNCTIONS_all_info()
    for i in [1..#MAGMAFUNCTIONS_registery[1]] do
        MAGMAFUNCTIONS_info(i);
    end for;
end procedure;
