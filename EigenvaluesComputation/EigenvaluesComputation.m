load "MAGMAFUNCTIONS.m";

procedure InitBruteForceRankDistributionMat()
    FILE := Open("DATASET_SUMMARY_KnownRankDistributions.txt", "w");
    FILE := Open("DATASET_KnownRankDistributions.txt", "w");
end procedure;

procedure InitBruteForceRank01DistributionTensors()
    FILE := Open("DATASET_SUMMARY_KnownRank01Distributions.txt", "w");
    FILE := Open("DATASET_KnownRank01Distributions.txt", "w");
end procedure;

function BruteForceRankDistributionMat(q,N1,N2)
    printf "Enumerating the rank distribution of subspaces of F_%o^%o x F_%o^%o by brute force...\n", q, N1, q, N2;
    nb_subspaces := Floor(&+ [MF_GBC(q,N1*N2, dim) : dim in [0..N1*N2]]);
    printf "Number of subspaces: %o\n", nb_subspaces;
    printf "Number of matrices: %o\n", &+ [q^dim * MF_GBC(q,N1*N2, dim) : dim in [1..N1*N2]];
    /*
        Given q a prime power, and N1, N2 a tuple of positive integers, 
        returns a tuple of tuples counting the number of subspaces of Fq^N1 x Fq^N2 of each rank distribution.

        The output is of the form 
           < .... , < <dim , RankDistribution , COUNT>, ... >
        where RankDistribution is a sequence of length Min(N1,N2)+1 of the form 
           [1,#of rank 1, #of rank 2, ... , #of rank Min(N1,N2)]>
        where dim = sum(RankDistribution) is the dimension a subspace of Fq^N1 x  Fq^N2 with that rank distribution would have, and COUNT is the number of subspaces of Fq^N1 x  Fq^N2 with that rank distribution.
    */

    AmbSpaces := KMatrixSpace(GF(q),N1*N2,N1*N2);
    Amb := KMatrixSpace(GF(q), N1, N2);

    res := < 
        < 0 , [Max(0,1-r) : r in [0..Min(N1,N2)]] , 1 >
    >;

    nb_done := 1;
    step := Max(1, nb_subspaces div 25);


    //Let us enumerate the RREF of AmbSpaces. That counts exactly once and only once all subspaces of Fq^N1 x Fq^N2 (rowspaces). 

    //Choice of the dimension of the subspace
    for dim in [1..N1*N2] do
        
        //Choice of the pivot columns 
        for pivotcols in Subsets({1..N1*N2}, dim) do
            pivotcols_seq := Sort(Setseq(pivotcols));

            //Choice of the pivot rows
            for pivotrows in Subsets({1..dim}, dim) do
                pivotrows_seq := Sort(Setseq(pivotrows));

                /*Choice of the entries in the non-pivot columns
                    For each pivot (i,j) before (i',j') there are ix(j'-j-1) entries to fill in the non-pivot columns between the pivot columns j and j'. There are also dim x (N1*N2 - j_last) after the last pivot column j_last. 
                */
                nb_entries_to_fill := 0;
                for t in [1..dim-1] do
                    nb_entries_to_fill +:= pivotrows_seq[t]*(pivotcols_seq[t+1] - pivotcols_seq[t] - 1);
                end for;
                nb_entries_to_fill +:= dim*(N1*N2 - pivotcols_seq[dim]);

                for entries in VectorSpace(GF(q), nb_entries_to_fill) do
                    nb_done +:= 1;
                    if nb_done mod step eq 0 then
                        printf "%o%%|", Round(100*nb_done/nb_subspaces);
                    end if;

                    //  =FILL THE MATRIX=

                    //Put the 1s in the pivot positions
                    M := ZeroMatrix(GF(q),dim, N1*N2);
                    for i in [1..dim] do
                        M[pivotrows_seq[i], pivotcols_seq[i]] := 1;
                    end for;

                    //Fill the non-pivot entries
                    entry_index := 1;
                    for t in [1..dim-1] do
                        //Fill the entries between pivotcols_seq[t] and pivotcols_seq[t+1]
                        for i in [1..pivotrows_seq[t]], j in [pivotcols_seq[t]+1..pivotcols_seq[t+1]-1] do
                            M[i, j] := entries[entry_index];
                            entry_index +:= 1;
                        end for;
                    end for;
                    //Fill the entries after the last pivot column
                    for i in [1..dim], j in [pivotcols_seq[dim]+1..N1*N2] do
                        M[i, j] := entries[entry_index];
                        entry_index +:= 1;
                    end for;


                    //  ===== THE RREF IS NOW DONE =====
                    // We can now compute the rank distribution of the row space of M, which is a subspace of Fq^N1 x Fq^N2.
                    TheVectorSpace := sub< Amb | [Amb!Eltseq(M[i]) : i in [1..dim]] >;
                    TheRankDistribution := [0 : r in [0..Min(N1,N2)]];
                    for M in TheVectorSpace do
                        TheRankDistribution[Rank(M)+1] +:= 1;
                    end for;

                    if exists(i){i : i in [1..#res] | res[i][2] eq TheRankDistribution} then
                        res[i][3] +:= 1;
                    else
                        res := res cat < <dim, TheRankDistribution, 1 > >;
                    end if;
                end for;                
            end for;
        end for;      
    end for;

    return res;
end function;

procedure TestBruteForceRankDistributionMat(q,N1,N2)
    res := BruteForceRankDistributionMat(q,N1,N2);
    for i in [1..#res] do
        printf "Found %o subspaces of F_%o^%o x F_%o^%o with rank distribution:\n%o\n\n", res[i][3], q, N1, q, N2, res[i][2];
    end for;
    //Test1: rank distribution makes sense with dimension in every item
    if exists(i){i : i in [1..#res] | q^res[i][1] ne &+res[i][2]} then
        error "Error: the dimension does not match the rank distribution in some item.";
    end if;
    //Test2: the total number of subspaces is correct
    if &+ [res[i][3] : i in [1..#res]] ne &+ [MF_GBC(q,N1*N2, dim) : dim in [0..N1*N2]] then
        error "Error: the total number of subspaces is not correct.";
    end if;
    printf "All tests passed for q=%o, N1=%o, N2=%o.\n", q, N1, N2;
end procedure;

procedure PrintBruteForceRankDistributionMat(q,N1,N2)
    res := BruteForceRankDistributionMat(q,N1,N2);
    printf "<\n";
    for x in res do
        x;
    end for;
    printf ">;\n";
end procedure;

function BruteForceRankOneDistribution(q,N)
    /*
        Given q a prime power, and N = [N1,..,Nk] a tuple of positive integers, 
        returns a tuple of tuples counting the number of subspaces of the space of tensors of size N1 x ... x Nk over Fq of each number of rank ones.

        The output is of the form 
           < .... , < <dim , nbrankones , COUNT>, ... >
        where nbrankones is an integer in [0..|Seg(Fq^N1 x ... x Fq^Nk)|] counting the number of rank one tensors in the subspace, where dim = sum(RankDistribution) is the dimension a subspace of Fq^(N1 x ... x Nk), and COUNT is the number of subspaces of Fq^(N1 x ... x Nk) with that number of rank one tensors and this dimension.
    */

    printf "Enumerating the number of rank ones of subspaces of F_%o^(%o) by brute force...\n", q, Sprint(N);
    nb_subspaces := Floor(&+ [MF_GBC(q,&*N, dim) : dim in [0..&*N]]);
    printf "Number of subspaces: %o\n", nb_subspaces;
    printf "Number of tensors: %o\n", &+ [q^dim * MF_GBC(q,&*N, dim) : dim in [1..&*N]];

    Nsize := &*N;
    AmbSpaces := KMatrixSpace(GF(q),Nsize,Nsize);
    Amb := VectorSpace(GF(q), Nsize);
    res := < 
        < 0 , 0 , 1 >
    >;

    nb_done := 1;
    step := Max(1, nb_subspaces div 100);

    //Choice of the dimension of the subspace
    for dim in [1..Nsize] do
        //Choice of the pivot columns 
        for pivotcols in Subsets({1..Nsize}, dim) do
            pivotcols_seq := Sort(Setseq(pivotcols));

            //Choice of the pivot rows
            for pivotrows in Subsets({1..dim}, dim) do
                pivotrows_seq := Sort(Setseq(pivotrows));

                /*Choice of the entries in the non-pivot columns
                    For each pivot (i,j) before (i',j') there are ix(j'-j-1) entries to fill in the non-pivot columns between the pivot columns j and j'. There are also dim x (Nsize - j_last) after the last pivot column j_last. 
                */
                nb_entries_to_fill := 0;
                for t in [1..dim-1] do
                    nb_entries_to_fill +:= pivotrows_seq[t]*(pivotcols_seq[t+1] - pivotcols_seq[t] - 1);
                end for;
                nb_entries_to_fill +:= dim*(Nsize - pivotcols_seq[dim]);

                for entries in VectorSpace(GF(q), nb_entries_to_fill) do
                    nb_done +:= 1;
                    if nb_done mod step eq 0 then
                        printf "%o%%|", Round(100*nb_done/nb_subspaces);
                    end if;

                    //  =FILL THE MATRIX=

                    //Put the 1s in the pivot positions
                    M := ZeroMatrix(GF(q),dim, Nsize);
                    for i in [1..dim] do
                        M[pivotrows_seq[i], pivotcols_seq[i]] := 1;
                    end for;

                    //Fill the non-pivot entries
                    entry_index := 1;
                    for t in [1..dim-1] do
                        //Fill the entries between pivotcols_seq[t] and pivotcols_seq[t+1]
                        for i in [1..pivotrows_seq[t]], j in [pivotcols_seq[t]+1..pivotcols_seq[t+1]-1] do
                            M[i, j] := entries[entry_index];
                            entry_index +:= 1;
                        end for;
                    end for;
                    //Fill the entries after the last pivot column
                    for i in [1..dim], j in [pivotcols_seq[dim]+1..Nsize] do
                        M[i, j] := entries[entry_index];
                        entry_index +:= 1;
                    end for;

                    //  ===== THE RREF IS NOW DONE =====
                    // We can now compute the number of rank ones tensors in the row space of M seen as a subspace of the space of tensors.

                    TheVectorSpace := sub< Amb | [Amb!Eltseq(M[i]) : i in [1..dim]] >;
                    nbrankones := 0;
                    for T in TheVectorSpace do
                        if MF_TensorIsTensorRankOne(T, N) then
                            nbrankones +:= 1;
                        end if;
                    end for;

                    if exists(i){i : i in [1..#res] | res[i][2] eq nbrankones and res[i][1] eq dim} then
                        res[i][3] +:= 1;
                    else
                        res := res cat < <dim, nbrankones, 1 > >;
                    end if;
                end for;
            end for;
        end for;
    end for;
    return res;
end function;

procedure StoreBruteForceRankDistributionMat(q,N1,N2)
    /*
        Each line of the dataset summary corresponds to a unique line to execute in the dataset file.

        A line in the dataset summary has the form 
           <q,[N1,N2]>
        A line in the dataset file has the form 
           return < ... , < <dim, RankDistribution , COUNT>, ... >;
        
        No line shoud be skiped in one entry.
    */

    /*
        Only stores N1 >= N2 > 1 !!
    */
    if N1 lt N2 or N2 le 1 then
        print "Only stores N1 >= N2 > 1 !!";
    end if;

    SummaryItem := "<" cat IntegerToString(q) cat ",[" cat IntegerToString(N1) cat "," cat IntegerToString(N2) cat "]>";
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRankDistributions.txt", "r");

    SUMMARY_str := Read(SUMMARY_FILE);
    if IsEof(SUMMARY_str) then
        SUMMARY_lines := [];
    else
        SUMMARY_lines := Split(SUMMARY_str, "\n");
    end if;

    if exists(i){i : i in [1..#SUMMARY_lines] | SUMMARY_lines[i] eq SummaryItem} then
        print "This item is already in the summary file.";
    else
        //Lets us add the item to the summary file and the data to the dataset file.
        res := BruteForceRankDistributionMat(q,N1,N2);
        SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRankDistributions.txt", "a");
        DATASET_FILE := Open("DATASET_KnownRankDistributions.txt", "a");

        //Construct the string to add to the dataset file.
        DatasetItem := "return <";

        for idx_x in [1..#res] do
            x := res[idx_x];
            //Contructing the string for the item x = <dim, RankDistribution, COUNT>
            x_string := "<" cat IntegerToString(x[1]) cat " , [" ;
            seq := x[2];
            seq_string := &cat [ IntegerToString(r) cat "," : r in seq[1..#seq-1] ];
            seq_string := seq_string cat IntegerToString(seq[#seq]);
            x_string := x_string cat seq_string cat "] , " cat IntegerToString(x[3]) cat ">";
            if idx_x ne #res then
                DatasetItem := DatasetItem cat x_string cat ",";
            else
                DatasetItem := DatasetItem cat x_string cat ">;";
            end if;
        end for;
        //Add the item to the files.
        Write(SUMMARY_FILE, SummaryItem cat "\n");
        Write(DATASET_FILE, DatasetItem cat "\n");
        printf "The item %o has been added to the summary file and the dataset file.\n", SummaryItem;
    end if;
end procedure;

procedure StoreBruteForceRank01DistributionTensors(q,N)
    /*
        Each line of the dataset summary corresponds to a unique line to execute in the dataset file.

        A line in the dataset summary has the form 
           <q,[N1,...,Nk]>
        A line in the dataset file has the form 
           return < ... , < <dim, nbrankones , COUNT>, ... >;
        
        No line shoud be skiped in one entry.
    */

    /*
        Only stores N1 >= N2 >= ... >= Nk > 1 !!
    */
    if N ne Reverse(Sort(N)) or N[#N] le 1 then
        print "Only stores N1 >= N2 >= ... >= Nk > 1 !!";
    end if;

    SummaryItem := "<" cat IntegerToString(q) cat ",[" cat &cat [ IntegerToString(N[i]) cat (i ne #N select "," else "") : i in [1..#N] ] cat "]>";
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRank01Distributions.txt", "r");
    SUMMARY_str := Read(SUMMARY_FILE);
    if IsEof(SUMMARY_str) then
        SUMMARY_lines := [];
    else
        SUMMARY_lines := Split(SUMMARY_str, "\n");
    end if;
    if exists(i){i : i in [1..#SUMMARY_lines] | SUMMARY_lines[i] eq SummaryItem} then
        print "This item is already in the summary file.";
    else
        //Lets us add the item to the summary file and the data to the dataset file.
        res := BruteForceRankOneDistribution(q,N);
        SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRank01Distributions.txt", "a");
        DATASET_FILE := Open("DATASET_KnownRank01Distributions.txt", "a");

        //Construct the string to add to the dataset file.
        DatasetItem := "return <";

        for idx_x in [1..#res] do
            x := res[idx_x];
            //Contructing the string for the item x = <dim, RankDistribution, COUNT>
            x_string := "<" cat IntegerToString(x[1]) cat " , " cat IntegerToString(x[2]) cat " , " cat IntegerToString(x[3]) cat ">";
            seq := x[2];
            if idx_x ne #res then
                DatasetItem := DatasetItem cat x_string cat ",";
            else
                DatasetItem := DatasetItem cat x_string cat ">;";
            end if;
        end for;
        //Add the item to the files.
        Write(SUMMARY_FILE, SummaryItem cat "\n");
        Write(DATASET_FILE, DatasetItem cat "\n");
        printf "The item %o has been added to the summary file and the dataset file.\n", SummaryItem;
    end if;
end procedure;

function KnownRankDistributions(q,N)
    /*
        Given q a prime power, and N = [N1,...,Nm] a tuple of positive integers, 
        returns a tuple of tuples counting the number of subspaces of Fq^N1 x ... x Fq^Nm of each rank distribution.

        The output is of the form 
           < .... , < <dim,RankDistribution> , COUNT>, ... >
        where RankDistribution is a sequence of length Min(N1,...,Nm)+1 of the form 
           [1,#of rank 1, #of rank 2, ... , #of rank Min(N1,...,Nm)]>
        where dim = sum(RankDistribution) is the dimension a subspace of Fq^N1 x ... x Fq^Nm with that rank distribution would have, and COUNT is the number of subspaces of Fq^N1 x ... x Fq^Nm with that rank distribution.
    */
    if N ne Reverse(Sort(N)) then
        KnownRankDistributions(q, Reverse(Sort(N)));
    end if;
    if #N eq 0 then
        error "N must be a non-empty tuple of positive integers.";
    end if;

    
    if #N eq 1 then
        res := <>;
        for dim in [0..N[1]] do
            x := <dim, [1,q^(N[1]) -1 ] , MF_GBC(q,N[1], dim)>;
            res := res cat <x>;
        end for;
        return res;
    end if;
    
    
    SummaryItem := "<" cat IntegerToString(q) cat ",[" cat &cat [ IntegerToString(N[i]) cat (i ne #N select "," else "") : i in [1..#N] ] cat "]>";
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRankDistributions.txt", "r");
    SUMMARY_str := Read(SUMMARY_FILE);
    SUMMARY_lines := Split(SUMMARY_str, "\n");
    if exists(i){i : i in [1..#SUMMARY_lines] | SUMMARY_lines[i] eq SummaryItem} then
        DATASET_FILE := Open("DATASET_KnownRankDistributions.txt", "r");
        DATASET_str := Read(DATASET_FILE);
        DATASET_lines := Split(DATASET_str, "\n");
        res := eval DATASET_lines[i];
        return res;
    end if;

    error "KnownRankDistributions is not implemented yet for this input.";
end function;

procedure DoTheComputations()
    for q in [2,3,4,5,7,8,9,11,13,17,19,23,29,31,37,41,43,47] do
        for N2 in [2..6] do
            for N1 in [N2..6] do
                nb_matrices_to_enumerate := &+ [q^dim * MF_GBC(q,N1*N2, dim) : dim in [0..N1*N2]];
                if nb_matrices_to_enumerate le 10^11 then
                    StoreBruteForceRankDistributionMat(q,N1,N2);
                end if;
            end for;
        end for;
    end for;
end procedure;

function SpectrumFromRank01Distribution(q,N)
    /*
        Given q a prime power, and N = [N1,...,Nm] a tuple of positive integers, 
        returns the spectrum of the Cayley graph of Fq^N1 x ... x Fq^Nm generated by the set of rank one tensors.

        Assuming that we know the number of rank ones in all subspaces of Fq^N1 x ... x Fq^N(j-1) x Fq^N(j+1) x ... x Fq^Nm for a certain j, then 
            lambda_A = A q^(Nj)/(q-1)  - |Seg(Fq^N1 x ... x Fq^N(j-1) x Fq^N(j+1) x ... x Fq^Nm )|/(q-1)
        is an eigenvalue of multiplicity
            mult_a = \sum_{k = 0}^{n_j} #( M \in \F_q^{(n_j,k)}  : rank(M) = k ) #(C <= \F_q^N_minus_j  :  codim C = k , A_1(C) = a ). 
    */
    m := #N;
    if m eq 0 then
        error "N must be a non-empty tuple of positive integers.";
    end if;
    if m eq 1 then
        return < <-1,q^(N[1]) -1> , <q^(N[1]) -1, 1> >;
    end if;
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRank01Distributions.txt", "r");
    SUMMARY_str := Read(SUMMARY_FILE);
    SUMMARY_lines := Split(SUMMARY_str, "\n");
    thej := -1;
    for j in [1..m] do
        SummaryItem := "<" cat IntegerToString(q) cat ",[";
        N_minus_j := Reverse(Sort([ N[t] : t in [1..m] | t ne j ]));
        SummaryItem := SummaryItem cat &cat [ IntegerToString(N_minus_j[i]) cat "," : i in [1..#N_minus_j-1] ] cat IntegerToString(N_minus_j[#N_minus_j]) cat "]>";
        if exists(i){i : i in [1..#SUMMARY_lines] | SUMMARY_lines[i] eq SummaryItem} then
            thej := j;
            break;
        end if;
    end for;
    if thej eq -1 then
        error "The rank distribution of subspaces of Fq^N1 x ... x Fq^N(j-1) x Fq^N(j+1) x ... x Fq^Nm is not known for any j.";
    end if;

    N_minus_thej := Reverse(Sort([ N[t] : t in [1..m] | t ne thej ]));
    Seg_N_minus_thej := &*([1] cat [ (q^(N_minus_thej[i]) - 1) : i in [1..#N_minus_thej] ])/(q-1)^(#N_minus_thej-1);
    printf "Spec of F_%o^N where N=%o, from rank distribution of subspaces of F_%o^%o.\n",q, Sprint(N), q, Sprint(N_minus_thej);

    DATASET_FILE := Open("DATASET_KnownRank01Distributions.txt", "r");
    DATASET_str := Read(DATASET_FILE);
    DATASET_lines := Split(DATASET_str, "\n");
    res := eval DATASET_lines[i];
    TheEigenvalues := <>;

    for x in res do
        k_x := &*([1] cat N_minus_thej) - x[1];
        if k_x le N[thej] then
            count_x := x[3];
            A1_x := x[2]; 
            lambda_x := (A1_x * q^(N[thej]) - Seg_N_minus_thej)/(q-1);
            mult_x := count_x * MF_NbMatrixRank(N[thej], k_x, q, k_x);
            TheEigenvalues := TheEigenvalues cat < <lambda_x, mult_x> >;
        end if;
    end for;

    //Makesure that the eigenvalues are distinct and sorted by increasing value.
    EVs := Reverse(Sort(Setseq({TheEigenvalues[i][1] : i in [1..#TheEigenvalues]})));
    TheEigenvaluesBis := <>;
    for lambda in EVs do
        TheMult := &+ [ TheEigenvalues[i][2] : i in [1..#TheEigenvalues] | TheEigenvalues[i][1] eq lambda ];
        TheEigenvaluesBis := TheEigenvaluesBis cat < <Floor(lambda), Floor(TheMult)> >;
    end for;

    return TheEigenvaluesBis;
end function;

function SpectrumFromRankDistribution(q,N)
    /*
        Given q a prime power, and N = [N1,...,Nm] a tuple of positive integers, 
        returns the spectrum of the Cayley graph of Fq^N1 x ... x Fq^Nm generated by
        the set of rank one tensors.

        Assuming that we know the rank distribution of subspaces of Fq^N1 x ... x Fq^N(j-1) x Fq^N(j+1) x ... x Fq^Nm for a certain j, then 
            lambda_A = A q^(Nj)/(q-1)  - |Seg(Fq^N1 x ... x Fq^N(j-1) x Fq^N(j+1) x ... x Fq^Nm )|/(q-1)
        is an eigenvalue of multiplicity
            mult_a = \sum_{k = 0}^{n_j} #( M \in \F_q^{(n_j,k)}  : rank(M) = k ) #(C <= \F_q^N_minus_j  :  codim C = k , A_1(C) = a ). 
    */

    m := #N;
    if m eq 0 then
        error "N must be a non-empty tuple of positive integers.";
    end if;
    if m eq 1 then
        return < <-1,q^(N[1]) -1> , <q^(N[1]) -1, 1> >;
    end if;

    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRankDistributions.txt", "r");
    SUMMARY_str := Read(SUMMARY_FILE);
    SUMMARY_lines := Split(SUMMARY_str, "\n");
    
    thej := -1;
    for j in [1..m] do
        SummaryItem := "<" cat IntegerToString(q) cat ",[";
        N_minus_j := Reverse(Sort([ N[t] : t in [1..m] | t ne j ]));
        SummaryItem := SummaryItem cat &cat [ IntegerToString(N_minus_j[i]) cat "," : i in [1..#N_minus_j-1] ] cat IntegerToString(N_minus_j[#N_minus_j]) cat "]>";
        if exists(i){i : i in [1..#SUMMARY_lines] | SUMMARY_lines[i] eq SummaryItem} then
            thej := j;
            break;
        end if;
    end for;
    if thej eq -1 then
        error "The rank distribution of subspaces of Fq^N1 x ... x Fq^N(j-1) x Fq^N(j+1) x ... x Fq^Nm is not known for any j.";
    end if;

    N_minus_thej := Reverse(Sort([ N[t] : t in [1..m] | t ne thej ]));
    Seg_N_minus_thej := &*([1] cat [ (q^(N_minus_thej[i]) - 1) : i in [1..#N_minus_thej] ])/(q-1)^(#N_minus_thej-1);
    printf "Spec of F_%o^N where N=%o, from rank distribution of subspaces of F_%o^%o.\n",q, "[" cat &cat([" "] cat [ IntegerToString(N[i]) cat " " : i in [1..m] ]) cat "]", q , "[" cat &cat([" "] cat [ IntegerToString(N_minus_thej[i]) cat " " : i in [1..#N_minus_thej] ]) cat "]";

    DATASET_FILE := Open("DATASET_KnownRankDistributions.txt", "r");
    DATASET_str := Read(DATASET_FILE);
    DATASET_lines := Split(DATASET_str, "\n");
    res := eval DATASET_lines[i];
    TheEigenvalues := <>;

    for x in res do
        k_x := &*([1] cat N_minus_thej) - x[1];
        if k_x le N[thej] then
            count_x := x[3];
            A1_x := x[2][2]; 
            lambda_x := (A1_x * q^(N[thej]) - Seg_N_minus_thej)/(q-1);
            mult_x := count_x * MF_NbMatrixRank(N[thej], k_x, q, k_x);
            TheEigenvalues := TheEigenvalues cat < <lambda_x, mult_x> >;
            printf "Found eigenvalue %o with multiplicity %o corresponding to subspaces of F_%o^%o with dimension %o and %o rank one tensors.\n", lambda_x, mult_x, q, Sprint(N_minus_thej), x[1], A1_x;
            printf "   a = %o \n   q^nj = %o \n   Seg = %o\n", A1_x, q^(N[thej]), Seg_N_minus_thej;
        end if;
    end for;

    //Makesure that the eigenvalues are distinct and sorted by increasing value.
    EVs := Reverse(Sort(Setseq({TheEigenvalues[i][1] : i in [1..#TheEigenvalues]})));
    TheEigenvaluesBis := <>;
    for lambda in EVs do
        TheMult := &+ [ TheEigenvalues[i][2] : i in [1..#TheEigenvalues] | TheEigenvalues[i][1] eq lambda ];
        TheEigenvaluesBis := TheEigenvaluesBis cat < <Floor(lambda), Floor(TheMult)> >;
    end for;

    return TheEigenvaluesBis;
end function;

procedure EigenvaluesStorageInit()
    /*
        Storage mode: 
            each line of the dataset summary corresponds to a unique line of data in both the eigenvalue dataset and the multiplicity dataset.
            namely, the i-th line of the dataset summary is of the form
                q, N1,...,Nm
            and the i-th line of the eigenvalue dataset like the  multiplicitites dataset are both of same length nb_distinct_eigenvalues and of the form
                lambda_1, lambda_2, ... , lambda_nb_distinct_eigenvalues
                mult_1, mult_2, ... , mult_nb_distinct_eigenvalues
            where lambda_i are sorted in strictly decreasing order.
    */
    FILE := Open("DATASET_SUMMARY_KnownSpectra.csv", "w");
    FILE := Open("DATASET_EIGENVALUES_KnownSpectra.csv", "w");
    FILE := Open("DATASET_MULTIPLICITIES_KnownSpectra.csv", "w");
end procedure;

procedure StoreSpectrum(q,N)
    if N ne Reverse(Sort(N)) then
        StoreSpectrum(q, Reverse(Sort(N)));
        return;
    end if;
    
    SummaryItem := IntegerToString(q) cat "," cat &cat [ IntegerToString(N[i]) cat "," : i in [1..#N-1] ] cat IntegerToString(N[#N]);
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownSpectra.csv", "r");
    SUMMARY_str := Read(SUMMARY_FILE);
    if IsEof(SUMMARY_str) then
        SUMMARY_lines := [];
    else
        SUMMARY_lines := Split(SUMMARY_str, "\n");
    end if;

    //Check if the item is already in the summary file.
    if exists(i){i : i in [1..#SUMMARY_lines] | SUMMARY_lines[i] eq SummaryItem} then
        printf "The spectrum of F_%o^%o is already in the dataset.\n", q, "[" cat &cat([" "] cat [ IntegerToString(N[i]) cat " " : i in [1..#N] ]) cat "]";
        return;
    end if;

    //If not, we compute it and add it to the files.
    TheSpectrum := <>;
    try
        TheSpectrum := SpectrumFromRankDistribution(q,N);
    catch e
        printf "SpectrumFromRankDistribution failed: %o\n", e`Object;
        try
            TheSpectrum := SpectrumFromRank01Distribution(q,N);
        catch e2
            error "Both ways failed";
        end try;
    end try;
    
    
    //Fill the files.
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownSpectra.csv", "a");
    
    EIGENVALUE_FILE := Open("DATASET_EIGENVALUES_KnownSpectra.csv", "a");
    EigenvaluesString := &cat [ IntegerToString(TheSpectrum[i][1]) cat "," : i in [1..#TheSpectrum-1] ] cat IntegerToString(TheSpectrum[#TheSpectrum][1]) cat "\n"; 
    
    MULTIPLICITY_FILE := Open("DATASET_MULTIPLICITIES_KnownSpectra.csv", "a");
    MultiplicitiesString := &cat [ IntegerToString(TheSpectrum[i][2]) cat "," : i in [1..#TheSpectrum-1] ] cat IntegerToString(TheSpectrum[#TheSpectrum][2]) cat "\n"; 
    

    Write(SUMMARY_FILE, SummaryItem cat "\n");
    Write(EIGENVALUE_FILE, EigenvaluesString);
    Write(MULTIPLICITY_FILE, MultiplicitiesString);
    printf "The spectrum of F_%o^%o has been computed and stored in the dataset.\n", q, "[" cat &cat([" "] cat [ IntegerToString(N[i]) cat " " : i in [1..#N] ]) cat "]";
end procedure;

procedure AllStoreSpectrum()
    SUMMARY_FILE := Open("DATASET_SUMMARY_KnownRankDistributions.txt", "r");
    SUMMARY_str := Read(SUMMARY_FILE);
    SUMMARY_lines := Split(SUMMARY_str, "\n");
    
    for x_str in SUMMARY_lines do
        x := eval "return " cat x_str cat ";";
        q := x[1];
        NN := x[2];
        lastNub := &*NN;
        for lastN in [2..lastNub] do
            StoreSpectrum(q, NN cat [lastN]);
        end for;
    end for;
end procedure;
    
