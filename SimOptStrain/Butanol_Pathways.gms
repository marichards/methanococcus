$onmulti
Sets
i augment new metabolites to the list
/
1butanol
/
j augment new reactions to the list
/
EX_1butanol_e
/
Parameters
LowerLimits(j)
/
EX_1butanol_e 0
/
UpperLimits(j)
/
EX_1butanol_e 1000
/
S(i,j)
/
1butanol.EX_1butanol_e -1
/;
$offmulti
Sets
k reactions in UnivS (p)
/
ferm
keto

/;
Parameters
LowerLimits_k(k) universal metabolic reaction maximum fluxes
/
ferm 0
keto 0
/
UpperLimits_k(k) universal metabolic reaction maximum fluxes
/
ferm 1000
keto 1000
/
U(i,k) contains the UnivS matrix
/
accoa.ferm      -2
nadh.ferm       -4
h.ferm          -4
1butanol.ferm   1
coa.ferm        2
nad.ferm        4
h2o.ferm        1

oaa.keto        -1
glu-L.keto      -1
accoa.keto      -1
atp.keto        -2
nadph.keto      -2
h2o.keto        -2
h.keto          -2
1butanol.keto   1
akg.keto        1
coa.keto        1
co2.keto        2
adp.keto        2
pi.keto         2
nadp.keto       2
nh4.keto        1
/;
