********************************************************************************
* SimOptStrain example code
* Kim J, Reed JL, Maravelias CT (2011) Large-Scale Bi-Level Strain Design
* Approaches and Mixed-Integer Programming Solution Techniques. PLoS ONE 6(9):
* e24162. doi:10.1371/journal.pone.0024162
*
* Joohnoon Kim, Jennifer L. Reed, and Christos T. Maravelias
* Department of Chemical and Biological Engineering
* University of Wisconsin-Madison, Madison, Wisconsin, USA
********************************************************************************

$Offlisting
$include Core_Metabolism.gms
$include Butanol_Pathways.gms

$onecho > cplex.o10
indic vbound(jk)$d(jk) 1
indic u2bound(jk)$d(jk) 0
indic u3bound(jk)$d(jk) 1
indic abound(k)$a(k) 0
indic dualbalance_k(k)$a(k) 1
eprhs 1e-8
epopt 1e-8
epint 0
$offecho

Sets
soln             possible solutions in the solution pool /s1*s5/
store(soln)      identified solutions in the solution pool
exchg(j)         used to define exchange reactions /EX_ac_e,EX_acald_e,EX_akg_e,EX_co2_e,EX_etoh_e,EX_for_e,EX_fru_e,EX_fum_e,EX_glc_e,EX_gln_L_e,EX_glu_L_e,EX_h_e,EX_h2o_e,EX_lac_D_e,EX_mal_L_e,EX_nh4_e,EX_o2_e,EX_pi_e,EX_pyr_e,EX_succ_e,EX_1butanol_e/
essen(j)         used to define essential reactions /ACONTa,ACONTb,Biomass_Ecoli_core_w_GAM,CS,ENO,EX_glc_e,EX_h_e,EX_nh4_e,EX_pi_e,GAPD,GLCpts,GLNS,ICDHyr,NH4t,PGK,PGM,PIt2r,RPI/
nouse(j)         used to define unused reactions /EX_fru_e,EX_fum_e,EX_gln_L_e,EX_mal_L_e,FRUpts2,FUMt2_2,GLNabc,MALt2_2/
trans(j)         used to define transport reactions /ACALDt,ACt2r,AKGt2r,CO2t,D_LACt2,ETOHt2r,FORt2,FORti,FRUpts2,FUMt2_2,GLCpts,GLNabc,GLUt2r,H2Ot,MALt2_2,NH4t,O2t,PIt2r,PYRt2r,SUCCt2_2,SUCCt3/
exclude(j)       used to define excluded reactions /ATPM,ATPS4r,ADK1,TKT1,TKT2/
bnd(j)           used to define bounded reactions
bnd_k(k)         used to define bounded non-native reactions
jk(j)            used to define reactions considered for deletion
subi(i)          used to define subset of metabolites
subj(j)          used to define subset of reactions;

jk(j)=not (exchg(j) or essen(j) or nouse(j) or trans(j) or exclude(j));
subj(j)=not nouse(j);
subi(i)$sum(subj$S(i,subj),1)=yes;

Scalar
alpha            minimum increase of the objective per deletion
beta             minimum increase of the objective per addition
Dmax             maximum number of deletions
Amax             maximum number of additions
Umax             bounds on the dual variables /1/
tolerance /1e-12/;

Parameters
c(j)             used to define the objective function for outer problem
pri(j)           used to define the objective function for inner problem
storeKnock(soln) store number of reaction deletions for each solution
storeAdded(soln) store number of reaction additions for each solution
storeD(j,soln)   store reaction deletions for each solution
storeA(k,soln)   store reaction additions for each solution
storeV(j,soln)   store reaction flux value for each solution
maxobj(soln)     maximum objective value for each solution
minobj(soln)     minimum objective value for each solution;

Variables
v(j)             flux values through reaction in network
v_k(k)           flux values through non-native reaction in network
u1(i)            dual variable for S*v = 0
u2(j)            dual variable for   v = 0
Obj              the objective function;

Positive Variables
u3(j)            dual variable for -v =l= -vmin
u3_k(k)          dual variable for -v_k =l= -vmin_k;

Binary Variables
d(j)             indicator for reactions
a(k)             indicator for non-native reactions;

Equation
calcobj          objective function
optobj           SimOptStrain objective function
maxknock         maximum number of allowed deletions
maxadded         maximum number of allowed additions
mingrowth        minimum growth requirement
primaldual       primal objective = dual objective
massbalance(i)   steady state mass balance
dualbalance(j)   dual constraint for primal variable v
dualbalance_k(k) dual constraint for primal variable v_k
vbound(j)        set bound for v from reaction indicators
u2bound(j)       set bound for u2 from reaction indicators
u3bound(j)       set bound for u3 from reaction indicators
abound(k)        set bound for v_k from reaction indicators
intcut(soln)     integer cut to generate diverse solutions;

optobj..         Obj=e=sum( j$subj(j),c(j)*v(j) ) - alpha*sum(jk, d(jk)) - beta*sum(k, a(k));
calcobj..        Obj=e=sum( j$subj(j),c(j)*v(j) );
maxknock..       sum( jk,d(jk) )=l=Dmax;
maxadded..       sum( k,a(k) )=l=Amax;
mingrowth..      sum( j$subj(j),pri(j)*v(j) )=g=0.1;

primaldual..     sum( j$subj(j),pri(j)*v(j) )=e= -sum( j$(subj(j) and bnd(j)),v.lo(j)*u3(j) )-sum( k$bnd_k(k),v_k.lo(k)*u3_k(k) );
massbalance(i)$subi(i).. sum( j$subj(j),S(i,j)*v(j) ) + sum( k,U(i,k)*v_k(k) )=e=0;
dualbalance(j)$subj(j).. sum( i$subi(i),S(i,j)*u1(i) )+u2(j)$jk(j)-u3(j)$bnd(j)=e=pri(j);
dualbalance_k(k).. sum( i$subi(i),U(i,k)*u1(i) )-u3_k(k)$bnd_k(k)=e=0;

vbound(jk)..     v(jk)=e=0;
u2bound(jk)..    u2(jk)=e=0;
u3bound(jk)..    u3(jk)=e=0;
abound(k)..      v_k(k)=e=0;

intcut(soln)$store(soln).. sum(jk$storeD(jk,soln),d(jk))-sum(jk$(not storeD(jk,soln)),d(jk))
                          +sum(k$storeA(k,soln),a(k))-sum(k$(not storeA(k,soln)),a(k))
                            =l=storeKnock(soln)-1+storeAdded(soln)-1;

Model FBA    /calcobj,massbalance/;
Model SimOptStrain /optobj,maxknock,maxadded,mingrowth,primaldual,massbalance,
                    dualbalance,dualbalance_k,vbound,u2bound,u3bound,abound,intcut/;

SimOptStrain.optfile=10;
SimOptStrain.optcr=0;
SimOptStrain.holdfixed=1;
SimOptStrain.reslim=1000000;

* Set the lower and upper bounds for variables
LowerLimits('ATPM')=8.39;
LowerLimits('EX_o2_e')=-Vmax;
LowerLimits('EX_co2_e')=-Vmax;
LowerLimits('EX_h_e')=-Vmax;
LowerLimits('EX_h2o_e')=-Vmax;
LowerLimits('EX_nh4_e')=-Vmax;
LowerLimits('EX_pi_e')=-Vmax;
LowerLimits('EX_glc_e')=-10;

bnd(j)$(LowerLimits(j) gt -Vmax)=yes;
bnd_k(k)$(LowerLimits_k(k) gt -Vmax)=yes;
LowerLimits(j)$(LowerLimits(j) = -Vmax)=-inf;
UpperLimits(j)=inf;
LowerLimits_k(k)$(LowerLimits_k(k) = -Vmax)=-inf;
UpperLimits_k(k)=inf;

v.lo(j)=LowerLimits(j);
v.up(j)=UpperLimits(j);
v_k.lo(k)=LowerLimits_k(k);
v_k.up(k)=UpperLimits_k(k);

* Define the outer and inner objectives
c('EX_1butanol_e')=1;
pri('Biomass_Ecoli_core_w_GAM')=1;
pri('EX_1butanol_e')=-1e-6;

* Calculate the theoretical maximum production
solve FBA using lp maximizing Obj;
abort$(Obj.l=0) "Maximum production rate is zero";
* Set the minimum increase of the objective per deletion and addition
alpha=0.01;
beta=0.01;

* Set the lower and upper bounds for variables
u2.lo(j)=-Umax; u2.up(j)=Umax; u3_k.up(k)=Umax;

* Define the allowed number of deletions and additions
Dmax=3;
Amax=1;

* Initialize parameter values
storeKnock(soln)=Dmax;   storeAdded(soln)=Amax;
StoreA(k,soln)=0;        storeD(j,soln)=0;        storeV(j,soln)=0;
maxobj(soln)=0;          minobj(soln)=0;

scalar iter;
file result /SimOptStrain_CoreNetwork.txt/;
result.pw=1000; result.pc=5; result.ps=130;
put result;
put 'Number';
put "Deletions";
iter=1;
while(iter lt Dmax, put ''; iter=iter+1;);
put 'Number';
put "Additions";
iter=1;
while(iter lt Amax, put ''; iter=iter+1;);
put 'Growth';
put 'MaxProd';
put 'MinProd';
put /;
putclose result;
result.ap=1;

Alias(soln,temp);
store(soln)=no;

loop(temp,
         v.lo(j)=LowerLimits(j);
         v.up(j)=UpperLimits(j);
         v_k.lo(k)=LowerLimits_k(k);
         v_k.up(k)=UpperLimits_k(k);

         solve SimOptStrain using mip maximizing Obj;

         store(temp)$((SimOptStrain.modelstat = 1 or 8) and Obj.l gt 1e-3)=yes;
         if (store(temp),
         storeKnock(temp)=sum(jk, d.l(jk));
         storeAdded(temp)=sum(k, a.l(k));
         storeD(j,temp)=d.l(j);
         storeA(k,temp)=a.l(k);
         storeV(j,temp)=v.l(j);
         maxobj(temp)=sum(j$c(j), v.l(j));

         v.fx(jk)$d.l(jk)=0;
         v.lo(j)$(abs(v.lo(j)) lt tolerance) = 0;
         v.up(j)$(abs(v.up(j)) lt tolerance) = 0;
         v.fx(j)$(pri(j)=1)=storeV(j,temp);
         v_k.fx(k)$(not a.l(k))=0;

         solve FBA using lp minimizing Obj;
         minobj(temp)=obj.l;
         v.lo(j)$(pri(j)=1)=0;
         v.up(j)$(pri(j)=1)=UpperLimits(j);

         iter=1;
         put result;
         put storeKnock(temp);
         loop(jk$(storeD(jk,temp) = 1), put jk.tl; );
         iter=storeKnock(temp);
         while(iter lt Dmax, put ''; iter=iter+1;);
         put storeAdded(temp);
         loop(k$(storeA(k,temp) = 1), put k.tl; );
         iter=storeAdded(temp);
         while(iter lt Amax, put ''; iter=iter+1;);
         put storeV('Biomass_Ecoli_core_w_GAM',temp):8:4;
         put maxobj(temp):8:4;
         put minobj(temp):8:4;
         put /;
         putclose result;
         );
);
