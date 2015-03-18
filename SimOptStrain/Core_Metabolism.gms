********************************************************************************
* Core E. coli metabolic model
* EcoSal Chapter 10.2.1 Reconstruction and Use of Microbial Metabolic Networks:
* the Core Escherichia coli Metabolic Model as an Educational Guide by Orth,
* Fleming, and Palsson (2010)
*
* Joohnoon Kim, Jennifer L. Reed, and Christos T. Maravelias
* Department of Chemical and Biological Engineering
* University of Wisconsin-Madison, Madison, Wisconsin, USA
********************************************************************************

$offlisting
$offdigit

Sets
i metabolites in S (m)
/13dpg,2pg,3pg,6pgc,6pgl,ac,ac_e,acald,acald_e,accoa,acon-C,actp,adp,akg,akg_e,amp,atp,cit,co2,co2_e,coa,dhap,e4p,etoh,etoh_e,f6p,fdp,formate,for_e,fru_e,fum,fum_e,g3p,g6p,glc-D_e,gln-L,gln-L_e,glu-L,glu-L_e,glx,h2o,h2o_e,h,h_e,icit,lac-D,lac-D_e,mal-L,mal-L_e,nad,nadh,nadp,nadph,nh4,nh4_e,o2,o2_e,oaa,pep,pi,pi_e,pyr,pyr_e,q8,q8h2,r5p,ru5p-D,s7p,succ,succ_e,succoa,xu5p-D
/
j reactions in S (n)
/ACALD,ACALDt,ACKr,ACONTa,ACONTb,ACt2r,ADK1,AKGDH,AKGt2r,ALCD2x,ATPM,ATPS4r,Biomass_Ecoli_core_w_GAM,CO2t,CS,CYTBD,D_LACt2,ENO,ETOHt2r,EX_ac_e,EX_acald_e,EX_akg_e,EX_co2_e,EX_etoh_e,EX_for_e,EX_fru_e,EX_fum_e,EX_glc_e,EX_gln_L_e,EX_glu_L_e,EX_h_e,EX_h2o_e,EX_lac_D_e,EX_mal_L_e,EX_nh4_e,EX_o2_e,EX_pi_e,EX_pyr_e,EX_succ_e,FBA,FBP,FORt2,FORti,FRD7,FRUpts2,FUM,FUMt2_2,G6PDH2r,GAPD,GLCpts,GLNS,GLNabc,GLUDy,GLUN,GLUSy,GLUt2r,GND,H2Ot,ICDHyr,ICL,LDH_D,MALS,MALt2_2,MDH,ME1,ME2,NADH16,NADTRHD,NH4t,O2t,PDH,PFK,PFL,PGI,PGK,PGL,PGM,PIt2r,PPC,PPCK,PPS,PTAr,PYK,PYRt2r,RPE,RPI,SUCCt2_2,SUCCt3,SUCDi,SUCOAS,TALA,THD2,TKT1,TKT2,TPI
/
Parameters
Vmax /1000/
UpperLimits(j) maximum value flux can take
LowerLimits(j) minimum value flux can take
/(
ACALD,ACALDt,ACKr,ACONTa,ACONTb,ACt2r,ADK1,AKGt2r,ALCD2x,ATPS4r,CO2t,D_LACt2,ENO,ETOHt2r,FBA,FUM,G6PDH2r,GAPD,GLUDy,GLUt2r,H2Ot,ICDHyr,LDH_D,MDH,NH4t,O2t,PGI,PGK,PGM,PIt2r,PTAr,PYRt2r,RPE,RPI,SUCOAS,TALA,TKT1,TKT2,TPI
)-1000/
S(i,j) contains the S matrix
/
acald.ACALD     -1.000000
accoa.ACALD     1.000000
coa.ACALD       -1.000000
h.ACALD 1.000000
nad.ACALD       -1.000000
nadh.ACALD      1.000000
acald.ACALDt    1.000000
acald_e.ACALDt  -1.000000
ac.ACKr -1.000000
actp.ACKr       1.000000
adp.ACKr        1.000000
atp.ACKr        -1.000000
acon-C.ACONTa   1.000000
cit.ACONTa      -1.000000
h2o.ACONTa      1.000000
acon-C.ACONTb   -1.000000
h2o.ACONTb      -1.000000
icit.ACONTb     1.000000
ac.ACt2r        1.000000
ac_e.ACt2r      -1.000000
h.ACt2r 1.000000
h_e.ACt2r       -1.000000
adp.ADK1        2.000000
amp.ADK1        -1.000000
atp.ADK1        -1.000000
akg.AKGDH       -1.000000
co2.AKGDH       1.000000
coa.AKGDH       -1.000000
nad.AKGDH       -1.000000
nadh.AKGDH      1.000000
succoa.AKGDH    1.000000
akg.AKGt2r      1.000000
akg_e.AKGt2r    -1.000000
h.AKGt2r        1.000000
h_e.AKGt2r      -1.000000
acald.ALCD2x    1.000000
etoh.ALCD2x     -1.000000
h.ALCD2x        1.000000
nad.ALCD2x      -1.000000
nadh.ALCD2x     1.000000
adp.ATPM        1.000000
atp.ATPM        -1.000000
h2o.ATPM        -1.000000
h.ATPM  1.000000
pi.ATPM 1.000000
adp.ATPS4r      -1.000000
atp.ATPS4r      1.000000
h2o.ATPS4r      1.000000
h.ATPS4r        3.000000
h_e.ATPS4r      -4.000000
pi.ATPS4r       -1.000000
3pg.Biomass_Ecoli_core_w_GAM    -1.496000
accoa.Biomass_Ecoli_core_w_GAM  -3.747800
adp.Biomass_Ecoli_core_w_GAM    59.810000
akg.Biomass_Ecoli_core_w_GAM    4.118200
atp.Biomass_Ecoli_core_w_GAM    -59.810000
coa.Biomass_Ecoli_core_w_GAM    3.747800
e4p.Biomass_Ecoli_core_w_GAM    -0.361000
f6p.Biomass_Ecoli_core_w_GAM    -0.070900
g3p.Biomass_Ecoli_core_w_GAM    -0.129000
g6p.Biomass_Ecoli_core_w_GAM    -0.205000
gln-L.Biomass_Ecoli_core_w_GAM  -0.255700
glu-L.Biomass_Ecoli_core_w_GAM  -4.941400
h2o.Biomass_Ecoli_core_w_GAM    -59.810000
h.Biomass_Ecoli_core_w_GAM      59.810000
nad.Biomass_Ecoli_core_w_GAM    -3.547000
nadh.Biomass_Ecoli_core_w_GAM   3.547000
nadp.Biomass_Ecoli_core_w_GAM   13.027900
nadph.Biomass_Ecoli_core_w_GAM  -13.027900
oaa.Biomass_Ecoli_core_w_GAM    -1.786700
pep.Biomass_Ecoli_core_w_GAM    -0.519100
pi.Biomass_Ecoli_core_w_GAM     59.810000
pyr.Biomass_Ecoli_core_w_GAM    -2.832800
r5p.Biomass_Ecoli_core_w_GAM    -0.897700
co2.CO2t        1.000000
co2_e.CO2t      -1.000000
accoa.CS        -1.000000
cit.CS  1.000000
coa.CS  1.000000
h2o.CS  -1.000000
h.CS    1.000000
oaa.CS  -1.000000
h2o.CYTBD       1.000000
h.CYTBD -2.000000
h_e.CYTBD       2.000000
o2.CYTBD        -0.500000
q8.CYTBD        1.000000
q8h2.CYTBD      -1.000000
h.D_LACt2       1.000000
h_e.D_LACt2     -1.000000
lac-D.D_LACt2   1.000000
lac-D_e.D_LACt2 -1.000000
2pg.ENO -1.000000
h2o.ENO 1.000000
pep.ENO 1.000000
etoh.ETOHt2r    1.000000
etoh_e.ETOHt2r  -1.000000
h.ETOHt2r       1.000000
h_e.ETOHt2r     -1.000000
ac_e.EX_ac_e    -1.000000
acald_e.EX_acald_e      -1.000000
akg_e.EX_akg_e  -1.000000
co2_e.EX_co2_e  -1.000000
etoh_e.EX_etoh_e        -1.000000
for_e.EX_for_e  -1.000000
fru_e.EX_fru_e  -1.000000
fum_e.EX_fum_e  -1.000000
glc-D_e.EX_glc_e        -1.000000
gln-L_e.EX_gln_L_e      -1.000000
glu-L_e.EX_glu_L_e      -1.000000
h_e.EX_h_e      -1.000000
h2o_e.EX_h2o_e  -1.000000
lac-D_e.EX_lac_D_e      -1.000000
mal-L_e.EX_mal_L_e      -1.000000
nh4_e.EX_nh4_e  -1.000000
o2_e.EX_o2_e    -1.000000
pi_e.EX_pi_e    -1.000000
pyr_e.EX_pyr_e  -1.000000
succ_e.EX_succ_e        -1.000000
dhap.FBA        1.000000
fdp.FBA -1.000000
g3p.FBA 1.000000
f6p.FBP 1.000000
fdp.FBP -1.000000
h2o.FBP -1.000000
pi.FBP  1.000000
formate.FORt2       1.000000
for_e.FORt2     -1.000000
h.FORt2 1.000000
h_e.FORt2       -1.000000
formate.FORti       -1.000000
for_e.FORti     1.000000
fum.FRD7        -1.000000
q8.FRD7 1.000000
q8h2.FRD7       -1.000000
succ.FRD7       1.000000
f6p.FRUpts2     1.000000
fru_e.FRUpts2   -1.000000
pep.FRUpts2     -1.000000
pyr.FRUpts2     1.000000
fum.FUM -1.000000
h2o.FUM -1.000000
mal-L.FUM       1.000000
fum.FUMt2_2     1.000000
fum_e.FUMt2_2   -1.000000
h.FUMt2_2       2.000000
h_e.FUMt2_2     -2.000000
6pgl.G6PDH2r    1.000000
g6p.G6PDH2r     -1.000000
h.G6PDH2r       1.000000
nadp.G6PDH2r    -1.000000
nadph.G6PDH2r   1.000000
13dpg.GAPD      1.000000
g3p.GAPD        -1.000000
h.GAPD  1.000000
nad.GAPD        -1.000000
nadh.GAPD       1.000000
pi.GAPD -1.000000
g6p.GLCpts      1.000000
glc-D_e.GLCpts  -1.000000
pep.GLCpts      -1.000000
pyr.GLCpts      1.000000
adp.GLNS        1.000000
atp.GLNS        -1.000000
gln-L.GLNS      1.000000
glu-L.GLNS      -1.000000
h.GLNS  1.000000
nh4.GLNS        -1.000000
pi.GLNS 1.000000
adp.GLNabc      1.000000
atp.GLNabc      -1.000000
gln-L.GLNabc    1.000000
gln-L_e.GLNabc  -1.000000
h2o.GLNabc      -1.000000
h.GLNabc        1.000000
pi.GLNabc       1.000000
akg.GLUDy       1.000000
glu-L.GLUDy     -1.000000
h2o.GLUDy       -1.000000
h.GLUDy 1.000000
nadp.GLUDy      -1.000000
nadph.GLUDy     1.000000
nh4.GLUDy       1.000000
gln-L.GLUN      -1.000000
glu-L.GLUN      1.000000
h2o.GLUN        -1.000000
nh4.GLUN        1.000000
akg.GLUSy       -1.000000
gln-L.GLUSy     -1.000000
glu-L.GLUSy     2.000000
h.GLUSy -1.000000
nadp.GLUSy      1.000000
nadph.GLUSy     -1.000000
glu-L.GLUt2r    1.000000
glu-L_e.GLUt2r  -1.000000
h.GLUt2r        1.000000
h_e.GLUt2r      -1.000000
6pgc.GND        -1.000000
co2.GND 1.000000
nadp.GND        -1.000000
nadph.GND       1.000000
ru5p-D.GND      1.000000
h2o.H2Ot        1.000000
h2o_e.H2Ot      -1.000000
akg.ICDHyr      1.000000
co2.ICDHyr      1.000000
icit.ICDHyr     -1.000000
nadp.ICDHyr     -1.000000
nadph.ICDHyr    1.000000
glx.ICL 1.000000
icit.ICL        -1.000000
succ.ICL        1.000000
h.LDH_D 1.000000
lac-D.LDH_D     -1.000000
nad.LDH_D       -1.000000
nadh.LDH_D      1.000000
pyr.LDH_D       1.000000
accoa.MALS      -1.000000
coa.MALS        1.000000
glx.MALS        -1.000000
h2o.MALS        -1.000000
h.MALS  1.000000
mal-L.MALS      1.000000
h.MALt2_2       2.000000
h_e.MALt2_2     -2.000000
mal-L.MALt2_2   1.000000
mal-L_e.MALt2_2 -1.000000
h.MDH   1.000000
mal-L.MDH       -1.000000
nad.MDH -1.000000
nadh.MDH        1.000000
oaa.MDH 1.000000
co2.ME1 1.000000
mal-L.ME1       -1.000000
nad.ME1 -1.000000
nadh.ME1        1.000000
pyr.ME1 1.000000
co2.ME2 1.000000
mal-L.ME2       -1.000000
nadp.ME2        -1.000000
nadph.ME2       1.000000
pyr.ME2 1.000000
h.NADH16        -4.000000
h_e.NADH16      3.000000
nad.NADH16      1.000000
nadh.NADH16     -1.000000
q8.NADH16       -1.000000
q8h2.NADH16     1.000000
nad.NADTRHD     -1.000000
nadh.NADTRHD    1.000000
nadp.NADTRHD    1.000000
nadph.NADTRHD   -1.000000
nh4.NH4t        1.000000
nh4_e.NH4t      -1.000000
o2.O2t  1.000000
o2_e.O2t        -1.000000
accoa.PDH       1.000000
co2.PDH 1.000000
coa.PDH -1.000000
nad.PDH -1.000000
nadh.PDH        1.000000
pyr.PDH -1.000000
adp.PFK 1.000000
atp.PFK -1.000000
f6p.PFK -1.000000
fdp.PFK 1.000000
h.PFK   1.000000
accoa.PFL       1.000000
coa.PFL -1.000000
formate.PFL 1.000000
pyr.PFL -1.000000
f6p.PGI 1.000000
g6p.PGI -1.000000
13dpg.PGK       1.000000
3pg.PGK -1.000000
adp.PGK 1.000000
atp.PGK -1.000000
6pgc.PGL        1.000000
6pgl.PGL        -1.000000
h2o.PGL -1.000000
h.PGL   1.000000
2pg.PGM -1.000000
3pg.PGM 1.000000
h.PIt2r 1.000000
h_e.PIt2r       -1.000000
pi.PIt2r        1.000000
pi_e.PIt2r      -1.000000
co2.PPC -1.000000
h2o.PPC -1.000000
h.PPC   1.000000
oaa.PPC 1.000000
pep.PPC -1.000000
pi.PPC  1.000000
adp.PPCK        1.000000
atp.PPCK        -1.000000
co2.PPCK        1.000000
oaa.PPCK        -1.000000
pep.PPCK        1.000000
amp.PPS 1.000000
atp.PPS -1.000000
h2o.PPS -1.000000
h.PPS   2.000000
pep.PPS 1.000000
pi.PPS  1.000000
pyr.PPS -1.000000
accoa.PTAr      -1.000000
actp.PTAr       1.000000
coa.PTAr        1.000000
pi.PTAr -1.000000
adp.PYK -1.000000
atp.PYK 1.000000
h.PYK   -1.000000
pep.PYK -1.000000
pyr.PYK 1.000000
h.PYRt2r        1.000000
h_e.PYRt2r      -1.000000
pyr.PYRt2r      1.000000
pyr_e.PYRt2r    -1.000000
ru5p-D.RPE      -1.000000
xu5p-D.RPE      1.000000
r5p.RPI -1.000000
ru5p-D.RPI      1.000000
h.SUCCt2_2      2.000000
h_e.SUCCt2_2    -2.000000
succ.SUCCt2_2   1.000000
succ_e.SUCCt2_2 -1.000000
h.SUCCt3        1.000000
h_e.SUCCt3      -1.000000
succ.SUCCt3     -1.000000
succ_e.SUCCt3   1.000000
fum.SUCDi       1.000000
q8.SUCDi        -1.000000
q8h2.SUCDi      1.000000
succ.SUCDi      -1.000000
adp.SUCOAS      1.000000
atp.SUCOAS      -1.000000
coa.SUCOAS      -1.000000
pi.SUCOAS       1.000000
succ.SUCOAS     -1.000000
succoa.SUCOAS   1.000000
e4p.TALA        1.000000
f6p.TALA        1.000000
g3p.TALA        -1.000000
s7p.TALA        -1.000000
h.THD2  2.000000
h_e.THD2        -2.000000
nad.THD2        1.000000
nadh.THD2       -1.000000
nadp.THD2       -1.000000
nadph.THD2      1.000000
g3p.TKT1        1.000000
r5p.TKT1        -1.000000
s7p.TKT1        1.000000
xu5p-D.TKT1     -1.000000
e4p.TKT2        -1.000000
f6p.TKT2        1.000000
g3p.TKT2        1.000000
xu5p-D.TKT2     -1.000000
dhap.TPI        -1.000000
g3p.TPI 1.000000
/;
