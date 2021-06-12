%NProcShared=4
%Mem=2GB
%oldchk=ethylene-anion.chk
%chk=anion
# B3LYP/6-31+G(d,p) Freq=SaveNM Geom=Check

anion

-1 2


--Link1--
%oldchk=ethylene-neutral.chk
%chk=neutral
# B3LYP/6-31+G(d,p) Freq=SaveNM Geom=Check

anion

0 1

--Link1--
%oldchk=anion
%chk=fc
# B3LYP/6-31+G(d,p) Freq=(FC,ReadFCHT) Geom=Check

photoionization

-1,2

PRTMAT=12 MAXC1=1000 MAXBANDS=1 SPECHWHM=0.5

neutral.chk
