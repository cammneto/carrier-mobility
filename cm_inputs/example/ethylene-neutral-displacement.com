%NProcShared=4
%Mem=2GB
%oldchk=ethylene-neutral.chk
%chk=neutral
# B3LYP/6-31+G(d,p) Freq=(SaveNM,noraman) Geom=Checkpoint guess=read

neutral form

0 1


--Link1--
%oldchk=ethylene-anion.chk
%chk=anion
# B3LYP/6-31+G(d,p) Freq=(SaveNM,noraman) Geom=Checkpoint guess=read

anion

-1 2

--Link1--
%oldchk=neutral
%chk=fc
# B3LYP/6-31+G(d,p) Freq=(FC,ReadFCHT) Geom=Check

photoionization

0,1

PRTMAT=12 MAXC1=1000 MAXBANDS=1 SPECHWHM=0.5

anion.chk

