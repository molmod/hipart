%mem=500MB
%chk=gaussian.chk
%nproc=1
#p B3LYP/6-31G(d) sp scf(fermi,tight) density(current)

Alanine

0 1
 O  1.657642604  0.744201101 -0.489836837
 O  1.102098787 -1.089074323  0.697793837
 N -1.404760606 -1.091053462 -0.458202548
 C -0.655372082  0.161523560 -0.415839794
 C -1.204817386  1.230334427  0.557894814
 C  0.774816114 -0.164020328 -0.012893818
 H -0.627559422  0.594521751 -1.422939548
 H -0.609640466  2.149682405  0.525945464
 H -2.235463406  1.483421899  0.284034915
 H -1.205447577  0.848510089  1.585349081
 H -2.403964854 -0.897766300 -0.440132543
 H -1.180805308 -1.635048404  0.373826645
 H  2.530514260  0.486012629 -0.137289368

