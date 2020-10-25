type	CV
fit	ON
fitdata	Fc-ME1_5_mVps.txt	ALB
sweeps	2
rate	0.05
steps	0.002
pots	-0.1	0.4	-0.1
temp	300.01
*area	0.00543	0.002	0.2
ircdrop	ON
*resist	100	0	400
*dlcap	0	0	1e-05
species	4
symbols	Ox	Red	A	B
conc	0.00835	0	0	0
diff	1e-05	1e-05	0	0
charge	0	1	0	1
adsco	0.1	0.1	10	0.1	0	0	0	0
qmax	5e-08
model	BV
r	1	-1	1	0	0	0.2	0.5	1000	0	10	1e-16	0.3	1000	1000
*r	2	0	-1	1	0	1	2	0	3	2.5	6.7
r	3	0	-2	0	2	1	2	3	4	5	6
beta	0.2
Dm	10
nlerror	1e-08
maxiter	50000
nlflag	0
