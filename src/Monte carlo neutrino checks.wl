(* ::Package:: *)

(* example entry in a notebook in the same directory as this file to run this code

Get[FileNameJoin[{NotebookDirectory[],"Monte carlo neutrino checks.wl"}]];
kappan[100 * G15TOFM, 5 / HBARC, 0.04870282860614022, 0.16 * 10^-3, 0.25, 20 / HBARC, 0.6, 0]
kappanmc[100 * G15TOFM, 5 / HBARC, 0.04870282860614022`, -0.09603718718498111`, -0.05445418660894468`, 20  / HBARC, 0.6, 0, 10]

*)

(*This package calculates Monte Carlo integrations of neutrino opacities and scattering to check accuracy of approximation methods*)
(*All results are given in units of GF^2 or GF^2 cos^2 theta_c*)

(*constants*)
SINTW = Sqrt[0.231]; (*Weinberg angle*)
HBARC = 197.3269718; (*all calculations will be done in fm*)
NUCMASS = 2 * 939.5653 * 938.272 / (939.5653 + 938.272) / HBARC; (*average nucleon mass*)
MSPLIT = (939.5653 - 938.272) / HBARC; (* nucleon mass splitting *)
ELECMASS = 0.511 / HBARC; (* electron mass, only for use in kinematics *)
GA = 1.27; (*no g-t discrepancy for electroweak*)
G15TOFM = 19.5 * Sqrt[4 * Pi / 137] / HBARC^2; (*convert B in 10^15 G to eB in fm^2 with the convention that alpha = e^2/4pi*)
GP = 5.5858; GN = -3.8263; (*gyromagnetic ratios for nucleons, note that GE=2.002 so the g-2 is negligible for leptons*)

(* functions for getting chemical potentials and distributions *)
nfd[e_, mu_, t_]:= 1 / (Exp[(e - mu) / t] + 1);(*If[(e - mu) / t > 50, 0, If[(e - mu) / t < -50, 1, 1 / (Exp[(e - mu) / t] + 1)]]*)
munofn[n_, eb_, t_]:= FindRoot[-n + NIntegrate[k^2 / (2 * Pi^2) * (nfd[k^2 / (2 * NUCMASS) + GN * eb / (4 * NUCMASS), mu, t] 
	+ nfd[k^2 / (2 * NUCMASS) - GN * eb / (4 * NUCMASS), mu, t]), {k, 0, Infinity}], {mu, t}];
mueofn[n_, eb_, t_]:= If[t^2 / eb > 100, FindRoot[-n + NIntegrate[k^2 / Pi^2 * nfd[k, mu, t], {k, 0, Infinity}], {mu, t}],
	FindRoot[-n + NIntegrate[eb / (4 * Pi^2) * nfd[Sqrt[k^2 + ELECMASS^2], mu, t] + Sum[eb / (2 * Pi^2) * nfd[Sqrt[k^2 + 2 * ne * eb + ELECMASS^2], mu, t], 
	{ne, 1, Ceiling[t^2 / eb * 5]}], {k, -Infinity, Infinity}], {mu, t}]];
mupofn[n_, eb_, t_]:= If[NUCMASS * t / eb > 100, FindRoot[-n + NIntegrate[k^2 / (2 * Pi^2) * (nfd[k^2 / (2 * NUCMASS) + GP * eb / (4 * NUCMASS), mu, t] 
	+ nfd[k^2 / (2 * NUCMASS) - GP * eb / (4 * NUCMASS), mu, t]), {k, 0, Infinity}], {mu, t}],
	FindRoot[-n + NIntegrate[eb / (4 * Pi^2) * nfd[k^2 / (2 * NUCMASS) - (GP - 2) * eb / (4 * NUCMASS), mu, t] 
	+ Sum[eb / (4 * Pi^2) * (nfd[k^2 / (2 * NUCMASS) + np * eb / NUCMASS + (GP - 2) * eb / (4 * NUCMASS), mu, t] + nfd[k^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * eb / (4 * NUCMASS), mu, t]), 
	{np, 1, Ceiling[NUCMASS * t * 5 / eb]}], {k, -Infinity, Infinity}], {mu, t}]];
(*lepton kinematics for cc interactions*)
thetapm[x_, sp_, sn_, eb_, ui_, pm_]:=HeavisideTheta[MSPLIT + pm * ELECMASS + ui - eb / (4 * NUCMASS) * (GN * sn - (GP - 2) * sp) - x];
e0pm[knu_, sp_, sn_, eb_, ui_, pm_]:= Module[{delta},
delta = MSPLIT + ui - eb / (4 * NUCMASS) * (GN * sn - (GP - 2) * sp);
If[Abs[knu + pm * delta] > ELECMASS, Sqrt[(knu + pm * delta)^2 - ELECMASS^2], 0]];
	
(*functions for degenerate neutron regime*)
ufunc[eb_, temp_]:= (NUCMASS * temp) / eb * (1 - Exp[-eb / (NUCMASS * temp)]);
f1[k0_, u_]:= Sqrt[Pi]/(2*u*Sqrt[1-u]) * (Sqrt[1-u] * Erf[k0] - Exp[-k0^2 * u] * Erf[k0 * Sqrt[1-u]]);
f2[k0_, u_]:= 1 / (4 * u^2 * Sqrt[1-u]) * (Sqrt[Pi] * Sqrt[1 - u] * (2 + u) * Erf[k0] - 2 * Sqrt[Pi] * Exp[-k0^2 * u] * (1 + k0^2 * u) * Erf[k0 * Sqrt[1-u]] 
	- 2 * k0 * u * Exp[-k0^2] * Sqrt[1-u]);
f3[k0_, u_]:= 1 / (4 * u * (1 - u)^(3/2)) * (Sqrt[Pi] * (1 - u)^(3/2) * Erf[k0] - Exp[-k0^2 * u] * Sqrt[Pi] * Erf[k0 * Sqrt[1 - u]] + 2 * u * Sqrt[1 - u] * k0 * Exp[-k0^2]);
g1[k0_, u_]:= Pi^2 * Exp[-k0^2] / (12 * k0 * Sqrt[1 - u]) * (Sqrt[1 - u] - k0 * u * Sqrt[Pi] * Erf[k0 * Sqrt[1 - u]]);
g2[k0_, u_]:= Pi^2 * Exp[-k0^2] / (12 * Sqrt[1 - u]) * (k0 * Sqrt[1 - u] - Exp[-k0^2 * (u - 1)] * (k0^2 * u - 1) * Sqrt[Pi] * Erf[k0 * Sqrt[1 - u]]);
g3[k0_, u_]:= Pi^2 * Exp[-k0^2] / (24 * (1 - u)^(3/2)) * (2 * k0 * Sqrt[1 - u]  - Exp[-k0^2 * (u - 1)] * u * Sqrt[Pi] * Erf[k0 * Sqrt[1 - u]]);
bigb[kfn_, u_, w_, eb_, t_, knu_, cost_]:= Module[{mt, kfnt, knut},
mt = NUCMASS * t;
kfnt = kfn / Sqrt[2 * mt];
knut = knu / Sqrt[2 * mt];
(2 * mt)^(3 / 2) / (2 * Pi)^2 * (f1[kfnt, u] + w * g1[kfnt, u] + knut^2 * ((f1[kfnt, u] + w * g1[kfnt, u]) * (u + cost^2 - u * cost^2)
	+ (f2[kfnt, u] + w * g2[kfnt, u]) * (u^2 * cost^2 - u^2) + (f3[kfnt, u] + w * g3[kfnt, u]) * (u^2 + 2 * cost^2 + u^2 * cost^2)))];

(*functions for nc scattering kernels*)
ffunc[q0_,qz_, qp_, temp_, eb_, s_, spr_]:= Module[{t, q, qpt, expterm, prefactor, sfunc, sarg},
q = Sqrt[qz^2 + qp^2];
t = Exp[-eb / (NUCMASS * temp)];
qpt = qp * Sqrt[eb * Sqrt[t] / (NUCMASS * temp * (1 - t))];
sarg = q * eb / (temp * Sqrt[2]) * (qz^2 + qpt^2) / (qz^2 * qpt^2);
sfunc = If[q0 > q, Erfc[sarg], 1 + Erf[sarg]];
expterm = Exp[-qpt^2 * NUCMASS * temp * Sqrt[t] / eb^2 - (q^2 + 2 * NUCMASS * q0 - (GP - 2) * eb * (spr - s) / 2) / (8 * NUCMASS * temp * (qz^2 + qpt^2))];
prefactor = Sqrt[Pi * NUCMASS / (2 * temp)] * 1 / (1 - t) * qpt / (qz^2 + qpt^2);
prefactor * expterm * sfunc];
npnorm[eb_, t_]:=Module[{mt},
mt = t * NUCMASS;
Sinh[eb / (2 * mt)] / Cosh[GP * eb / (4 * mt)] * (2 * Pi / mt)^(3 / 2) * mt / eb];

protkern[eb_, t_, rhop_, q0_, q_, cosq_, s_, spr_]:= Module[{mt, ebmt, npnorm},
mt = NUCMASS * t;
ebmt = eb / mt;
npnorm = Sinh[ebmt / 2] / Cosh[GP * ebmt / 4] * (2 * Pi)^(3/2) / (eb * Sqrt[mt]);
eb / (2 * Pi) * rhop * npnorm * (ffunc[eb, t, q0, q, cosq, s, spr] - rhop * npnorm * ffunc[eb, t / 2, q0, q, cosq, s, spr] * Exp[q0 / t])];

neutkern[mun_, t_, q0_, q_]:= If[q0 == 0, 0, If[q0 > 0, NUCMASS^2 * t / (2 * Pi * q) * (Log[1 + Exp[mun / t]] - Log[1 + Exp[(mun - q0) / t]]) / (Exp[q0 / t] - 1),
	NUCMASS^2 * t / (2 * Pi * q) * (Log[1 + Exp[(mun + q0) / t]] - Log[1 + Exp[mun / t]]) / (Exp[q0 / t] - 1)]];
	
eleckern[mue_, eb_, t_, q0_, q_, cosq_]:= If[q0 <= q, 0, eb / (2 * Pi) * t * Exp[-q^2 * (1 - cosq^2) / (2 * eb)] * (Log[1 + Exp[mue / t]] - Log[1 + Exp[(mue - q0) / t]]) / (Exp[q0 / t] - 1)];

(*pretty sure these aren't used anymore -MLK 12-17-25 *)
alpha[delta_, knu_]:=If[delta > 0, 4 * knu^2 + 12 * knu * delta + 6 * delta^2, If[-knu < delta, 4 * (knu + delta)^3 / knu, 0]];
beta[delta_, knu_]:= If[delta > 0, knu^2 + 4 * knu * delta, If[-knu < delta, (knu + delta)^4 / knu^2, 0]];
asspr[s_, spr_, cost_]:=If[s == spr, 1 / 2 * (1 + GA^2 + 2 * GA * s * (1 - 4 * SINTW^2) * cost), GA^2];
bsspr[s_, spr_, cost_]:=If[s == spr, 1 / 2 * (1 - GA^2 + 2 * GA^2 * cost + 2 * GA * s * (1 - 4 * SINTW^2) * cost), GA^2 * (-1 + 2 * cost)];

(*Matrix elements for charged current*)
(*note these do not match the normalization convention in the paper but that factor is cancelled in the final kappa expression*)
mredcc[sp_, sn_, nezero_, cost_]:=If[nezero, If[sn == 1, If[sp == 1, 4 * (1 + GA)^2 * (1 + cost), 0], If[sp == 1, 16 * GA^2 * (1 - cost), 4 * (1 - GA)^2 * (1 + cost)]],
	If[sn == 1, If[sp == 1, 2 * (1 + GA)^2 * (1 + cost) + 2 * (1 - GA)^2 * (1 - cost), 8 * GA^2 * (1 + cost)], If[sp == 1, 8 * GA^2 * (1 - cost),
	2 * (1 + GA)^2 * (1 - cost) + 2 * (1 - GA)^2 * (1 + cost)]]];
vtildecc[x_, sp_, sn_, cost_]:=If[x <= 0, 0, If[x < 1 / 2, mredcc[sp, sn, True, cost], mredcc[sp, sn, True, cost] + (2 * x - 1) * mredcc[sp, sn, False, cost]]];

(*Matrix elements for neutral current*)
mredncp[s_, spr_, cost_, costpr_, phi_]:= Module[{sint, sintpr},
sint = Sqrt[1 - cost^2];
sintpr = Sqrt[1 - costpr^2];
If[s == spr, 1 / 2 * ((1 - 4 * SINTW^2)^2 * (1 + cost * costpr + sint * sintpr * Cos[phi]) + GA^2 * (1 + cost * costpr - sint * sintpr * Cos[phi]) - 2 * GA * s * (1 - 4 * SINTW^2) * (cost + costpr)),
	GA^2 * (1 - cost * costpr)]];
mredncpint[s_, spr_, cost_]:=If[s == spr, 1/2 * ((1 - 4 * SINTW^2)^2 + GA^2 
	- 2 (1 - 4 * SINTW^2)* GA * s * cost), GA^2];
mredncn[s_, spr_, cost_, costpr_, phi_]:= Module[{sint, sintpr},
sint = Sqrt[1 - cost^2];
sintpr = Sqrt[1 - costpr^2];
If[s == spr, 1 / 2 * ((1 + cost * costpr + sint * sintpr * Cos[phi]) + GA^2 * (1 + cost * costpr - sint * sintpr * Cos[phi]) - 2 * GA * s * (cost + costpr)),
	GA^2 * (1 - cost * costpr)]];
mredncnint[s_, spr_, cost_]:=If[s == spr, 1/2 * (1 + GA^2 - 2 * GA * s * cost), GA^2];
mredncintexp[s_, spr_, cost_, knu_, t_]:= Module[{mt, u, iu},
mt = NUCMASS * t;
u = knu^2 / (2 * mt);
iu = Coth[u] - 1 / u;
1 / (4 * u) * (1 - Exp[-u]) * If[s == spr, (1 - iu + GA^2 * (1 + (2 * cost^2 - 1) * iu)) - 2 * GA * s * cost * (1 + iu),
	GA^2 * (2 - 2 * cost * iu)]];

(* mrednce is for electron type neutrinos, mredncex is for muon and tau neutrinos *)	
mrednce[h_, hpr_, cost_, costpr_]:=If[h != hpr, 0, If[h == 1, 4 * SINTW^4 * (1 - cost) * (1 - costpr), (1 + 2 * SINTW^2)^2 * (1 + cost) * (1 + costpr)]];
mredncex[h_, hpr_, cost_, costpr_]:=If[h != hpr, 0, If[h == 1, 4 * SINTW^4 * (1 - cost) * (1 - costpr), (1 - 2 * SINTW^2)^2 * (1 + cost) * (1 + costpr)]];

(*just the helicity weights for electron scattering since 1-hc can be absorbed into the kinematics*)
mredncestrip[h_, hpr_]:=If[h != hpr, 0, If[h == 1, 4 * SINTW^4, (1 + 2 * SINTW^2)^2]];
mredncexstrip[h_, hpr_]:=If[h != hpr, 0, If[h == 1, 4 * SINTW^4, (1 - 2 * SINTW^2)^2]];

(*Opacities - analytical approx*)
(*neutron charged current - nondegenerate*)
kappan[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sp, sn, pdenom},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb / (8 * Pi * Cosh[GN * ebmt / 4]);
pdenom = ((1 + Exp[-eb / mt]) * Cosh[GP * eb / (4 * mt)] + (1 - Exp[-eb / mt]) * Sinh[(GP - 2) * eb / (4 * mt)]);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t]) * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * Exp[GN * sn * ebmt / 4]
	 * (1 - nb * yp * Exp[(GP - 2) * sp * ebmt / 4] / pdenom * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, 1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, 1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

(*neutron charged current - degenerate*)
kappandegen[eb_, t_, mue_, mun_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor, pdenom},
mt = NUCMASS * t;
kfn = Sqrt[2 * mun * NUCMASS];
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb / (8 * Pi);
spinsum = Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t]) * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * (1
	- Exp[(GP - 2) * sp * ebmt / 4]/(2 * Cosh[GP * ebmt / 4]) * (1 - Exp[-ebmt]) * bigb[Sqrt[(kfn^2 - GN * sn * eb / 2)], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8], eb, t, knu, cost]),
	{sp, {-1, 1}}, {sn, {-1, 1}}];
prefactor * spinsum];

(*proton charged current - nondegenerate neutrons*)
kappap[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sn, sp, pdenom},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * yp * eb * Exp[eb / (2 * NUCMASS * t)] / (8 * Pi * Cosh[GP * eb / (4 * mt)]);
spinsum = Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], -mue, t]) * (1 - thetapm[knu, sp, sn, eb, ui, 1]) * If[sp == sn, 
	vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost],
	vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), -sp, -sn, cost]] * Exp[(GP - 2) * sp * ebmt / 4] 
	 * (1 - nb * (1 - yp) * Exp[GN * sn * ebmt / 4] / Cosh[GN * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, -1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, -1]^2)/(4*mt)]), {sp, {-1, 1}}, {sn, {-1, 1}}];
prefactor * spinsum];

(*proton charged current - degenerate neutrons *)
kappapdegen[eb_, t_, mue_, mun_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor, pdenom},
mt = NUCMASS * t;
kfn = Sqrt[2 * mun * NUCMASS];
ebmt = eb / mt;
prefactor = nb * yp * eb * Exp[ebmt / 2] / (8 * Pi * Cosh[GP * ebmt / 4]);
spinsum = Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], -mue, t]) * (1 - thetapm[knu, sp, sn, eb, ui, 1]) * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * Exp[(GP - 2) * sp * ebmt / 4] * (1
	- (1 - Exp[-ebmt]) * bigb[Sqrt[kfn^2 - GN * sn * eb / 2], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8], eb, t, knu, cost]),
	{sp, {-1, 1}}, {sn, {-1, 1}}];
prefactor * spinsum];

(*Differential emissivities - analytic approx*)
epsn[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:=Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sp, sn},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb * knu^3 / (64 * Pi^4 * Cosh[GN * ebmt / 4]);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], mue, t]) * thetapm[knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * Exp[GN * sn * ebmt / 4]
	 * (1 - nb * yp * Exp[(GP - 2) * sp * ebmt / 4] / Cosh[GP * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, -1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, -1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

epsp[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:=Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sp, sn},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * yp * eb * knu^3 / (64 * Pi^4 * Cosh[GP * ebmt / 4]);
spinsum = Sum[Sum[nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t] * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * Exp[(GP - 2) * sp * ebmt / 4]
	 * (1 - nb * (1 - yp) * Exp[GN * sn * ebmt / 4] / Cosh[GN * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, 1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, 1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

epsndegen[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor},
mt = NUCMASS * t;
kfn = (3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3)));
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb * knu^3 / (64 * Pi^4);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], mue, t]) * thetapm[knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * (1
	- Exp[(GP - 2) * sp * ebmt / 4]/(2 * Cosh[GP * ebmt / 4]) * (1 - Exp[-ebmt]) * bigb[Sqrt[kfn^2 - GN * sn * eb / 2], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8]]),
	{sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

epspdegen[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor},
mt = NUCMASS * t;
kfn = (3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3)));
ebmt = eb / mt;
prefactor = nb * yp * eb * knu^3 / (64 * Pi^4 * Cosh[GP * ebmt / 4]);
spinsum = Sum[Sum[nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t] * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * (1
	- (1 - Exp[-ebmt]) * bigb[Sqrt[kfn^2 - GN * sn * eb / 2], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8]]),
	{sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

(*Double diff synchroton rate - analytical approx*)
synchp[eb_, t_, nb_, yp_, knu_, knupr_, cost_, costpr_, phi_]:=Module[{s, spr, q, cosq},
cosq = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 + 2 * knu * knupr * cosq];
knu^2 * knupr^2 / (2 * Pi)^6 * Sum[Sum[mredncp[s, spr, cost, costpr, phi] * protkern[eb, t, nb * yp, knu + knupr, q, cosq, s, spr],{s, {-1, 1}}], {spr, {-1, 1}}]];

synchn[eb_, t_, nb_, yp_, knu_, knupr_, cost_, costpr_, phi_]:=Module[{q, cosq, mun, q0},
cosq = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 + 2 * knu * knupr * cosq];
mun = ((3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3))))^2 / (2 * NUCMASS);
q0 = -GN * eb / (2 * NUCMASS);
If[q0 > knu + knupr, knu^2 * knupr^2 / (2 * Pi)^6 * mredncn[-1, 1, cost, costpr, phi] * neutkern[mun, t, q0, q], 0]];
	
(*Elastic scattering opacity - analytic approx*)

deltasspr[eb_, g_, s_, spr_]:= - g * (s - spr) * eb / (4 * NUCMASS);
k0func[eb_, g_, s_, spr_, q0_, q_]:= Abs[NUCMASS / q * (q0 - deltasspr[eb, g, s, spr])];

kappanncdiff[eb_, t_, n_, knu_, cost_, knupr_, costpr_, phi_]:= Module[{prefactor, mt, q, q0, cosnunupr, spinsum},
cosnunupr = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosnunupr];
q0 = knupr - knu;
mt = NUCMASS * t;
prefactor = knupr^2 * n / (4 * (2 * Pi)^(5/2) * Cosh[GN * eb / (4 * mt)]) * Sqrt[NUCMASS / (t * q^2)];
spinsum = Sum[mredncn[s, spr, cost, costpr, phi] * Exp[GN * s * eb / (4 * mt)] * (Exp[-k0func[eb, GN, s, spr, q0, q]^2 / (2 * mt)]
	- n * Exp[GN * spr * eb / (4 * mt)] / (Sqrt[2] * Cosh[GN * eb / (4 * mt)]) * Pi^2 / (mt * q) * Exp[-q^2 / (4 * mt)]
	* (Erf[(q + 2 * k0func[eb, GN, s, spr, q0, q]) / (2 * Sqrt[mt])] + Erf[(q - 2 * k0func[eb, GN, s, spr, q0, q]) / (2 * Sqrt[mt])])),
	{s, {-1, 1}}, {spr, {-1, 1}}];
prefactor * spinsum];

kappannc[eb_, t_, n_, knu_, cost_]:= Module[{prefactor, spinsum, mt, dss},
mt = NUCMASS * t;
prefactor = n / (4 * Pi * Cosh[GN * eb / (4 * mt)]); 
spinsum = Sum[dss = deltasspr[eb, GN, s, spr];
	(knu + dss)^2 * Exp[GN * s * eb / (4 * mt)] * 
	(mredncnint[s, spr, cost] - n * Exp[GN * spr * eb / (4 * mt)] / Cosh[GN * eb / (4 * mt)]
	* (Pi / mt)^(3/2) * mredncintexp[s, spr, cost, knu, t]), {s, {-1, 1}}, {spr, {-1, 1}}];
prefactor * spinsum];

kappanncdiffdegen[eb_, t_, mu_, knu_, cost_, knupr_, costpr_, phi_]:=Module[{q0, q, cosnunupr, e0, prefactor, spinsum},
cosnunupr = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosnunupr];
q0 = knupr - knu;
prefactor = knupr^2 * NUCMASS^2 / (64 * Pi^5);
spinsum = Sum[e0 = k0func[eb, GN, s, spr, q0, q]^2 / (2 * NUCMASS) - GN * eb * s / (4 * NUCMASS);
mredncn[s, spr, cost, costpr, phi] * (q0 + t * 
Log[(1 + Exp[(e0 + q0 - mu) / t]) / (1 + Exp[(e0 - mu) / t])]) / (Exp[q0 / t] - 1),{s, {-1, 1}}, {spr, {-1, 1}}];
prefactor * spinsum];

kappanncdegen[eb_, t_, mu_, knu_, cost_]:=Module[{prefactor, spinsum},
prefactor = NUCMASS / (16 * Pi^4) * (2 * Pi^4 * t^3 / 15 + Pi^2 * t * knu^2 / 3 + 8 * t^2 * knu * Zeta[3]);
spinsum = Sum[mredncnint[s, spr, cost], {s, {-1, 1}}, {spr, {-1, 1}}];
prefactor * spinsum];
	
kappapncdiff[eb_, t_, n_, knu_, cost_, knupr_, costpr_, phi_]:= Module[
{prefactor, mt, q, q0, cosq, qz, qperp, cosnunupr, spinsum, exponent, zpr, z, 
polaravg, polaravgpr, texp, blockingpref, abar},
cosnunupr = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosnunupr];
q0 = knupr - knu;
qz = knupr * costpr - knu * cost;
cosq = qz / q;
qperp = Sqrt[q^2 - qz^2];
mt = NUCMASS * t;
texp = Exp[-eb / mt];
z = qperp^2 * Sqrt[texp] / (eb * (1 - texp));
zpr = qperp^2 * texp^(1/4) / (eb * (1 - Sqrt[texp]));
prefactor = (knupr^2 * n) / (4 * (2 * Pi)^(5 / 2) * Cosh[GP * eb / (4 * mt)] * Sqrt[texp]) * Sqrt[NUCMASS / (t * qz^2)];
spinsum = Sum[abar = Round[NUCMASS * Abs[q0 - deltasspr[eb, GP - 2, s, spr]] / eb];
exponent = - (NUCMASS * q0 + (GP - 2) * (s - spr) * eb / 4)^2 / (2 * qz^2 * mt);
polaravg = cosq^2 * Exp[exponent] + (1 - cosq^2) * BesselI[abar, z] * Exp[-qperp^2 / (2 * eb) * (1 + t) / (1 - t)];
polaravgpr = cosq^2 * Exp[2 * exponent] + (1 - cosq^2) * BesselI[abar, zpr] * Exp[-qperp^2 / (2 * eb) * (1 + Sqrt[t]) / (1 - Sqrt[t])];
blockingpref = n * Sinh[eb / (2 * mt)] * Exp[q0 / t] / (Cosh[GP * eb / (4 * mt)]) * 2 * Pi^(3/2) / (eb * Sqrt[mt]);
mredncp[s, spr, cost, costpr, phi] * (Exp[(GP - 2) * s * eb / (4 * mt)] * polaravg
- Exp[(GP - 2) * s * eb / (2 * mt)] * blockingpref * polaravgpr),
{s, {-1, 1}}, {spr, {-1, 1}}];
prefactor * spinsum];

kappapnc[eb_, t_, n_, knu_, cost_]:= Module[{mt, prefactor, spinsum, texp, sumpref},
mt = NUCMASS * t;
texp = Exp[-eb / mt];
prefactor = n / (4 * Pi * Cosh[GP * eb / (4 * mt)] * Sqrt[texp]);
spinsum = Sum[sumpref = mredncpint[s, spr, cost] * (knu + deltasspr[eb, GP - 2, s, spr])^2 * Exp[(GP - 2) * eb / (4 * mt)];
sumpref * (1 - n * Sinh[eb / (2 * mt)] / Cosh[GP * eb / (4 * mt)] * (2 * Pi)^(3/2) 
/ (eb * Sqrt[mt]) * Exp[(GP - 2) * spr * eb / (4 * mt)]), {s, {-1, 1}}, {spr, {-1, 1}}];
prefactor * spinsum];

jfunc[knu_, mue_, t_, hc_, eb_]:=Module[{kbar, mubar, tbar, ifunc, q0},
kbar = knu / Sqrt[eb];
mubar = mue / Sqrt[eb];
tbar = t / Sqrt[eb];
ifunc[qbar_]:= If[qbar > 0, Log[(1 + Exp[mubar / tbar]) / (1 + Exp[(mubar - qbar) / tbar])],
	Log[(1 + Exp[(mubar + 2 * qbar) / tbar]) / (1 + Exp[(mubar + qbar) / tbar])]] / (Exp[qbar / tbar] - 1);
NIntegrate[q0 = (kperp^2 - kbar^2 * (1 - hc^2)) / (2 * kbar * (1 + hc)); 
	kperp * (1 - hc) * (1 - (hc * kbar - q0) / Sqrt[kperp^2 + (hc * kbar - q0)^2])
	* Exp[-kperp^2 / 2] * BesselI[0, kperp * kbar * Sqrt[1 - hc^2]] 
	* ifunc[q0], {kperp, 0, Infinity}]
];

kappaenc[eb_, t_, mue_, knu_, cost_, ise_]:=Module[{prefactor, rest, h, sgn},
sgn = If[ise, 1, -1];
prefactor = eb^2 * t / (8 * Pi^3) * Exp[-knu^2 * (1 - cost^2) / (2 * eb)];
rest = 4 * SINTW^4 * jfunc[knu, mue, t, cost, eb] + (1 + 2 * sgn * SINTW^2)^2 * jfunc[knu, mue, t, -cost, eb];
prefactor * rest];


(*Opacities - Monte Carlo*)

(*generated precompiled inr*)
inrgenerator[{n1_, n2_}]:=Compile[
	{{x, _Real}},
	Module[{in, inmin, innew, inminnew, rtemp, n, r},
	n = Max[{n1, n2}];
	r = Min[{n1, n2}];
	in = Sqrt[x^n / Factorial[n]] * Exp[-x / 2];
	inmin = Sqrt[x^(n - 1) / Factorial[n - 1]] * Exp[-x / 2];
	For[rtemp = 1, rtemp <= r, rtemp++,
	innew = - Sqrt[x / rtemp] * in + Sqrt[n / rtemp] * inmin;
	inminnew = Sqrt[n / x] * (innew - Sqrt[rtemp / n] * inmin);
	in = innew;
	inmin = inminnew;];
	in],
	CompilationTarget -> "C",
	RuntimeOptions -> "Speed"
];

(*does monte carlo sum over landau levels*)
llmc[integrandfunc_, eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, nsamples_, isn_]:=Module[{vals, nemax, npmax, flatweights, flatpairs, probs, samples, subfunc, weightfunc},
nemax = Floor[If[isn, (Max[{2 * t, knu, mue + 2 * t}] + ui + MSPLIT + (GP - 2 - GN) * eb / (4 * NUCMASS))^2 / (2 * eb), (Max[{2 * t, knu, mue + 2 * t}] - ui - MSPLIT + (GP - 2 - GN) * eb / (4 * NUCMASS))^2 / (2 * eb)]];
npmax = Floor[If[mup > 0, (mup + 5 * t) * NUCMASS / eb, 5 * t * NUCMASS / eb]];
weightfunc[np_, ne_]:=Exp[-np * eb / (t * NUCMASS)] / (nemax + 1) * (1 - Exp[- eb / (t *  NUCMASS)]);
subfunc[np_, ne_]:= integrandfunc[eb, t, mue, mun, mup, knu, cost, ui, ne, np]  / weightfunc[np, ne];
flatweights = Flatten[Table[weightfunc[np, ne], {np, 0, npmax}, {ne, 0, nemax}]];
flatpairs = Flatten[Table[{np, ne}, {np, 0, npmax}, {ne, 0, nemax}], 1];
probs = flatweights / Total[flatweights];
samples = RandomChoice[probs -> flatpairs, nsamples];
vals = (subfunc @@@ samples);
Total[vals] / nsamples];

llmcncdiff[integrandfunc_, eb_, t_, mu_, knu_, cost_, knupr_, costpr_, phi_, nsamples_, isprot_]:=Module[{nmax, rangen, vals, flatweights, flatpairs, probs, samples, subfunc, weightfunc, weightnorm, inrsamples},
nmax = Floor[If[mu > 0, If[isprot,(mu + 5 * t) * NUCMASS / (2 * eb), (mu + 5 * t)^2 / (2 * eb) + 1], If[isprot, (5 * NUCMASS * t) / (2 * eb), (5 * t)^2 / (2 * eb) + 1]]];
rangen = Range[0, nmax];
weightfunc[n_, npr_]:=Exp[- (n + Abs[n - npr]) * eb / (t * NUCMASS)];
weightnorm =  (1 - Exp[- eb / (t *  NUCMASS)]) * (1 - Exp[- 2 * eb / (t * NUCMASS)]) / 
	(1 + Exp[-eb / (t * NUCMASS)] + Exp[- 2 * eb / (t * NUCMASS)]);
subfunc[{n_, npr_, inr_}]:= integrandfunc[eb, t, mu, knu, cost, knupr, costpr, phi, n, npr, inr] / (weightfunc[n, npr] * weightnorm);
flatweights = Flatten[Table[weightfunc[n, npr], {n, rangen}, {npr, rangen}]];
flatpairs = Flatten[Table[{n, npr}, {n, rangen}, {npr, rangen}], 1];
(*probs = flatweights / Total[flatweights];*)
samples = RandomChoice[flatweights -> flatpairs, nsamples];
inrsamples = ParallelMap[inrgenerator, samples];
vals = ParallelMap[subfunc, Join[samples, List /@ inrsamples, 2]];
{Total[vals] / nsamples, StandardDeviation[vals]}];

llmcnc[integrandfunc_, eb_, t_, mu_, knu_, cost_, nsamples_, isprot_]:=Module[{nmax, rangen, vals, flatweights, flatpairs, probs, samples, subfunc, weightfunc, weightnorm, inrsamples},
nmax = Floor[If[mu > 0, (mu + 5 * t)^2 / (2 * eb), If[isprot, (5 * NUCMASS * t) / (2 * eb), (5 * t)^2 / (2 * eb)]]];
rangen = Range[0, nmax];
weightfunc[n_, npr_]:=Exp[- (n + Abs[n - npr]) * eb / (t * NUCMASS)];
weightnorm = (1 - Exp[- eb / (t *  NUCMASS)]) * (1 - Exp[- 2 * eb / (t * NUCMASS)]) / 
	(1 + Exp[-eb / (t * NUCMASS)] + Exp[- 2 * eb / (t * NUCMASS)]);
subfunc[{n_, npr_, inr_}]:= integrandfunc[eb, t, mu, knu, cost, n, npr, inr] / (weightfunc[n, npr] * weightnorm);
flatweights = Flatten[Table[weightfunc[n, npr], {n, rangen}, {npr, rangen}]];
flatpairs = Flatten[Table[{n, npr}, {n, rangen}, {npr, rangen}], 1];
(*probs = flatweights / Total[flatweights];*)
samples = RandomChoice[flatweights -> flatpairs, nsamples];
inrsamples = ParallelMap[inrgenerator, samples];
vals = ParallelMap[subfunc, Join[samples, List /@ inrsamples, 2]];
Total[vals] / nsamples];

(*Calculates inn^2 for monte carlo integration*)
nmcme[mun_, mup_, mue_, kn_, cosn_, phin_, knu_, cost_, sn_, sp_, np_, ne_, eb_, t_, ui_, cose_, elec_]:= Module[{wperpsq, nezero, inn, innpmin, innemin, innepmin},
wperpsq = kn^2 * (1 - cosn^2) + knu^2 * (1 - cost^2) - 2 * kn * knu * Sqrt[1 - cosn^2] * Sqrt[1 - cost^2] * Cos[phin];
nezero = If[ne == 0, True, False];
inn = If[ne <= np, (wperpsq / (2 * eb))^(np - ne) * Exp[-wperpsq / (2 * eb)] * Factorial[ne] / Factorial[np] * LaguerreL[ne, np - ne, wperpsq / (2 * eb)]^2,
	(wperpsq / (2 * eb))^(ne - np) * Exp[-wperpsq / (2 * eb)] * Factorial[np] / Factorial[ne] * LaguerreL[np, ne - np, wperpsq / (2 * eb)]^2];
innpmin = If[np == 0, 0, If[ne <= np - 1, (wperpsq / (2 * eb))^(np - 1 - ne) * Exp[-wperpsq / (2 * eb)] * Factorial[ne] / Factorial[np - 1] * LaguerreL[ne, np - 1 - ne, wperpsq / (2 * eb)]^2,
	(wperpsq / (2 * eb))^(ne - np + 1) * Exp[-wperpsq / (2 * eb)] * Factorial[np - 1] / Factorial[ne] * LaguerreL[np - 1, ne - np + 1, wperpsq / (2 * eb)]^2]];
innemin = If[ne == 0, 0, If[ne <= np + 1, (wperpsq / (2 * eb))^(np + 1 - ne) * Exp[-wperpsq / (2 * eb)] * Factorial[ne - 1] / Factorial[np] * LaguerreL[ne - 1, np + 1 - ne, wperpsq / (2 * eb)]^2,
	(wperpsq / (2 * eb))^(ne - np - 1) * Exp[-wperpsq / (2 * eb)] * Factorial[np] / Factorial[ne - 1] * LaguerreL[np, ne - np - 1, wperpsq / (2 * eb)]^2]];
innepmin = If[ne == 0 || np == 0, 0, If[ne <= np, (wperpsq / (2 * eb))^(np - ne) * Exp[-wperpsq / (2 * eb)] * Factorial[ne - 1] / Factorial[np - 1] * LaguerreL[ne - 1, np - ne, wperpsq / (2 * eb)]^2,
	(wperpsq / (2 * eb))^(ne - np) * Exp[-wperpsq / (2 * eb)] * Factorial[np - 1] / Factorial[ne - 1] * LaguerreL[np - 1, ne - np, wperpsq / (2 * eb)]^2]];
If[sn == 1, If[sp == 1, 2 * (1 + elec * GA)^2 * (1 + cose) * (1 + cost) * inn + 2 * (1 - elec * GA)^2 * (1 - cose) * (1 - cost) * innemin, 8 * GA^2 * (1 - cose) * (1 + cost) * innepmin],
	If[sp == 1, 8 * GA^2 * (1 + cose) * (1 - cost) * inn, 2 * (1 + elec * GA)^2 * (1 - cose) * (1 - cost) * innepmin + 2 * (1 - elec * GA)^2 * (1 + cose) * (1 + cost) * innpmin]]];

(*does a spin sum with charged current matrix elements*)
nmcspinsum[mun_, mup_, mue_, kn_, cosn_, phin_, knu_, cost_, np_, ne_, eb_, t_, ui_, sp_, sn_]:= Module[{kzsoltemp, ep, en, ee, eppl, epmin, kze},
ep = (kn * cosn)^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * sp * eb / (4 * NUCMASS);
en = kn^2 / (2 * NUCMASS) - GN * sn * eb / (4 * NUCMASS);
ee = en + knu - ep + ui + MSPLIT;
kze = Sqrt[ee^2 - 2 * ne * eb - ELECMASS^2];
eppl = (kn * cosn + knu * cost + kze)^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * sp * eb / (4 * NUCMASS);
epmin = (kn * cosn + knu * cost - kze)^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * sp * eb / (4 * NUCMASS);
If[ee < Sqrt[2 * ne * eb], 0, 1 / 2 * ee * kn^2 / Sqrt[ee^2 - 2 * ne * eb]
	 * eb / (2 * Pi)^4 * nfd[en, mun, t] * (1 - nfd[ee, mue, t]) * ((1 - nfd[epmin, mup, t]) * nmcme[mun, mup, mue, kn, cosn, phin, knu, cost, sn, sp, np, ne, eb, t, ui, kze / ee, 1]
	 + (1 - nfd[eppl, mup, t]) * nmcme[mun, mup, mue, kn, cosn, phin, knu, cost, sn, sp, np, ne, eb, t, ui, -kze / ee, 1])]];

(*opacity for capture on neutrons with specified landau levels*)
kappanmcnenp[eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, ne_, np_]:= Module[{knmax},
knmax = If[mun > 0, Sqrt[2 * (mun + 5 * t) * NUCMASS], Sqrt[2 * (5 * t) * NUCMASS]];
Quiet[NIntegrate[Sum[nmcspinsum[mun, mup, mue, kn, cosn, phin, knu, cost, np, ne, eb, t, ui, sp, sn], {sp, {-1, 1}}, {sn, {-1, 1}}], 
	 {kn, 0, knmax}, {cosn, -1, 1}, {phin, 0, 2 * Pi}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}]]];
(*total opacity for capture on neutrons*)	
kappanmc[{eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, nsamples_}]:= llmc[kappanmcnenp, eb, t, mue, mun, mup, knu, cost, ui, nsamples, True];

(*integrand for capture on protons*)
pmcspinsum[mun_, mup_, mue_, kn_, cosn_, phin_, knu_, cost_, np_, ne_, eb_, t_, ui_, sp_, sn_]:= Module[{helfix, kzsoltemp, ep, en, ee, kzp, kze, eppl, epmin},
kzp = kn * cosn (*kzpx/.kzsoltemp;*);
ep = Sqrt[NUCMASS^2 + kzp^2 + 2 * np * eb - (GP - 2) * sp * eb / 2] - NUCMASS;
en = Sqrt[NUCMASS^2 + kn^2] - GN * sn * eb / (4 * NUCMASS) - NUCMASS;
ee = ep + knu - en - ui - MSPLIT;
kze = Sqrt[ee^2 - 2 * ne * eb];
eppl = (kn * cosn + knu * cost + kze)^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * sp * eb / (4 * NUCMASS);
epmin = (kn * cosn + knu * cost - kze)^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * sp * eb / (4 * NUCMASS);
helfix = If[sp == sn, 1, -1];
(*note that since the positron has helicity + in the weak interaction, the sign for kze / ee in the matrix element is switched*)
If[ee < Sqrt[2 * ne * eb], 0, 1 / 2 * ee * kn^2 / Sqrt[ee^2 - 2 * ne * eb]
	 * eb / (2 * Pi)^4 * (1 - nfd[en, mun, t]) * (1 - nfd[ee, -mue, t]) * (nfd[epmin, mup, t] * nmcme[mun, mup, mue, kn, cosn, phin, knu, cost, sn * helfix, sp * helfix, np, ne, eb, t, ui, kze / ee, -1]
	 + nfd[eppl, mup, t] * nmcme[mun, mup, mue, kn, cosn, phin, knu, cost, sn * helfix, sp * helfix, np, ne, eb, t, ui, -kze / ee, -1])]];

(*opacity for capture on protons with specified landau levels*)
kappapmcnenp[eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, ne_, np_]:= Module[{knmax},
knmax = If[mun > 0, Sqrt[2 * (mun + 20 * t) * NUCMASS], Sqrt[2 * (20 * t) * NUCMASS]];
Quiet[NIntegrate[Sum[pmcspinsum[mun, mup, mue, kn, cosn, phin, knu, cost, np, ne, eb, t, ui, sp, sn], {sp, {-1, 1}}, {sn, {-1, 1}}], 
	 {kn, 0, knmax}, {cosn, -1, 1}, {phin, 0, 2 * Pi}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}]]];
(*total opacity for capture on protons*)	
kappapmc[{eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, nsamples_}]:= llmc[kappapmcnenp, eb, t, mue, mun, mup, knu, cost, ui, nsamples, False];
	
(*elastic opacities*)
(*integrand for scattering on protons*)
kappapncintegrand[eb_, t_, mup_, knu_, cost_, knupr_, costpr_, phi_, s_, spr_, n_, npr_, inr_]:= Module[{qz, ep, eppr, kz, wperpsq},
qz = knupr * costpr - knu * cost;
kz = (knupr - knu + qz^2 / (2 * NUCMASS) + (npr - n) * eb / NUCMASS + (GP - 2) * eb / (4 * NUCMASS) * (s - spr)) * NUCMASS / qz;
ep = kz^2 / (2 * NUCMASS) + n * eb / NUCMASS - (GP - 2) * s * eb / (4 * NUCMASS);
eppr = (kz - qz)^2 / (2 * NUCMASS) + npr * eb / NUCMASS - (GP - 2) * spr * eb / (4 * NUCMASS);
wperpsq = knu^2 * (1 - cost^2) + knupr^2 * (1 - costpr^2) - 2 * knu * knupr * Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
eb * knupr^2 * NUCMASS / (Abs[qz] * (2 * Pi)^4) * mredncp[s, spr, cost, costpr, phi] * nfd[ep, mup, t] * (1 - nfd[eppr, mup, t]) * inr[wperpsq / (2 * eb)]^2];

(*opacity for in one LL for proton scattering*)
kappapncdiffnnpr[eb_, t_, mup_, knu_, cost_, knupr_, costpr_, phi_, n_, npr_, inr_]:=Sum[kappapncintegrand[eb, t, mup, knu, cost, knupr, costpr, phi, s, spr, n, npr, inr], 
	{s, {-1, 1}}, {spr, {-1, 1}}]
kappapncnnpr[eb_, t_, mup_, knu_, cost_, n_, npr_, inr_]:=NIntegrate[Sum[kappapncintegrand[eb, t, mup, knu, cost, knupr, costpr, phi, s, spr, n, npr, inr], 
	{s, {-1, 1}}, {spr, {-1, 1}}], {costpr, -1, 1}, {phi, 0, 2 * Pi}, {knupr, 0, 20 * t}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}];
	
(*full opacity for scattering on protons*)
kappapmcncdiff[{eb_, t_, mup_, knu_, cost_, knupr_, costpr_, phi_, nsamples_}]:=If[mup < 0, 
	llmcncdiff[kappapncdiffnnpr, eb, t, mup, knu, cost, knupr, costpr, phi, nsamples, True],
	-1];
kappapmcnc[{eb_, t_, mup_, knu_, cost_, nsamples_}]:= If[mup < 0,
	llmcnc[kappapncnnpr, eb, t, mup, knu, cost, nsamples, True],
	-1];

(*diff opacity for neutrons*)
kappanncintegrand[eb_, t_, mun_, knu_, cost_, knupr_, costpr_, phi_, en_, s_, spr_]:=Module[{q, enpr, kn, knpr, cosnunupr, tritest},
cosnunupr = cost * costpr + Sqrt[(1 - cost^2) * (1 - costpr^2)] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosnunupr];
kn = Sqrt[Max[{0, 2 * NUCMASS * (en + GN * eb * s / (4 * NUCMASS))}]];
enpr = en + knu - knupr;
knpr = Sqrt[Max[{0, 2 * NUCMASS * (enpr + GN * eb * spr / (4 * NUCMASS))}]];
tritest = If[kn + knpr > q && kn + q > knpr && knpr + q > kn, 1, 0];
tritest * knupr^2 * NUCMASS^2 / (2 * q * (2 * Pi)^4) * mredncn[s, spr, cost, costpr, phi] * nfd[en, mun, t] * (1 - nfd[enpr, mun, t])];

kappanmcncdiff[{eb_, t_, mun_, knu_, cost_, knupr_, costpr_, phi_}]:=NIntegrate[Sum[kappanncintegrand[eb, t, mun, knu, cost, knupr, costpr, phi, en, s, spr], {s, {-1, 1}}, {spr, {-1, 1}}],
	{en, Max[{0, knupr - knu}], Infinity}];
kappanmcnc[{eb_, t_, mun_, knu_, cost_}]:=NIntegrate[Sum[kappanncintegrand[eb, t, mun, knu, cost, knupr, costpr, phi, en, s, spr], {s, {-1, 1}}, {spr, {-1, 1}}],
	{en, 0, If[mun > 0, mun + 5 * t, 5 * t]}, {costpr, -1, 1}, {phi, 0, 2 * Pi}, {knupr, 0, en + knu}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}];
	
(* opacity for scattering on electrons *)
kappaencintegrand[eb_, t_, mue_, knu_, cost_, costpr_, phi_, h_, iselec_]:= Module[{knupr, qperp, q0, q, qz, cosnunupr, prefactor, helsum, ifunc},
prefactor = eb * t / (32 * Pi^4);
knupr = Max[{0, knu * (h + cost) / (h + costpr)}];
q0 = knupr - knu;
cosnunupr = cost * costpr + Sqrt[(1 - cost^2) * (1 - costpr^2)] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosnunupr];
qz = knupr * costpr - knu * cost;
qperp = Sqrt[q^2 - qz^2];
ifunc = If[q0 > 0, Log[(1 + Exp[mue / t]) / (1 + Exp[(mue - q0) / t])] / (Exp[q0 / t] - 1),
	Log[(1 + Exp[(mue + 2 * q0) / t])/ (1 + Exp[(mue + q0) / t])] / (Exp[q0 / t] - 1)];
helsum = knupr^2 / Abs[h + costpr] * ifunc * Exp[-qperp^2 / (2 * eb)] * If[iselec, mrednce[h, h, cost, costpr], mredncex[h, h, cost, costpr]];
prefactor * helsum];

kappaemcnc[{eb_, t_, mue_, knu_, cost_, iselec_}]:= NIntegrate[Sum[kappaencintegrand[eb, t, mue, knu, cost, costpr, phi, h, iselec], 
	{h, {-1, 1}}], {costpr, -1, 1}, {phi, 0, 2 * Pi}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}];
	
(* monte carlo tools *)
(* draw a random set of conditions and neutrino properties *)
randomcond[{t1_, t2_, n1_, n2_, yp1_, yp2_, eb1_, eb2_}]:={RandomReal[{t1, t2}], Exp[RandomReal[{Log[n1], Log[n2]}]], RandomReal[{yp1, yp2}], Exp[RandomReal[{Log[eb1], Log[eb2]}]]};
randomnu[{knu1_, knu2_, cost1_, cost2_, ui1_, ui2_}]:={RandomReal[{knu1, knu2}], RandomReal[{cost1, cost2}], RandomReal[{ui1, ui2}]};
randompair[{knu1_, knu2_, cost1_, cost2_, phi1_, phi2_}]:={RandomReal[{knu1, knu2}], RandomReal[{cost1, cost2}], 
	RandomReal[{knu1, knu2}], RandomReal[{cost1, cost2}], RandomReal[{phi1, phi2}]};

opacitywrapper[t1_, t2_, n1_, n2_, yp1_, yp2_, eb1_, eb2_, knu1_, knu2_, cost1_, cost2_, ui1_, ui2_, numcond_, numnu_, numkern_]:=Module[
{icond, mulist, inu, cond, nus, params, resultsn, 
	resultsp, finalresults, numsamples, counter, temp},
cond = Map[randomcond, ConstantArray[{t1, t2, n1, n2, yp1, yp2, eb1, eb2}, numcond]];
nus = Map[randomnu, ConstantArray[{knu1, knu2, cost1, cost2, ui1, ui2}, numnu]];
mulist = {};
Print["Starting mu"];
For[icond = 1, icond < numcond + 1, icond++,
AppendTo[mulist, {(mu/.Quiet[mueofn[cond[[icond]][[2]] * cond[[icond]][[3]], cond[[icond]][[4]], cond[[icond]][[1]]]]),
(mu/.Quiet[munofn[cond[[icond]][[2]] * (1 - cond[[icond]][[3]]), cond[[icond]][[4]], cond[[icond]][[1]]]]),
(mu/.Quiet[mupofn[cond[[icond]][[2]] * cond[[icond]][[3]], cond[[icond]][[4]], cond[[icond]][[1]]]])}]];
params = {};
For[inu = 1, inu < numnu + 1, inu++, For[icond = 1, icond < numcond + 1, icond++,
numsamples = Min[{Max[{Ceiling[10 * (cond[[icond]][[1]] / (20 / HBARC))^3 / (cond[[icond]][[4]] / (100 * G15TOFM))^2], 10}], 60}];
For[counter = 0, counter < 10, counter++, AppendTo[params, {cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], mulist[[icond]][[2]], mulist[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], numsamples}]]]];
resultsn = {};
resultsp = {};
CloseKernels[];
LaunchKernels[numkern];
Print["Starting neutron opacity"];
temp = ParallelMap[kappanmc, params];
AppendTo[resultsn, temp];
Print["Done with neutrons, starting proton opacity"];
temp = ParallelMap[kappapmc, params];
AppendTo[resultsp, temp];
CloseKernels[];
finalresults = {};
For[inu = 1, inu < numnu + 1, inu++, For[icond = 1, icond < numcond + 1, icond++,
AppendTo[finalresults, {Mean[resultsn[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]],
	StandardDeviation[resultsn[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	Mean[resultsp[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	StandardDeviation[resultsp[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]],
	kappan[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], cond[[icond]][[2]], cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]]],
	kappandegen[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], mulist[[icond]][[2]], cond[[icond]][[2]], cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]]],
	kappap[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], cond[[icond]][[2]], cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]]],
	kappapdegen[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], mulist[[icond]][[2]], cond[[icond]][[2]], cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]]],
	cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], cond[[icond]][[2]], cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]]}]]];
finalresults];

scatteringwrapper[t1_, t2_, n1_, n2_, yp1_, yp2_, eb1_, eb2_, knu1_, knu2_, cost1_, cost2_, phi1_, phi2_, numcond_, numnu_, numkern_]:=Module[
{icond, mulist, inu, cond, nus, paramsdn, paramsn, paramsp, paramsdp, paramse, 
	resultsdn, resultsn, resultsdp, resultsp, resultse, finalresults, nump, counter},
cond = Map[randomcond, ConstantArray[{t1, t2, n1, n2, yp1, yp2, eb1, eb2}, numcond]];
nus = Map[randompair, ConstantArray[{knu1, knu2, cost1, cost2, phi1, phi2}, numnu]];
mulist = {};
Print["Starting mu"];
For[icond = 1, icond < numcond + 1, icond++,
AppendTo[mulist, {(mu/.Quiet[mueofn[cond[[icond]][[2]] * cond[[icond]][[3]], cond[[icond]][[4]], cond[[icond]][[1]]]]),
(mu/.Quiet[munofn[cond[[icond]][[2]] * (1 - cond[[icond]][[3]]), cond[[icond]][[4]], cond[[icond]][[1]]]]),
(mu/.Quiet[mupofn[cond[[icond]][[2]] * cond[[icond]][[3]], cond[[icond]][[4]], cond[[icond]][[1]]]])}]];
paramsdn = {};
paramsn = {};
paramsdp = {};
paramsp = {};
paramse = {};
For[inu = 1, inu < numnu + 1, inu++, For[icond = 1, icond < numcond + 1, icond++,
nump = Min[{Max[{10, Ceiling[10 * ((cond[[icond]][[1]] + nus[[inu]][[1]]) / (20 / HBARC))^2 / (cond[[icond]][[4]] / (100 * G15TOFM))^2]}], 30}];
Print[nump];
For[counter = 0, counter < 10, counter++, AppendTo[paramsdn, {cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[2]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], nus[[inu]][[4]], nus[[inu]][[5]]}]];
For[counter = 0, counter < 10, counter++, AppendTo[paramsn, {cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[2]], nus[[inu]][[1]], nus[[inu]][[2]]}]];
For[counter = 0, counter < 10, counter++, AppendTo[paramsdp, {cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], nus[[inu]][[4]], nus[[inu]][[5]], nump}]];
For[counter = 0, counter < 10, counter++, AppendTo[paramsp, {cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nump}]];
For[counter = 0, counter < 10, counter++, AppendTo[paramse, {cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], nus[[inu]][[1]], nus[[inu]][[2]], True}]]];
];
resultsdn = {};
resultsn = {};
resultsdp = {};
resultsp = {};
resultse = {};
CloseKernels[];
LaunchKernels[numkern];
Print["Starting proton opacity"];
AppendTo[resultsp, ParallelMap[kappapmcnc, paramsp]];
Print["Starting electron opacity"];
AppendTo[resultse, ParallelMap[kappaemcnc, paramse]];
Print["Starting diff neutron opacity"];
AppendTo[resultsdn, ParallelMap[kappanmcncdiff, paramsdn]];
Print["Starting neutron opacity"];
AppendTo[resultsn, ParallelMap[kappanmcnc, paramsn]];
Print["Starting diff proton opacity"];
AppendTo[resultsdp, ParallelMap[kappapmcncdiff, paramsdp]];
CloseKernels[];
Print[resultsp];
Print[resultsdp];
finalresults = {};
For[inu = 1, inu < numnu + 1, inu++, For[icond = 1, icond < numcond + 1, icond++,
AppendTo[finalresults, {Mean[resultsdn[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	StandardDeviation[resultsdn[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]],
	Mean[resultsn[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	StandardDeviation[resultsn[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	Mean[resultsdp[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	StandardDeviation[resultsdp[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]],
	Mean[resultsp[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	StandardDeviation[resultsp[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	Mean[resultse[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	StandardDeviation[resultse[[1,((inu - 1) * numcond + icond) * 10 - 9 ;; ((inu - 1) * numcond + icond) * 10]]], 
	kappanncdiff[cond[[icond]][[4]], cond[[icond]][[1]], cond[[icond]][[2]] * (1 - cond[[icond]][[3]]), nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], nus[[inu]][[4]], nus[[inu]][[5]]],
	kappannc[cond[[icond]][[4]], cond[[icond]][[1]], cond[[icond]][[2]] * (1 -  cond[[icond]][[3]]), nus[[inu]][[1]], nus[[inu]][[2]]],
	kappanncdiffdegen[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[2]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], nus[[inu]][[4]], nus[[inu]][[5]]],
	kappanncdegen[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[2]], nus[[inu]][[1]], nus[[inu]][[2]]],
	kappapncdiff[cond[[icond]][[4]], cond[[icond]][[1]], cond[[icond]][[2]] * cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], nus[[inu]][[4]], nus[[inu]][[5]]],
	kappapnc[cond[[icond]][[4]], cond[[icond]][[1]], cond[[icond]][[2]] * cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]]],
	kappaenc[cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], nus[[inu]][[1]], nus[[inu]][[2]], True],
	cond[[icond]][[4]], cond[[icond]][[1]], mulist[[icond]][[1]], cond[[icond]][[2]], cond[[icond]][[3]], nus[[inu]][[1]], nus[[inu]][[2]], nus[[inu]][[3]], nus[[inu]][[4]], nus[[inu]][[5]]}];]];
finalresults];
	
