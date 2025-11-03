(* ::Package:: *)

(* example entry in a notebook to run this code

Get[FileNameJoin[{NotebookDirectory[],"Monte carlo neutrino checks.wl"}]]
kappan[100 * G15TOFM, 5 / HBARC, 0.04870282860614022, 0.16 * 10^-3, 0.25, 20 / HBARC, 0.6, 0]
kappanmc[100 * G15TOFM, 5 / HBARC, 0.04870282860614022`, -0.09603718718498111`, -0.05445418660894468`, 20  / HBARC, 0.6, 0, 10]

*)

(*This package calculates Monte Carlo integrations of neutrino opacities and scattering to check accuracy of approximation methods*)

(*constants*)
SINTW = Sqrt[0.231]; (*Weinberg angle*)
HBARC = 197.3269718; (*all calculations will be done in fm*)
NUCMASS = 2 * 939.5653 * 938.272 / (939.5653 + 938.272) / HBARC; (*average nucleon mass*)
MSPLIT = (939.5653 - 938.272) / HBARC;
ELECMASS = 0.511 / HBARC;
GA = 1.27; (*no g-t discrepancy for electroweak*)
G15TOFM = 19.5 * Sqrt[4 * Pi / 137] / HBARC^2; (*convention that alpha = e^2/4pi*)
GP = 5.5858; GN = -3.8263; (*for anomalous magnetic moments, note that GE=2.002 so the g-2 is negligible for leptons*)

(* functions for getting chemical potentials and distributions *)
nfd[e_, mu_, t_]:= 1 / (Exp[(e - mu) / t] + 1);(*If[(e - mu) / t > 50, 0, If[(e - mu) / t < -50, 1, 1 / (Exp[(e - mu) / t] + 1)]]*)
munofn[n_, eb_, t_]:= FindRoot[-n + NIntegrate[k^2 / (2 * Pi^2) * (nfd[k^2 / (2 * NUCMASS) + GN * eb / (4 * NUCMASS), mu, t] 
	+ nfd[k^2 / (2 * NUCMASS) - GN * eb / (4 * NUCMASS), mu, t]), {k, 0, Infinity}], {mu, t}];
mueofn[n_, eb_, t_]:= FindRoot[-n + NIntegrate[eb / (4 * Pi^2) * nfd[Sqrt[k^2 + ELECMASS^2], mu, t] + Sum[eb / (2 * Pi^2) * nfd[Sqrt[k^2 + 2 * ne * eb + ELECMASS^2], mu, t], 
	{ne, 1, Ceiling[t^2 / eb * 3]}], {k, -Infinity, Infinity}], {mu, t}];
mupofn[n_, eb_, t_]:= FindRoot[-n + NIntegrate[eb / (4 * Pi^2) * nfd[k^2 / (2 * NUCMASS) - (GP - 2) * eb / (4 * NUCMASS), mu, t] 
	+ Sum[eb / (4 * Pi^2) * (nfd[k^2 / (2 * NUCMASS) + np * eb / NUCMASS + (GP - 2) * eb / (4 * NUCMASS), mu, t] + nfd[k^2 / (2 * NUCMASS) + np * eb / NUCMASS - (GP - 2) * eb / (4 * NUCMASS), mu, t]), 
	{np, 1, t^2 / eb * 5}], {k, -Infinity, Infinity}], {mu, t}];
thetapm[x_, sp_, sn_, eb_, ui_, pm_]:=HeavisideTheta[MSPLIT + pm * ELECMASS + ui - eb / (4 * NUCMASS) * (GN * sn - GP * sp) - x];
e0pm[knu_, sp_, sn_, eb_, ui_, pm_]:= Module[{delta},
delta = MSPLIT + ui - eb / (4 * NUCMASS) * (GN * sn - GP * sp);
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
bigb[kfn_, u_, w_, eb_, t_, knu_, cost_]:=Module[{mt, kfnt, knut},
mt = NUCMASS * t;
kfnt = kfn / Sqrt[2 * mt];
knut = knu / Sqrt[2 * mt];
(2 * mt)^(3 / 2) / (2 * Pi)^2 * (f1[kfnt, u] + w * g1[kfnt, u] + knut^2 * ((f1[kfnt, u] + w * g1[kfnt, u]) * (u + cost^2 - u * cost^2)
	+ (f2[kfnt, u] + w * g2[kfnt, u]) * (u^2 * cost^2 - u^2) + (f3[kfnt, u] + w * g3[kfnt, u]) * (u^2 + 2 * cost^2 + u^2 * cost^2)))];

(*functions for elastic scattering kernels*)
ffunc[q0_,qz_, qp_, temp_, eb_, s_, spr_]:= Module[{t, q, qpt, expterm, prefactor, sfunc, sarg},
q = Sqrt[qz^2 + qp^2];
t = Exp[-eb / (NUCMASS * temp)];
qpt = qp * Sqrt[eb * Sqrt[t] / (NUCMASS * temp * (1 - t))];
sarg = q * eb / (temp * Sqrt[2]) * (qz^2 + qpt^2) / (qz^2 * qpt^2);
sfunc = If[q0 > q, Erfc[sarg], 1 + Erf[sarg]];
expterm = Exp[-qpt^2 * NUCMASS * temp * Sqrt[t] / eb^2 - (q^2 + 2 * NUCMASS * q0 - GP * eb * (spr - s) / 2) / (8 * NUCMASS * temp * (qz^2 + qpt^2))];
prefactor = Sqrt[Pi * NUCMASS / (2 * temp)] * 1 / (1 - t) * qpt / (qz^2 + qpt^2);
prefactor * expterm * sfunc];
npnorm[eb_, t_]:=Module[{mt},
mt = t * NUCMASS;
Sinh[eb / (2 * mt)] / Cosh[GP * eb / (4 * mt)] * (2 * Pi / mt)^(3 / 2) * mt / eb];

protkern[eb_, t_, rhop_, q0_, q_, cosq_, s_, spr_]:= Module[{mt, ebmt, npnorm},
mt = NUCMASS * t;
ebmt = eb / mt;
npnorm = Sinh[ebmt / 2] / Cosh[GP * ebmt / 4] * (2 * Pi)^(3/2) / (eb * Sqrt[mt]);
eb / (2 * pi) * rhop * npnorm * (ffunc[eb, t, q0, q, cosq, s, spr] - rhop * npnorm * ffunc[eb, t / 2, q0, q, cosq, s, spr] * Exp[q0 / t])];

neutkern[mun_, t_, q0_, q_]:= If[q0 == 0, 0, If[q0 > 0, NUCMASS^2 * t / (2 * Pi * q) * (Log[1 + Exp[mun / t]] - Log[1 + Exp[(mun - q0) / t]]) / (Exp[q0 / t] - 1),
	NUCMASS^2 * t / (2 * Pi * q) * (Log[1 + Exp[(mun + q0) / t]] - Log[1 + Exp[mun / t]]) / (Exp[q0 / t] - 1)]];
	
eleckern[mue_, eb_, t_, q0_, q_, cosq_]:= If[q0 <= q, 0, eb / (2 * Pi) * t * Exp[-q^2 * (1 - cosq^2) / (2 * eb)] * (Log[1 + Exp[mue / t]] - Log[1 + Exp[(mue - q0) / t]]) / (Exp[q0 / t] - 1)];

alpha[delta_, knu_]:=If[delta > 0, 4 * knu^2 + 12 * knu * delta + 6 * delta^2, If[-knu < delta, 4 * (knu + delta)^3 / knu, 0]];
beta[delta_, knu_]:= If[delta > 0, knu^2 + 4 * knu * delta, If[-knu < delta, (knu + delta)^4 / knu^2, 0]];
asspr[s_, spr_, cost_]:=If[s == spr, 1 / 2 * (1 + GA^2 + 2 * GA * s * (1 - 4 * SINTW^2) * cost), GA^2];
bsspr[s_, spr_, cost_]:=If[s == spr, 1 / 2 * (1 - GA^2 + 2 * GA^2 * cost + 2 * GA * s * (1 - 4 * SINTW^2) * cost), GA^2 * (-1 + 2 * cost)];

(*Matrix elements for charged current*)
mredcc[sp_, sn_, nezero_, cost_]:=If[nezero, If[sn == 1, If[sp == 1, 2 * (1 + GA)^2 * (1 + cost), 0], If[sp == 1, 8 * GA^2 * (1 - cost), 2 * (1 - GA)^2 * (1 + cost)]],
	If[sn == 1, If[sp == 1, 2 * (1 + GA)^2 * (1 + cost) + 2 * (1 - GA)^2 * (1 - cost), 8 * GA^2 * (1 + cost)], If[sp == 1, 8 * GA^2 * (1 - cost),
	2 * (1 + GA)^2 * (1 - cost) + 2 * (1 - GA)^2 * (1 + cost)]]];
vtildecc[x_, sp_, sn_, cost_]:=If[x <= 0, 0, If[x < 1 / 2, mredcc[sp, sn, True, cost], mredcc[sp, sn, True, cost] + (2 * x - 1) * mred[sp, sn, False, cost]]];

(*Matrix elements for neutral current*)
mredncp[s_, spr_, cost_, costpr_, phi_]:= Module[{sint, sintpr},
sint = Sqrt[1 - cost^2];
sintpr = Sqrt[1 - costpr^2];
If[s == spr, 1 / 2 * ((1 - 4 * SINTW^2)^2 * (1 + cost + costpr + sint * sintpr * Cos[phi]) + GA^2 * (1 + cost * costpr - sint * sintpr * Cos[phi]) + 2 * GA * s * (1 - 4 * SINTW^2) * (cost + costpr)),
	GA^2 * (1 + cost * costpr - sint * sintpr * Cos[phi])]];
mredncn[s_, spr_, cost_, costpr_, phi_]:= Module[{sint, sintpr},
sint = Sqrt[1 - cost^2];
sintpr = Sqrt[1 - costpr^2];
If[s == spr, 1 / 2 * ((1 + cost + costpr + sint * sintpr * Cos[phi]) + GA^2 * (1 + cost * costpr - sint * sintpr * Cos[phi]) + 2 * GA * s * (cost + costpr)),
	GA^2 * (1 + cost * costpr - sint * sintpr * Cos[phi])]];
mrednce[h_, hpr_, cost_, costpr_]:=If[h != hpr, 0, If[h == 1, 32 * SINTW^4 * (1 - cost) * (1 - costpr), 8 * (1 - 2 * SINTW^2)^2 * (1 + cost) * (1 + costpr)]];
mredncex[h_, hpr_, cost_, costpr_]:=If[h != hpr, 0, If[h == 1, 32 * SINTW^4 * (1 - cost) * (1 - costpr), 8 * (1 + 2 * SINTW^2)^2 * (1 + cost) * (1 + costpr)]];

(*Opacities - analytical approx*)
kappan[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sp, sn},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb / (8 * Pi * Cosh[GN * ebmt / 4]);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t]) * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * Exp[GN * sn * ebmt / 4]
	 * (1 - nb * yp * Exp[GP * sp * ebmt / 4] / Cosh[GP * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, 1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, 1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

kappandegen[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor},
mt = NUCMASS * t;
kfn = (3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3)));
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb / (8 * Pi);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t]) * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * (1
	- Exp[GP * sp * ebmt / 4]/(2 * Cosh[GP * ebmt / 4]) * (1 - Exp[-ebmt]) * bigb[Sqrt[kfn^2 - GN * sn * eb / 2], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8]]),
	{sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

kappap[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sn, sp},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * yp * eb / (8 * Pi * Cosh[GP * ebmt / 4]);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], -mue, t]) * (1 - thetapm[knu, sp, sn, eb, ui, 1]) * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * Exp[GP * sp * ebmt / 4]
	 * (1 - nb * (1 - yp) * Exp[GN * sn * ebmt / 4] / Cosh[GN * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, -1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, -1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

kappapdegen[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor},
mt = NUCMASS * t;
kfn = (3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3)));
ebmt = eb / mt;
prefactor = nb * yp * eb / (8 * Pi * Cosh[GP * ebmt / 4]);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], -mue, t]) * (1 - thetapm[knu, sp, sn, eb, ui, 1]) * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * (1
	- (1 - Exp[-ebmt]) * bigb[Sqrt[kfn^2 - GN * sn * eb / 2], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8]]),
	{sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

(*Differential emissivities - analytic approx*)
epsn[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:=Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sp, sn},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb * knu^3 / (64 * Pi^4 * Cosh[GN * ebmt / 4]);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], mue, t]) * thetapm[knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * Exp[GN * sn * ebmt / 4]
	 * (1 - nb * yp * Exp[GP * sp * ebmt / 4] / Cosh[GP * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, -1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, -1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

epsp[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:=Module[{kperp, kz, prefactor, mt, ebmt, spinsum, sp, sn},
kz = knu * cost;
kperp = knu * Sqrt[1 - cost^2];
mt = NUCMASS * t;
ebmt = eb / mt;
prefactor = nb * yp * eb * knu^3 / (64 * Pi^4 * Cosh[GP * ebmt / 4]);
spinsum = Sum[Sum[nfd[e0pm[knu, sp, sn, eb, ui, 1], mue, t] * thetapm[-knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, 1]^2 / (2 * eb), sp, sn, cost] * Exp[GP * sp * ebmt / 4]
	 * (1 - nb * (1 - yp) * Exp[GN * sn * ebmt / 4] / Cosh[GN * ebmt / 4] * (Pi / mt)^(3/2) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) * Cosh[kz * e0pm[knu, sp, sn, eb, ui, 1] / (2 * mt)]
	 * Exp[-kperp^2 / (2 * mt) * (1 - Exp[-ebmt]) / (ebmt + 1 - Exp[-ebmt]) - (kz^2 + e0pm[knu, sp, sn, eb, ui, 1]^2)/(4*mt)]), {sp, {-1, 1}}], {sn, {-1, 1}}];
prefactor * spinsum];

epsndegen[eb_, t_, mue_, nb_, yp_, knu_, cost_, ui_]:= Module[{kfn, sp, sn, mt, ebmt, spinsum, prefactor},
mt = NUCMASS * t;
kfn = (3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3)));
ebmt = eb / mt;
prefactor = nb * (1 - yp) * eb * knu^3 / (64 * Pi^4);
spinsum = Sum[Sum[(1 - nfd[e0pm[knu, sp, sn, eb, ui, -1], mue, t]) * thetapm[knu, sp, sn, eb, ui, -1] * vtildecc[e0pm[knu, sp, sn, eb, ui, -1]^2 / (2 * eb), sp, sn, cost] * (1
	- Exp[GP * sp * ebmt / 4]/(2 * Cosh[GP * ebmt / 4]) * (1 - Exp[-ebmt]) * bigb[Sqrt[kfn^2 - GN * sn * eb / 2], (1 - Exp[-ebmt]) / ebmt, Exp[GN * sn * ebmt / 8] / Cosh[GN * ebmt / 8]]),
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

(*Differential elastic scattering rate - analytic approx*)
dkappancp[eb_, t_, nb_, yp_, knu_, knupr_, cost_, costpr_, phi_]:=Module[{s, spr, q, cosq},
cosq = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosq];
knu^2 * knupr^2 / (2 * Pi)^6 * Sum[Sum[mredncp[s, spr, cost, costpr, phi] * protkern[eb, t, nb * yp, knu - knupr, q, cosq, s, spr],{s, {-1, 1}}], {spr, {-1, 1}}]];

dkappancn[eb_, t_, nb_, yp_, knu_, knupr_, cost_, costpr_, phi_]:=Module[{s, spr, q, cosq, mun},
cosq = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosq];
mun = ((3 * Pi^2 * nb * (1 - yp))^(1 / 3) * (1 + Pi^2 * t^2 * NUCMASS^2 / (6 * (3 * Pi^2 * n)^(4 / 3))))^2 / (2 * NUCMASS);
knu^2 * knupr^2 / (2 * Pi)^6 * Sum[Sum[If[-GN * eb * (s - spr) / (4 * NUCMASS) + knu - knupr > 0, 
	mredncn[s, spr, cost, costpr, phi] * neutkern[mun, t, -GN * eb * (s - spr) / (4 * NUCMASS), q], 0],{s, {-1, 1}}], {spr, {-1, 1}}]];
	
dkappance[eb_, t_, mue_, knu_, knupr_, cost_, costpr_, phi_]:=Module[{q, cosq, h, hpr},
cosq = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosq];
knu^2 * knupr^2 / (2 * Pi)^6 * Sum[Sum[mrednce[h, hpr, cost, costpr] * eleckern[mue, eb, t, knupr - knu, q, cosq], {h, {-1, 1}}], {hpr, {-1, 1}}]];

dkappancex[eb_, t_, mue_, knu_, knupr_, cost_, costpr_, phi_]:=Module[{q, cosq, h, hpr},
cosq = cost * costpr + Sqrt[1 - cost^2] * Sqrt[1 - costpr^2] * Cos[phi];
q = Sqrt[knu^2 + knupr^2 - 2 * knu * knupr * cosq];
knu^2 * knupr^2 / (2 * Pi)^6 * Sum[Sum[mredncex[h, hpr, cost, costpr] * eleckern[mue, eb, t, knupr - knu, q, cosq], {h, {-1, 1}}], {hpr, {-1, 1}}]];
	
(*Elastic scattering opacity - analytic approx*)

kappancn[eb_, t_, knu_, cost_]:= NUCMASS^2 * t / (6 * (2 * Pi)^3) * Sum[Sum[asspr[s, spr, cost] * alpha[-GN * eb * (s - spr) / (4 * NUCMASS), knu] 
	+ bsspr[s, spr, cost] * beta[-GN * eb * (s - spr) / (4 * NUCMASS), knu], {s, {-1, 1}}], {spr, {-1, 1}}];
	
kappance[eb_, t_, mue_, knu_, cost_]:=Module[{prefactor, rest},
prefactor = 32 * eb^2 * knu / Pi * (1 + cost) / (1 - Exp[-knu * (1 + cost) / t]) * Exp[-knu^2 * (1 - cost^2) / (2 * eb + 4 * t * knu * (1 - cost))];
rest = 2 * SINTW^4 * t * Sqrt[eb] * (1 - cost) / (eb + 2 * t * knu * (1 - cost)) + (1 - 2 * SINTW^2)^2 * t^2 * eb * (1 + cost) / (eb + 2 * t * knu * (1 - cost))^2 * 
	(1 + t * knu^3 * (1 - cost - cost^2 + cost^3) / (eb * (eb + 2 * t * knu * (1 - cost))));
prefactor * rest];

kappancex[eb_, t_, mue_, knu_, cost_]:=Module[{prefactor, rest},
prefactor = 32 * eb^2 * knu / Pi * (1 + cost) / (1 - Exp[-knu * (1 + cost) / t]) * Exp[-knu^2 * (1 - cost^2) / (2 * eb + 4 * t * knu * (1 - cost))];
rest = 2 * SINTW^4 * t * Sqrt[eb] * (1 - cost) / (eb + 2 * t * knu * (1 - cost)) + (1 + 2 * SINTW^2)^2 * t^2 * eb * (1 + cost) / (eb + 2 * t * knu * (1 - cost))^2 * 
	(1 + t * knu^3 * (1 - cost - cost^2 + cost^3) / (eb * (eb + 2 * t * knu * (1 - cost))));
prefactor * rest];



(*Opacities - Monte Carlo*)
llmc[integrandfunc_, eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, nsamples_]:=Module[{nemax, npmax, rangee, rangep, weights, flatweights, flatpairs, probs, samples, ne, np, subfunc, weightfunc},
nemax = Floor[If[mue > 0, (mue + 5 * t)^2 / (2 * eb), (5 * t)^2 / (2 * eb)]];
npmax = Floor[If[mup > 0, (mup + 5 * t) * NUCMASS / eb, 5 * t * NUCMASS / eb]];
rangee = Range[0, nemax];
rangep = Range[0, npmax];
weightfunc[np_, ne_]:=Exp[-np * eb / (t * NUCMASS)];
subfunc[np_, ne_]:= integrandfunc[eb, t, mue, mun, mup, knu, cost, ui, ne, np] / weightfunc[np, ne];
weights = Table[weightfunc[np, ne], {np, rangep}, {ne, rangee}];
flatweights = Flatten[weights];
flatpairs = Flatten[Table[{np, ne}, {np, rangep}, {ne, rangee}], 1];
probs = flatweights / Total[flatweights];
samples = RandomChoice[probs -> flatpairs, nsamples];
Total[(subfunc @@@ samples) * Total[flatweights]] / nsamples];

nmcme[mun_, mup_, mue_, kn_, cosn_, phin_, knu_, cost_, sn_, sp_, np_, ne_, eb_, t_, ui_]:= Module[{wperpsq},
wperpsq = kn^2 * (1 - cosn^2) + knu^2 * (1 - cost^2) - 2 * kn * knu * Sqrt[1 - cosn^2] * Sqrt[1 - cost^2] * Cos[phin];
If[ne <= np, (wperpsq / (2 * eb))^(np - ne) * Exp[-wperpsq / (2 * eb)] * Factorial[ne] / Factorial[np] * LaguerreL[ne, np - ne, wperpsq / (2 * eb)]^2,
	(wperpsq / (2 * eb))^(ne - np) * Exp[-wperpsq / (2 * eb)] * Factorial[np] / Factorial[ne] * LaguerreL[np, ne - np, wperpsq / (2 * eb)]^2]
];

nmcspinsum[mun_, mup_, mue_, kn_, cosn_, phin_, knu_, cost_, np_, ne_, eb_, t_, ui_, sp_, sn_]:= Module[{ep, en, ee, kzp},
kzp = kn * cosn + knu * cost;
ep = Sqrt[NUCMASS^2 + kzp^2 + 2 * np * eb - (GP - 2) * sp * eb / 2] - NUCMASS;
en = Sqrt[NUCMASS^2 + kn^2] - GN * sn * eb / (4 * NUCMASS) - NUCMASS;
ee = en + knu - ep + ui + MSPLIT;
If[ee < Sqrt[2 * ne * eb], 0, 1 / 2 * ee * kn^2 / Sqrt[ee^2 - 2 * ne * eb]
	 * eb / (2 * Pi)^4 * nfd[en, mun, t] * (1 - nfd[ee, mue, t]) * (1 - nfd[ep, mup, t]) * If[ne == 0, mredcc[sp, sn, True, cost], mredcc[sp, sn, False, cost]]]]

kappanmcnenp[eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, ne_, np_]:= Module[{knmax, sn, sp, kn, cosn, phin},
knmax = If[mun > 0, Sqrt[2 * (mun + 5 * t) * NUCMASS], Sqrt[2 * (5 * t) * NUCMASS]];
NIntegrate[Sum[Sum[nmcspinsum[mun, mup, mue, kn, cosn, phin, knu, cost, np, ne, eb, t, ui, sp, sn], {sp, {-1, 1}}], {sn, {-1, 1}}] * nmcme[mun, mup, mue, kn, cosn, phin, knu, cost, sn, sp, np, ne, eb, t, ui], 
	 {kn, 0, knmax}, {cosn, -1, 1}, {phin, 0, 2 * Pi}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}]];
	
kappanmc[eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, nsamples_]:= llmc[kappanmcnenep, eb, t, mue, mun, mup, knu, cost, ui, nsamples];
	
pmcintegrand[mun_, mup_, mue_, kn_, cosn_, phin_, knu_, cost_, sn_, sp_, np_, ne_, eb_, t_, ui_]:= Module[{wperpsq, matelt, kzp, ep, en, ee, kze},
wperpsq = kn^2 * (1 - cosn^2) + knu^2 * (1 - cost^2) - 2 * kn * knu * Sqrt[1 - cosn^2] * Sqrt[1 - cost^2] * Cos[phin];
kzp = kn * cosn + knu * cost;
ep = Sqrt[NUCMASS^2 + kzp^2 + 2 * np * eb - (GP - 2) * sp * eb / 2] - NUCMASS;
en = Sqrt[NUCMASS^2 + kn^2] - GN * sn * eb / (4 * NUCMASS) - NUCMASS;
ee = ep + knu - en - ui - MSPLIT;
matelt = If[ne == 0, mredcc[sp, sn, True, cost], mredcc[sp, sn, False, cost]] * If[
	ne <= np, (wperpsq / (2 * eb))^(np - ne) * Exp[-wperpsq / (2 * eb)] * Factorial[ne] / Factorial[np] * LaguerreL[ne, np - ne, wperpsq / (2 * eb)]^2,
	(wperpsq / (2 * eb))^(ne - np) * Exp[-wperpsq / (2 * eb)] * Factorial[np] / Factorial[ne] * LaguerreL[np, ne - np, wperpsq / (2 * eb)]^2];
If[ee < Sqrt[2 * ne * eb], 0, 1 / 2 * ee * kn^2 / Sqrt[ee^2 - 2 * ne * eb] * eb / (2 * Pi)^4 * matelt * (1 - nfd[en, mun, t]) * (1 - nfd[ee, -mue, t]) * nfd[ep, mup, t]]];

kappapmcnenp[eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, ne_, np_]:= Module[{knmax, sn, sp, kn, cosn, phin},
knmax = If[mun > 0, Sqrt[2 * (mun + 5 * t) * NUCMASS], Sqrt[2 * (5 * t) * NUCMASS]];
NIntegrate[Sum[Sum[pmcintegrand[mun, mup, mue, kn, cosn, phin, knu, cost, sn, sp, np, ne, eb, t, ui], 
	{sp, {-1, 1}}], {sn, {-1, 1}}], {kn, 0, knmax}, {cosn, -1, 1}, {phin, 0, 2 * Pi}, Method -> {"AdaptiveMonteCarlo", MaxPoints -> 1000000}]];
	
kappapmc[eb_, t_, mue_, mun_, mup_, knu_, cost_, ui_, nsamples_]:=llmc[kappapmcnenep, eb, t, mue, mun, mup, knu, cost, ui, nsamples];
	

