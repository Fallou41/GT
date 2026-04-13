cd "C:\Users\USER\Desktop\ISE3_2025\GT\base"
clear all
set more off
use "base"


*********************************Description des données**********************************
**** Création des variables logarithme
gen lnagri = ln(agriculture)
gen lnman =ln(manufacture)
gen lncom = ln(commerce)
gen lntran = ln(transport)
gen lnheb = ln(hebergement)
gen lnimmo = ln(immobilier)
gen lnautre = ln(autre)

tsset Annee


********* Cointégration
local vars "lnagri lnman lncom lntran lnheb lnimmo lnautre"
varsoc `vars'
vecrank `vars', lag(4) trace
vecrank `vars', lag(4) max
*ttttt


*************** Calcul des persistances globale


gen lninfor = ln(informel)
** On regresse la variable différenciée
gen dlninfor = D.lninfor



dfuller lninfor

dfuller dlninfor



ac dlninfor

pac dlninfor

twoway (line dlninfor Annee)

twoway (line lninfor Annee)




*************** Calcul des persistances globale univariée
***La persistance avec un AR
** Sélection de lag optimal
varsoc dlninfor, maxlag(5) // Le lag optimal est 5. 
*** On a un AR(5)
regress dlninfor L(1/5).dlninfor
scalar somme_phi = _b[L1.dlninfor] + _b[L2.dlninfor] + _b[L3.dlninfor] + _b[L4.dlninfor] + _b[L5.dlninfor]

* Mesure de persistance P = 1 / (1 - somme_phi)
scalar psi1 = 1 / (1 - somme_phi)
scalar P    = abs(psi1)
di P // Le choc sur le secteur informel est persistant


***** Mesure de la persistance par la densité
*On récupère les résidus de la régréssion
predict resid_dlninfor, residuals
* Étape 2 : Variance de court terme σ² (variance des innovations)
quietly summarize resid_dlninfor
scalar sigma2 = r(Var)
di "  Variance des innovations σ² : " %10.6f sigma2
* Étape 3 : Densité spectrale à la fréquence nulle via Newey-West
* f(0) = (1/2π) · [γ(0) + 2·Σ w(j)·γ(j)]
* avec w(j) = 1 - j/(q+1)  (noyau de Bartlett)
* et q = nombre de lags retenus (règle : q = int(T^(1/3)))
* Étape 3 : Densité spectrale à la fréquence nulle via Newey-West
* f(0) = (1/2π) · [γ(0) + 2·Σ w(j)·γ(j)]
* avec w(j) = 1 - j/(q+1)  (noyau de Bartlett)
* et q = nombre de lags retenus (règle : q = int(T^(1/3)))

quietly summarize resid_dlninfor
scalar T_eff = r(N)
scalar q_nw  = int(T_eff^(1/3))   /* bandwidth de Bartlett */
di "  Nombre d'obs. (T)           : " T_eff
di "  Bandwidth de Bartlett (q)   : " q_nw

* Calcul de γ(0) : autocovariance d'ordre 0
quietly corr resid_dlninfor, covariance
scalar gamma0 = r(Var_1)

* Calcul des autocovariances γ(j) et somme pondérée
scalar sum_weighted_gamma = 0

forvalues j = 1/`=q_nw' {
    * Poids de Bartlett : w(j) = 1 - j/(q+1)
    scalar w_j = 1 - `j' / (q_nw + 1)
    
    * Autocovariance γ(j) = Cov(resid_t, resid_{t-j})
    quietly corr resid_dlninfor L`j'.resid_dlninfor, covariance
    scalar gamma_j = r(cov_12)
    
    scalar sum_weighted_gamma = sum_weighted_gamma + w_j * gamma_j
}

* Densité spectrale à ω = 0 :
* f(0) = (1/2π) · [γ(0) + 2 · Σ w(j)·γ(j)]
scalar f0_hat = (1 / (2 * _pi)) * (gamma0 + 2 * sum_weighted_gamma)

* Étape 4 : Mesure de persistance
* P² = 2π · f(0) / σ²
scalar P_sq  = (2 * _pi * f0_hat) / sigma2
scalar P_hat = sqrt(abs(P_sq))


********* Persistance multisectorielle
local vars "D.lnagri D.lnman D.lncom D.lntran D.lnheb D.lnimmo D.lnautre"
varsoc `vars' //Il y a 4 retards d'après les critères d'information


/*===========================================================
  PERSISTANCE MULTISECTORIELLE VIA DENSITÉ SPECTRALE
  Pesaran, Pierse & Lee (1993) — Noyau de Bartlett (NW)
  7 secteurs | VAR(4)

  Formule DIRECTE sur les résidus :
    P²ᵢ = f̂(0)ᵢᵢ / σ²ᵢ

  où f̂(0)ᵢᵢ = densité spectrale à ω=0 du résidu du secteur i
             estimée directement par noyau de Bartlett
     σ²ᵢ    = variance de court terme du résidu du secteur i
===========================================================*/
capture drop d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre
gen d_lnagri = D.lnagri 
gen d_lnman = D.lnman
gen d_lncom = D.lncom
gen d_lntran = D.lntran
gen d_lnheb = D.lnheb
gen d_lnimmo = D.lnimmo
gen d_lnautre = D.lnautre

local secteurs "d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre"
local noms     "Agri Manufact Commerce Transport Hébergement Immo Autres"
local m = 7
local p = 4

/*-----------------------------------------------------------
  ÉTAPE 1 : Estimation du VAR(4) et récupération des résidus
-----------------------------------------------------------*/
var `secteurs', lags(1/`p')

scalar T_eff = e(N)
scalar q_nw  = int(T_eff^(1/3))

di "  Observations          : " T_eff
di "  Bandwidth Bartlett q  : " q_nw

* Récupérer les résidus de chaque équation
foreach var of local secteurs {
    quietly predict resid_`var' if e(sample), equation(`var') residuals
}

/*-----------------------------------------------------------
  ÉTAPE 2 : Pour chaque secteur i, estimer f̂(0)ᵢᵢ

  f̂(0)ᵢᵢ = γᵢ(0) + 2 · Σⱼ₌₁ᵍ w(j) · γᵢ(j)

  où γᵢ(j) = autocovariance d'ordre j du résidu du secteur i
     w(j)   = 1 - j/(q+1)   (poids de Bartlett)
     σ²ᵢ    = γᵢ(0)
-----------------------------------------------------------*/

di _newline(1)
di "============================================================"
di "  PERSISTANCE — DENSITÉ SPECTRALE DIRECTE (Bartlett/NW)"
di "  VAR(4) | 7 secteurs | q = " q_nw
di "============================================================"
di "  " %15s "Secteur" %14s "σ²ᵢ = γ(0)" %14s "f̂(0)ᵢᵢ" %8s "P²" %8s "P" "   Interprétation"
di "  " "{hline 75}"

local i = 1
foreach var of local secteurs {

    * --- σ²ᵢ = γᵢ(0) : variance du résidu ---
    quietly summarize resid_`var'
    scalar sigma2_`var' = r(Var)

    * --- f̂(0)ᵢᵢ : densité spectrale à fréquence nulle ---
    * On part de γ(0) puis on ajoute les autocovariances pondérées
    scalar f0_`var' = sigma2_`var'

    forvalues j = 1/`=q_nw' {

        * Poids de Bartlett
        scalar w_j = 1 - `j' / (q_nw + 1)

        * Autocovariance γᵢ(j) = Cov(resid_t, resid_{t-j})
        quietly corr resid_`var' L`j'.resid_`var', covariance
        scalar gamma_j = r(cov_12)

        * f̂(0) += 2 · w(j) · γ(j)
        scalar f0_`var' = f0_`var' + 2 * w_j * gamma_j

    }

    * --- P²ᵢ = f̂(0)ᵢᵢ / σ²ᵢ ---
    scalar P_sq_`var' = f0_`var' / sigma2_`var'
    scalar P_`var'    = sqrt(abs(P_sq_`var'))

    local nom : word `i' of `noms'
    local interp ""
    if      abs(P_`var' - 1) < 0.05  local interp "Choc permanent (≈RW)"
    else if P_`var' > 1               local interp "Amplification"
    else                              local interp "Atténuation partielle"

    di "  " %15s "`nom'" %14.6f sigma2_`var' %14.6f f0_`var' ///
       %8.4f P_sq_`var' %8.4f P_`var' "   `interp'"

    local i = `i' + 1
}

di "  " "{hline 75}"
di "============================================================"


/*===========================================================
  PERSISTANCE MULTISECTORIELLE — MODÈLE CONTRAINT
  Pesaran, Pierse & Lee (1993) — Noyau de Bartlett (NW)
  7 secteurs | p = 4

  Contrainte imposée : pour chaque secteur i, l'équation
  ne contient que :
    - ses propres valeurs passées        : L(1/4).d_lni
    - la somme des autres secteurs       : L(1/4).reste_i

  Soit pour le secteur i :
  d_lni_t = μ + Σₖ αₖ·d_lni_{t-k} + Σₖ βₖ·reste_i_{t-k} + εᵢₜ

  avec reste_i = Σⱼ≠ᵢ d_lnj   (agrégat non pondéré)

  → Modèle estimé équation par équation par OLS (suridentifié)
  → Persistance via densité spectrale directe sur les résidus
===========================================================*/

local secteurs "d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre"
local noms     "Agri Manufact Commerce Transport Hébergement Immo Autres"
local m  = 7
local p  = 4

/*-----------------------------------------------------------
  ÉTAPE 1 : Construire la variable "reste_i" pour chaque secteur
  reste_i = somme de tous les secteurs sauf i
-----------------------------------------------------------*/

foreach var of local secteurs {

    * Somme de tous les secteurs
	capture drop total_all
    egen total_all = rowtotal(`secteurs')

    * reste_i = total - secteur i
	capture drop reste_`var'
    gen reste_`var' = total_all - `var'

    drop total_all
}

di _newline(1)
di "  Variables 'reste' construites pour chaque secteur"
di "  (somme simple des 6 autres secteurs)"

/*-----------------------------------------------------------
  ÉTAPE 2 : Estimation équation par équation (OLS contraint)
  Pour chaque secteur i :
    d_lni = μ + α1·L.d_lni + ... + α4·L4.d_lni
              + β1·L.reste_i + ... + β4·L4.reste_i + ε
  → 2×4 + 1 = 9 paramètres par équation
    (vs 7×4 + 1 = 29 dans le VAR non contraint)
-----------------------------------------------------------*/

di _newline(1)
di "============================================================"
di "  ÉTAPE 2 : ESTIMATIONS OLS CONTRAINTES (équation par équation)"
di "============================================================"

local i = 1
foreach var of local secteurs {

    local nom : word `i' of `noms'
    di _newline(1) "  >>> Secteur : `nom' (`var')"

    * Régression contrainte
    regress `var' L(1/`p').`var' L(1/`p').reste_`var'

    * Stocker les résidus
	capture drop resid_`var'
    quietly predict resid_`var' if e(sample), residuals

    * Afficher AIC et BIC pour diagnostic
    scalar aic_`var' = e(N) * ln(e(rss)/e(N)) + 2 * e(df_m)
    scalar bic_`var' = e(N) * ln(e(rss)/e(N)) + ln(e(N)) * e(df_m)
    di "    AIC = " %8.2f aic_`var' "  |  BIC = " %8.2f bic_`var'
    di "    Nombre de paramètres : " e(df_m) " (vs 29 dans VAR non contraint)"

    local i = `i' + 1
}

/*-----------------------------------------------------------
  ÉTAPE 3 : Persistance via densité spectrale directe
  sur les résidus du modèle contraint

  f̂(0)ᵢᵢ = γᵢ(0) + 2·Σⱼ₌₁ᵍ w(j)·γᵢ(j)
  w(j)    = 1 - j/(q+1)   (noyau de Bartlett)
  q       = int(T^(1/3))

  P²ᵢ = f̂(0)ᵢᵢ / σ²ᵢ
-----------------------------------------------------------*/

* Nombre d'observations et bandwidth
quietly regress `= word("`secteurs'", 1)' L(1/`p').`= word("`secteurs'", 1)' ///
    L(1/`p').reste_`= word("`secteurs'", 1)'
scalar T_eff = e(N)
scalar q_nw  = int(T_eff^(1/3))

di _newline(1)
di "============================================================"
di "  ÉTAPE 3 : PERSISTANCE — DENSITÉ SPECTRALE (Bartlett/NW)"
di "  Modèle contraint | p = 4 | q = " q_nw
di "============================================================"
di "  " %15s "Secteur" %14s "σ²ᵢ = γ(0)" %14s "f̂(0)ᵢᵢ" ///
   %8s "P²" %8s "P" "   Interprétation"
di "  " "{hline 78}"

local i = 1
foreach var of local secteurs {

    * --- σ²ᵢ = γᵢ(0) ---
    quietly summarize resid_`var'
    scalar sigma2_`var' = r(Var)

    * --- f̂(0)ᵢᵢ par Bartlett ---
    scalar f0_`var' = sigma2_`var'

    forvalues j = 1/`=q_nw' {
        scalar w_j = 1 - `j' / (q_nw + 1)
        quietly corr resid_`var' L`j'.resid_`var', covariance
        scalar gamma_j = r(cov_12)
        scalar f0_`var' = f0_`var' + 2 * w_j * gamma_j
    }
	* --- P²ᵢ = f̂(0)ᵢᵢ / σ²ᵢ ---
    scalar P_sq_`var' = f0_`var' / sigma2_`var'
    scalar P_`var'    = sqrt(abs(P_sq_`var'))

    local nom : word `i' of `noms'
    local interp ""
    if      abs(P_`var' - 1) < 0.05  local interp "Choc permanent (≈RW)"
    else if P_`var' > 1               local interp "Amplification"
    else                              local interp "Atténuation partielle"

    di "  " %15s "`nom'" %14.6f sigma2_`var' %14.6f f0_`var' ///
       %8.4f P_sq_`var' %8.4f P_`var' "   `interp'"

    local i = `i' + 1

} 


di "  " "{hline 78}"

/*-----------------------------------------------------------
  ÉTAPE 4 : Comparaison VAR non contraint vs modèle contraint
-----------------------------------------------------------*/

di _newline(1)
di "============================================================"
di "  COMPARAISON : Modèle contraint vs VAR(4) non contraint"
di "============================================================"
di "  " %15s "Secteur" %20s "P contraint" %20s "Nb paramètres"
di "  " "{hline 58}"

local i = 1
foreach var of local secteurs {
    local nom : word `i' of `noms'
    di "  " %15s "`nom'" %20.4f P_`var' %20.0f 9
    local i = `i' + 1
}

di "  " "{hline 58}"
di "  VAR(4) non contraint : 29 paramètres par équation"
di "  Modèle contraint     :  9 paramètres par équation"
di "  Gain                 : 20 paramètres en moins (-69%)"
di "============================================================"


/*-----------------------------------------------------------
  TEST DU RAPPORT DE VRAISEMBLANCE (LR)
  H0 : les contraintes sont valides
  LR = 2·(lnL_NC - lnL_C) ~ χ²(r)
  r  = K_NC - K_C = 7×29 - 7×9 = 140
-----------------------------------------------------------*/

* --- lnL_NC : VAR(4) non contraint ---
quietly var `secteurs', lags(1/`p')
scalar ll_NC = e(ll)

* --- lnL_C : modèle contraint via |Sigma_C| ---
* Sigma_C construite à partir des résidus déjà disponibles
matrix Sigma_C = J(`m', `m', 0)
local i = 1
foreach vi of local secteurs {
    local j = 1
    foreach vj of local secteurs {
        quietly corr resid_`vi' resid_`vj', covariance
        matrix Sigma_C[`i', `j'] = r(cov_12)
        local j = `j' + 1
    }
    local i = `i' + 1
}

scalar ll_C  = -(T_eff / 2) * (`m' * ln(2 * _pi) + ln(det(Sigma_C)) + `m')

* --- Statistique LR ---
scalar LR_stat = 2 * (ll_NC - ll_C)
scalar r_ddl   = `m' * (2*`p'+1+1) - `m' * (2*`p'+1)  /* = 7×(29-9) = 140 */
scalar p_val   = chi2tail(r_ddl, LR_stat)

di _newline(1)
di "============================================================"
di "  TEST LR : Modèle contraint vs VAR(4) non contraint"
di "============================================================"
di "  lnL non contraint  : " %10.4f ll_NC
di "  lnL contraint      : " %10.4f ll_C
di "  LR = 2·(lnL_NC - lnL_C) = " %10.4f LR_stat
di "  Degrés de liberté  : " r_ddl
di "  p-valeur           : " %10.4f p_val
di "  Seuil 5%  χ²(" r_ddl ") : " %10.4f invchi2tail(r_ddl, 0.05)
di "  Seuil 1%  χ²(" r_ddl ") : " %10.4f invchi2tail(r_ddl, 0.01)
di "  {hline 55}"
if p_val >= 0.05 {
    di "  → Non-rejet H0 : contraintes VALIDES"
    di "    Le modèle contraint est préférable"
}
else {
    di "  → Rejet H0 : contraintes INVALIDES"
    di "    Préférer le VAR(4) non contraint"
}
di "============================================================"

/*-----------------------------------------------------------
  NETTOYAGE
-----------------------------------------------------------*/
foreach var of local secteurs {
    drop reste_`var' resid_`var'
}
