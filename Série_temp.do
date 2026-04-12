cd "D:\Users\hp\Documents\bibliothèque de travail\2 Elève ingénieur\ISE\Cours\ISE3\Tronc commun\UE15 Recherche\Groupe de travail\data\raw"
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
********************************************************************************
* Estimation de la régression: pour tester la significativité de la tendance et de la constante
********************************************************************************
** Création des différence
local vars "lnagri lnman lncom lntran lnheb lnimmo lnautre"

foreach var of local vars {
    gen d_`var' = `var' - L.`var'
}
************* Estimation des équations 

foreach var in lnagri lnman lncom lntran lnheb lnimmo lnautre {

    di "====================================="
    di "Résultats pour `var'"
    di "====================================="

    foreach p in 0 1 2 3 4 {

        di "----- Lag = `p' -----"

        if `p' == 0 {
            reg d_`var' L.`var' Annee
        }
        else {
            reg d_`var' L.`var' L(1/`p').d_`var' Annee
        }

    }
}
******** Résultats
********lnagri
**lag0:tendance et constante significative/
**lag1:tendance et constante significative/
**lag2:tendance et constante significative/
**lag3:tendance et constante non significative
**lag4:tendance et constante non significative
*****lnman: 
** lag 0: tendance et constante significative/
**lag1: tendance et constante significative/
**lag2: tendance et constante significative/
**lag3:tendance et constante significative/
**lag4:tendance et constante non significative
*****lncom
**lag0:tendance et constante significative/
**lag1:tendance et constante non significative
**lag2:tendance et constante significative
**lag3:tendance et constante non significative
**lag4:tendance et constante significative
*****lntran
**lag0:tendance et constante non significative
**lag1:tendance et constante significative
**lag2:tendance et constante significative
**lag3:tendance et constante significative
**lag4:tendance et constante significative
*****lnheb
**lag0:tendance et constante significative
**lag1:tendance et constante significative
**lag2:tendance et constante significative
**lag3:tendance et constante significative
**lag4:tendance et constante significative
*****lnimo
**lag0:tendance et constante non significative
**lag1:tendance et constante significative
**lag2:tendance et constante significative
**lag3:tendance et constante significative
**lag4:tendance et constante significative
*****lnautre 
**lag0:tendance et constante significative
**lag1:tendance et constante significative
**lag2:tendance et constante significative
**lag3:tendance et constante significative
**lag4:tendance et constante significative
********************************************************************************
*					Test ADF
********************************************************************************
* Boucle sur toutes les variables et leurs lags avec choix automatique du test ADF
* lnagri
dfuller lnagri, lags(0) trend            // lag0 : tendance + constante significative
dfuller lnagri, lags(1) trend            // lag1 : tendance + constante significative
dfuller lnagri, lags(2) trend            // lag2 : tendance + constante significative
dfuller lnagri, noconstant      // lag3 : constante non significative
dfuller lnagri, noconstant      // lag4 : constante non significative

* lnman
dfuller lnman, lags(0) trend
dfuller lnman, lags(1) trend
dfuller lnman, lags(2) trend
dfuller lnman, lags(3) trend
dfuller lnman, lags(4)  noconstant

* lncom
dfuller lncom, lags(0) trend
dfuller lncom, noconstant
dfuller lncom, lags(2) trend
dfuller lncom, noconstant
dfuller lncom, lags(4) trend

* lntran
dfuller lntran,  noconstant
dfuller lntran, lags(1) trend
dfuller lntran, lags(2) trend
dfuller lntran, lags(3) trend
dfuller lntran, lags(4) trend

* lnheb
dfuller lnheb, lags(0) trend
dfuller lnheb, lags(1) trend
dfuller lnheb, lags(2) trend
dfuller lnheb, lags(3) trend
dfuller lnheb, lags(4) trend

* lnimo
dfuller lnimmo, noconstant
dfuller lnimmo, lags(1) trend
dfuller lnimmo, lags(2) trend
dfuller lnimmo, lags(3) trend
dfuller lnimmo, lags(4) trend

* lnautre
dfuller lnautre, lags(0) trend
dfuller lnautre, lags(1) trend
dfuller lnautre, lags(2) trend
dfuller lnautre, lags(3) trend
dfuller lnautre, lags(4) trend




*** Détermination de p optimal 


* =====================================================================
*affichage direct varsoc variable par variable
* =====================================================================
foreach var in d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre {
    di _newline as txt ">> `var'"
    varsoc `var', maxlag(`maxlag')
}



local maxlag 4

* -- Ouvrir le fichier Excel et créer les en-têtes --------------------
putexcel set "resultats_lag_AIC_BIC.xlsx", sheet("Résultats") replace

putexcel A1 = "Variable"  B1 = "Lag"  C1 = "LL"  D1 = "LR"  E1 = "df" ///
         F1 = "p-value"   G1 = "AIC"  H1 = "BIC" I1 = "HQIC" J1= "SBIC"

local row = 2

foreach var in  d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre {

    di _newline as txt ">> `var'"
    varsoc `var', maxlag(`maxlag')
    matrix V = r(stats)

    * -- Écrire le nom de la variable sur la première ligne du bloc ----
    putexcel A`row' = "`var'"

    * -- Écrire chaque lag --------------------------------------------
    forvalues p =0/`maxlag' {
        putexcel A`row' = "`var'"       ///
                 B`row' = `p'           ///
                 C`row' = (V[`p'+1, 2])   ///
                 D`row' = (V[`p'+1, 3])   ///
                 E`row' = (V[`p'+1, 4])   ///
                 F`row' = (V[`p'+1, 5])   ///
                 G`row' = (V[`p'+1, 6])   ///
                 H`row' = (V[`p'+1, 7])   ///
                 I`row' = (V[`p'+1, 8]) ///
				 J`row' = (V[`p'+1, 9])
        local ++row
    }

    * -- Ligne vide entre les variables --------------------------------
    local ++row
}

di as txt _newline "Fichier Excel exporté : resultats_lag_AIC_BIC.xlsx"




















********************************************************************************
* TEST KPSS
********************************************************************************

* Créer une matrice pour stocker les résultats

matrix kpss_results = J(7, 5, .)

* Liste des variables
local vars "lnagri lnman lncom lntran lnheb lnimmo lnautre"
local names `" "Agriculture" "Manufacture" "Commerce" "Transport" "Hébergement" "Immobilier" "Autre" "'

* Boucle pour exécuter les tests et stocker les résultats
local row = 1
foreach var of local vars {
    forvalues lag = 0/4 {
        quietly kpss `var', maxlag(`lag') trend
        matrix kpss_results[`row', `lag'+1] = r(stat)
    }
    local row = `row' + 1
}

* Exporter vers Excel
putexcel set "resultats_kpss.xlsx", replace

* En-têtes
putexcel A1 = "Secteurs", bold
putexcel B1 = "KPSS(0)", bold
putexcel C1 = "KPSS(1)", bold
putexcel D1 = "KPSS(2)", bold
putexcel E1 = "KPSS(3)", bold
putexcel F1 = "KPSS(4)", bold

* Remplir les données
local row = 2
forvalues i = 1/7 {
    local name: word `i' of `names'
    putexcel A`row' = "`name'"
    putexcel B`row' = matrix(kpss_results[`i',1]), nformat(number_d2)
    putexcel C`row' = matrix(kpss_results[`i',2]), nformat(number_d2)
    putexcel D`row' = matrix(kpss_results[`i',3]), nformat(number_d2)
    putexcel E`row' = matrix(kpss_results[`i',4]), nformat(number_d2)
    putexcel F`row' = matrix(kpss_results[`i',5]), nformat(number_d2)
    local row = `row' + 1
}

* Ajouter une note
putexcel A9 = "Note: Valeurs critiques KPSS (avec tendance):"
putexcel A10 = "1% = 0.216, 5% = 0.146, 10% = 0.119"
putexcel A11 = "H0: Série stationnaire (rejet si stat > valeur critique)"
varsoc d_lnheb, maxlag(4)


****************************************************************************
* Test DFA sur la différence
*****************************************************************************

** Test de DFA

** Créer une matrice pour stocker les résultats
matrix results = J(7, 5, .)

* Liste des variables
local vars "d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre"
local names `" "Agriculture" "Manufacture" "Commerce" "Transport" "Hebergement" "Immobilier" "Autre" "'

* Boucle pour exécuter les tests et stocker les résultats
local row = 1
foreach var of local vars {
    forvalues lag = 0/3 {
        quietly dfuller `var', lags(`lag') trend
        matrix results[`row', `lag'+1] = r(Zt)
    }
    local row = `row' + 1
}

* Exporter vers Excel
putexcel set "resultats2_dfa.xlsx", replace

* En-têtes
putexcel A1 = "Secteurs", bold
putexcel B1 = "DFA(0)", bold
putexcel C1 = "DFA(1)", bold
putexcel D1 = "DFA(2)", bold
putexcel E1 = "DFA(3)", bold
putexcel F1 = "DFA(4)", bold

* Remplir les données
local row = 2
forvalues i = 1/7 {
    local name: word `i' of `names'
    putexcel A`row' = "`name'"
    putexcel B`row' = matrix(results[`i',1]), nformat(number_d2)
    putexcel C`row' = matrix(results[`i',2]), nformat(number_d2)
    putexcel D`row' = matrix(results[`i',3]), nformat(number_d2)
    putexcel E`row' = matrix(results[`i',4]), nformat(number_d2)
    putexcel F`row' = matrix(results[`i',5]), nformat(number_d2)
    local row = `row' + 1
}

* Ajouter une note
putexcel A9 = "Note: Valeurs critiques de MacKinnon (avec tendance):"
putexcel A10 = "1% = -4.22, 5% = -3.53, 10% = -3.20"

********* Cointégration
local vars "lnagri lnman lncom lntran lnheb lnimmo lnautre"
varsoc `vars'
vecrank `vars', lag(4) trace
vecrank `vars', lag(4) max
*ttttt

<<<<<<< HEAD
*************** Calcul des persistances globale


gen lninfor = ln(informel)
** On regresse la variable différenciée
gen dlninfor = D.lninfor
** Sélection de lag optimal
varsoc dlninfor // Le lag optimal est 0. On a une marche aléatoire


dfuller lninfor

dfuller dlninfor



ac dlninfor

pac dlninfor

twoway (line dlninfor Annee)

twoway (line lninfor Annee)


* Tester la robustesse
foreach p in 0 1 2 3 4 5 6 7{
    if `p' == 0 {
        di "AR(0) : P = 1"
=======
*************** Calcul des persistances globale univariée
***La persistance avec un AR
gen lninfor = ln(informel) //Comme la variable n'est pas stationnaire
** On regresse la variable différenciée
gen dlninfor = D.lninfor
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
>>>>>>> ba09e9fcdde3491ae62cd51df508979baccebc0a
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
<<<<<<< HEAD
} 
=======

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






>>>>>>> ba09e9fcdde3491ae62cd51df508979baccebc0a



