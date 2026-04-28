cd "C:\Users\USER\Desktop\ISE3_2025\GT\base"

*cd "D:\Users\hp\Documents\bibliothèque de travail\2 Elève ingénieur\ISE\Cours\ISE3\Tronc commun\UE15 Recherche\Groupe de travail\data\raw"

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


/*===========================================================
  DELTA METHOD — se(P_hat) pour P_hat = sqrt(f̂(0)/σ²)
  Variable : dlninfor | AR(5)
===========================================================*/

local p = 5


quietly summarize resid_dlninfor
scalar T_eff  = r(N)
scalar sigma2 = r(Var)
scalar q_nw   = int(T_eff^(1/3))
local  q      = q_nw

/*-----------------------------------------------------------
  Autocovariances γ̂(0)...γ̂(2q) stockées dans GAM
-----------------------------------------------------------*/

local qmax = 2 * `q'
matrix GAM = J(`qmax'+1, 1, 0)

forvalues j = 0/`qmax' {
    if `j' == 0 {
        matrix GAM[1, 1] = sigma2
    }
    else {
        quietly corr resid_dlninfor L`j'.resid_dlninfor, covariance
        if r(N) > 0 matrix GAM[`j'+1, 1] = r(cov_12)
    }
}

/*-----------------------------------------------------------
  f̂(0) = γ(0) + 2·Σⱼ w(j)·γ(j)
-----------------------------------------------------------*/

scalar f0_hat = sigma2
forvalues j = 1/`q' {
    scalar w_j  = 1 - `j' / (q_nw + 1)
    scalar f0_hat = f0_hat + 2 * w_j * GAM[`j'+1, 1]
}
scalar P_hat = sqrt(abs(f0_hat / sigma2))

/*-----------------------------------------------------------
  Matrice Cov(γ̂(j), γ̂(k)) = (1/T)·Σₗ γ̂(|l+j|)·γ̂(|l+k|)
-----------------------------------------------------------*/

matrix Cov_GAM = J(`q'+1, `q'+1, 0)

forvalues j = 0/`q' {
    forvalues k = 0/`q' {
        scalar cov_jk = 0
        forvalues l = -`q'/`q' {
            local lj = abs(`l' + `j')
            local lk = abs(`l' + `k')
            if `lj' <= `qmax' & `lk' <= `qmax' {
                scalar cov_jk = cov_jk + GAM[`lj'+1,1] * GAM[`lk'+1,1]
            }
        }
        matrix Cov_GAM[`j'+1, `k'+1] = cov_jk / T_eff
    }
}

/*-----------------------------------------------------------
  Gradient ∂f̂(0)/∂γ(j) :
    j = 0 → 1
    j ≥ 1 → 2·w(j)
-----------------------------------------------------------*/

matrix grad_f0 = J(`q'+1, 1, 0)
matrix grad_f0[1, 1] = 1
forvalues j = 1/`q' {
    scalar w_j = 1 - `j' / (q_nw + 1)
    matrix grad_f0[`j'+1, 1] = 2 * w_j
}

/*-----------------------------------------------------------
  Var(f̂(0)) = grad' · Cov_GAM · grad
  se(P_hat) = se(f̂(0)) / (2 · P_hat · σ²)
-----------------------------------------------------------*/

matrix Var_f0        = grad_f0' * Cov_GAM * grad_f0
scalar se_f0         = sqrt(Var_f0[1,1])
scalar se_Phat       = se_f0 / (2 * P_hat * sigma2)
scalar IC_lo         = P_hat - 1.96 * se_Phat
scalar IC_hi         = P_hat + 1.96 * se_Phat

di _newline(1)
di "============================================"
di "  DELTA METHOD — P_hat = sqrt(f̂(0)/σ²)"
di "============================================"
di "  P_hat      : " %8.4f P_hat
di "  se(P_hat)  : " %8.4f se_Phat
di "  IC 95% bas : " %8.4f IC_lo
di "  IC 95% haut: " %8.4f IC_hi
di "============================================"


/*===========================================================
  BOOTSTRAP — ÉCART-TYPE DE LA PERSISTANCE UNIVARIÉE
  Variable : dlninfor | AR(5)
  Les deux mesures : P (somme des phi) et P_hat (densité)
===========================================================*/

/*===========================================================
  BOOTSTRAP CORRIGÉ — ÉCART-TYPE DE LA PERSISTANCE UNIVARIÉE
  Variable : dlninfor | AR(5)
===========================================================*/

local B   = 500
local p   = 5

/*-----------------------------------------------------------
  ÉTAPE 1 : Identifier les observations complètes
  L'AR(5) utilise les obs à partir de t = p+1
  → on rééchantillonne uniquement sur ces obs
-----------------------------------------------------------*/

* Valeurs ajustées sur le sample effectif
capture drop fitted_dlninfor
gen fitted_dlninfor = dlninfor - resid_dlninfor

* Identifier le range des obs non manquantes
quietly summarize resid_dlninfor
scalar T_eff = r(N)
scalar q_nw  = int(T_eff^(1/3))

* Trouver première et dernière obs non manquante
quietly count if resid_dlninfor != .
scalar n_complete = r(N)

* Stocker les résidus dans un vecteur (obs complètes seulement)
* pour le tirage avec remise
quietly gen obs_id = _n if resid_dlninfor != .
quietly summarize obs_id
scalar first_obs = r(min)
scalar last_obs  = r(max)

di "  Obs. complètes pour AR(`p') : " n_complete
di "  Première obs. : " first_obs "  |  Dernière : " last_obs
di "  Bandwidth q   : " q_nw
di "  P observé  (somme phi) : " %8.4f P
di "  P observé  (densité)   : " %8.4f P_hat

/*-----------------------------------------------------------
  ÉTAPE 2 : BOUCLE BOOTSTRAP
-----------------------------------------------------------*/

matrix BOOT_P     = J(`B', 1, .)
matrix BOOT_P_hat = J(`B', 1, .)

di _newline(1) "  Bootstrap en cours (B = `B')..."

preserve

forvalues b = 1/`B' {
    quietly {

        restore, preserve

        * --- a. Tirer n_complete indices dans [first_obs, last_obs] ---
        * avec remise, uniquement parmi les obs complètes
        gen resid_star = .
        forvalues t = 1/`=n_complete' {
            * indice aléatoire dans les obs complètes
            scalar idx = first_obs + floor(runiform() * n_complete)
            * si idx dépasse last_obs, recadrer
            if idx > last_obs scalar idx = last_obs
            replace resid_star = resid_dlninfor[idx] ///
                in `=first_obs + `t' - 1'
        }

        * --- b. Reconstruire dlninfor* sur les obs complètes ---
        replace dlninfor = fitted_dlninfor + resid_star ///
            if obs_id >= first_obs & obs_id <= last_obs

        * --- c. Ré-estimer AR(5) ---
        capture drop resid_b
        regress dlninfor L(1/`p').dlninfor ///
            if obs_id >= first_obs & obs_id <= last_obs

        * Vérifier qu'on a assez d'obs
        if e(N) < `p' + 2 {
            drop resid_star
            continue
        }

        * P* via somme des phi
        scalar sp = 0
        forvalues k = 1/`p' {
            scalar sp = sp + _b[L`k'.dlninfor]
        }

        * Eviter division par zéro ou explosion
        if abs(1 - sp) < 1e-6 {
            matrix BOOT_P[`b', 1] = .
        }
        else {
            matrix BOOT_P[`b', 1] = abs(1 / (1 - sp))
        }

        * Résidus pour densité spectrale
        predict resid_b if obs_id >= first_obs & obs_id <= last_obs, residuals

        * P*_hat via densité spectrale
        summarize resid_b
        if r(N) < q_nw + 2 {
            matrix BOOT_P_hat[`b', 1] = .
            drop resid_star
            continue
        }

        scalar sig2_b = r(Var)
        scalar f0_b   = sig2_b

        forvalues j = 1/`=q_nw' {
            scalar w_j = 1 - `j' / (q_nw + 1)
            capture corr resid_b L`j'.resid_b, covariance
            if _rc == 0 {
                scalar f0_b = f0_b + 2 * w_j * r(cov_12)
            }
        }

        matrix BOOT_P_hat[`b', 1] = sqrt(abs(f0_b / sig2_b))

        drop resid_star
    }

    if mod(`b', 100) == 0 di "  ... `b'/`B'"
}

restore

/*-----------------------------------------------------------
  ÉTAPE 3 : Calcul écarts-types et IC 95%
  On ignore les réplications manquantes (.)
-----------------------------------------------------------*/

svmat BOOT_P,     name(boot_P_tmp)
svmat BOOT_P_hat, name(boot_Phat_tmp)

* --- P via somme phi ---
quietly summarize boot_P_tmp1
scalar se_P     = r(sd)
scalar mean_P_b = r(mean)
scalar B_valid_P = r(N)
quietly centile boot_P_tmp1, centile(2.5 97.5)
scalar ic_lo_P  = r(c_1)
scalar ic_hi_P  = r(c_2)

* --- P via densité ---
quietly summarize boot_Phat_tmp1
scalar se_Phat     = r(sd)
scalar mean_Phat_b = r(mean)
scalar B_valid_Phat = r(N)
quietly centile boot_Phat_tmp1, centile(2.5 97.5)
scalar ic_lo_Phat  = r(c_1)
scalar ic_hi_Phat  = r(c_2)

drop boot_P_tmp1 boot_Phat_tmp1

/*-----------------------------------------------------------
  ÉTAPE 4 : Affichage
-----------------------------------------------------------*/

di _newline(1)
di "=========================================================="
di "  RÉSULTATS BOOTSTRAP — Persistance univariée dlninfor"
di "  AR(`p') | B = `B' | q = " q_nw
di "=========================================================="
di "  " %25s " " %8s "P obs." %10s "Moy.boot" ///
   %8s "se(P)" %22s "IC 95%"  %8s "B valides"
di "  " "{hline 85}"
di "  " %25s "Somme phi" ///
   %8.4f P %10.4f mean_P_b %8.4f se_P ///
   "   [" %6.4f ic_lo_P " ; " %6.4f ic_hi_P "]" ///
   %8.0f B_valid_P
di "  " %25s "Densité spectrale" ///
   %8.4f P_hat %10.4f mean_Phat_b %8.4f se_Phat ///
   "   [" %6.4f ic_lo_Phat " ; " %6.4f ic_hi_Phat "]" ///
   %8.0f B_valid_Phat
di "  " "{hline 85}"
di "  Biais (somme phi)  : " %8.4f mean_P_b   - P
di "  Biais (densité)    : " %8.4f mean_Phat_b - P_hat
di "=========================================================="

capture drop fitted_dlninfor obs_id




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
  DELTA METHOD — se(P_i) pour chaque secteur
  P_i = sqrt(f̂(0)ᵢᵢ / σ²ᵢ)
  se(P_i) = se(f̂(0)ᵢᵢ) / (2 · P_i · σ²ᵢ)
===========================================================*/

local q = q_nw

di _newline(1)
di "============================================================"
di "  DELTA METHOD — se(P_i) | VAR(4) | 7 secteurs | q=" q_nw
di "============================================================"
di "  " %15s "Secteur" %8s "P_i" %10s "se(P_i)" %22s "IC 95%" 
di "  " "{hline 60}"

foreach var of local secteurs {

    local qmax = 2 * `q'

    /*-------------------------------------------------------
      Autocovariances γ̂(0)...γ̂(2q) du résidu du secteur
    -------------------------------------------------------*/
    matrix GAM_`var' = J(`qmax'+1, 1, 0)
    matrix GAM_`var'[1, 1] = sigma2_`var'

    forvalues j = 1/`qmax' {
        quietly corr resid_`var' L`j'.resid_`var', covariance
        if r(N) > 0 matrix GAM_`var'[`j'+1, 1] = r(cov_12)
    }

    /*-------------------------------------------------------
      Matrice Cov(γ̂(j), γ̂(k)) = (1/T)·Σₗ γ̂(|l+j|)·γ̂(|l+k|)
    -------------------------------------------------------*/
    matrix Cov_GAM_`var' = J(`q'+1, `q'+1, 0)

    forvalues j = 0/`q' {
        forvalues k = 0/`q' {
            scalar cov_jk = 0
            forvalues l = -`q'/`q' {
                local lj = abs(`l' + `j')
                local lk = abs(`l' + `k')
                if `lj' <= `qmax' & `lk' <= `qmax' {
                    scalar cov_jk = cov_jk + ///
                        GAM_`var'[`lj'+1,1] * GAM_`var'[`lk'+1,1]
                }
            }
            matrix Cov_GAM_`var'[`j'+1, `k'+1] = cov_jk / T_eff
        }
    }

    /*-------------------------------------------------------
      Gradient ∂f̂(0)/∂γ(j) :
        j = 0 → 1
        j ≥ 1 → 2·w(j)
    -------------------------------------------------------*/
    matrix grad_`var' = J(`q'+1, 1, 0)
    matrix grad_`var'[1, 1] = 1
    forvalues j = 1/`q' {
        scalar w_j = 1 - `j' / (q_nw + 1)
        matrix grad_`var'[`j'+1, 1] = 2 * w_j
    }

    /*-------------------------------------------------------
      Var(f̂(0)) = grad' · Cov_GAM · grad
      se(P_i)   = se(f̂(0)) / (2 · P_i · σ²ᵢ)
    -------------------------------------------------------*/
    matrix Var_f0_`var' = grad_`var'' * Cov_GAM_`var' * grad_`var'
    scalar se_f0_`var'  = sqrt(abs(Var_f0_`var'[1,1]))
    scalar se_P_`var'   = se_f0_`var' / (2 * P_`var' * sigma2_`var')
    scalar IC_lo_`var'  = P_`var' - 1.96 * se_P_`var'
    scalar IC_hi_`var'  = P_`var' + 1.96 * se_P_`var'

    local nom : word `i' of `noms'
    di "  " %15s "`var'" %8.4f P_`var' %10.4f se_P_`var' ///
       "   [" %6.4f IC_lo_`var' " ; " %6.4f IC_hi_`var' "]"
}

di "  " "{hline 60}"
di "  IC 95% = P_i ± 1.96·se(P_i)"
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


/*===========================================================
  DELTA METHOD — se(P_i) pour le modèle contraint
  Même logique que le VAR : on travaille sur les résidus
  déjà disponibles (resid_`var')
===========================================================*/

local q = q_nw

di _newline(1)
di "============================================================"
di "  DELTA METHOD — se(P_i) | Modèle contraint | q=" q_nw
di "============================================================"
di "  " %15s "Secteur" %8s "P_i" %10s "se(P_i)" %22s "IC 95%"
di "  " "{hline 60}"

local i = 1
foreach var of local secteurs {

    local qmax = 2 * `q'

    /*-------------------------------------------------------
      Autocovariances γ̂(0)...γ̂(2q)
    -------------------------------------------------------*/
    matrix GAM_`var' = J(`qmax'+1, 1, 0)
    matrix GAM_`var'[1, 1] = sigma2_`var'

    forvalues j = 1/`qmax' {
        quietly corr resid_`var' L`j'.resid_`var', covariance
        if r(N) > 0 matrix GAM_`var'[`j'+1, 1] = r(cov_12)
    }

    /*-------------------------------------------------------
      Cov(γ̂(j), γ̂(k)) = (1/T)·Σₗ γ̂(|l+j|)·γ̂(|l+k|)
    -------------------------------------------------------*/
    matrix Cov_GAM_`var' = J(`q'+1, `q'+1, 0)

    forvalues j = 0/`q' {
        forvalues k = 0/`q' {
            scalar cov_jk = 0
            forvalues l = -`q'/`q' {
                local lj = abs(`l' + `j')
                local lk = abs(`l' + `k')
                if `lj' <= `qmax' & `lk' <= `qmax' {
                    scalar cov_jk = cov_jk + ///
                        GAM_`var'[`lj'+1,1] * GAM_`var'[`lk'+1,1]
                }
            }
            matrix Cov_GAM_`var'[`j'+1, `k'+1] = cov_jk / T_eff
        }
    }

    /*-------------------------------------------------------
      Gradient : ∂f̂(0)/∂γ(0) = 1  |  ∂f̂(0)/∂γ(j) = 2w(j)
    -------------------------------------------------------*/
    matrix grad_`var' = J(`q'+1, 1, 0)
    matrix grad_`var'[1, 1] = 1
    forvalues j = 1/`q' {
        scalar w_j = 1 - `j' / (q_nw + 1)
        matrix grad_`var'[`j'+1, 1] = 2 * w_j
    }

    /*-------------------------------------------------------
      Var(f̂(0)) = grad' · Cov_GAM · grad
      se(P_i)   = se(f̂(0)) / (2 · P_i · σ²ᵢ)
    -------------------------------------------------------*/
    matrix Var_f0_`var'  = grad_`var'' * Cov_GAM_`var' * grad_`var'
    scalar se_f0_`var'   = sqrt(abs(Var_f0_`var'[1,1]))
    scalar se_P_`var'    = se_f0_`var' / (2 * P_`var' * sigma2_`var')
    scalar IC_lo_`var'   = P_`var' - 1.96 * se_P_`var'
    scalar IC_hi_`var'   = P_`var' + 1.96 * se_P_`var'

    local nom : word `i' of `noms'
    di "  " %15s "`nom'" %8.4f P_`var' %10.4f se_P_`var' ///
       "   [" %6.4f IC_lo_`var' " ; " %6.4f IC_hi_`var' "]"

    local i = `i' + 1
}

di "  " "{hline 60}"
di "  IC 95% = P_i ± 1.96·se(P_i)"
di "============================================================"

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
*scalar ll_NC = e(ll_dfk)

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


/*===========================================================
  MODÈLE CONTRAINT SÉLECTIF — SEUIL |t| > 1.5
  Pour chaque secteur i :
    1. Régression complète : L(1/4).vi + L(1/4).vj pour tout j≠i
    2. Pour chaque j≠i : tester si AU MOINS UN lag a |t| > 1.5
    3. Construire reste_sig_i = somme des vj retenus
    4. Ré-estimer : L(1/4).vi + L(1/4).reste_sig_i
    5. Persistance via densité spectrale
===========================================================*/

local secteurs "d_lnagri d_lnman d_lncom d_lntran d_lnheb d_lnimmo d_lnautre"
local noms     "Agri Manufact Commerce Transport Hébergement Immo Autres"
local m   = 7
local p   = 4
local seuil_t = 1.5

/*-----------------------------------------------------------
  ÉTAPE 1 : Pour chaque secteur i, régression complète
  puis vérifier pour chaque j≠i si au moins un lag
  de j a |t| > 1.5
-----------------------------------------------------------*/

di _newline(1)
di "============================================================"
di "  SÉLECTION — Au moins un lag |t| > `seuil_t' par variable"
di "============================================================"

local i = 1
foreach vi of local secteurs {

    local nom_i : word `i' of `noms'
    di _newline(1) "  >>> Équation : `nom_i' (`vi')"
    di "  " "{hline 55}"

    * Régression complète : propres lags + lags de chaque j≠i
    *quietly regress `vi' L(1/`p').`vi' ///
	 regress `vi' L(1/`p').`vi' ///
        L(1/`p').d_lnagri L(1/`p').d_lnman  L(1/`p').d_lncom  ///
        L(1/`p').d_lntran L(1/`p').d_lnheb  L(1/`p').d_lnimmo ///
        L(1/`p').d_lnautre

    * Pour chaque j≠i : vérifier si au moins un lag |t| > seuil
    local sig_vars_`vi' ""
    local j = 1
    foreach vj of local secteurs {

        if "`vi'" != "`vj'" {

            local nom_j : word `j' of `noms'

            * Tester chaque lag de vj
            local at_least_one = 0
            local t_vals ""

            forvalues k = 1/`p' {
                scalar t_k = _b[L`k'.`vj'] / _se[L`k'.`vj']
				di "nn"
				
				local t_vals "`t_vals' `=string(abs(t_k), "%5.2f")'"
                *local t_vals "`t_vals' " %5.2f abs(t_k)

                if abs(t_k) > `seuil_t' {
                    local at_least_one = 1
                }
            }
			

            if `at_least_one' == 1 {
                local sig_vars_`vi' "`sig_vars_`vi'' `vj'"
                di "  RETENU  : `nom_j' — |t| max > `seuil_t'"
            }
            else {
                di "  exclu   : `nom_j' — aucun lag |t| > `seuil_t'"
            }
        }
        local j = `j' + 1
    }

    if "`sig_vars_`vi''" == "" {
        di "  → Aucune variable externe retenue : AR(`p') pur"
    }
    else {
        di "  → Variables retenues : `sig_vars_`vi''"
    }

    local i = `i' + 1
}


/*-----------------------------------------------------------
  ÉTAPE 2 : Construire reste_sig_i = somme des vj retenus
-----------------------------------------------------------*/

di _newline(1)
di "============================================================"
di "  ÉTAPE 2 : CONSTRUCTION DES AGRÉGATS SÉLECTIFS"
di "============================================================"

local i = 1
foreach vi of local secteurs {

    local nom_i : word `i' of `noms'
    capture drop reste_sig_`vi'

    if "`sig_vars_`vi''" == "" {
        gen reste_sig_`vi' = 0
        di "  `nom_i' : reste_sig = 0"
    }
    else {
        gen reste_sig_`vi' = 0
        foreach vj of local sig_vars_`vi' {
            replace reste_sig_`vi' = reste_sig_`vi' + `vj'
        }
        di "  `nom_i' : reste_sig = Σ(`sig_vars_`vi'')"
    }
    local i = `i' + 1
}

/*-----------------------------------------------------------
  ÉTAPE 3 : Ré-estimation contrainte
  d_lni = μ + Σα·L(1/4).d_lni + Σβ·L(1/4).reste_sig_i + ε
-----------------------------------------------------------*/

di _newline(1)
di "============================================================"
di "  ÉTAPE 3 : ESTIMATION CONTRAINTE SÉLECTIVE"
di "============================================================"

local i = 1
foreach vi of local secteurs {

    local nom_i : word `i' of `noms'
    di _newline(1) "  >>> `nom_i' (`vi')"

    capture drop resid_sel_`vi'

    if "`sig_vars_`vi''" == "" {
        regress `vi' L(1/`p').`vi'
        di "    Modèle : AR(`p') pur"
    }
    else {
        regress `vi' L(1/`p').`vi' L(1/`p').reste_sig_`vi'
        di "    Modèle : AR(`p') + L(1/`p').reste_sig"
        di "    Composantes de reste_sig : `sig_vars_`vi''"
    }

    quietly predict resid_sel_`vi' if e(sample), residuals
    scalar k_sel_`vi' = e(df_m)

    di "    R²         = " %6.4f e(r2)
    di "    Paramètres = " e(df_m)
    di "    AIC        = " %8.2f e(N)*ln(e(rss)/e(N)) + 2*e(df_m)

    local i = `i' + 1
}

/*-----------------------------------------------------------
  ÉTAPE 4 : Persistance via densité spectrale (Bartlett)
  P²ᵢ = f̂(0)ᵢᵢ / σ²ᵢ
-----------------------------------------------------------*/

quietly regress `= word("`secteurs'",1)' L(1/`p').`= word("`secteurs'",1)'
scalar T_eff = e(N)
scalar q_nw  = int(T_eff^(1/3))

di _newline(1)
di "============================================================"
di "  ÉTAPE 4 : PERSISTANCE — DENSITÉ SPECTRALE (Bartlett/NW)"
di "  Seuil |t| > `seuil_t' | q = " q_nw
di "============================================================"
di "  " %15s "Secteur" %14s "σ²ᵢ" %14s "f̂(0)ᵢᵢ" ///
   %8s "P²" %8s "P" "   Interprétation"
di "  " "{hline 75}"

local i = 1
foreach vi of local secteurs {

    quietly summarize resid_sel_`vi'
    scalar sigma2_`vi' = r(Var)
    scalar f0_`vi'     = sigma2_`vi'

    forvalues j = 1/`=q_nw' {
        scalar w_j = 1 - `j' / (q_nw + 1)
        quietly corr resid_sel_`vi' L`j'.resid_sel_`vi', covariance
        scalar f0_`vi' = f0_`vi' + 2 * w_j * r(cov_12)
    }

    scalar P_sq_`vi' = f0_`vi' / sigma2_`vi'
    scalar P_`vi'    = sqrt(abs(P_sq_`vi'))

    local nom_i : word `i' of `noms'
    local interp ""
    if      abs(P_`vi' - 1) < 0.05  local interp "Choc permanent (≈RW)"
    else if P_`vi' > 1               local interp "Amplification"
    else                             local interp "Atténuation partielle"

    di "  " %15s "`nom_i'" %14.6f sigma2_`vi' %14.6f f0_`vi' ///
       %8.4f P_sq_`vi' %8.4f P_`vi' "   `interp'"

    local i = `i' + 1
}

di "  " "{hline 75}"


/*===========================================================
  DELTA METHOD — se(P_i) | Modèle sélectif |t| > 1.5
===========================================================*/

local q = q_nw

di _newline(1)
di "============================================================"
di "  DELTA METHOD — se(P_i) | Modèle sélectif | q=" q_nw
di "============================================================"
di "  " %15s "Secteur" %8s "P_i" %10s "se(P_i)" %22s "IC 95%"
di "  " "{hline 60}"

local i = 1
foreach vi of local secteurs {

    local qmax = 2 * `q'

    /*-------------------------------------------------------
      Autocovariances γ̂(0)...γ̂(2q) sur resid_sel_`vi'
    -------------------------------------------------------*/
    matrix GAM_`vi' = J(`qmax'+1, 1, 0)
    matrix GAM_`vi'[1, 1] = sigma2_`vi'

    forvalues j = 1/`qmax' {
        quietly corr resid_sel_`vi' L`j'.resid_sel_`vi', covariance
        if r(N) > 0 matrix GAM_`vi'[`j'+1, 1] = r(cov_12)
    }

    /*-------------------------------------------------------
      Cov(γ̂(j), γ̂(k)) = (1/T)·Σₗ γ̂(|l+j|)·γ̂(|l+k|)
    -------------------------------------------------------*/
    matrix Cov_GAM_`vi' = J(`q'+1, `q'+1, 0)

    forvalues j = 0/`q' {
        forvalues k = 0/`q' {
            scalar cov_jk = 0
            forvalues l = -`q'/`q' {
                local lj = abs(`l' + `j')
                local lk = abs(`l' + `k')
                if `lj' <= `qmax' & `lk' <= `qmax' {
                    scalar cov_jk = cov_jk + ///
                        GAM_`vi'[`lj'+1,1] * GAM_`vi'[`lk'+1,1]
                }
            }
            matrix Cov_GAM_`vi'[`j'+1, `k'+1] = cov_jk / T_eff
        }
    }

    /*-------------------------------------------------------
      Gradient : ∂f̂(0)/∂γ(0) = 1  |  ∂f̂(0)/∂γ(j) = 2w(j)
    -------------------------------------------------------*/
    matrix grad_`vi' = J(`q'+1, 1, 0)
    matrix grad_`vi'[1, 1] = 1
    forvalues j = 1/`q' {
        scalar w_j = 1 - `j' / (q_nw + 1)
        matrix grad_`vi'[`j'+1, 1] = 2 * w_j
    }

    /*-------------------------------------------------------
      Var(f̂(0)) = grad' · Cov_GAM · grad
      se(P_i)   = se(f̂(0)) / (2 · P_i · σ²ᵢ)
    -------------------------------------------------------*/
    matrix Var_f0_`vi' = grad_`vi'' * Cov_GAM_`vi' * grad_`vi'
    scalar se_f0_`vi'  = sqrt(abs(Var_f0_`vi'[1,1]))
    scalar se_P_`vi'   = se_f0_`vi' / (2 * P_`vi' * sigma2_`vi')
    scalar IC_lo_`vi'  = P_`vi' - 1.96 * se_P_`vi'
    scalar IC_hi_`vi'  = P_`vi' + 1.96 * se_P_`vi'

    local nom_i : word `i' of `noms'
    di "  " %15s "`nom_i'" %8.4f P_`vi' %10.4f se_P_`vi' ///
       "   [" %6.4f IC_lo_`vi' " ; " %6.4f IC_hi_`vi' "]"

    local i = `i' + 1
}

di "  " "{hline 60}"
di "  IC 95% = P_i ± 1.96·se(P_i)"
di "============================================================"

/*-----------------------------------------------------------
  ÉTAPE 5 : TEST LR vs VAR(4) non contraint
-----------------------------------------------------------*/

quietly var `secteurs', lags(1/`p')
scalar ll_NC = e(ll)
*scalar ll_NC = e(ll_dfk)
scalar K_NC  = `m' * (`m'*`p' + 1)

matrix Sigma_sel = J(`m', `m', 0)
local i = 1
foreach vi of local secteurs {
    local j = 1
    foreach vj of local secteurs {
        quietly corr resid_sel_`vi' resid_sel_`vj', covariance
        matrix Sigma_sel[`i', `j'] = r(cov_12)
        local j = `j' + 1
    }
    local i = `i' + 1
}

scalar ll_sel = -(T_eff/2)*(`m'*ln(2*_pi) + ln(det(Sigma_sel)) + `m')

scalar K_sel = 0
foreach vi of local secteurs {
    scalar K_sel = K_sel + k_sel_`vi' + 1
}
scalar r_ddl = K_NC - K_sel
scalar LR    = 2 * (ll_NC - ll_sel)
scalar pv    = chi2tail(r_ddl, LR)

di _newline(1)
di "============================================================"
di "  TEST LR : Modèle sélectif vs VAR(4) non contraint"
di "============================================================"
di "  lnL NC       : " %10.4f ll_NC
di "  lnL sélectif : " %10.4f ll_sel
di "  LR stat.     : " %10.4f LR
di "  DDL          : " r_ddl
di "  p-valeur     : " %10.4f pv
di "  Seuil 5%     : " %10.4f invchi2tail(r_ddl, 0.05)
di "  " "{hline 45}"
if pv >= 0.05 {
    di "  → Non-rejet H0 : modèle sélectif VALIDE"
}
else {
    di "  → Rejet H0 : modèle sélectif non validé"
}
di "============================================================"

* Nettoyage
foreach vi of local secteurs {
    capture drop reste_sig_`vi' resid_sel_`vi'
}




