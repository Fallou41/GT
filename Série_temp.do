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