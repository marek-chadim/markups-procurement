* Install missing SSC packages from network (non-interactive)
set checksum off
foreach pkg in ddml pdslasso lassopack ivreghdfe psmatch2 nnmatch oster sensemakr_scc {
    cap ssc install `pkg', replace
    if _rc == 0 dis "[installed] `pkg'"
    else dis "[fail ssc] `pkg' rc=" _rc
}

* Some live on GitHub only — try SJ-archive/GitHub installers
cap net install sensemakr, from("https://raw.githubusercontent.com/chadhazlett/sensemakr/master/stata") replace
cap net install did_multiplegt, from("https://raw.githubusercontent.com/chaisemartinPackages/did_multiplegt/master") replace

* Verify
foreach pkg in ddml pdslasso rlasso ivreghdfe psmatch2 nnmatch oster sensemakr did_multiplegt {
    cap which `pkg'
    if _rc == 0 dis "[OK] `pkg'"
    else dis "[MISS] `pkg'"
}
