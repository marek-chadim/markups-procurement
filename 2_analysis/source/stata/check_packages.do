* Check which SSC packages needed by Phase 1 are available
foreach pkg in ddml pdslasso rlasso reghdfe ivreghdfe csdid sunab did_multiplegt teffects psmatch2 nnmatch oster bacondecomp xtabond2 estout binscatter sensemakr {
    cap which `pkg'
    if _rc == 0 {
        dis "[OK]    `pkg'"
    }
    else {
        dis "[MISS]  `pkg'"
    }
}
