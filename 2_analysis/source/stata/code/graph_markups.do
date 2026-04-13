*===============================================================================
* graph_markups.do — shared Stata graph style for markups-procurement figures
*
* Stata analog of theme_markups.R and style_markups.py. Defines the Paul Tol
* bright palette as global macros, provides consistent plot-region styling,
* and exposes default export options. Include this at the top of any .do file
* producing figures so R, Python, and Stata all share the same palette.
*
* Usage:
*     do "$code/graph_markups.do"
*     twoway (line y x, lcolor("${markups_blue}")) ///
*            (line z x, lcolor("${markups_pink}") lpattern(dash)) ///
*         , ${markups_gropts} ytitle("...") xtitle("...")
*     graph export "$output/fig.pdf", replace
*
* Reference: Healy (2026), Data Visualization: A Practical Introduction,
* Princeton UP, 2nd ed., Chapter 8 (Refine your Plots).
*===============================================================================

*-------------------------------------------------------------------------------
* Paul Tol bright palette (colorblind-safe, up to 7 categories)
* Source: https://personal.sron.nl/~pault/
* Matches tol_bright in theme_markups.R and MARKUPS_COLORS in style_markups.py
*
* Stata 18+ accepts hex strings directly via lcolor("#4477AA").
* Older versions need space-separated RGB via lcolor("68 119 170").
*-------------------------------------------------------------------------------

* Hex (Stata 18+)
global tol_blue_hex    "#4477AA"
global tol_red_hex     "#EE6677"
global tol_green_hex   "#228833"
global tol_yellow_hex  "#CCBB44"
global tol_cyan_hex    "#66CCEE"
global tol_purple_hex  "#AA3377"
global tol_grey_hex    "#BBBBBB"

* RGB (Stata 15-17 compatible)
global tol_blue_rgb    "68 119 170"
global tol_red_rgb     "238 102 119"
global tol_green_rgb   "34 136 51"
global tol_yellow_rgb  "204 187 68"
global tol_cyan_rgb    "102 204 238"
global tol_purple_rgb  "170 51 119"
global tol_grey_rgb    "187 187 187"

*-------------------------------------------------------------------------------
* Named shortcuts (use these in figure scripts for readability).
* Defaults to hex; swap to ${tol_blue_rgb} etc. if on older Stata.
*-------------------------------------------------------------------------------
global markups_blue   "${tol_blue_hex}"
global markups_pink   "${tol_purple_hex}"
global markups_red    "${tol_red_hex}"
global markups_green  "${tol_green_hex}"
global markups_yellow "${tol_yellow_hex}"
global markups_grey   "${tol_grey_hex}"

*-------------------------------------------------------------------------------
* Consistent plot region / graph region styling.
* Clean white background, modest margins, no heavy borders. Pass via ${markups_gropts}
* to any twoway, scatter, histogram, binscatter, or coefplot call.
*-------------------------------------------------------------------------------
global markups_gropts "graphregion(color(white) margin(medsmall)) plotregion(color(white) margin(medsmall)) bgcolor(white)"

*-------------------------------------------------------------------------------
* Export options
* - PDF is vector: DPI is nominal; Stata handles font embedding automatically.
* - PNG/TIF needs explicit width to approximate 300 DPI for a 6.5in figure:
*   300 dpi * 6.5 in = 1950 px. Round to 2000.
*-------------------------------------------------------------------------------
global markups_export_raster "width(2000) height(1400)"

di _newline "graph_markups.do: loaded Paul Tol palette + markups_gropts"
