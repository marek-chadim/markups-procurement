use "../../output/stata/analysis_panel.dta", clear
dis "=== Searching for profit/profitability variables ==="
ds *profit* *ebit* *income* *earn* *return* *roe* *roa* *margin*
dis "=== Full variable list (first 60) ==="
ds, detail
