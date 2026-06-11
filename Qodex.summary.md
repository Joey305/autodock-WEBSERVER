# Qodex.summary

## Task
Expand `7_Graph.py` into a richer docking analytics and PDF reporting tool.

## Original Goal
The user wants to expand graph generation beyond basic histograms and boxplots so docking results can show robust ligand rankings, hit rates, tautomer/protomer/conformer state performance, and a wrapped PDF report.

## Assumptions
- More negative Vina affinities are always better, so ranking metrics treat lower affinity values as stronger performance.
- A consensus rank should balance best-score performance with robustness signals such as median score, top-k means, hit rates, and outlier gap.
- Existing users may already depend on the older `Plots/*.csv` outputs, so new `summary_tables/` outputs should be additive rather than replacing the old locations.
- Plotting and PDF generation should degrade cleanly when matplotlib is unavailable instead of making summary CSV generation fail.

## Files Inspected
- `/Users/jxs794/Documents/autodock-WEBSERVER/7_Graph.py`: inspected the existing summary/statistics flow, per-ligand plots, state summaries, and CLI.
- `/Users/jxs794/Documents/autodock-WEBSERVER/ligand_naming.py`: confirmed the existing ligand/state normalization helpers.
- `/Users/jxs794/Documents/autodock-WEBSERVER/ligand_manifest.py`: confirmed manifest merge helpers and chemical metadata columns used by downstream summaries.
- `/Users/jxs794/Documents/autodock-WEBSERVER/tests/test_ligand_provenance.py`: inspected the existing graph/state summary tests and extended them.
- `/Users/jxs794/Documents/autodock-WEBSERVER/runner_templates.py`: updated generated portable package docs.
- `/Users/jxs794/Documents/autodock-WEBSERVER/README.md`: updated user-facing documentation for the expanded analytics workflow.

## Files Changed
- `/Users/jxs794/Documents/autodock-WEBSERVER/7_Graph.py`: substantially expanded the analytics pipeline, new metrics, ranking logic, heatmaps, manifest JSON output, optional PDF reporting, and compatibility fixes.
- `/Users/jxs794/Documents/autodock-WEBSERVER/tests/test_ligand_provenance.py`: added graph summary, consensus ranking, manifest, and output-generation regression coverage.
- `/Users/jxs794/Documents/autodock-WEBSERVER/runner_templates.py`: documented the richer graph/report workflow in generated package docs.
- `/Users/jxs794/Documents/autodock-WEBSERVER/README.md`: documented new ranking metrics, hit-rate logic, consensus ranking, and the optional PDF workflow.

## Files Created
- No new source files were created in this pass.

## Implementation Summary
`7_Graph.py` now keeps the original histogram and boxplot workflow but adds a richer analytics layer around the docking score CSV. The script now computes stronger per-ligand ranking metrics, writes expanded ligand/variant/state summary CSVs, creates ranking and hit-rate visualization assets when matplotlib is available, writes a machine-readable `Plots/graph_report_manifest.json`, and can optionally assemble a human-readable PDF report with summary tables and plots using `PdfPages`.

The summary logic now distinguishes between “best single pose” and “consistently good ligand” performance. Ligand-level rows include best, worst, mean, median, standard deviation, quantiles, interquartile range, score range, top-3/5/10 mean and median, hit counts/rates across configurable thresholds, outlier gap, and a consensus rank. Variant, tautomer, protomer, and state summaries reuse the same core robustness metrics and keep provenance fields such as `best_ligand_variant`, `best_state_tag`, `best_vina_pose`, `best_outfile`, and chemical metadata where present.

## Key Decisions
- Selected summary statistics: added quantiles, IQR, top-k means/medians, hit counts/rates, outlier gap, and robustness score so users can compare both peak scores and stability.
- Consensus ranking formula: ranks best affinity, median affinity, top-5 mean, top-10 mean, hit rate at `<= -9`, and outlier gap within each receptor, then averages those ranks into `consensus_rank_score`.
- Hit-rate thresholds: default thresholds are `-8`, `-8.5`, `-9`, `-9.5`, and `-10`, matching common “how many genuinely strong poses did this ligand produce?” review questions.
- State/tautomer heatmap design: rows are `LigandBase`, columns are `StateTag`, and the displayed value is configurable as `best`, `median`, or `top5_mean`.
- PDF report structure: cover page, executive summary table, receptor plot sections, per-ligand appendix, and output paths page using `matplotlib.backends.backend_pdf.PdfPages`.
- Output folder organization: new authoritative summary CSVs live under `Plots/summary_tables/`, but compatibility copies are still written to the top-level `Plots/` directory.
- Backward compatibility: original CLI usage, original histogram/boxplot outputs, and original summary locations were preserved.

## Commands Run
- `git status --short`
- `sed -n '1,260p' 7_Graph.py`
- `sed -n '260,620p' 7_Graph.py`
- `sed -n '620,1040p' 7_Graph.py`
- `rg -n "build_summary|_summarize_group|write_state_summary_tables|plot_|boxplot|histogram|summary|variant|state|Binding_Affinity|LigandBase|LigandVariant|TautomerTag|StateTag|argparse|PdfPages|FutureWarning|labels" 7_Graph.py tests README.md runner_templates.py`
- `sed -n '280,390p' tests/test_ligand_provenance.py`
- `python3 -m py_compile 7_Graph.py ligand_naming.py ligand_manifest.py runner_templates.py`
- `python3 -m unittest discover -s tests -v`
- `find "/Users/jxs794/Downloads/job 3" -maxdepth 3 -type f -name "*vina_docking_scores_sorted.csv" | sort`
- `python3 7_Graph.py --csv "/Users/jxs794/Downloads/job 3/Docking_Results_Receptors_Ligands_Ligands_PDBQT_8Poses_20260611_1718_10Poses_2026-06-11_17-24-06_2026-06-11_18-23-06_vina_docking_scores_sorted.csv" --top-n 4 --write-state-summaries --write-ranking-plots --write-hit-rate-plots --write-state-heatmaps --write-pdf-report --outdir "$TMPDIR/Plots"`
- `python3 - <<'PY' ...` to inspect the generated `graph_report_manifest.json`

## Validation Results
- Syntax checks passed for:
  - `/Users/jxs794/Documents/autodock-WEBSERVER/7_Graph.py`
  - `/Users/jxs794/Documents/autodock-WEBSERVER/ligand_naming.py`
  - `/Users/jxs794/Documents/autodock-WEBSERVER/ligand_manifest.py`
  - `/Users/jxs794/Documents/autodock-WEBSERVER/runner_templates.py`
- Unit tests passed:
  - `48` tests total
  - `7` skipped because RDKit-dependent tests are intentionally skipped when RDKit is unavailable
- Added regression coverage for:
  - richer ligand summary metrics
  - consensus ranking behavior
  - expanded state summary columns
  - run-manifest and summary-table generation without plot dependencies
- Sample graph generation against the real score CSV in `/Users/jxs794/Downloads/job 3` succeeded for summaries and manifest output.
- `Plots/graph_report_manifest.json` was written and included the expected keys (`input_csv`, `output_dir`, `generated_at`, `top_n`, `hit_thresholds`, `summary_csvs`, `plot_files`, `pdf_report`, `warnings`).
- The previous pandas `groupby.apply` future-warning path is no longer used in the graph summary flow.
- The matplotlib `labels` deprecation warning was addressed with a compatibility wrapper that uses `tick_labels` when supported and falls back to `labels` for older versions.

## Known Issues
- This environment does not currently have matplotlib available, so the sample validation run could not generate actual PNG plots or the PDF report.
- The sample run therefore validated summary CSVs and manifest JSON generation, but not rendered plot files.
- Plot/PDF creation will still work only in environments where matplotlib is installed and importable.

## Manual Verification
1. Run:
   `python 7_Graph.py --csv auto --top-n 10 --write-state-summaries --write-ranking-plots --write-hit-rate-plots --write-state-heatmaps --write-pdf-report`
2. Confirm `Plots/summary_tables/` contains:
   - ligand summary
   - variant summary
   - ligand-base summary
   - tautomer summary
   - protomer summary
   - state summary
3. Open the ligand summary CSV and confirm it contains:
   - `best_affinity`
   - `median_affinity`
   - `top5_mean_affinity`
   - `hit_rate_below_-9`
   - `outlier_gap`
   - `consensus_rank`
4. Confirm receptor plot folders contain:
   - `*_ranking_bar.png`
   - `*_consensus_rank.png`
   - `*_hit_rate_heatmap.png`
   - `*_state_heatmap.png` when state metadata exists
5. Confirm `Plots/docking_report_<stem>.pdf` is created when matplotlib is available and `--write-pdf-report` is enabled.
6. Confirm `Plots/graph_report_manifest.json` lists the generated summary CSVs, plot files, PDF path, and warnings.

## Suggested Next Prompt
Add a web-server results dashboard that displays the PDF report, summary tables, plots, and direct links to PyMOL/SDF/PDBQT top-hit files.
