#!/usr/bin/env python3
"""
generate_docking_viewer.py  ·  VinaScope
══════════════════════════════════════════════════════════════════════
Generates a fully self-contained, serverless HTML docking viewer from
AutoDock Vina output PDBQT files.

Features
--------
  · 3-D interactive viewer  (3Dmol.js, WebGL)
  · All binding poses ranked by affinity with energy bars
  · Interaction panel — H-bond / Hydrophobic / π-stacking /
    Salt bridge / Van der Waals, grouped by residue
  · Interactions CSV exported automatically on first load
  · "Save Complex as PNG" button
  · VinaScope branding — unique viewer identity

Usage
-----
  python generate_docking_viewer.py                         # defaults
  python generate_docking_viewer.py receptor.pdbqt pose.pdbqt
  python generate_docking_viewer.py -r REC.pdbqt -l POSE.pdbqt -o out.html

Requirements
------------
  Python ≥ 3.6 · standard library only (base64, argparse, pathlib, sys)
  Internet on first open (3Dmol.js + Google Fonts via CDN)
"""

import argparse
import base64
import sys
from pathlib import Path


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="VinaScope — generate a self-contained HTML docking viewer.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("positional", nargs="*",
                   help="receptor.pdbqt  vina_pose.pdbqt  (positional shorthand)")
    p.add_argument("-r", "--receptor", default="receptor.pdbqt",
                   help="Receptor PDBQT  (default: receptor.pdbqt)")
    p.add_argument("-l", "--ligand",   default="vina_pose.pdbqt",
                   help="Vina poses PDBQT (default: vina_pose.pdbqt)")
    p.add_argument("-o", "--output",   default="VinaScope-DockingViewer.html",
                   help="Output HTML file (default: VinaScope-DockingViewer.html)")
    args = p.parse_args()
    if len(args.positional) >= 1:
        args.receptor = args.positional[0]
    if len(args.positional) >= 2:
        args.ligand   = args.positional[1]
    return args


def load_b64(path: str) -> str:
    """Read a file and return its Base-64–encoded content as an ASCII string."""
    p = Path(path)
    if not p.exists():
        print(f"ERROR: file not found — {path}", file=sys.stderr)
        sys.exit(1)
    return base64.b64encode(p.read_bytes()).decode("ascii")


# ─── HTML TEMPLATE ─────────────────────────────────────────────────────────────
# Uses %%RECEPTOR_B64%% / %%POSE_B64%% as injection tokens (no f-string escaping
# needed — every JS brace is written as-is).

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>VinaScope · Molecular Docking Viewer</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;500;600;700&family=Space+Mono:wght@400;700&display=swap" rel="stylesheet">
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.0/3Dmol-min.js"></script>
<style>
/* ── TOKENS ── */
:root{
  --bg:#f4f6f9;--panel:#ffffff;--surface:#edf0f5;--surface2:#e2e7ef;
  --border:#d0d7e2;--border-hi:#b6c0d0;
  --accent:#0d9488;--accent-dim:rgba(13,148,136,.08);
  --txt:#0f172a;--txt-muted:#475569;--txt-dim:#94a3b8;
  --g1:#059669;--g2:#0d9488;
  --amber:#d97706;--purple:#7c3aed;--red:#dc2626;
  --radius:5px;
  --font:'Space Grotesk',system-ui,sans-serif;
  --mono:'Space Mono','Fira Code',monospace;
}

*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html,body{width:100%;height:100%;background:var(--bg);color:var(--txt);font-family:var(--font);overflow:hidden}

/* ── LAYOUT ── */
#app{display:flex;width:100vw;height:100vh}
#viewer-wrap{flex:1;position:relative;overflow:hidden}
#viewer{width:100%;height:100%}

/* ── SIDEBAR ── */
#sidebar{
  width:288px;min-width:288px;height:100vh;
  background:var(--panel);
  border-right:1px solid var(--border);
  display:flex;flex-direction:column;
  overflow:hidden;z-index:10;
  box-shadow:4px 0 24px rgba(0,0,0,.08);
}

/* header */
#sidebar-header{
  padding:16px 16px 14px;
  border-bottom:1px solid var(--border);
  flex-shrink:0;
  background:linear-gradient(160deg,#dce3ed 0%,var(--panel) 100%);
}
.logo-row{display:flex;align-items:center;gap:10px;margin-bottom:9px}
.logo-icon{
  width:30px;height:30px;flex-shrink:0;
  background:var(--accent-dim);border:1px solid rgba(13,148,136,.3);
  border-radius:7px;display:flex;align-items:center;justify-content:center;
}
.logo-icon svg{width:18px;height:18px;color:var(--accent)}
.logo-text{}
.logo-name{font-size:15px;font-weight:700;letter-spacing:-.01em;
  background:linear-gradient(90deg,var(--g1),var(--accent));
  -webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;}
.logo-sub{font-size:8px;font-weight:500;letter-spacing:.15em;text-transform:uppercase;color:var(--txt-muted);margin-top:1px}
.ligand-chip{
  display:inline-flex;align-items:center;gap:6px;
  font-family:var(--mono);font-size:10px;color:var(--txt-muted);
  background:var(--surface);border:1px solid var(--border);
  padding:4px 8px;border-radius:var(--radius);
}
.ligand-chip svg{width:10px;height:10px;opacity:.5}

/* scrollable body */
#sidebar-body{
  flex:1;overflow-y:auto;padding:12px 0;
  scrollbar-width:thin;scrollbar-color:var(--border-hi) transparent;
}
#sidebar-body::-webkit-scrollbar{width:3px}
#sidebar-body::-webkit-scrollbar-thumb{background:var(--border-hi);border-radius:2px}

.section{padding:0 13px;margin-bottom:14px}
.sec-label{
  font-size:8px;font-weight:700;letter-spacing:.15em;text-transform:uppercase;
  color:var(--txt-dim);margin-bottom:8px;padding:0 1px;
  display:flex;align-items:center;gap:6px;
}
.sec-label::after{content:'';flex:1;height:1px;background:var(--border)}
.divider{height:1px;background:var(--border);margin:2px 13px 14px}

/* ── POSE LIST ── */
.pose-item{
  display:flex;align-items:center;gap:7px;
  padding:7px 8px;border-radius:var(--radius);
  cursor:pointer;margin-bottom:2px;
  border:1px solid transparent;
  transition:background .12s,border-color .12s;
}
.pose-item:hover{background:var(--surface)}
.pose-item.active{background:var(--surface);border-color:var(--border-hi)}
.pose-dot{width:7px;height:7px;border-radius:50%;flex-shrink:0}
.pose-info{flex:1;min-width:0}
.pose-rank{font-size:10px;font-weight:600;color:var(--txt-muted);margin-bottom:2px;display:flex;align-items:center;gap:4px}
.pose-bar-wrap{height:2px;background:var(--border);border-radius:1px;overflow:hidden}
.pose-bar{height:100%;border-radius:1px}
.pose-rmsd{font-family:var(--mono);font-size:7.5px;color:var(--txt-dim);margin-top:2px}
.pose-score{font-family:var(--mono);font-size:11px;font-weight:700;flex-shrink:0;text-align:right}
.pose-unit{font-family:var(--mono);font-size:7px;color:var(--txt-dim);text-align:right;margin-top:1px}
.badge-best{
  font-size:6px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;
  color:var(--g1);background:rgba(5,150,105,.1);
  border:1px solid rgba(5,150,105,.3);padding:1px 4px;border-radius:3px;
}
.pose-eye{
  width:17px;height:17px;display:flex;align-items:center;justify-content:center;
  border-radius:3px;flex-shrink:0;opacity:.3;transition:opacity .12s;cursor:pointer;
}
.pose-eye:hover{opacity:1}
.pose-eye.hidden{opacity:.15}
.pose-eye svg{width:12px;height:12px}

/* ── CONTROLS ── */
.ctrl-group{margin-bottom:10px}
.ctrl-label{font-size:9px;font-weight:500;color:var(--txt-muted);margin-bottom:4px}
.btn-row{display:flex;gap:3px;flex-wrap:wrap}
.btn{
  font-family:var(--font);font-size:9px;font-weight:600;
  padding:4px 8px;border-radius:4px;
  border:1px solid var(--border);background:transparent;
  color:var(--txt-muted);cursor:pointer;
  transition:all .12s;letter-spacing:.02em;
}
.btn:hover{border-color:var(--border-hi);color:var(--txt);background:var(--surface)}
.btn.active{background:var(--accent-dim);border-color:var(--accent);color:var(--accent)}
.btn-action{
  width:100%;padding:7px 10px;font-size:9.5px;font-weight:600;
  display:flex;align-items:center;justify-content:center;gap:6px;
  border-radius:var(--radius);border:1px solid var(--border);
  background:transparent;color:var(--txt-muted);cursor:pointer;
  transition:all .15s;font-family:var(--font);letter-spacing:.02em;
}
.btn-action:hover{border-color:var(--border-hi);color:var(--txt);background:var(--surface2)}
.btn-action.primary{border-color:rgba(45,212,191,.35);color:var(--accent);background:var(--accent-dim)}
.btn-action.primary:hover{border-color:var(--accent);background:rgba(45,212,191,.18)}
.btn-action svg{width:11px;height:11px;flex-shrink:0}
.toggle-row{display:flex;align-items:center;gap:7px;padding:3px 0;cursor:pointer}
.toggle-row input{accent-color:var(--accent);width:11px;height:11px;cursor:pointer}
.toggle-row span{font-size:10px;color:var(--txt-muted)}

/* ── INTERACTION TABLE ── */
.int-hdr{
  display:flex;align-items:center;padding:3px 1px 5px;
  border-bottom:1px solid var(--border);margin-bottom:3px;
}
.int-hdr span{font-size:7.5px;font-weight:700;letter-spacing:.1em;text-transform:uppercase;color:var(--txt-dim)}
.ihc-res{flex:0 0 92px}
.ihc-typ{flex:1}
.ihc-dst{flex:0 0 28px;text-align:right}

.int-row{display:flex;align-items:flex-start;padding:4px 1px;border-bottom:1px solid rgba(22,34,56,.7);gap:3px}
.int-row:last-child{border-bottom:none}
.int-res{font-family:var(--mono);font-size:9px;color:var(--txt);flex:0 0 92px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;padding-top:1px}
.int-types{flex:1;display:flex;flex-wrap:wrap;gap:2px}
.int-dist{font-family:var(--mono);font-size:8px;color:var(--txt-muted);flex:0 0 28px;text-align:right;padding-top:2px}
.int-badge{font-size:6.5px;font-weight:700;padding:1px 4px;border-radius:3px;text-transform:uppercase;letter-spacing:.04em;white-space:nowrap}
.ib-hbond      {background:rgba(13,148,136,.12);color:#0d9488;border:1px solid rgba(13,148,136,.25)}
.ib-hydrophobic{background:rgba(217,119,6,.12);color:#d97706;border:1px solid rgba(217,119,6,.25)}
.ib-pistacking {background:rgba(124,58,237,.12);color:#7c3aed;border:1px solid rgba(124,58,237,.25)}
.ib-saltbridge {background:rgba(220,38,38,.12);color:#dc2626;border:1px solid rgba(220,38,38,.25)}
.ib-vanderwaals{background:rgba(71,85,105,.12);color:#475569;border:1px solid rgba(71,85,105,.25)}
.int-empty{font-size:9.5px;color:var(--txt-dim);padding:8px 1px;font-style:italic}
.int-summary{font-family:var(--mono);font-size:8px;color:var(--txt-dim);margin-top:6px;padding:0 1px}

/* ── FOOTER ── */
#sidebar-footer{
  padding:10px 13px;border-top:1px solid var(--border);flex-shrink:0;
  background:var(--panel);
}
.stats-row{display:flex;justify-content:space-between}
.stat-item{text-align:center}
.stat-val{font-family:var(--mono);font-size:13px;font-weight:700;color:var(--accent)}
.stat-lbl{font-size:7px;text-transform:uppercase;letter-spacing:.1em;color:var(--txt-dim);margin-top:1px}

/* ── LOADING ── */
#loading{
  position:fixed;inset:0;background:var(--bg);
  display:flex;flex-direction:column;align-items:center;justify-content:center;
  z-index:100;transition:opacity .5s;
}
#loading.hidden{opacity:0;pointer-events:none}
.loading-ring{
  width:40px;height:40px;
  border:2px solid var(--border);border-top-color:var(--accent);
  border-radius:50%;animation:spin .75s linear infinite;margin-bottom:14px;
}
@keyframes spin{to{transform:rotate(360deg)}}
.loading-label{font-size:11px;color:var(--txt-muted);letter-spacing:.06em}
.loading-name{
  font-size:20px;font-weight:700;
  background:linear-gradient(90deg,var(--g1),var(--accent));
  -webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;
  margin-bottom:20px;letter-spacing:-.01em;
}

/* ── HUD ── */
#hud{position:absolute;top:14px;right:14px;display:flex;flex-direction:column;align-items:flex-end;gap:7px;pointer-events:none;z-index:5}
.hud-card{
  background:rgba(255,255,255,.9);border:1px solid var(--border);
  border-radius:var(--radius);padding:8px 12px;backdrop-filter:blur(12px);
  pointer-events:auto;min-width:120px;
}
.hud-lbl{font-size:7.5px;text-transform:uppercase;letter-spacing:.13em;color:var(--txt-muted);margin-bottom:2px}
.hud-val{font-family:var(--mono);font-size:21px;font-weight:700;line-height:1}
.hud-unit{font-size:8.5px;color:var(--txt-muted);margin-left:2px}
.hud-pose{font-family:var(--mono);font-size:12px;font-weight:700;color:var(--txt)}

/* ── TOAST ── */
#toast{
  position:fixed;bottom:18px;right:18px;
  background:rgba(255,255,255,.95);border:1px solid rgba(13,148,136,.4);
  color:var(--accent);padding:9px 14px;border-radius:var(--radius);
  font-size:10px;letter-spacing:.04em;
  z-index:200;opacity:0;transform:translateY(6px);
  transition:opacity .25s,transform .25s;pointer-events:none;
  display:flex;align-items:center;gap:8px;
  font-family:var(--font);font-weight:500;
}
#toast.show{opacity:1;transform:translateY(0)}
#toast svg{width:12px;height:12px;flex-shrink:0;color:var(--g1)}
</style>
</head>
<body>

<!-- Loading screen -->
<div id="loading">
  <div class="loading-name">VinaScope</div>
  <div class="loading-ring"></div>
  <div class="loading-label">Loading molecular data…</div>
</div>

<!-- Toast notification -->
<div id="toast">
  <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="2">
    <polyline points="2,9 6,13 14,4"/>
  </svg>
  <span id="toast-msg"></span>
</div>

<div id="app">

  <!-- ══ SIDEBAR ══ -->
  <aside id="sidebar">
    <div id="sidebar-header">
      <div class="logo-row">
        <div class="logo-icon">
          <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.4">
            <circle cx="12" cy="12" r="2.4" fill="currentColor" stroke="none"/>
            <ellipse cx="12" cy="12" rx="9.5" ry="3.6"/>
            <ellipse cx="12" cy="12" rx="9.5" ry="3.6" transform="rotate(60 12 12)"/>
            <ellipse cx="12" cy="12" rx="9.5" ry="3.6" transform="rotate(120 12 12)"/>
          </svg>
        </div>
        <div class="logo-text">
          <div class="logo-name">VinaScope</div>
          <div class="logo-sub">AutoDock Vina · Docking Viewer</div>
        </div>
      </div>
      <div class="ligand-chip" id="ligand-label">
        <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.5">
          <circle cx="4" cy="8" r="2"/><circle cx="12" cy="4" r="2"/><circle cx="12" cy="12" r="2"/>
          <line x1="6" y1="8" x2="10" y2="5"/><line x1="6" y1="8" x2="10" y2="11"/>
        </svg>
        <span>Ligand: —</span>
      </div>
    </div>

    <div id="sidebar-body">

      <!-- POSES -->
      <div class="section">
        <div class="sec-label">Binding Poses</div>
        <div id="pose-list"></div>
      </div>

      <div class="divider"></div>

      <!-- RECEPTOR STYLE -->
      <div class="section">
        <div class="sec-label">Receptor</div>
        <div class="ctrl-group">
          <div class="ctrl-label">Style</div>
          <div class="btn-row" id="rec-style-btns">
            <button class="btn active" data-v="cartoon">Cartoon</button>
            <button class="btn" data-v="surface">Surface</button>
            <button class="btn" data-v="line">Lines</button>
            <button class="btn" data-v="stick">Sticks</button>
          </div>
        </div>
        <div class="ctrl-group">
          <div class="ctrl-label">Colour</div>
          <div class="btn-row" id="rec-color-btns">
            <button class="btn active" data-v="spectrum">Spectrum</button>
            <button class="btn" data-v="chain">Chain</button>
            <button class="btn" data-v="ss">Struct</button>
          </div>
        </div>
        <label class="toggle-row">
          <input type="checkbox" id="surface-toggle">
          <span>VDW surface overlay</span>
        </label>
      </div>

      <div class="divider"></div>

      <!-- LIGAND STYLE -->
      <div class="section">
        <div class="sec-label">Ligand</div>
        <div class="ctrl-group">
          <div class="ctrl-label">Style</div>
          <div class="btn-row" id="lig-style-btns">
            <button class="btn active" data-v="stick">Sticks</button>
            <button class="btn" data-v="sphere">Spheres</button>
            <button class="btn" data-v="ball">Ball&amp;Stick</button>
          </div>
        </div>
        <div class="ctrl-group">
          <div class="ctrl-label">Colour by</div>
          <div class="btn-row" id="lig-color-btns">
            <button class="btn active" data-v="pose">Pose rank</button>
            <button class="btn" data-v="element">Element</button>
          </div>
        </div>
        <label class="toggle-row">
          <input type="checkbox" id="all-poses-toggle" checked>
          <span>Show all poses</span>
        </label>
      </div>

      <div class="divider"></div>

      <!-- INTERACTIONS -->
      <div class="section">
        <div class="sec-label">Interactions — Pose <span id="int-pose-num">1</span></div>
        <div class="int-hdr">
          <span class="ihc-res">Residue</span>
          <span class="ihc-typ">Type</span>
          <span class="ihc-dst">Å</span>
        </div>
        <div id="interaction-list"><div class="int-empty">Computing…</div></div>
        <div id="int-summary" class="int-summary"></div>
      </div>

      <div class="divider"></div>

      <!-- EXPORT & CAMERA -->
      <div class="section">
        <div class="sec-label">Export &amp; Camera</div>
        <div style="display:flex;flex-direction:column;gap:5px">
          <button class="btn-action primary" id="btn-save-png">
            <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6">
              <rect x="1" y="1" width="14" height="14" rx="2"/>
              <circle cx="5.5" cy="5.5" r="1.4"/>
              <polyline points="1,11 5,7 8,10 11,6.5 15,11"/>
            </svg>
            Save Complex as PNG
          </button>
          <button class="btn-action" id="btn-export-csv">
            <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.6">
              <path d="M8 2v8M5 7l3 3 3-3"/><path d="M2 12v1a1 1 0 001 1h10a1 1 0 001-1v-1"/>
            </svg>
            Re-export Interactions CSV
          </button>
          <div style="display:flex;gap:4px;margin-top:1px">
            <button class="btn-action" id="btn-reset-view" style="flex:1">Reset View</button>
            <button class="btn-action" id="btn-focus-lig" style="flex:1">Focus Ligand</button>
          </div>
        </div>
      </div>

    </div><!-- /sidebar-body -->

    <div id="sidebar-footer">
      <div class="stats-row">
        <div class="stat-item"><div class="stat-val" id="stat-poses">—</div><div class="stat-lbl">Poses</div></div>
        <div class="stat-item"><div class="stat-val" id="stat-best">—</div><div class="stat-lbl">Best kcal/mol</div></div>
        <div class="stat-item"><div class="stat-val" id="stat-atoms">—</div><div class="stat-lbl">Rec. Atoms</div></div>
      </div>
    </div>
  </aside>

  <!-- ══ VIEWER ══ -->
  <div id="viewer-wrap">
    <div id="viewer"></div>
    <div id="hud">
      <div class="hud-card">
        <div class="hud-lbl">Binding Affinity</div>
        <div><span class="hud-val" id="hud-score">—</span><span class="hud-unit">kcal/mol</span></div>
      </div>
      <div class="hud-card">
        <div class="hud-lbl">Active Pose</div>
        <div class="hud-pose" id="hud-pose">—</div>
      </div>
    </div>
  </div>

</div><!-- /#app -->

<script>
// ═══════════════════════════════════════════════════════════════════
//  EMBEDDED MOLECULAR DATA  (Base-64 PDBQT)
// ═══════════════════════════════════════════════════════════════════
const RECEPTOR_B64 = "%%RECEPTOR_B64%%";
const POSE_B64     = "%%POSE_B64%%";

function b64ToStr(b64) {
  const bytes = atob(b64), arr = new Uint8Array(bytes.length);
  for (let i = 0; i < bytes.length; i++) arr[i] = bytes.charCodeAt(i);
  return new TextDecoder("utf-8").decode(arr);
}

// ═══════════════════════════════════════════════════════════════════
//  POSE COLOUR PALETTE  (rank 0=best green → rank N=red)
// ═══════════════════════════════════════════════════════════════════
const PALETTE = [
  "#00ffa3","#00e8b8","#2dd4bf","#38bdf8",
  "#818cf8","#a78bfa","#e879f9","#fb923c","#f87171"
];
const poseHex = i => PALETTE[i % PALETTE.length];

// ═══════════════════════════════════════════════════════════════════
//  PDBQT PARSERS
// ═══════════════════════════════════════════════════════════════════
function parsePoses(pdbqt) {
  const models = []; let cur = null;
  for (const raw of pdbqt.split("\n")) {
    const line = raw.trimEnd();
    if (line.startsWith("MODEL")) {
      cur = { score:null, rmsd_lb:null, rmsd_ub:null, raw:[] };
      continue;
    }
    if (line.startsWith("ENDMDL") && cur) {
      cur.raw.push(line); models.push(cur); cur = null; continue;
    }
    if (cur) {
      if (line.startsWith("REMARK VINA RESULT:")) {
        const p = line.replace("REMARK VINA RESULT:","").trim().split(/\s+/);
        cur.score=parseFloat(p[0]); cur.rmsd_lb=parseFloat(p[1]); cur.rmsd_ub=parseFloat(p[2]);
      }
      cur.raw.push(line);
    }
  }
  return models;
}

const COORD_RX = /^\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;

function parseAtoms(pdbStr) {
  /* Returns full atom records including PDBQT charge + adtype fields. */
  const atoms = [];
  for (const raw of pdbStr.split("\n")) {
    const line = raw.trimEnd();
    if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) continue;
    if (line.length < 54) continue;
    const name    = line.substring(12,16).trim();
    const resname = line.substring(17,20).trim();
    const chain   = line.substring(21,22).trim() || "A";
    const resi    = parseInt(line.substring(22,26)) || 0;
    const m = COORD_RX.exec(line.substring(30));
    if (!m) continue;
    const x=parseFloat(m[1]), y=parseFloat(m[2]), z=parseFloat(m[3]);
    // PDBQT cols 54+: occupancy  bfactor  partialCharge  adType
    const parts  = line.substring(54).trim().split(/\s+/);
    const charge = parts.length>=3 ? parseFloat(parts[2])||0 : 0;
    const adtype = parts.length>=4 ? parts[3] : "";
    atoms.push({ name, resname, chain, resi, x, y, z, charge, adtype });
  }
  return atoms;
}

function poseToPDB(pose) {
  return pose.raw
    .filter(l=>l.startsWith("ATOM")||l.startsWith("HETATM")||l.startsWith("TER")||l.startsWith("ENDMDL"))
    .join("\n");
}

function receptorToPDB(pdbqt) {
  return pdbqt.split("\n")
    .filter(l=>l.startsWith("ATOM")||l.startsWith("HETATM")||l.startsWith("TER")||l.startsWith("END"))
    .join("\n");
}

function countAtoms(pdb) {
  return pdb.split("\n").filter(l=>l.startsWith("ATOM")||l.startsWith("HETATM")).length;
}

function extractLigandName(pdbqt) {
  const m = pdbqt.match(/REMARK\s+Name\s*=\s*(\S+)/);
  return m ? m[1] : "Ligand";
}

// ═══════════════════════════════════════════════════════════════════
//  INTERACTION DETECTION
//  Rules validated against AutoDock PDBQT atom types:
//    OA/NA/SA = H-bond acceptors
//    HD       = H-bond donor hydrogen
//    N        = backbone/charged nitrogen (donor)
//    A        = aromatic carbon
//    C        = aliphatic carbon
// ═══════════════════════════════════════════════════════════════════
const ACC  = new Set(["OA","NA","SA"]);
const HYDR = new Set(["C","A"]);
const AROM = new Set(["PHE","TYR","TRP","HIS"]);
const NEG  = { ASP:new Set(["OD1","OD2"]), GLU:new Set(["OE1","OE2"]) };
const POS  = { LYS:new Set(["NZ"]), ARG:new Set(["NH1","NH2","NE"]), HIS:new Set(["ND1","NE2"]) };

const TYPE_CSS = {
  "H-bond"      : "ib-hbond",
  "Hydrophobic" : "ib-hydrophobic",
  "π-stacking"  : "ib-pistacking",
  "Salt bridge" : "ib-saltbridge",
  "Van der Waals":"ib-vanderwaals"
};

function dist3(a,b) {
  const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
  return Math.sqrt(dx*dx+dy*dy+dz*dz);
}

function classifyPair(la, ra, d) {
  const lt=la.adtype, rt=ra.adtype;
  /* H-bond — heavy-atom donor/acceptor pairs (≤3.5 Å) */
  if (d<=3.5) {
    const laAcc=ACC.has(lt), raAcc=ACC.has(rt);
    const laDon=(lt==="N"||laAcc), raDon=(rt==="N"||raAcc);
    if ((laAcc&&raDon)||(raAcc&&laDon)) return "H-bond";
  }
  /* H-bond — HD hydrogen to acceptor (≤2.5 Å) */
  if (d<=2.5 && lt==="HD" && ACC.has(rt)) return "H-bond";
  if (d<=2.5 && rt==="HD" && ACC.has(lt)) return "H-bond";
  /* Salt bridge — charged side-chains (≤4.0 Å) */
  if (d<=4.0) {
    if (NEG[ra.resname]?.has(ra.name) && la.charge >  0.3) return "Salt bridge";
    if (POS[ra.resname]?.has(ra.name) && la.charge < -0.3) return "Salt bridge";
  }
  /* π-stacking — both aromatic-type atoms, aromatic receptor residue (≤5.5 Å) */
  if (d<=5.5 && lt==="A" && rt==="A" && AROM.has(ra.resname)) return "π-stacking";
  /* Hydrophobic — non-polar carbon/aromatic pairs (≤4.5 Å) */
  if (d<=4.5 && HYDR.has(lt) && HYDR.has(rt)) return "Hydrophobic";
  /* Van der Waals catch-all — exclude bare H contacts (≤3.8 Å) */
  if (d<=3.8 && lt!=="HD" && rt!=="HD") return "Van der Waals";
  return null;
}

function computeInteractions(ligAtoms, recAtoms) {
  /* Pre-filter: only receptor atoms inside ligand bounding sphere + 6 Å */
  const cx=ligAtoms.reduce((s,a)=>s+a.x,0)/ligAtoms.length;
  const cy=ligAtoms.reduce((s,a)=>s+a.y,0)/ligAtoms.length;
  const cz=ligAtoms.reduce((s,a)=>s+a.z,0)/ligAtoms.length;
  const rmax = Math.max(...ligAtoms.map(a=>{
    const dx=a.x-cx,dy=a.y-cy,dz=a.z-cz; return Math.sqrt(dx*dx+dy*dy+dz*dz);
  }));
  const nearby = recAtoms.filter(ra=>{
    const dx=ra.x-cx,dy=ra.y-cy,dz=ra.z-cz;
    return Math.sqrt(dx*dx+dy*dy+dz*dz)<=rmax+6;
  });

  const rows = [];
  for (const la of ligAtoms) {
    for (const ra of nearby) {
      const d = dist3(la,ra);
      if (d>5.5) continue;
      const type = classifyPair(la,ra,d);
      if (!type) continue;
      rows.push({
        ligand_atom    : la.name,
        ligand_adtype  : la.adtype,
        ligand_charge  : la.charge.toFixed(3),
        receptor_chain : ra.chain,
        receptor_resname:ra.resname,
        receptor_resi  : ra.resi,
        receptor_atom  : ra.name,
        receptor_adtype: ra.adtype,
        distance       : d,
        type
      });
    }
  }
  rows.sort((a,b)=>a.distance-b.distance);
  return rows;
}

function groupByResidue(rows) {
  const map = new Map();
  for (const r of rows) {
    const key = `${r.receptor_resname} ${r.receptor_resi}:${r.receptor_chain}`;
    if (!map.has(key)) map.set(key, { key, types:new Set(), minDist:Infinity });
    const g = map.get(key);
    g.types.add(r.type);
    g.minDist = Math.min(g.minDist, r.distance);
  }
  return [...map.values()].sort((a,b)=>a.minDist-b.minDist);
}

// ═══════════════════════════════════════════════════════════════════
//  RENDER INTERACTION TABLE
// ═══════════════════════════════════════════════════════════════════
function renderInteractions(rows, poseIdx) {
  document.getElementById("int-pose-num").textContent = poseIdx+1;
  const list = document.getElementById("interaction-list");
  const sumEl = document.getElementById("int-summary");
  if (!rows || !rows.length) {
    list.innerHTML = '<div class="int-empty">No interactions detected within 5.5 Å</div>';
    sumEl.textContent = ""; return;
  }
  const groups = groupByResidue(rows);
  list.innerHTML = groups.map(g => {
    const badges = [...g.types].map(t =>
      `<span class="int-badge ${TYPE_CSS[t]||"ib-vanderwaals"}">${t}</span>`
    ).join("");
    return `<div class="int-row">
      <div class="int-res" title="${g.key}">${g.key}</div>
      <div class="int-types">${badges}</div>
      <div class="int-dist">${g.minDist.toFixed(1)}</div>
    </div>`;
  }).join("");
  sumEl.textContent = `${rows.length} contacts · ${groups.length} residues`;
}

// ═══════════════════════════════════════════════════════════════════
//  CSV BUILDER
// ═══════════════════════════════════════════════════════════════════
function buildCSV(allInts, poses) {
  const header = [
    "pose","score_kcal_mol",
    "ligand_atom","ligand_adtype","ligand_charge",
    "receptor_chain","receptor_resname","receptor_resnum",
    "receptor_atom","receptor_adtype",
    "distance_A","interaction_type"
  ].join(",");
  const lines = [];
  allInts.forEach((rows,i) => rows.forEach(r => lines.push([
    i+1, poses[i].score,
    r.ligand_atom, r.ligand_adtype, r.ligand_charge,
    r.receptor_chain, r.receptor_resname, r.receptor_resi,
    r.receptor_atom, r.receptor_adtype,
    r.distance.toFixed(2), r.type
  ].join(","))));
  return header+"\n"+lines.join("\n");
}

// ═══════════════════════════════════════════════════════════════════
//  HELPERS: download blob · toast
// ═══════════════════════════════════════════════════════════════════
function downloadBlob(data, filename, mime) {
  const blob = new Blob([data],{type:mime});
  const url  = URL.createObjectURL(blob);
  const a    = Object.assign(document.createElement("a"),{href:url,download:filename});
  document.body.appendChild(a); a.click();
  document.body.removeChild(a); URL.revokeObjectURL(url);
}

let _toastT;
function toast(msg) {
  const el = document.getElementById("toast");
  document.getElementById("toast-msg").textContent = msg;
  el.classList.add("show"); clearTimeout(_toastT);
  _toastT = setTimeout(()=>el.classList.remove("show"), 3400);
}

// ═══════════════════════════════════════════════════════════════════
//  MAIN
// ═══════════════════════════════════════════════════════════════════
(function main() {

  // ── Decode & parse ──────────────────────────────────────────────
  const receptorData = b64ToStr(RECEPTOR_B64);
  const poseData     = b64ToStr(POSE_B64);
  const poses        = parsePoses(poseData);
  const ligandName   = extractLigandName(poseData);
  const receptorPDB  = receptorToPDB(receptorData);

  // Header stats
  document.getElementById("ligand-label").querySelector("span").textContent = "Ligand: "+ligandName;
  document.getElementById("stat-poses").textContent = poses.length;
  document.getElementById("stat-best").textContent  = poses[0]?.score?.toFixed(1) ?? "—";
  document.getElementById("stat-atoms").textContent = countAtoms(receptorPDB);

  // ── 3Dmol viewer ────────────────────────────────────────────────
  const viewer = $3Dmol.createViewer(document.getElementById("viewer"),{
    backgroundColor:"#f4f6f9", antialias:true
  });
  const recModel = viewer.addModel(receptorPDB,"pdb");

  let recStyle="cartoon", recColor="spectrum", recSurfObj=null, surfObj=null, surfOn=false;

  function applyRecStyle() {
    viewer.setStyle({model:recModel.getID()},{});
    const cs = recColor==="chain"?"chain": recColor==="ss"?"ssJmol": undefined;
    if (recSurfObj) { viewer.removeSurface(recSurfObj); recSurfObj = null; }
    if (recStyle==="cartoon")
      viewer.setStyle({model:recModel.getID()},{cartoon:{color:recColor==="spectrum"?"spectrum":undefined,colorscheme:cs}});
    else if (recStyle==="surface") {
      viewer.setStyle({model:recModel.getID()},{line:{colorscheme:cs??"element",linewidth:0.3}});
      recSurfObj = viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.88,colorscheme:$3Dmol.elementColors.rasmol},{model:recModel.getID()});
    }
    else if (recStyle==="line")
      viewer.setStyle({model:recModel.getID()},{line:{colorscheme:cs??"element",linewidth:1.5}});
    else
      viewer.setStyle({model:recModel.getID()},{stick:{colorscheme:cs??"element",radius:0.12}});
    if (surfObj) { viewer.removeSurface(surfObj); surfObj=null; }
    if (surfOn)
      surfObj = viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.35,color:"#0a2040"},{model:recModel.getID()});
  }

  const poseMdls   = poses.map(p=>viewer.addModel(poseToPDB(p),"pdb"));
  const poseVis    = poses.map(()=>true);
  let curPose=0, ligStyle="stick", ligColor="pose", showAll=true;

  function applyLigStyles() {
    poseMdls.forEach((mdl,i)=>{
      const vis=poseVis[i]&&(showAll||i===curPose);
      const active=i===curPose, hex=poseHex(i), op=vis?(active?1:.55):0;
      if (!vis) { viewer.setStyle({model:mdl.getID()},{}); return; }
      const cs=ligColor==="element"?{colorscheme:"element"}:{color:hex};
      const r=active?.22:.13;
      if (ligStyle==="sphere")
        viewer.setStyle({model:mdl.getID()},{sphere:{...cs,radius:active?.5:.35,opacity:op}});
      else if (ligStyle==="ball")
        viewer.setStyle({model:mdl.getID()},{stick:{...cs,radius:r,opacity:op},sphere:{...cs,radius:.3,opacity:op}});
      else
        viewer.setStyle({model:mdl.getID()},{stick:{...cs,radius:r,opacity:op}});
    });
  }

  // ── Pre-compute interactions for ALL poses ───────────────────────
  const recAtoms      = parseAtoms(receptorPDB);
  const allInteractions = poses.map((_,i)=>
    computeInteractions(parseAtoms(poseToPDB(poses[i])), recAtoms)
  );

  // ── Pose list ────────────────────────────────────────────────────
  const bestScore  = Math.min(...poses.map(p=>p.score));
  const worstScore = Math.max(...poses.map(p=>p.score));

  function buildPoseList() {
    const c = document.getElementById("pose-list"); c.innerHTML="";
    poses.forEach((pose,i)=>{
      const frac = poses.length>1?(pose.score-worstScore)/(bestScore-worstScore):1;
      const barW = Math.round(16+frac*78), hex=poseHex(i);
      const item = document.createElement("div");
      item.className = "pose-item"+(i===0?" active":"");
      item.innerHTML = `
        <div class="pose-dot" style="background:${hex};box-shadow:0 0 5px ${hex}55"></div>
        <div class="pose-info">
          <div class="pose-rank">
            Pose ${i+1}
            ${i===0?'<span class="badge-best">BEST</span>':''}
          </div>
          <div class="pose-bar-wrap"><div class="pose-bar" style="width:${barW}%;background:${hex}"></div></div>
          <div class="pose-rmsd">lb ${pose.rmsd_lb.toFixed(2)} / ub ${pose.rmsd_ub.toFixed(2)} Å</div>
        </div>
        <div>
          <div class="pose-score" style="color:${hex}">${pose.score.toFixed(1)}</div>
          <div class="pose-unit">kcal/mol</div>
        </div>
        <div class="pose-eye" title="Toggle visibility">
          <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.4">
            <ellipse cx="8" cy="8" rx="5" ry="3.5"/>
            <circle cx="8" cy="8" r="1.5" fill="currentColor" stroke="none"/>
          </svg>
        </div>`;
      item.addEventListener("click", e=>{
        if (e.target.closest(".pose-eye")) return;
        selectPose(i);
      });
      item.querySelector(".pose-eye").addEventListener("click", e=>{
        e.stopPropagation();
        poseVis[i]=!poseVis[i];
        e.currentTarget.classList.toggle("hidden",!poseVis[i]);
        applyLigStyles(); viewer.render();
      });
      c.appendChild(item);
    });
  }

  function selectPose(i) {
    curPose=i;
    document.querySelectorAll(".pose-item").forEach((el,j)=>el.classList.toggle("active",j===i));
    const hex=poseHex(i);
    document.getElementById("hud-score").textContent = poses[i].score.toFixed(1);
    document.getElementById("hud-score").style.color = hex;
    document.getElementById("hud-pose").textContent  = `Pose ${i+1} of ${poses.length}`;
    applyLigStyles();
    renderInteractions(allInteractions[i], i);
    viewer.render();
  }

  // ── Button wiring ────────────────────────────────────────────────
  function wireGroup(sel, key, cb) {
    document.querySelectorAll(sel).forEach(btn=>btn.addEventListener("click",()=>{
      document.querySelectorAll(sel).forEach(b=>b.classList.remove("active"));
      btn.classList.add("active"); cb(btn.dataset[key]);
    }));
  }
  wireGroup("#rec-style-btns .btn","v", v=>{ recStyle=v; applyRecStyle(); viewer.render(); });
  wireGroup("#rec-color-btns .btn","v", v=>{ recColor=v; applyRecStyle(); viewer.render(); });
  wireGroup("#lig-style-btns .btn","v", v=>{ ligStyle=v; applyLigStyles(); viewer.render(); });
  wireGroup("#lig-color-btns .btn","v", v=>{ ligColor=v; applyLigStyles(); viewer.render(); });

  document.getElementById("surface-toggle").addEventListener("change",e=>{
    surfOn=e.target.checked; applyRecStyle(); viewer.render();
  });
  document.getElementById("all-poses-toggle").addEventListener("change",e=>{
    showAll=e.target.checked; applyLigStyles(); viewer.render();
  });
  document.getElementById("btn-reset-view").addEventListener("click",()=>{
    viewer.zoomTo(); viewer.render();
  });
  document.getElementById("btn-focus-lig").addEventListener("click",()=>{
    viewer.zoomTo({model:poseMdls[curPose].getID()}); viewer.render();
  });

  // ── Save PNG ─────────────────────────────────────────────────────
  document.getElementById("btn-save-png").addEventListener("click",()=>{
    viewer.render();
    let uri;
    try { uri = viewer.pngURI(); }
    catch(e) {
      const canvas = document.querySelector("#viewer canvas");
      uri = canvas ? canvas.toDataURL("image/png") : null;
    }
    if (!uri) { toast("PNG export not available in this browser"); return; }
    const comma = uri.indexOf(',');
    const mime = uri.substring(5, comma).split(';')[0];
    const bstr = atob(uri.substring(comma + 1));
    const buf = new Uint8Array(bstr.length);
    for (let i = 0; i < bstr.length; i++) buf[i] = bstr.charCodeAt(i);
    downloadBlob(new Blob([buf], {type: mime}), `vinascope_pose${curPose+1}.png`, mime);
    toast(`Pose ${curPose+1} saved as PNG`);
  });

  // ── CSV — manual re-export button ────────────────────────────────
  document.getElementById("btn-export-csv").addEventListener("click",()=>{
    downloadBlob(buildCSV(allInteractions,poses),"vina_interactions.csv","text/csv;charset=utf-8;");
    toast("Interactions CSV exported");
  });

  // ── First render ─────────────────────────────────────────────────
  buildPoseList();
  applyRecStyle();
  applyLigStyles();
  selectPose(0);
  viewer.zoomTo();
  viewer.render();

  // Dismiss loading overlay, then auto-export CSV
  setTimeout(()=>{
    const el = document.getElementById("loading");
    el.classList.add("hidden");
    setTimeout(()=>{
      el.remove();
      downloadBlob(buildCSV(allInteractions,poses),"vina_interactions.csv","text/csv;charset=utf-8;");
      toast("Interactions CSV saved automatically");
    }, 500);
  }, 1000);

})();
</script>
</body>
</html>"""


def build_html(receptor_b64: str, pose_b64: str) -> str:
    """Inject Base-64 data into the HTML template via token replacement."""
    return (HTML_TEMPLATE
            .replace("%%RECEPTOR_B64%%", receptor_b64)
            .replace("%%POSE_B64%%",     pose_b64))


# ─── ENTRY POINT ──────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    print(f"  Receptor : {args.receptor}")
    print(f"  Poses    : {args.ligand}")

    receptor_b64 = load_b64(args.receptor)
    pose_b64     = load_b64(args.ligand)

    print(f"  Receptor : {len(receptor_b64):>8,} chars (Base-64)")
    print(f"  Poses    : {len(pose_b64):>8,} chars (Base-64)")

    html = build_html(receptor_b64, pose_b64)

    out = Path(args.output)
    out.write_text(html, encoding="utf-8")
    print(f"  Output   : {out}  ({out.stat().st_size:,} bytes)")
    print()
    print("  Open the HTML in any modern browser.")
    print("  Interaction CSV exports automatically on first load.")

if __name__ == "__main__":
    print("\n  VinaScope — Molecular Docking Viewer Generator")
    print("  " + "─"*46)
    main()
    print("  " + "─"*46 + "\n")
