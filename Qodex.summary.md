# Qodex.summary

## Task
Optimize AutoDock-Vina PrepServer for AI-search, OpenAI crawler discovery, and machine-readable public documentation.

## Original Goal
The user wants ChatGPT/OpenAI bots and assistant-style systems to find and understand the site when users ask about docking, AutoDock Vina, docking preparation, headless APIs, and related workflows. They also asked for robots.txt and other discovery improvements.

## Assumptions
- Public educational and documentation pages should be discoverable and indexable.
- Generated workspaces, downloads, uploads, auth routes, and mutation-oriented API surfaces should not be promoted as crawl targets.
- `OAI-SearchBot`, `GPTBot`, and `ChatGPT-User` should be allowed for this site unless the owner later chooses a different crawl policy.
- `llms.txt` and `llms-full.txt` are helpful machine-readable conventions, not guaranteed ranking or citation mechanisms.
- The Build page remains internally linked and documented, but it is omitted from the sitemap because its browser JavaScript auto-creates a workspace on page load.
- Scientific claims should remain preparation-focused and should not imply docking accuracy, clinical validation, or replacement of expert review.

## Files Inspected
- `app.py` — reviewed public routes, generated/download/API/auth surfaces, and route placement for root discovery files.
- `templates/base.html` — reviewed canonical URL, Open Graph metadata, JSON-LD, and head links.
- `templates/api.html` — reviewed headless API copy, examples, endpoint summary, and shell snippets.
- `templates/home.html` — reviewed top-level AI/search landing-page wording.
- `templates/workflow.html` — reviewed workflow explanation and API links.
- `templates/documentation.html` — reviewed documentation hub structure.
- `templates/troubleshooting.html` — reviewed debugging summaries and safety guidance.
- `templates/open_source.html` — reviewed source/developer/deployment context.
- `templates/partials/footer.html` — reviewed crawlable internal navigation.
- `static/css/app.css` — reviewed existing card/section styling before adding AI summary styles.
- `README.md` — reviewed public project description and route/API context.
- `Qodex.summary.md` — reviewed previous SEO/content summary before replacing it with this pass.

## Files Changed
- `app.py` — added root routes for `/robots.txt`, `/sitemap.xml`, `/llms.txt`, and `/llms-full.txt` using static source files.
- `templates/base.html` — changed canonical and Open Graph URLs to stable query-free `request.base_url`, advertised sitemap and LLM summary links, and added conservative JSON-LD creator/sameAs fields.
- `templates/home.html` — added a visible AI-friendly page summary near the top.
- `templates/workflow.html` — added a visible AI-friendly workflow summary.
- `templates/documentation.html` — added a visible AI-friendly documentation summary.
- `templates/api.html` — added a visible AI-friendly API summary and fixed broken `curl` line continuations/JSON quoting in API examples.
- `templates/troubleshooting.html` — added a visible AI-friendly troubleshooting summary.
- `templates/open_source.html` — added a visible AI-friendly developer/source summary.
- `static/css/app.css` — added reusable `.ai-summary-card` styles.
- `Qodex.summary.md` — updated with this task summary.

## Files Created
- `static/robots.txt` — crawler policy allowing public discovery while disallowing auth, generated downloads, workspace/upload/tmp/admin paths, and mutation-oriented API surfaces.
- `static/sitemap.xml` — sitemap for stable public documentation pages.
- `static/llms.txt` — concise public Markdown summary for AI assistants and LLM-oriented discovery.
- `static/llms-full.txt` — longer public Markdown reference covering routes, workflow, API groups, accepted inputs, package modes, troubleshooting, and safety notes.
- `templates/partials/ai_summary.html` — reusable visible page-summary partial.

## Implementation Summary
The site now exposes root-level crawler and AI-discovery files, advertises stable canonical/sitemap/LLM links from the base template, provides concise visible summaries on high-value pages, and fixes API shell examples so assistant-style systems can understand and cite the public documentation more reliably. No browser workflow behavior, API contract, workspace creation logic, package generation logic, or route names were changed.

## Key Decisions
- Used static files plus explicit Flask root routes so discovery content is easy to inspect in the repo and available at canonical URLs.
- Allowed `OAI-SearchBot`, `GPTBot`, and `ChatGPT-User` while disallowing private/generated surfaces.
- Omitted `/build` from `sitemap.xml` because the page auto-creates workspaces via JavaScript, but left it available through internal links and `llms.txt`.
- Created `llms-full.txt` because the content could be safely assembled from public page summaries without exposing secrets, workspace paths, generated files, logs, or credentials.
- Added visible AI-summary blocks rather than hidden bot-only text.
- Skipped breadcrumb JSON-LD in this pass because the existing related-page/internal-link pattern is accurate and adding route-specific breadcrumb data cleanly would require more template plumbing than the task needs.

## Commands Run
- `git status --short` — inspected working tree before and after edits.
- `git diff --stat`, `git diff`, `git diff --cached --stat`, `git diff --cached` — inspected current unstaged/staged changes before editing.
- `find . -maxdepth 4 \( -name "robots.txt" -o -name "sitemap.xml" -o -name "llms.txt" -o -name "llms-full.txt" \) -print` — confirmed discovery files did not already exist.
- `rg "robots|sitemap|llms|canonical|og:|application/ld\\+json|FAQPage|WebApplication|BreadcrumbList" -n .` — audited current metadata/discovery signals.
- `rg "def (home|about|workflow|documentation|docking_problems|open_source|api_docs|troubleshooting|contact|modules|build)" app.py` — confirmed public route definitions.
- `rg "download|workspace|upload|login|logout|admin|api/v1/workspaces" app.py templates -n` — reviewed routes that should not be promoted to crawlers.
- `sed` inspections of `templates/base.html`, `templates/api.html`, public pages, `README.md`, and `Qodex.summary.md` — reviewed content and metadata before editing.
- `rg "request.url" templates app.py` — confirmed query-string canonical usage was removed.
- `rg 'curl .*" $|curl .*upload" $|^-H |^-F |^-d ' templates/api.html` — checked for broken API shell continuation patterns.
- `rg "PORTAL_SECRET|api key|token|/tmp/autodock|/download\\?|/api/v1/workspaces/.+download" static/robots.txt static/sitemap.xml static/llms.txt static/llms-full.txt` — checked public discovery files for sensitive/private signals; only intentional robots disallow and public API documentation references matched.
- `python -m py_compile app.py` — passed.
- `python -m flask --app app routes` — failed because `flask_login` is missing in the active Python environment.
- `test -f static/robots.txt`, `test -f static/sitemap.xml`, `test -f static/llms.txt`, `test -f static/llms-full.txt` — confirmed files exist.
- XML parse script with `xml.etree.ElementTree` — parsed `static/sitemap.xml` successfully.
- Public text-file safety script — checked `static/robots.txt`, `static/llms.txt`, and `static/llms-full.txt`.
- Jinja parse script — parsed 21 templates successfully.
- `rg '```' templates` — no stray Markdown fences found.
- `rg "url_for\\(" templates` — reviewed template endpoint usage.
- `rg "Disallow: /$" static/robots.txt robots.txt 2>/dev/null || true` — confirmed root is not disallowed.
- `git diff --check` — passed.
- `python -m unittest discover -s tests -v` — ran the existing suite; see validation results.
- `git show HEAD:__pycache__/app.cpython-313.pyc > __pycache__/app.cpython-313.pyc` — restored generated bytecode changed by validation.

## Validation Results
- `python -m py_compile app.py` passed.
- `static/sitemap.xml` parsed successfully as XML.
- `static/robots.txt`, `static/llms.txt`, and `static/llms-full.txt` passed the public text-file safety assertions.
- 21 Jinja templates parsed successfully.
- `rg "request.url"` returned no remaining matches after canonical cleanup.
- No stray Markdown fences were found in templates.
- `git diff --check` passed.
- `python -m flask --app app routes` could not run because `flask_login` is not installed in the active Python environment.
- `python -m unittest discover -s tests -v` did not fully pass:
  - App-import tests failed because `flask_login` is missing.
  - Existing PyMOL-builder tests failed with missing attributes, CLI argument mismatch, and existing expectation mismatches unrelated to this AI-discovery pass.
  - Center resolver, ligand provenance, ligand upload normalization, parser workflow, and several docking-runner/summary tests passed before the unrelated failures.

## Known Issues
- Local dynamic route testing is blocked until the active Python environment has `flask_login` and the rest of the app dependencies installed.
- The broad unittest suite still has unrelated existing PyMOL-builder failures.
- `llms.txt` and `llms-full.txt` are helpful conventions, but they do not guarantee ChatGPT/OpenAI ranking or citation.
- Root-file serving should be manually verified in the production deployment after release: `/robots.txt`, `/sitemap.xml`, `/llms.txt`, and `/llms-full.txt`.

## Manual Verification
1. Open `https://autodockvina.com/robots.txt` and confirm `OAI-SearchBot`, `GPTBot`, and `ChatGPT-User` are allowed while generated/private surfaces are disallowed.
2. Open `https://autodockvina.com/sitemap.xml` and confirm only stable public pages are listed.
3. Open `https://autodockvina.com/llms.txt` and `https://autodockvina.com/llms-full.txt` and confirm they contain only public, safe, factual summaries.
4. Open `/`, `/workflow`, `/documentation`, `/api`, `/troubleshooting`, and `/open-source` and confirm the visible Page summary blocks render correctly.
5. Open `/api` and visually confirm the `curl` examples preserve line continuations and JSON quoting.
6. Open `/build` and confirm the workflow still loads normally.

## Suggested Next Prompt
Activate the full project dependency environment, run Flask route/client smoke tests for `/robots.txt`, `/sitemap.xml`, `/llms.txt`, `/llms-full.txt`, and all public pages, then manually verify the production deployment URLs after release.
