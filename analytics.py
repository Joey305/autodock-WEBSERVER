"""Best-effort, aggregate-safe telemetry for AutoDock-Vina PrepServer."""
from __future__ import annotations

import json
import os
import secrets
import threading
import re
import time
from http.cookiejar import CookieJar
from urllib.request import Request, build_opener, HTTPCookieProcessor

from flask import Response, g, request

EVENTS = {
    "page_view": {""}, "build_opened": {""},
    "workspace_created": {"browser", "headless"}, "receptor_staged": {"upload", "fetch"}, "box_defined": {"interactive", "api"}, "receptor_prepared": {"browser", "headless"}, "ligand_staged": {"upload", "curated"}, "package_mode_selected": {"portable", "hpc"},
    "workflow_selected": {"interactive", "bound_ligand", "batch", "headless"},
    "package_build_started": {"interactive", "bound_ligand", "batch", "headless"},
    "package_build_completed": {"interactive", "bound_ligand", "batch", "headless"},
    "package_build_failed": {"interactive", "bound_ligand", "batch", "headless"},
    "package_downloaded": {""}, "results_opened": {""}, "viewer_opened": {""},
    "api_docs_opened": {""}, "examples_opened": {""},
    "headless_client_downloaded": {"interactive", "bound_ligand"},
    "api_example_copied": {"interactive", "bound_ligand"},
    "workflow_failed": {"validation", "timeout", "dependency", "server"},
    "package_build_duration": {"under_10s", "10_to_30s", "30_to_120s", "over_120s"},
}


def enabled() -> bool:
    return os.environ.get("AUTODOCK_ANALYTICS_ENABLED", "").lower() in {"1", "true", "yes", "on"}


def _receiver() -> tuple[str, str] | None:
    url = os.environ.get("AUTODOCK_ANALYTICS_RANDY_URL", "").strip().rstrip("/")
    token = os.environ.get("AUTODOCK_ANALYTICS_RANDY_TOKEN", "").strip()
    return (url, token) if enabled() and url and token else None


def _identifier(name: str) -> str:
    cached = getattr(g, name, None)
    if cached:
        return cached
    value = request.cookies.get(name, "")
    value = value if len(value) >= 12 and value.replace("_", "").replace("-", "").isalnum() else secrets.token_urlsafe(18)
    setattr(g, name, value)
    return value


def _device(user_agent: str) -> str:
    value = user_agent.lower()
    return "mobile" if any(x in value for x in ("mobile", "android", "iphone")) else "tablet" if "ipad" in value else "desktop"


def _browser(user_agent: str) -> str:
    value = user_agent.lower()
    return "firefox" if "firefox" in value else "safari" if "safari" in value and "chrome" not in value else "chrome" if "chrome" in value else "other"


def _platform(user_agent: str) -> str:
    """Return a broad OS family only; never preserve the user-agent string."""
    value = user_agent.lower()
    if "android" in value: return "android"
    if any(item in value for item in ("iphone", "ipad", "ios")): return "ios"
    if "windows" in value: return "windows"
    if "mac os" in value or "macintosh" in value: return "macos"
    if "linux" in value or "x11" in value: return "linux"
    return "other"


def _referrer(value: str) -> str:
    # Store a broad acquisition category only; never retain a referring domain.
    from urllib.parse import urlparse
    host = urlparse(value).hostname or ""
    public_host = urlparse(os.environ.get("AUTODOCK_PUBLIC_BASE_URL", "")).hostname or ""
    if not host: return "direct"
    if public_host and host == public_host: return "internal"
    if "youtube" in host: return "youtube"
    if "github" in host: return "github"
    if any(item in host for item in ("google.", "bing.", "duckduckgo.", "yahoo.")): return "search"
    if any(item in host for item in ("x.com", "twitter.", "linkedin.", "facebook.", "reddit.")): return "social"
    return "other_referral"


def _duration_bucket(seconds: float) -> str:
    return "under_10s" if seconds < 10 else "10_to_30s" if seconds < 30 else "30_to_120s" if seconds < 120 else "over_120s"


def begin_request() -> None:
    """Capture an in-process elapsed time only for coarse build-duration buckets."""
    g.autodock_analytics_started = time.monotonic()


def _send(payload: dict[str, str]) -> None:
    receiver = _receiver()
    if not receiver:
        return
    url, token = receiver
    try:
        req = Request(url + "/event", data=json.dumps(payload).encode(), headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json", "User-Agent": "autodock-analytics/1.0"}, method="POST")
        build_opener(HTTPCookieProcessor(CookieJar())).open(req, timeout=2.0).close()
    except Exception:
        pass


def _campaign() -> str:
    value = request.args.get("utm_campaign", "").strip().lower()
    return value if re.fullmatch(r"youtube-[a-z0-9-]{1,60}", value) else ""

def track(event_type: str, feature: str = "", status: str = "", path: str | None = None, campaign: str | None = None) -> tuple[str, str]:
    """Queue a validated event without reading request bodies or scientific data."""
    visitor, session = _identifier("autodock_visitor"), _identifier("autodock_session")
    if event_type not in EVENTS or feature not in EVENTS[event_type]:
        return visitor, session
    safe_path = path or request.path
    if not safe_path.startswith("/") or "?" in safe_path or len(safe_path) > 240:
        return visitor, session
    payload = {"event_id": secrets.token_urlsafe(18), "event_type": event_type, "visitor_id": visitor, "session_id": session, "path": safe_path, "referrer": _referrer(request.referrer or ""), "device": _device(request.user_agent.string), "browser": _browser(request.user_agent.string), "platform": _platform(request.user_agent.string), "feature": feature, "status": status, "campaign": campaign if campaign is not None else _campaign()}
    threading.Thread(target=_send, args=(payload,), daemon=True).start()
    return visitor, session


def attach(response: Response) -> Response:
    if not enabled():
        return response
    visitor, session = _identifier("autodock_visitor"), _identifier("autodock_session")
    if not request.cookies.get("autodock_visitor"):
        response.set_cookie("autodock_visitor", visitor, max_age=31536000, secure=True, httponly=True, samesite="Lax")
    if not request.cookies.get("autodock_session"):
        response.set_cookie("autodock_session", session, max_age=1800, secure=True, httponly=True, samesite="Lax")
    if request.method == "GET" and response.status_code < 400 and not request.path.startswith(("/static/", "/api/", "/download")):
        track("page_view")
        page_events = {"/build": "build_opened", "/results": "results_opened", "/viz/project": "viewer_opened", "/api": "api_docs_opened", "/documentation": "api_docs_opened", "/workflow": "examples_opened"}
        if request.path in page_events:
            track(page_events[request.path])
    build_routes = {"/api/build": "interactive", "/api/v1/headless/package": "headless"}
    milestones = {"/api/workspace": ("workspace_created", "browser"), "/api/receptors/upload": ("receptor_staged", "upload"), "/api/receptors/fetch": ("receptor_staged", "fetch"), "/api/receptor/center": ("box_defined", "interactive"), "/api/prep/start": ("receptor_prepared", "browser"), "/api/prep/start_one": ("receptor_prepared", "browser"), "/api/ligands/upload": ("ligand_staged", "upload"), "/api/ligands/curated": ("ligand_staged", "curated")}
    if request.method == "POST" and response.status_code < 400:
        if request.path in milestones: track(*milestones[request.path])
        elif request.path == "/api/v1/workspaces": track("workspace_created", "headless")
        elif request.path.endswith("/receptors/upload"): track("receptor_staged", "upload")
        elif request.path.endswith("/receptors/fetch"): track("receptor_staged", "fetch")
        elif request.path.endswith(("/centers/resolve", "/centers/save")): track("box_defined", "api")
        elif request.path.endswith("/prep/start"): track("receptor_prepared", "headless")
        elif request.path.endswith("/ligands/upload"): track("ligand_staged", "upload")
        elif request.path.endswith(("/ligands/curated", "/ligands/extract")): track("ligand_staged", "curated")
    if request.method == "POST" and (request.path in build_routes or (request.path.startswith("/api/v1/workspaces/") and request.path.endswith("/build"))):
        workflow = build_routes.get(request.path, "headless")
        track("workflow_selected", workflow)
        track("package_build_started", workflow)
        track("package_build_completed" if response.status_code < 400 else "package_build_failed", workflow, "success" if response.status_code < 400 else "failure")
        # Browser builds measure elapsed time in the page so navigation and
        # response delivery are included.  This server-side fallback is only
        # needed for the headless API.
        started = getattr(g, "autodock_analytics_started", None)
        if workflow == "headless" and started is not None:
            track("package_build_duration", _duration_bucket(time.monotonic() - started), "success" if response.status_code < 400 else "failure")
    if request.method == "GET" and response.status_code < 400 and (request.path == "/download" or request.path.endswith("/download")):
        track("package_downloaded")
    if request.method == "GET" and response.status_code < 400 and request.path.startswith("/api/v1/clients/"):
        workflow = "bound_ligand" if "bound_ligand" in request.path else "interactive"
        track("workflow_selected", workflow)
        track("headless_client_downloaded", workflow)
    return response


def browser_event(data: object) -> bool:
    if not isinstance(data, dict):
        return False
    event_type, feature = data.get("event_type"), data.get("feature", "")
    if event_type not in EVENTS or feature not in EVENTS[event_type]:
        return False
    campaign = str(data.get("campaign") or "")
    track(str(event_type), str(feature), campaign=campaign if re.fullmatch(r"youtube-[a-z0-9-]{1,60}", campaign) else "")
    return True
