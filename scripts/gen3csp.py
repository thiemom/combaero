"""Measured material properties from the Georgia Tech Gen3 CSP database.

https://gen3csp.gatech.edu/compare/ plots thermal conductivity, specific
heat and thermal diffusivity against temperature, each with its 95 %
confidence band, for alloys, ceramics and a few fluids measured for the
DOE Gen3 concentrating-solar-power program. The page has no API: a GET of

    /compare/?materials=<id>&property=y1|y2|y3&errors=on

returns the chart script with the data inline as `dat1 = [[T_C, y], ...]`
and the band as `bdat1 = [[a, b], ...]` (upper and lower, in either
order). One material per request; a comma list is ignored. This script
reads those arrays (plain number parsing, nothing executed) and writes
one CSV per material.

    uv run scripts/gen3csp.py list                     # material ids
    uv run scripts/gen3csp.py files                    # spreadsheets (not every material has one)
    uv run scripts/gen3csp.py fetch ss310 haynes230 --out <dir>

The chart rounds some materials to 2-3 digits (SS310) where the
spreadsheet carries full precision; both agree to that rounding. The
data are Georgia Tech's measurements: downloads belong under sources/
(git-ignored, #846) and each material used gets a validation/sources.yaml
entry. Density (y4) is not measured on most materials (all zeros) and is
not fetched.
"""

from __future__ import annotations

import argparse
import csv
import html
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from pathlib import Path

BASE = "https://gen3csp.gatech.edu"
PROPERTIES = {
    "k": ("y1", "Thermal Conductivity", "W/m-K"),
    "cp": ("y2", "Specific Heat", "J/g-K"),
    "alpha": ("y3", "Thermal Diffusivity", "mm2/s"),
}
_ARRAY = r"\[\s*((?:\[\s*[-+0-9.eE]+\s*,\s*[-+0-9.eE]+\s*\]\s*,?\s*)*)\]"
_DAT = re.compile(r"^[ \t]*dat1[ \t]*=[ \t]*" + _ARRAY, re.MULTILINE)
_BDAT = re.compile(r"^[ \t]*bdat1[ \t]*=[ \t]*" + _ARRAY, re.MULTILINE)
_PAIR = re.compile(r"\[\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\]")
_SERIES_LABEL = re.compile(r"^[ \t]*label:[ \t]*'([^']*)'", re.MULTILINE)
_OPTION = re.compile(r"<option value='([^']+)'\s*>([^<]+)</option>")
_TIMEOUT = 30


def _html_complete(body: str) -> bool:
    return "</html>" in body.lower()


def _json_complete(body: str) -> bool:
    try:
        json.loads(body)
    except ValueError:
        return False
    return True


def _get(url: str, attempts: int = 6, complete: Callable[[str], bool] = _html_complete) -> str:
    # Replies sometimes arrive cut short (IncompleteRead on a chunked
    # HTTP/1.1 body, or a silently short HTTP/1.0 one) or time out, seen
    # through the sandbox proxy. Only a complete reply (a page that reaches
    # </html>, or JSON that parses) is accepted; retry. An HTTP error status
    # is an answer, not a transport fault, and is raised at once.
    http.client.HTTPConnection._http_vsn = 10
    http.client.HTTPConnection._http_vsn_str = "HTTP/1.0"
    req = urllib.request.Request(url, headers={"User-Agent": "combaero-gen3csp/1"})
    for i in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:  # noqa: S310 -- fixed https host
                body = r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError:
            raise
        except (http.client.HTTPException, OSError):
            body = ""
        if complete(body):
            return body
        time.sleep(0.5 * (i + 1))
    raise RuntimeError(f"{url}: no complete reply after {attempts} attempts")


def _pairs(body: str) -> list[tuple[float, float]]:
    return [(float(a), float(b)) for a, b in _PAIR.findall(body)]


def list_materials() -> list[tuple[str, str]]:
    """(id, display name) for every material in the compare form."""
    page = _get(BASE + "/compare/")
    form = page[page.find("id='input_2_1'") :]
    form = form[: form.find("</select>")]
    return [(i, html.unescape(n).strip()) for i, n in _OPTION.findall(form)]


def list_files() -> list[str]:
    """Spreadsheet and PDF uploads, via the WordPress media API."""
    out: list[str] = []
    page = 1
    while True:
        q = urllib.parse.urlencode({"per_page": 100, "page": page, "_fields": "source_url"})
        try:
            items = json.loads(_get(f"{BASE}/wp-json/wp/v2/media?{q}", complete=_json_complete))
        except urllib.error.HTTPError:
            break
        if not items:
            break
        out += [
            m["source_url"]
            for m in items
            if m["source_url"].lower().endswith((".xlsx", ".pdf", ".csv"))
        ]
        page += 1
    return out


def fetch_property(material: str, prop: str) -> tuple[str, list[tuple[float, float, float, float]]]:
    """(series label, rows of (T_C, value, lower, upper)) for one property."""
    code = PROPERTIES[prop][0]
    q = urllib.parse.urlencode({"materials": material, "property": code, "errors": "on"})
    page = _get(f"{BASE}/compare/?{q}")
    dat, bdat = _DAT.search(page), _BDAT.search(page)
    if dat is None:
        raise ValueError(f"{material}/{prop}: no data series on the page (unknown id?)")
    values = _pairs(dat.group(1))
    band = _pairs(bdat.group(1)) if bdat else []
    if band and len(band) != len(values):
        raise ValueError(f"{material}/{prop}: {len(values)} points but {len(band)} band rows")
    labels = _SERIES_LABEL.findall(page)
    rows = []
    for i, (t, y) in enumerate(values):
        lo, hi = (min(band[i]), max(band[i])) if band else (float("nan"), float("nan"))
        rows.append((t, y, lo, hi))
    return (labels[0] if labels else material), rows


def fetch(material: str, out_dir: Path) -> Path:
    by_t: dict[float, dict[str, float]] = {}
    label = material
    for prop in PROPERTIES:
        label, rows = fetch_property(material, prop)
        for t, y, lo, hi in rows:
            if y == 0.0 and lo == 0.0 and hi == 0.0:
                continue  # not measured
            by_t.setdefault(t, {}).update({prop: y, prop + "_lo": lo, prop + "_hi": hi})
    cols = [c for p in PROPERTIES for c in (p, p + "_lo", p + "_hi")]
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{material}.csv"
    with path.open("w", newline="") as f:
        f.write(
            f"# Georgia Tech Gen3 CSP database, {label} ({BASE}/compare/?materials={material})\n"
        )
        f.write(
            "# units: T_C degC; "
            + "; ".join(f"{p} {u}" for p, (_, _, u) in PROPERTIES.items())
            + "\n"
        )
        f.write("# _lo/_hi: 95% confidence band as plotted; blank = not measured at that T\n")
        w = csv.writer(f)
        w.writerow(["T_C", *cols])
        for t in sorted(by_t):
            w.writerow([f"{t:g}", *[f"{by_t[t][c]:.10g}" if c in by_t[t] else "" for c in cols]])
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="material ids and names")
    sub.add_parser("files", help="downloadable spreadsheets and PDFs")
    f = sub.add_parser("fetch", help="write <out>/<id>.csv for each material")
    f.add_argument("materials", nargs="+")
    f.add_argument("--out", type=Path, default=Path("sources/materials/gatech_gen3csp"))
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for i, n in list_materials():
            print(f"{i:16s} {n}")
    elif a.cmd == "files":
        print("\n".join(list_files()))
    else:
        for m in a.materials:
            print(fetch(m, a.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
