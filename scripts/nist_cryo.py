"""Cryogenic material property fits from NIST's cryogenics materials pages.

https://trc.nist.gov/cryogenics/materials/materialproperties.htm indexes
about 40 materials (stainless steels, aluminium alloys, copper, Inconel
718, Ti-6Al-4V, G-10, polymers, ...). Each page prints curve-fit
coefficients in HTML tables: rows a..i (or a..f), one column per property,
with units, data range, equation range and the fit error relative to the
data. The fit's FORM is stated in text after the table, and it differs
between properties and pages. Most use log10 y = sum_n c_n (log10 T)^n.
Linear expansion is a polynomial in T with a constant below T_low. OFHC
copper conductivity is a rational form in T^0.5.

    uv run scripts/nist_cryo.py list
    uv run scripts/nist_cryo.py fetch "310 Stainless" --out <dir>   # -> <dir>/<slug>.json
    uv run scripts/nist_cryo.py eval "310 Stainless" "Thermal Conductivity" 77 300

`fetch` stores every table with the verbatim form text that follows it,
plus a `form` tag. `eval` evaluates only the tags it recognises (log10poly,
poly_tlow) and refuses the rest, so read those pages. Coefficients are
NIST's, so files belong under sources/ (git-ignored, #846) with a
validation/sources.yaml entry citing the page and access date.
"""

from __future__ import annotations

import argparse
import contextlib
import html
import http.client
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

BASE = "https://trc.nist.gov"
INDEX = BASE + "/cryogenics/materials/materialproperties.htm"
_TIMEOUT = 30
_COEF_ROWS = "abcdefghi"


def _get(url: str, attempts: int = 6) -> str:
    # Replies sometimes arrive cut short (IncompleteRead on a chunked
    # HTTP/1.1 body, or a silently short HTTP/1.0 one) or time out, seen through the
    # sandbox proxy. Only a page that reaches </html> is accepted; retry.
    http.client.HTTPConnection._http_vsn = 10
    http.client.HTTPConnection._http_vsn_str = "HTTP/1.0"
    req = urllib.request.Request(url, headers={"User-Agent": "combaero-nist-cryo/1"})
    for i in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:  # noqa: S310 -- fixed https host
                body = r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError:
            raise
        except (http.client.HTTPException, OSError):
            body = ""
        if "</html>" in body.lower():
            return body
        time.sleep(0.5 * (i + 1))
    raise RuntimeError(f"{url}: no complete reply after {attempts} attempts")


def _text(fragment: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", fragment)).split())


class _Tables(HTMLParser):
    """Top-level tables as lists of rows of cell text, plus their end offsets."""

    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self.row: list[str] | None = None
        self.cell: list[str] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "table":
            self.tables.append([])
        elif tag == "tr" and self.tables:
            self.row = []
            self.tables[-1].append(self.row)
        elif tag in ("td", "th") and self.row is not None:
            self.cell = []

    def handle_endtag(self, tag: str) -> None:
        if tag in ("td", "th") and self.cell is not None and self.row is not None:
            self.row.append(" ".join("".join(self.cell).split()))
            self.cell = None

    def handle_data(self, data: str) -> None:
        if self.cell is not None:
            self.cell.append(data)


def list_materials() -> list[tuple[str, str]]:
    page = _get(INDEX)
    out = []
    for href, label in re.findall(
        r'href="(/cryogenics/materials/[^"]+)"[^>]*>(.*?)</a>', page, re.S
    ):
        name = _text(label) or href.rsplit("/", 2)[-2]
        out.append((name, BASE + urllib.parse.quote(href, safe="/%()_.-")))
    return out


def _classify(form_text: str) -> str:
    t = form_text.replace(" ", "")
    if "T0.5" in t or "T^0.5" in t:
        return "rational_sqrtT"
    if "log10y=a+b(log10T)" in t:
        return "log10poly"
    if "y=a+bT+cT2" in t:
        return "poly_tlow"
    return "unknown"


def parse_page(page: str) -> list[dict]:
    """Every coefficient table on a page, with the form text that follows it."""
    p = _Tables()
    p.feed(page)
    # Split the raw HTML at table starts to recover the prose between tables.
    chunks = re.split(r"(?i)<table", page)[1:]
    is_coef = [bool(t) and len(t[0]) >= 2 and any(r and r[0] == "a" for r in t) for t in p.tables]
    out: list[dict] = []
    for i, table in enumerate(p.tables):
        if not is_coef[i]:
            continue
        props = table[0][1:]
        for col, prop in enumerate(props, start=1):
            entry: dict = {"property": prop, "coefficients": {}}
            for row in table[1:]:
                if len(row) <= col:
                    continue
                raw_key, val = row[0].strip(". ").lower(), row[col]
                key = raw_key.rstrip(">").strip()
                if raw_key.endswith(">"):  # "f>": the value below T_low
                    entry["below_T_low"] = float(val) if val.strip(".") else None
                elif key in _COEF_ROWS and len(key) == 1:
                    with contextlib.suppress(ValueError):
                        entry["coefficients"][key] = float(val)
                elif key.startswith("unit"):
                    entry["units"] = val
                elif key.startswith("tlow"):
                    entry["T_low_K"] = float(val) if val.strip(".") else None
                else:
                    entry[key] = val
            out.append(entry)
        # The form text: prose and small tables up to the next coefficient table.
        parts = [chunks[i].split("</table>", 1)[-1]] if i < len(chunks) else []
        j = i + 1
        while j < len(p.tables) and j < len(chunks) and not is_coef[j]:
            parts.append(chunks[j])
            j += 1
        form = _text("<table".join(parts))[:600]
        for e in out[-len(props) :]:
            e["form_text"] = form
    # A form stated once applies to the tables before the next statement.
    last = ""
    for e in out:
        if "of the form" in e["form_text"]:
            last = e["form_text"]
        e["form_text"] = e["form_text"] if "of the form" in e["form_text"] else last
        e["form"] = _classify(e["form_text"])
    return out


def _resolve(name: str) -> tuple[str, str]:
    q = name.lower()
    hits = [
        (n, u)
        for n, u in list_materials()
        if q in n.lower() or q in urllib.parse.unquote(u).lower()
    ]
    if not hits:
        raise SystemExit(f"no material matching {name!r}; see `list`")
    if len(hits) > 1:
        exact = [h for h in hits if h[0].lower() == name.lower()]
        if len(exact) != 1:
            raise SystemExit(f"{name!r} matches several: " + ", ".join(h[0] for h in hits))
        hits = exact
    return hits[0]


def evaluate(entry: dict, T: float) -> float:
    c = entry["coefficients"]
    if entry["form"] == "log10poly":
        x = math.log10(T)
        return 10.0 ** sum(c.get(k, 0.0) * x**i for i, k in enumerate(_COEF_ROWS))
    if entry["form"] == "rational_sqrtT":  # OFHC copper conductivity
        r = math.sqrt(T)
        num = c["a"] + c["c"] * r + c["e"] * T + c["g"] * T * r + c["i"] * T * T
        den = 1.0 + c["b"] * r + c["d"] * T + c["f"] * T * r + c["h"] * T * T
        return 10.0 ** (num / den)
    if entry["form"] == "poly_tlow":
        t_low = entry.get("T_low_K")
        if t_low is not None and t_low > T:
            return entry["below_T_low"]
        return sum(c.get(k, 0.0) * T**i for i, k in enumerate("abcde"))
    raise SystemExit(
        f"form {entry['form']!r} not evaluated here; read the page: {entry['form_text'][:200]}"
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    f = sub.add_parser("fetch")
    f.add_argument("material")
    f.add_argument("--out", type=Path, default=Path("sources/materials/nist_cryogenics"))
    e = sub.add_parser("eval")
    e.add_argument("material")
    e.add_argument("property")
    e.add_argument("T", type=float, nargs="+", help="temperature in K")
    a = ap.parse_args(argv)
    if a.cmd == "list":
        for n, u in list_materials():
            print(f"{n:40s} {u}")
        return 0
    name, url = _resolve(a.material)
    entries = parse_page(_get(url))
    if a.cmd == "fetch":
        a.out.mkdir(parents=True, exist_ok=True)
        path = a.out / (re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_") + ".json")
        path.write_text(
            json.dumps({"material": name, "url": url, "fits": entries}, indent=1) + "\n"
        )
        print(path)
        return 0
    matches = [x for x in entries if x["property"].lower() == a.property.lower()]
    if not matches:
        raise SystemExit(f"no {a.property!r}; have: " + ", ".join(x["property"] for x in entries))
    for T in a.T:
        for m in matches:
            rng = m.get("equation range") or m.get("equation range (k)")
            if rng is None and "low range" in m:
                rng = m["low range"].rstrip("Kk") + "-" + m.get("high range", "?").rstrip("Kk")
            rng = rng or "?"
            lo, _, hi = rng.partition("-")
            inside = lo and hi and float(lo) <= T <= float(hi.rstrip("Kk"))
            if len(matches) > 1 and not inside:
                continue
            flag = "" if inside else f"  OUTSIDE equation range {rng} K"
            print(f"{T:g} K  {evaluate(m, T):.6g} {m.get('units', '')}{flag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
