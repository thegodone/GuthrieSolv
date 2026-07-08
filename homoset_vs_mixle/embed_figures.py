#!/usr/bin/env python3
"""Inline PNG figures into report.html as base64 data URIs (the artifact CSP blocks
external images). Idempotent: replaces the figures section if already present."""
import base64, re
from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
HTML = HERE / "report.html"

FIGS = [
    ("unit_calibration_trajectory.png",
     "Active-learning calibration trajectory — FreeSolv MAE (held-out) stays flat while coverage "
     "grows across annealed outer iterations, then the loop self-terminates."),
    ("guthrie_homoset_vs_mixle.png",
     "GuthrieSolv reconciliation vs FreeSolv truth: Homoset and mixle both hug the diagonal; the "
     "third panel colours each molecule green where mixle beats Homoset."),
]


def datauri(p: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def block() -> str:
    cards = []
    for name, cap in FIGS:
        f = OUT / name
        if not f.exists():
            continue
        cards.append(
            f'<figure style="margin:0 0 22px">'
            f'<img src="{datauri(f)}" alt="{cap}" '
            f'style="width:100%;max-width:100%;border:1px solid var(--hair);border-radius:10px"/>'
            f'<figcaption class="foot" style="margin-top:8px">{cap}</figcaption></figure>')
    return ('<section id="figures">\n<h2>Figures</h2>\n'
            '<p class="sub">Held-out calibration trajectory and the method-vs-truth scatter.</p>\n'
            + "\n".join(cards) + "\n</section>\n")


def main():
    html = HTML.read_text()
    new = block()
    if '<section id="figures">' in html:
        html = re.sub(r'<section id="figures">.*?</section>\n', new, html, flags=re.S)
    else:
        # insert before the Verdict section
        html = html.replace('<section>\n  <h2>Verdict</h2>', new + '\n<section>\n  <h2>Verdict</h2>', 1)
    HTML.write_text(html)
    print(f"embedded {len(FIGS)} figures; report.html now {len(html)//1024} KB")


if __name__ == "__main__":
    main()
