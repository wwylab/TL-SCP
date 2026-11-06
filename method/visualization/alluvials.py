"""Draw mutation-level alluvial diagrams from aggregate annotations."""

from pathlib import Path
import csv
import sys

COLORS = {"Early": "F8766D", "Mid": "00BA38", "Late": "619CFF", "Tied": "666666"}
W = 8.4
UNIT = 0.075


def lines(rows):
    # Approximately three names per row; never let a line exceed the text column.
    result, row, length = [], [], 0
    for r in sorted(rows, key=lambda r: r["mutation"]):
        n = len(r["mutation"]) + (2 if row else 0)
        if row and length + n > 19:
            result.append(row)
            row = []
            length = 0
            n = len(r["mutation"])
        row.append(r)
        length += n
    if row:
        result.append(row)
    return result


def geometry(rows, left_order, right_order):
    groups = [
        {g: [r for r in rows if r[side] == g] for g in order}
        for side, order in [("left", left_order), ("right", right_order)]
    ]
    # Spacing accommodates text independently of ribbon thickness; counts set bar height.
    spans = [
        {g: max(len(q) * UNIT, 0.62 + len(lines(q)) * 0.34) for g, q in group.items()}
        for group in groups
    ]
    h = max(sum(s.values()) + 0.24 * (len(s) - 1) for s in spans)
    return groups, spans, h


def panel(spec, geo, panel_id, H):
    rows, lo, ro, ln, rn, title = spec
    groups, spans, _ = geo
    centers = []
    for sizes in spans:
        gap = (H - sum(sizes.values())) / (len(sizes) - 1)
        y = H
        cc = {}
        for g, h in sizes.items():
            cc[g] = y - h / 2
            y -= h + gap
        centers.append(cc)
    left, right = groups
    lc, rc = centers
    curL = {g: lc[g] + len(q) * UNIT / 2 for g, q in left.items()}
    curR = {g: rc[g] + len(q) * UNIT / 2 for g, q in right.items()}
    tex = [
        r"\begin{tikzpicture}[x=1cm,y=1cm]",
        rf"\path[use as bounding box] (0,-.10) rectangle ({W},{H+1.30:.5f});",
    ]
    for l in lo:
        for r in ro:
            n = sum(z["left"] == l and z["right"] == r for z in rows)
            if not n:
                continue
            lt = curL[l]
            lb = lt - n * UNIT
            rt = curR[r]
            rb = rt - n * UNIT
            col = "panel" + str(panel_id) + r
            tex.append(
                r"\path[fill="
                + col
                + r",fill opacity=.57,draw=none] "
                + f"(2.95,{lt:.5f}) .. controls (3.72,{lt:.5f}) and (4.28,{rt:.5f}) .. (5.05,{rt:.5f}) -- "
                + f"(5.05,{rb:.5f}) .. controls (4.28,{rb:.5f}) and (3.72,{lb:.5f}) .. (2.95,{lb:.5f}) -- cycle;"
            )
            curL[l] = lb
            curR[r] = rb
    for side, gs, cs, cx, bounds in [
        ("left", left, lc, 1.36, (2.79, 2.95)),
        ("right", right, rc, 6.64, (5.05, 5.21)),
    ]:
        for g, q in gs.items():
            cy = cs[g]
            h = len(q) * UNIT
            fill = "black!3" if side == "left" else "panel" + str(panel_id) + g + "!9"
            tex.append(
                r"\draw[fill="
                + fill
                + r",draw=black!50,line width=.35pt] "
                + f"({bounds[0]},{cy-h/2:.5f}) rectangle ({bounds[1]},{cy+h/2:.5f});"
            )
            content = (
                r"{\fontsize{8.5}{10}\selectfont\bfseries "
                + g
                + f" ({len(q)})"
                + r"}\\[3pt]"
            )
            content += r"\\[1pt]".join(
                ", ".join(
                    r"\textcolor{panel"
                    + str(panel_id)
                    + z["right"]
                    + r"}{\texttt{"
                    + z["mutation"]
                    + "}}"
                    for z in row
                )
                for row in lines(q)
            )
            tex.append(
                r"\node[align=center,inner sep=0pt,font=\fontsize{7.4}{8.5}\selectfont] at "
                + f"({cx},{cy:.5f})"
                + r"{\begin{tabular}{@{}c@{}}"
                + content
                + r"\end{tabular}};"
            )
    for cx, text in [(1.36, ln), (6.64, rn)]:
        tex.append(
            r"\node[align=center,inner sep=0pt,font=\fontsize{9}{10.5}\selectfont\bfseries] at "
            + f"({cx},{H+.55:.5f})"
            + " {"
            + text
            + "};"
        )
    tex.append(r"\end{tikzpicture}")
    body = "\n".join(tex)
    for z in rows:
        assert body.count(r"\texttt{" + z["mutation"] + "}") == 2
    return body


def load(path):
    return list(csv.DictReader(path.open()))


def label(s, scp=False):
    if s == "Unknown":
        return "Group 4" if scp else "Unclassified"
    return s.replace("-onset", "")


def build(specs, name):
    geos = [geometry(*s[:3]) for s in specs]
    H = max(g[2] for g in geos)
    defs = "\n".join(
        r"\definecolor{panel" + str(i) + g + "}{HTML}{" + v + "}"
        for i in range(len(specs))
        for g, v in COLORS.items()
    )
    bodies = [panel(s, g, i, H) for i, (s, g) in enumerate(zip(specs, geos))]
    header = (
        r"\documentclass[11pt,border=0pt]{standalone}"
        + "\n"
        + r"\usepackage{mathpazo}\usepackage[T1]{fontenc}\usepackage{tikz}\usepackage{xcolor}\usepackage{graphicx}"
        + "\n"
        + r"\begin{document}"
        + "\n"
        + defs
        + "\n"
    )
    if len(bodies) == 2:
        content = (
            r"\begin{minipage}{6.5in}\centering"
            + "\n"
            + r"\begin{minipage}[t]{0.49\textwidth}\vspace{0pt}\centering (a) YTA versus TL-SCP\par\smallskip\resizebox{\linewidth}{!}{"
            + "\n"
            + bodies[0]
            + "\n"
            + r"}\end{minipage}\hfill"
            + "\n"
            + r"\begin{minipage}[t]{0.49\textwidth}\vspace{0pt}\centering (b) SCP versus TL-SCP\par\smallskip\resizebox{\linewidth}{!}{"
            + "\n"
            + bodies[1]
            + "\n"
            + r"}\end{minipage}\end{minipage}"
        )
    else:
        content = bodies[0]
    if len(specs) == 2:
        content = content.replace("(a) YTA versus TL-SCP", specs[0][5])
        content = content.replace("(b) SCP versus TL-SCP", specs[1][5])
    Path(name).write_text(header + content + "\n" + r"\end{document}" + "\n")


def main():
    rows = load(Path(sys.argv[1]))
    specs = []
    for panel_id in sorted({row["panel"] for row in rows}):
        subset = [row for row in rows if row["panel"] == panel_id]
        first = subset[0]
        specs.append(
            (
                subset,
                first["left_order"].split("|"),
                first["right_order"].split("|"),
                first["left_name"],
                first["right_name"],
                first["title"],
            )
        )
    build(specs, sys.argv[2])


if __name__ == "__main__":
    main()
