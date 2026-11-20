import io
import os
import re
import shutil

SRC = r"d:\atmos\npj_submission\sn-article_revised-5.tex"
DST = r"d:\atmos\aas_submission"
TPL = r"D:\atmos\aas_template\AsampleLATEX\AampleLATEX"

os.makedirs(DST, exist_ok=True)
for f in ("atmosphere.sty", "command.tex"):
    shutil.copy2(os.path.join(TPL, f), os.path.join(DST, f))
shutil.copy2(r"d:\atmos\npj_submission\sn-bibliography.bib", DST)
for f in ("fig3_spectrum.png", "fig_signflip_finite.png", "fig3b_cities.png",
          "fig_cities_24win.png", "fig_trunc_cities.png", "fig_cities_lvg.png",
          "fig_alloc.png", "fig_optimal.png"):
    shutil.copy2(os.path.join(r"d:\atmos\npj_submission", f), DST)

s = io.open(SRC, encoding="utf-8").read()

ab = re.search(r"\\abstract\{(.*?)\}\s*\n", s, re.S).group(1).strip()
body = s[s.index(r"\section{Introduction}"):s.index(r"\backmatter")].rstrip() + "\n"
body = body.replace(r"\cite{Kang2024}", r"\citep{Kang2024}")
body = body.replace(r"\citep{Jansma2025}", r"\citep{Jansma2025}")

head = r"""\documentclass[12pt]{article}
\usepackage{atmosphere}
\input{command}
\usepackage{amsmath}
\usepackage{booktabs}
\usepackage[authoryear,round]{natbib}
\bibliographystyle{apalike}

\begin{document}

\title{\vspace*{-20mm}\fontsize{12}{12}\bf
Exact response tables reveal high-order interactions in time-resolved ozone control
\vspace*{-3mm}}

\date{}

\author{\\[-17mm]\fontsize{12}{12}
First AUTHOR\footnote{\setlength{\baselineskip}{12pt}Corresponding author: First AUTHOR\newline\hspace*{6.2mm}Email: iauthor@gmail.com}$^{*1}$,
Second AUTHOR$^{2}$, and Third AUTHOR$^{3}$}

\maketitle

\add{$^1$College of Environmental and Resource Sciences, Fujian Normal University, No. 8 Xuefu South Road, Fuzhou 350117, China\\[-3mm]
$^2$Department, Organization, Street, City 10587, China\\[-3mm]
$^3$Department, Organization, Street, City 610101, China}

\abs{""" + ab + r"""}

\vspace*{-10mm}

\keywords{Ozone sensitivity, MDA8, emission control, M\"obius decomposition, interaction spectrum, cooperative games}

\vspace*{-5mm}

\noindent
{\bf Article Highlights:}

\noindent$\bullet$
All 65,536 discrete NO$_x$ and VOC shutdown packages are enumerated exactly in a CB05 box model.

\noindent$\bullet$
Regulatory MDA8 is dominated by high-order interactions; low-order truncation fails where terms cancel.

\noindent$\bullet$
In NO$_x$-saturated Dallas, first-order predictions reverse the sign of 13.4\% of finite packages.

\noindent$\bullet$
Pairwise interactions recover the exact budget-constrained optimum in every game and budget tested.

"""

tail = r"""
\noindent\textbf{Supplementary information.} The Supplementary Information contains additional spectral diagnostics (Tables S1--S2), the value-level integration checks (Table S3), and a description of the released response tables (Section S3).

\noindent\textbf{Data availability.} The complete response tables for the four $n=16$ games (65{,}536 coalitions per game, columns as described in Supplementary Section~S3) are provided as Supplementary Data and are available at [repository DOI / accession to be inserted].

\noindent\textbf{Code availability.} The box-model driver, the KPP chemical mechanism, and the analysis scripts that generate the spectra and the tables in this manuscript are available at [repository URL to be inserted].

\noindent\textbf{Competing interests.} [Insert declaration.]

\noindent\textbf{Author contributions.} [Insert statement.]

\acknowledgmentsname{[Insert funding and acknowledgements.]}

\bibliography{sn-bibliography}

\end{document}
"""

out = head + body + tail
io.open(os.path.join(DST, "aas-article.tex"), "w", encoding="utf-8").write(out)
print("wrote", os.path.join(DST, "aas-article.tex"), len(out), "chars")
