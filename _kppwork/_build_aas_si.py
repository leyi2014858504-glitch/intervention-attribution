import io
import os

SRC = r"d:\atmos\npj_submission\sn-si.tex"
DST = r"d:\atmos\aas_submission\aas-si.tex"

s = io.open(SRC, encoding="utf-8").read()
body = s[s.index(r"\section{Additional spectral diagnostics}"):s.index(r"\end{document}")]
body = body.replace("Sec.~3 of the main text", "Sec.~2 of the main text")

head = r"""%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary Information (AAS submission)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\documentclass[12pt]{article}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{booktabs}
\usepackage{graphicx}
\usepackage[margin=1in]{geometry}
\renewcommand{\thetable}{S\arabic{table}}
\setlength{\emergencystretch}{2em}

\begin{document}

\begin{center}
{\large\bf Supplementary Information}\\[3mm]
{\large for: Exact response tables reveal high-order interactions in time-resolved ozone control}
\end{center}

\vspace{4mm}

"""

tail = "\n\\end{document}\n"

io.open(DST, "w", encoding="utf-8").write(head + body + tail)
print("wrote", DST)
