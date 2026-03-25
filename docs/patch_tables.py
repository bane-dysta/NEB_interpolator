from pathlib import Path
import re
p=Path('out/tmp.tex')
text=p.read_text()
# Improve wrapping in tables by converting default ll/lll longtables to p-columns.
text=text.replace(r'\begin{longtable}[]{@{}ll@{}}', r'\begin{longtable}[]{@{}>{\raggedright\arraybackslash}p{0.24\linewidth}>{\raggedright\arraybackslash}p{0.70\linewidth}@{}}')
text=text.replace(r'\begin{longtable}[]{@{}lll@{}}', r'\begin{longtable}[]{@{}>{\raggedright\arraybackslash}p{0.25\linewidth}>{\raggedright\arraybackslash}p{0.16\linewidth}>{\raggedright\arraybackslash}p{0.51\linewidth}@{}}')
# slightly smaller tables to avoid edge collisions
text=text.replace(r'\begin{longtable}', r'\begingroup\small\setlength{\tabcolsep}{5pt}\begin{longtable}')
text=text.replace(r'\end{longtable}', r'\end{longtable}\endgroup')
# keep title page clean if table starts right after title
p.write_text(text)
print('[OK] patched', p)
