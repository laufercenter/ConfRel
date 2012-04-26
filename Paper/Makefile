LC = pdflatex
BC = bibtex

PAPER_ROOT = paper

TEXFILES = \
	paper.tex \
	paper.bib 

all: ${PAPER_ROOT}.pdf ${PAPER_ROOT}.view

${PAPER_ROOT}.view: ${PAPER_ROOT}.pdf 
	/usr/bin/open ${PAPER_ROOT}.pdf

${PAPER_ROOT}.pdf: ${TEXFILES}
	${LC} ${PAPER_ROOT}
	${BC} ${PAPER_ROOT}
	${LC} ${PAPER_ROOT}
	${LC} ${PAPER_ROOT}

clean:
	rm -f *.aux *.log *.bbl *.blg *.dvi *.tbx *.fgx *.fff *.lof *.lot *.ttt

