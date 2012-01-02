cpv.pdf: cpv.tex
	rm -f cpv.aux cpv.pdf cpv.log
	pdflatex cpv >pdflatex.err
	pdflatex cpv >>pdflatex.err

all:
	cpv.pdf
