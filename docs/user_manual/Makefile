# Sarat: Simple makefile to make pdf
# If references are added, uncomment bibtex etc.

FILE=quick_guide
all: $(FILE)
$(FILE):
	pdflatex $(FILE)
	pdflatex $(FILE)
	# bibtex $(FILE)
	# bibtex $(FILE)
	# pdflatex $(FILE)
	# pdflatex $(FILE)
	# dvips $(FILE).dvi
	# dvipdf $(FILE).dvi
	# ps2pdf $(FILE).ps
clean:
	rm -f *~ *.log *.aux *.lof *.lot *.toc *.bbl *.blg *.loa *.out
allclean:
	rm -f *~ *.log *.aux *.lof *.lot *.toc *.bbl *.blg *.loa *.out *.dvi *.ps *.pdf 
