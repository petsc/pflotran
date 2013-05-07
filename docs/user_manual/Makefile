# Sarat: Simple makefile to make pdf
# If references are added, uncomment bibtex etc.


USER_MANUAL = \
	user_manual.tex \
	benchmark_problems.tex \
	code_structure.tex \
	example_input_files.tex \
	fortran_coding_standard.tex \
	governing_eqns.tex \
	installation.tex \
	intro.tex \
	keywords.tex \
	laur.tex \
	local.tex \
	pflotran_flow_diag.tex \
	quick_start.tex \
	references.tex \
	regression-tests.tex \
	title.tex \
	visualization.tex 

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

user_manual : $(USER_MANUAL)
	pdflatex user_manual
	pdflatex user_manual
	pdflatex user_manual

clean:
	rm -f *~ *.log *.aux *.lof *.lot *.toc *.bbl *.blg *.loa *.out
allclean:
	rm -f *~ *.log *.aux *.lof *.lot *.toc *.bbl *.blg *.loa *.out *.dvi *.ps *.pdf 
