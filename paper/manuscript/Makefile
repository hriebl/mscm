manuscript.pdf: manuscript.tex
	sed -i "s/collaborator/author/g" references.bib
	pdflatex -halt-on-error $< || { touch $@ && false; }
	bibtex $(<:.tex=.aux)
	pdflatex -halt-on-error $<
	pdflatex -halt-on-error $<
