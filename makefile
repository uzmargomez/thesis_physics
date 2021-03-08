DOCNAME=thesis_uzmar

report:
	pdflatex ${DOCNAME}.tex
	bibtex ${DOCNAME}.aux
	pdflatex ${DOCNAME}.tex
	pdflatex ${DOCNAME}.tex

clean:
	rm *.blg *.bbl *.aux *.log *.toc *.out