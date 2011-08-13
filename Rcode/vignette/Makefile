
both: run_code
	pdflatex vignette.tex

run_code: vignette.Rnw
	echo 'Sweave("vignette.Rnw", eval=TRUE)' | R --vanilla --quiet

clean :
	rm -f vignette.aux vignette.log vignette.pdf vignette.tex vignette.toc fig/*.pdf

