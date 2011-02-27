TEXI2DVI = R CMD texi2dvi --pdf
SWEAVE = R CMD Sweave
RM = rm -rf

all: intro2fitdistrplus.pdf clean

noclean: intro2fitdistrplus.pdf

intro2fitdistrplus.pdf : intro2fitdistrplus.Rnw
	$(SWEAVE) intro2fitdistrplus.Rnw
	$(TEXI2DVI) intro2fitdistrplus.tex

clean:
	$(RM) Rplots*.eps Rplots*.pdf Rplots*.jpg intro2fitdistrplus-*.eps intro2fitdistrplus-*.pdf intro2fitdistrplus-*.jpg *.aux *.log *.blg *.bbl *.brf *.out *.tex *.rel *.toc *~ auto

allclean: 
	$(RM) Rplots*.eps Rplots*.pdf Rplots*.jpg intro2fitdistrplus-*.eps intro2fitdistrplus-*.pdf intro2fitdistrplus-*.jpg *.aux *.log *.blg *.bbl *.brf *.out *.tex *.rel *.toc *~ auto intro2fitdistrplus.pdf
