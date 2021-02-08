WP=2
NUMBER=1
TITLE=main

NAME=D$(WP).$(NUMBER)-$(TITLE)
TEX=$(NAME).tex
PDF=$(NAME).pdf

default: $(TEX) references.bib
	pdflatex $(TEX)
	pdflatex $(TEX)
	bibtex $(TEX)
	pdflatex $(TEX)

create:
	@test -f $(TEX) && \
	echo  "File $(TEX) already exists" ||\
        cp main.tex $(TEX)


