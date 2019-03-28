SRC=main

all:
	pdflatex $(SRC)
	bibtex   $(SRC)
	pdflatex $(SRC)
	pdflatex $(SRC)

clean:
	rm -f *~ $(SRC).{log,aux,bbl,blg,toc}
