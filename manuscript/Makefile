# Compila el artículo. Requiere latexmk y una distribución LaTeX con
# siunitx, booktabs. numbers.tex y figures/ los genera el repositorio
# (scripts/network_benchmark/10_* y 11_*); no editarlos a mano.
all: FSR_TMR.pdf

FSR_TMR.pdf: FSR_TMR.tex numbers.tex figures/*.pdf
	latexmk -pdf -interaction=nonstopmode FSR_TMR.tex   # o: tectonic FSR_TMR.tex

clean:
	latexmk -C FSR_TMR.tex

.PHONY: all clean
