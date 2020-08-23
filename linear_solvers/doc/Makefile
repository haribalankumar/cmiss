# -*- Makefile -*-
# 
# For use with GNU make.

#-----------------------------------------------------------------------------
# Makefile for building the CMISS solver manual.
#
# Original by Stuart Norris 2002-05-15
# Changes:
#-----------------------------------------------------------------------------

PS   = solvers.ps
DVI  = solvers.dvi
TEX  = solvers.tex solvers_01.tex solvers_02.tex solvers_03.tex\
       solvers_04.tex solvers_05.tex

$(DVI): $(TEX)
	latex solvers
	latex solvers

$(PS): $(DVI)
	dvips solvers -o solvers.ps

clean:
	rm -f *~ solvers.log texput.log *.bak *.aux *.dvi *.toc

