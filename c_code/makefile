#$Id: makefile,v 1.1 2013/08/20 11:49:54 zjcao Exp $

include makefile.inc

.SUFFIXES: .o .c

.c.o:
	${CC} -c $< $(filein) -o $@

# projects
SEOBNRE: Panyi_elip.o Panyimain_elip.o Panyicomm_elip.o PNwaveformPRD544813.o
	$(CC) -o $@ Panyi_elip.o Panyimain_elip.o \
	Panyicomm_elip.o PNwaveformPRD544813.o $(LDLIBS) $(GSLIB)
