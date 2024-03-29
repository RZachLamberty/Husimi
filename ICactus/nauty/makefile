# makefile for nauty 2.4, including gtools
# makefile.  Generated from makefile.in by configure.

SHELL=/bin/sh
CC=g++
CFLAGS= -O 
SAFECFLAGS= -O
LDFLAGS=
LOK=1         # 0 if no 64-bit integers
TESTPROGS=dreadtest dreadtestS dreadtestS1 dreadtest4K dreadtest1 dreadtestW1 dreadtestL1 dreadtestL
EXEEXT=
BIGTEST=1         # set to 0 if less than 512MB is available

SMALL=-DMAXN=WORDSIZE
L1=-DMAXN=WORDSIZE -DWORDSIZE=64
W1=-DMAXN=WORDSIZE -DWORDSIZE=32
S1=-DMAXN=WORDSIZE -DWORDSIZE=16
S=-DWORDSIZE=16
L=-DWORDSIZE=64
FOURK=-DMAXN=4096

all : nauty gtools ;

nauty : dreadnaut ;

gtools : copyg listg labelg dretog amtog geng complg shortg showg NRswitchg \
  biplabg addedgeg deledgeg countg pickg genrang newedgeg catg genbg directg \
  multig planarg gentourng ;

rng.o: rng.c
	${CC} -c ${CFLAGS} rng.c

dreadnaut: dreadnaut.c naututil.o nauty.o nautil.o nautinv.o naugraph.o rng.o
	${CC} -o dreadnaut ${SAFECFLAGS} dreadnaut.c \
	    naututil.o nauty.o nautil.o nautinv.o naugraph.o rng.o ${LDFLAGS}
naututil.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} naututil.c
nautil.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} nautil.c
nauty.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} nauty.c
nautinv.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} nautinv.c
nautaux.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} nautaux.c
naugraph.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} naugraph.c
nausparse.o : nauty.h nausparse.h nausparse.c
	${CC} -c ${CFLAGS} nausparse.c

dreadnaut1: dreadnaut.c naututil1.o nauty1.o nautil1.o \
            nautinv1.o naugraph1.o rng.o
	${CC} -o dreadnaut1 ${SAFECFLAGS} ${SMALL} dreadnaut.c naututil1.o \
	    nauty1.o nautil1.o nautinv1.o naugraph1.o rng.o ${LDFLAGS}
naututil1.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${SMALL} -o naututil1.o naututil.c
nautil1.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${SMALL} -o nautil1.o nautil.c
nauty1.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${SMALL} -o nauty1.o nauty.c
nautinv1.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${SMALL} -o nautinv1.o nautinv.c
nautaux1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${SMALL} -o nautaux1.o nautaux.c
naugraph1.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${SMALL} -o naugraph1.o naugraph.c

dreadnaut4K: dreadnaut.c naututil4K.o nauty4K.o nautil4K.o \
             nautinv4K.o naugraph4K.o rng.o
	${CC} -o dreadnaut4K ${SAFECFLAGS} ${FOURK} dreadnaut.c naututil4K.o \
	     nauty4K.o nautil4K.o nautinv4K.o naugraph4K.o rng.o ${LDFLAGS}
naututil4K.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${FOURK} -o naututil4K.o naututil.c
nautil4K.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${FOURK} -o nautil4K.o nautil.c
nauty4K.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${FOURK} -o nauty4K.o nauty.c
nautinv4K.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${FOURK} -o nautinv4K.o nautinv.c
nautaux4K.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${FOURK} -o nautaux4K.o nautaux.c
naugraph4K.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${FOURK} -o naugraph4K.o naugraph.c

dreadnautL1: dreadnaut.c naututilL1.o nautyL1.o nautilL1.o nautinvL1.o \
	     naugraphL1.o rng.o
	${CC} -o dreadnautL1 ${SAFECFLAGS} ${L1} dreadnaut.c naututilL1.o \
	     nautyL1.o nautilL1.o nautinvL1.o naugraphL1.o rng.o ${LDFLAGS}
naututilL1.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${L1} -o naututilL1.o naututil.c
nautilL1.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${L1} -o nautilL1.o nautil.c
nautyL1.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${L1} -o nautyL1.o nauty.c
nautinvL1.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${L1} -o nautinvL1.o nautinv.c
nautauxL1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${L1} -o nautauxL1.o nautaux.c
naugraphL1.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${L1} -o naugraphL1.o naugraph.c

dreadnautW1: dreadnaut.c naututilW1.o nautyW1.o nautilW1.o nautinvW1.o \
	     naugraphW1.o rng.o
	${CC} -o dreadnautW1 ${SAFECFLAGS} ${W1} dreadnaut.c naututilW1.o \
	     nautyW1.o nautilW1.o nautinvW1.o naugraphW1.o rng.o ${LDFLAGS}
naututilW1.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${W1} -o naututilW1.o naututil.c
nautilW1.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${W1} -o nautilW1.o nautil.c
nautyW1.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${W1} -o nautyW1.o nauty.c
nautinvW1.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${W1} -o nautinvW1.o nautinv.c
nautauxW1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${W1} -o nautauxW1.o nautaux.c
naugraphW1.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${W1} -o naugraphW1.o naugraph.c
gtoolsW1.o : gtools.h gtools.c
	${CC} -c ${CFLAGS} ${W1} -o gtoolsW1.o gtools.c 


dreadnautS1: dreadnaut.c naututilS1.o nautyS1.o nautilS1.o nautinvS1.o \
	     naugraphS1.o rng.o
	${CC} -o dreadnautS1 ${SAFECFLAGS} ${S1} dreadnaut.c naututilS1.o \
	     nautyS1.o nautilS1.o nautinvS1.o naugraphS1.o rng.o ${LDFLAGS}
naututilS1.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${S1} -o naututilS1.o naututil.c
nautilS1.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${S1} -o nautilS1.o nautil.c
nautyS1.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${S1} -o nautyS1.o nauty.c
nautinvS1.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${S1} -o nautinvS1.o nautinv.c
nautauxS1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${S1} -o nautauxS1.o nautaux.c
naugraphS1.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${S1} -o naugraphS1.o naugraph.c

dreadnautL: dreadnaut.c naututilL.o nautyL.o nautilL.o nautinvL.o \
	     naugraphL.o rng.o
	${CC} -o dreadnautL ${SAFECFLAGS} ${L} dreadnaut.c naututilL.o \
	     nautyL.o nautilL.o nautinvL.o naugraphL.o rng.o ${LDFLAGS}
naututilL.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${L} -o naututilL.o naututil.c
nautilL.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${L} -o nautilL.o nautil.c
nautyL.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${L} -o nautyL.o nauty.c
nautinvL.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${L} -o nautinvL.o nautinv.c
nautauxL.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${L} -o nautauxL.o nautaux.c
naugraphL.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${L} -o naugraphL.o naugraph.c

dreadnautS: dreadnaut.c naututilS.o nautyS.o nautilS.o nautinvS.o \
	     naugraphS.o rng.o
	${CC} -o dreadnautS ${SAFECFLAGS} ${S} dreadnaut.c naututilS.o \
	     nautyS.o nautilS.o nautinvS.o naugraphS.o rng.o ${LDFLAGS}
naututilS.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${S} -o naututilS.o naututil.c
nautilS.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${S} -o nautilS.o nautil.c
nautyS.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${S} -o nautyS.o nauty.c
nautinvS.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${S} -o nautinvS.o nautinv.c
nautauxS.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${S} -o nautauxS.o nautaux.c
naugraphS.o: nauty.h naugraph.c
	${CC} -c ${CFLAGS} ${S} -o naugraphS.o naugraph.c

clean:
	rm -f *.o config.log config.cache config.status
	rm -f dreadtest${EXEEXT} dreadtestL${EXEEXT} \
		dreadtestS${EXEEXT} dreadtestS1${EXEEXT} \
	        dreadtest4K${EXEEXT} dreadtest1${EXEEXT} \
		dreadtestL1${EXEEXT} dreadtestW1${EXEEXT}

gtools.h : nauty.h naututil.h nausparse.h
	touch gtools.h

gtools.o : gtools.h gtools.c
	${CC} -c ${CFLAGS} gtools.c

gtnauty.o : gtools.h gtnauty.c
	${CC} -c ${CFLAGS} gtnauty.c

gutil1.o : gtools.h gutils.h gutil1.c
	${CC} -c ${CFLAGS} gutil1.c

gutil2.o : gtools.h gutils.h gutil2.c
	${CC} -c ${CFLAGS} gutil2.c

catg : gtools.h catg.c gtools.o
	${CC} -o catg ${CFLAGS} catg.c gtools.o ${LDFLAGS}

copyg : gtools.h copyg.c gtools.o
	${CC} -o copyg ${CFLAGS} copyg.c gtools.o ${LDFLAGS}

listg : gtools.h listg.c gtools.o nautil.o
	${CC} -o listg ${CFLAGS} listg.c gtools.o nautil.o ${LDFLAGS}

labelg : gtools.h naututil.h labelg.c gtools.o gtnauty.o nauty.o \
            nautil.o naugraph.o nautinv.o nausparse.h nausparse.o \
	    gutils.h gutil2.o
	${CC} -o labelg ${SAFECFLAGS} labelg.c gtools.o gtnauty.o \
	    nauty.o nautil.o naugraph.o nautinv.o nausparse.o \
	    gutil2.o ${LDFLAGS}

shortg : gtools.h shortg.c gtools.o gtnauty.o nauty.o nautil.o naugraph.o \
            nausparse.h nausparse.o gutils.h gutil2.o nautinv.h nautinv.o
	${CC} -o shortg ${CFLAGS} shortg.c gtools.o gtnauty.o nautinv.o \
            nauty.o nautil.o naugraph.o nausparse.o gutil2.o ${LDFLAGS}

dretog : gtools.h dretog.c gtools.o naututil.o nautil.o naugraph.o rng.o
	${CC} -o dretog ${CFLAGS} dretog.c \
		gtools.o naututil.o nautil.o rng.o naugraph.o ${LDFLAGS}

amtog : gtools.h amtog.c gtools.o
	${CC} -o amtog ${CFLAGS} amtog.c gtools.o ${LDFLAGS}

genbg : gtools.h genbg.c gtoolsW1.o nautyW1.o nautilW1.o naugraphW1.o
	${CC} -o genbg ${CFLAGS} ${W1} genbg.c gtoolsW1.o \
		nautyW1.o nautilW1.o naugraphW1.o ${LDFLAGS}

genbgL : gtools.h genbg.c gtoolsL1.o nautyL1.o nautilL1.o naugraphL1.o
	${CC} -o genbgL ${CFLAGS} ${L1} -DMAXN1=30 genbg.c gtoolsL1.o \
	    nautyL1.o nautilL1.o naugraphL1.o ${LDFLAGS}

geng : gtools.h geng.c gtools.o nauty1.o nautil1.o naugraph1.o
	${CC} -o geng ${CFLAGS} -DMAXN=32 geng.c gtools.o nauty1.o \
		nautil1.o naugraph1.o ${LDFLAGS}

geng16 : gtools.h geng.c gtools.o nauty1.o nautil1.o naugraph1.o
	${CC} -o geng16 ${CFLAGS} -DMAXN=16 geng.c gtools.o nauty1.o \
		nautil1.o naugraph1.o ${LDFLAGS}

geng24 : gtools.h geng.c gtools.o nauty1.o nautil1.o naugraph1.o
	${CC} -o geng24 ${CFLAGS} -DMAXN=24 geng.c gtools.o nauty1.o \
		nautil1.o naugraph1.o ${LDFLAGS}

gentourng : gtools.h gentourng.c gtools.o nauty1.o nautil1.o naugraph1.o
	${CC} -o gentourng ${CFLAGS} -DMAXN=24 gentourng.c gtools.o nauty1.o \
		nautil1.o naugraph1.o ${LDFLAGS}

genrang : gtools.h genrang.c gtools.o naututil.o nautil.o naugraph.o \
	    rng.o nausparse.o
	${CC} -o genrang ${CFLAGS} genrang.c nausparse.o \
		gtools.o nautil.o naututil.o naugraph.o rng.o ${LDFLAGS}

complg : gtools.h complg.c gtools.o gtnauty.o nauty.o nautil.o \
           naugraph.o nausparse.o
	${CC} -o complg ${CFLAGS} complg.c gtools.o gtnauty.o \
               nauty.o nautil.o naugraph.o nausparse.o ${LDFLAGS}

biplabg : gtools.h biplabg.c gtools.o nautil.o naugraph.o gutil1.o
	${CC} -o biplabg ${CFLAGS} biplabg.c \
		gtools.o nautil.o naugraph.o gutil1.o ${LDFLAGS}

NRswitchg : gtools.h NRswitchg.c gtools.o gtnauty.o nauty.o nautil.o \
               nausparse.o naugraph.o
	${CC} -o NRswitchg ${CFLAGS} NRswitchg.c gtools.o gtnauty.o \
		 nauty.o nautil.o naugraph.o nausparse.o ${LDFLAGS}

deledgeg : gtools.h deledgeg.c gtools.o gtnauty.o nauty.o nautil.o \
               naugraph.o nausparse.o
	${CC} -o deledgeg ${CFLAGS} deledgeg.c gtools.o gtnauty.o \
		 nauty.o nautil.o naugraph.o nausparse.o ${LDFLAGS}

addedgeg : gtools.h addedgeg.c gtools.o gtnauty.o nauty.o nautil.o \
               naugraph.o gutil1.o nausparse.o
	${CC} -o addedgeg ${CFLAGS} addedgeg.c gtools.o gtnauty.o \
		 nauty.o nautil.o naugraph.o gutil1.o nausparse.o ${LDFLAGS}

newedgeg : gtools.h newedgeg.c nausparse.o \
             gtools.o gtnauty.o nauty.o nautil.o naugraph.o
	${CC} -o newedgeg ${CFLAGS} newedgeg.c gtools.o gtnauty.o \
		 nauty.o nautil.o naugraph.o nausparse.o ${LDFLAGS}

pickg : gtools.h testg.c splay.c nausparse.o \
	     gtools.o gtnauty.o nauty.o nautil.o naugraph.o gutil1.o gutil2.o
	${CC} -o pickg ${SAFECFLAGS} testg.c gtools.o gtnauty.o gutil1.o \
	        gutil2.o nauty.o nautil.o naugraph.o nausparse.o ${LDFLAGS}

countg : gtools.h testg.c splay.c nausparse.o \
	     gtools.o gtnauty.o nauty.o nautil.o naugraph.o gutil1.o gutil2.o
	${CC} -o countg ${SAFECFLAGS} testg.c gtools.o gtnauty.o gutil1.o \
	         gutil2.o nauty.o nautil.o naugraph.o nausparse.o ${LDFLAGS}

directg : gtools.h naugroup.h directg.c \
             gtools.o nauty.o nautil.o naugraph.o naugroup.o
	${CC} -o directg ${CFLAGS} directg.c gtools.o \
                nauty.o nautil.o naugraph.o naugroup.o ${LDFLAGS}

multig : gtools.h naugroup.h multig.c \
             gtools.o nauty.o nautil.o naugraph.o naugroup.o
	${CC} -o multig ${CFLAGS} multig.c gtools.o \
                nauty.o nautil.o naugraph.o naugroup.o ${LDFLAGS}

planarg : gtools.h planarg.c gtools.o planarity.c
	${CC} -o planarg ${CFLAGS} \
		planarg.c planarity.c gtools.o ${LDFLAGS}

showg : showg.c
	${CC} -o showg ${CFLAGS} showg.c ${LDFLAGS}

readg :
	echo "*** readg is now called showg ***"
	make showg

gtoolsL1.o : gtools.h gtools.c
	${CC} -c ${CFLAGS} ${L1} -o gtoolsL1.o gtools.c

gtnautyL1.o : gtools.h gtnauty.c
	${CC} -c ${CFLAGS} ${L1} -o gtnautyL1.o gtnauty.c

sumlines : sumlines.c
	${CC} -o sumlines ${CFLAGS} sumlines.c ${LDFLAGS} -lgmp

nautyex1 : nauty.h nauty.o nautil.o naugraph.o nautyex1.c
	${CC} -o nautyex1 ${CFLAGS} nautyex1.c nauty.o nautil.o \
		naugraph.o ${LDFLAGS}

nautyex2 : nauty.h nauty.o nautil.o naugraph.o nautyex2.c
	${CC} -o nautyex2 ${CFLAGS} \
		nautyex2.c nauty.o nautil.o naugraph.o ${LDFLAGS}

naugroup.o : nauty.h naugroup.h naugroup.c
	${CC} -c ${CFLAGS} naugroup.c

nautyex3 : nauty.h nauty.o nautil.o naugraph.o naugroup.o nautyex3.c
	${CC} -o nautyex3 ${CFLAGS} \
		nautyex3.c nauty.o nautil.o naugraph.o naugroup.o ${LDFLAGS}

nautyex4 : nauty.h nausparse.h nauty.o nautil.o \
		nausparse.o naugroup.o nautyex4.c
	${CC} -o nautyex4 ${CFLAGS} \
		nautyex4.c nauty.o nautil.o nausparse.o naugroup.o ${LDFLAGS}

nautyex5 : nauty.h nausparse.h nauty.o nautil.o \
		nausparse.o naugroup.o nautyex5.c
	${CC} -o nautyex5 ${CFLAGS} \
		nautyex5.c nauty.o nautil.o nausparse.o naugroup.o ${LDFLAGS}

dreadtest: dreadtest.c naututil.o nauty.o nautil.o nautinv.o naugraph.o rng.o
	${CC} -o dreadtest ${SAFECFLAGS} dreadtest.c \
	    naututil.o nauty.o nautil.o nautinv.o naugraph.o rng.o ${LDFLAGS}
	./dreadtest <testS.dre
	./dreadtest <testW.dre
	./dreadtest <testL.dre
	./dreadtest <testM.dre
	./dreadtest <testa.dre
	test ${BIGTEST} -eq 0 || ./dreadtest < testB.dre
	rm dreadtest${EXEEXT}

dreadtestL: dreadtest.c \
		naututilL.o nautyL.o nautilL.o nautinvL.o naugraphL.o rng.o
	${CC} -o dreadtestL ${SAFECFLAGS} ${L} dreadtest.c naututilL.o \
	    nautyL.o nautilL.o nautinvL.o naugraphL.o rng.o ${LDFLAGS}
	./dreadtestL < testS.dre
	./dreadtestL < testW.dre
	./dreadtestL < testL.dre
	./dreadtestL < testM.dre
	./dreadtestL < testa.dre
	rm dreadtestL${EXEEXT}

dreadtestS: dreadtest.c \
		naututilS.o nautyS.o nautilS.o nautinvS.o naugraphS.o rng.o
	${CC} -o dreadtestS ${SAFECFLAGS} ${S} dreadtest.c naututilS.o \
	    nautyS.o nautilS.o nautinvS.o naugraphS.o rng.o ${LDFLAGS}
	./dreadtestS < testS.dre
	./dreadtestS < testW.dre
	./dreadtestS < testL.dre
	./dreadtestS < testM.dre
	./dreadtestS < testa.dre
	rm dreadtestS${EXEEXT}

dreadtest4K: dreadtest.c \
		naututil4K.o nauty4K.o nautil4K.o nautinv4K.o naugraph4K.o rng.o
	${CC} -o dreadtest4K ${SAFECFLAGS} ${FOURK} dreadtest.c naututil4K.o \
	    nauty4K.o nautil4K.o nautinv4K.o naugraph4K.o rng.o ${LDFLAGS}
	./dreadtest4K <testS.dre
	./dreadtest4K <testW.dre
	./dreadtest4K < testM.dre
	rm dreadtest4K${EXEEXT}

dreadtest1: dreadtest.c \
		naututil1.o nauty1.o nautil1.o nautinv1.o naugraph1.o rng.o
	${CC} -o dreadtest1 ${SAFECFLAGS} ${SMALL} dreadtest.c naututil1.o \
	    nauty1.o nautil1.o nautinv1.o naugraph1.o rng.o ${LDFLAGS}
	./dreadtest1 < testS.dre
	./dreadtest1 < testW.dre
	rm dreadtest1${EXEEXT}

dreadtestL1: dreadtest.c \
		naututilL1.o nautyL1.o nautilL1.o nautinvL1.o naugraphL1.o rng.o
	${CC} -o dreadtestL1 ${SAFECFLAGS} ${L1} dreadtest.c naututilL1.o \
	    nautyL1.o nautilL1.o nautinvL1.o naugraphL1.o rng.o ${LDFLAGS}
	./dreadtestL1 <testS.dre
	./dreadtestL1 <testW.dre
	./dreadtestL1 <testL.dre
	rm dreadtestL1${EXEEXT}

dreadtestW1: dreadtest.c \
		naututilW1.o nautyW1.o nautilW1.o nautinvW1.o naugraphW1.o rng.o
	${CC} -o dreadtestW1 ${SAFECFLAGS} ${W1} dreadtest.c naututilW1.o \
	    nautyW1.o nautilW1.o nautinvW1.o naugraphW1.o rng.o ${LDFLAGS}
	./dreadtestW1 <testS.dre
	./dreadtestW1 <testW.dre
	rm dreadtestW1${EXEEXT}

dreadtestS1: dreadtest.c \
		naututilS1.o nautyS1.o nautilS1.o nautinvS1.o naugraphS1.o rng.o
	${CC} -o dreadtestS1 ${SAFECFLAGS} ${S1} dreadtest.c naututilS1.o \
	    nautyS1.o nautilS1.o nautinvS1.o naugraphS1.o rng.o ${LDFLAGS}
	./dreadtestS1 <testS.dre

checks: ${TESTPROGS}

checks6: checks6.c gtools.h gtools.o
	${CC} -o checks6 ${SAFECFLAGS} checks6.c gtools.o

bigs: gtnauty.o gutil2.o nausparse.o nautil.o naututil.o gtools.o \
	naugraph.o nautaux.o nautinv.o nauty.o
	rm -f gtnautyB.o gutil2B.o nausparseB.o nautilB.o naututilB.o \
	  gtoolsB.o naugraphB.o nautauxB.o nautinvB.o nautyB.o
	cp gtnauty.o gtnautyB.o
	cp gutil2.o gutil2B.o
	cp nausparse.o nausparseB.o
	cp nautil.o nautilB.o
	cp naututil.o naututilB.o
	cp gtools.o gtoolsB.o
	cp naugraph.o naugraphB.o
	cp nautaux.o nautauxB.o
	cp nautinv.o nautinvB.o
	cp nauty.o nautyB.o

# @edit_msg@
