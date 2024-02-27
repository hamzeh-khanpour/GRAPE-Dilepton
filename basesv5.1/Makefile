# Makefile for library BASES ver 5.1 
#
SHELL		= /bin/csh
AR		= /usr/ccs/bin/ar
#
# for SUN.OS
RANLIB		= /usr/ccs/bin/ranlib
# for HP/SGI
#RANLIB		= /bin/ar ts
#
FC		= f77
FFLAGS		= -O
#
CC              = /opt/SUNWspro/bin/cc
#
LIBBASES	= ../lib/libbases.a
#
SRCS	= bases.f bhinit.f bhplot.f bhrset.f bhsave.f bschck.f \
bsdate.f bsetgu.f bsetgv.f bsinfo.f bsinit.f bsdims.f bsgrid.f \
bsparm.f bsintg.f bslist.f bsordr.f bsprnt.f bsread.f bstcnv.f \
bstime_sun.f bsutim.f bswrit.f dhfill.f dhinit.f dhplot.f drn.f \
drnset.f shcler.f shfill.f shinit.f shplot.f shrset.f shupdt.f \
spchck.f sphist.f spinfo.f sprgen.f spring.f xhchck.f xhfill.f \
xhinit.f xhordr.f xhplot.f xhrnge.f xhscle.f bshbok.f sphbok.f \
bsgetw.f bsputw.f imtime_.c ixtime_.c uxdate_.c

OBJS	= bases.o bhinit.o bhplot.o bhrset.o bhsave.o bschck.o \
bsdate.o bsetgu.o bsetgv.o bsinfo.o bsinit.o bsdims.o bsgrid.o \
bsparm.o bsintg.o bslist.o bsordr.o bsprnt.o bsread.o bstcnv.o \
bstime_sun.o bsutim.o bswrit.o dhfill.o dhinit.o dhplot.o drn.o \
drnset.o shcler.o shfill.o shinit.o shplot.o shrset.o shupdt.o \
spchck.o sphist.o spinfo.o sprgen.o spring.o xhchck.o xhfill.o \
xhinit.o xhordr.o xhplot.o xhrnge.o xhscle.o bshbok.o sphbok.o \
bsgetw.o bsputw.o imtime_.o ixtime_.o uxdate_.o

all:	${LIBBASES}

${LIBBASES}:	${OBJS}
	${AR} rc ${LIBBASES} $?
	${RANLIB} ${LIBBASES}

tags:${SRCS}
	${FTAGS} ${SRCS}

clean:
	\rm -f ${OBJS} core
