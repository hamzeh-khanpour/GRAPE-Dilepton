# Makefile for KINEMATICS library ver 1.0 
#
SHELL		= /bin/csh
AR		= /usr/ccs/bin/ar
#
# for SUN.OS
RANLIB	= /usr/ccs/bin/ranlib
# for HP/SGI
# RANLIB	= /bin/ar ts
#
FC		= f77
FFLAGS		= -O
#
LIBKINEM      = ../lib/libkinem.a
#
SRCS  =  arrclr.f arrcpy.f determ.f ktolab.f labtok.f memclc.f \
memclr.f minor.f  minvr2.f minvrt.f mmmult.f mtxclr.f mtxcpy.f \
mtxinv.f mvmult.f pboost.f rotmtx.f vvvmlt.f wtolab.f
OBJS  =  arrclr.o arrcpy.o determ.o ktolab.o labtok.o memclc.o \
memclr.o minor.o  minvr2.o minvrt.o mmmult.o mtxclr.o mtxcpy.o \
mtxinv.o mvmult.o pboost.o rotmtx.o vvvmlt.o wtolab.o

all:	${LIBKINEM}

${LIBKINEM}:	${OBJS}
	${AR} rc ${LIBKINEM} $?
	${RANLIB} ${LIBKINEM}

tags:${SRCS}
	${FTAGS} ${SRCS}

clean:
	\rm -f ${OBJS} core
