# Makefile for library CHANEL ver 2.1
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
LIBCHANEL      = ../lib/libchanel.a
#
SRCS   = sminit.f smconf.f smconv.f smcons.f smffv_new2.f \
smffs.f  smvvv.f  smsvv.f  smssv.f  smsss.f  smvvvv.f smssvv.f \
smssss.f smextv.f smintv.f sminv0.f smintf.f smextf.f smggg.f  smgggg.f \
smprpd.f pol.f    spltq.f  phaseq.f ffvmm1.f ffvmm2.f ffvmm3.f \
ffvmm4.f ffvmm5.f ffvmm6.f ffvmm7.f ffvmm8.f ffvmm9.f splt.f   \
vvv.f    vvvv.f   vvs.f    ssv.f    pola.f   ffs0.f   ffs.f    \
ffv.f    ffv0.f \
smvww.f wwva.f \
polp.f   smintp.f smintt.f smcpol.f

OBJS   = sminit.o smconf.o smconv.o smcons.o smffv_new2.o \
smffs.o  smvvv.o  smsvv.o  smssv.o  smsss.o  smvvvv.o smssvv.o \
smssss.o smextv.o smintv.o sminv0.o smintf.o smextf.o smggg.o  smgggg.o \
smprpd.o pol.o    spltq.o  phaseq.o ffvmm1.o ffvmm2.o ffvmm3.o \
ffvmm4.o ffvmm5.o ffvmm6.o ffvmm7.o ffvmm8.o ffvmm9.o splt.o   \
vvv.o    vvvv.o   vvs.o    ssv.o    pola.o   ffs0.o   ffs.o    \
ffv.o    ffv0.o \
smvww.o wwva.o \
polp.o   smintp.o smintt.o smcpol.o

all:	${LIBCHANEL}

${LIBCHANEL}:	${OBJS}
	${AR} rc ${LIBCHANEL} $?
	${RANLIB} ${LIBCHANEL}

tags:${SRCS}
	${FTAGS} ${SRCS}

clean:
	\rm -f ${OBJS} core
