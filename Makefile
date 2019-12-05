# SYSTEM     = x86-64_sles10_4.1
# use for CAEN 
SYSTEM = x86-64_linux

LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = $(CPLEXHOME)/cplex
CONCERTDIR    = $(CPLEXHOME)/concert
# use for CAEN
# CPLEXDIR      = /usr/caen/cplex-12.6/cplex
# CONCERTDIR    = /usr/caen/cplex-12.6/concert

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0
CC  = gcc -O0
JAVAC = javac 

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC -fno-strict-aliasing
JOPT  = -classpath $(CPLEXDIR)/lib/cplex.jar -O

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread
CLNFLAGS  = -lcplex -lm -lpthread
JAVA      = java  -d64 -Djava.library.path=$(CPLEXDIR)/bin/x86-64_linux -classpath $(CPLEXJARDIR):


all:
	make all_c
	make all_cpp
	make all_java

execute: all
	make execute_c
	make execute_cpp
	make execute_java
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data
EXSRCC        = $(EXDIR)/src/c
EXSRCCX       = $(EXDIR)/src/c_x
EXSRCCPP      = $(EXDIR)/src/cpp
EXSRCJAVA     = $(EXDIR)/src/java

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT)
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------


C_EX = lpex1 lpex2 lpex3 lpex4 lpex5 lpex6 lpex7 lpex8 \
       mipex1 mipex2 mipex3 mipex4 miqpex1 netex1 netex2 \
       qcpex1 qpex1 qpex2 indefqpex1 globalqpex1 globalmiqpex1 \
       steel diet fixnet foodmanu adpreex1 \
       admipex1 admipex2 admipex3 admipex4 admipex5 admipex6 admipex7 \
       populate tuneset bendersatsp socpex1 qcpdual

CX_EX = xlpex1 xlpex2 xlpex3 xlpex4 xlpex5 xlpex6 xlpex7 xlpex8 \
       xmipex1 xmipex2 xmipex3 xmipex4 xmiqpex1 xnetex1 xnetex2 \
       xqcpex1 xqpex1 xqpex2 xindefqpex1 xglobalqpex1 xglobalmiqpex1 \
       xsteel xdiet xfixnet xfoodmanu xadpreex1 \
       xadmipex1 xadmipex2 xadmipex3 xadmipex4 xadmipex5 xadmipex6 xadmipex7 \
       xpopulate xtuneset xbendersatsp xsocpex1 xqcpdual

CPP_EX = blend cutstock etsp facility fixcost1 foodmanufact \
         iloadmipex1 iloadmipex2 iloadmipex3 iloadmipex4 \
         iloadmipex5 iloadmipex6 iloindefqpex1 iloglobalqpex1 ilodiet \
         ilolpex1 ilolpex2 ilolpex3 ilolpex4 ilolpex6 ilolpex7 \
         ilomipex1 ilomipex2 ilomipex3 ilomipex4 ilomiqpex1 \
         ilogoalex1 ilogoalex2 ilogoalex3 iloqcpex1 iloqpex1 iloqpex2 \
         iloqpex3 inout1 inout3 mixblend rates transport ilosteel \
         warehouse ilopopulate ilotuneset ilobendersatsp ilosocpex1 iloqcpdual

JAVA_EX = Blend.class MixBlend.class CutStock.class Diet.class \
          Etsp.class Facility.class FixCost1.class \
          FoodManufact.class InOut1.class InOut3.class \
          Rates.class Steel.class Transport.class \
          LPex1.class LPex2.class LPex3.class LPex4.class \
          LPex6.class LPex7.class IndefQPex1.class GlobalQPex1.class \
          MIPex1.class MIPex2.class MIPex3.class MIPex4.class MIQPex1.class \
          AdMIPex1.class AdMIPex2.class AdMIPex3.class \
          AdMIPex4.class AdMIPex5.class AdMIPex6.class QCPex1.class \
          QPex1.class QPex2.class QPex3.class Goalex1.class Goalex2.class \
          TuneSet.class BendersATSP.class \
          Goalex3.class Populate.class Warehouse.class CplexServer.class \
          SocpEx1.class QCPDual.class

all_c: $(C_EX) $(CX_EX)

all_cpp: $(CPP_EX)

all_java: $(JAVA_EX)

execute_c: $(C_EX) $(CX_EX)
	 ./indefqpex1
	 ./globalqpex1 $(EXDATA)/nonconvexqp.lp g
	 ./globalmiqpex1 $(EXDATA)/nonconvexmiqp.lp
	 ./lpex1 -r
	 ./lpex2 $(EXDATA)/example.mps p
	 ./lpex3
	 ./lpex4
	 ./lpex5
	 ./lpex6
	 ./lpex7 $(EXDATA)/afiro.mps p
	 ./lpex8
	 ./mipex1
	 ./mipex2 $(EXDATA)/mexample.mps
	 ./mipex3
	 ./mipex4 $(EXDATA)/p0033.mps l
	 ./miqpex1
	 ./netex1
	 ./netex2 $(EXDATA)/infnet.net
	 ./qcpex1
	 ./qpex1
	 ./qpex2 $(EXDATA)/qpex.lp o
	 ./steel -r
	 ./steel -c
	 ./diet -r $(EXDATA)/diet.dat
	 ./fixnet
	 ./foodmanu
	 ./populate $(EXDATA)/location.lp
	 ./tuneset $(EXDATA)/p0033.mps
	 ./bendersatsp 1 $(EXDATA)/atsp.dat
	 ./socpex1
	 ./qcpdual
	 ./adpreex1
	 ./admipex1 $(EXDATA)/mexample.mps
	 ./admipex2 $(EXDATA)/p0033.mps
	 ./admipex3 $(EXDATA)/sosex3.lp
	 ./admipex4
	 ./admipex5
	 ./admipex6 $(EXDATA)/mexample.mps
	 ./admipex7 $(EXDATA)/mexample.mps
	 ./xindefqpex1
	 ./xglobalqpex1 $(EXDATA)/nonconvexqp.lp g
	 ./xglobalmiqpex1 $(EXDATA)/nonconvexmiqp.lp
	 ./xlpex1 -r
	 ./xlpex2 $(EXDATA)/example.mps p
	 ./xlpex3
	 ./xlpex4
	 ./xlpex5
	 ./xlpex6
	 ./xlpex7 $(EXDATA)/afiro.mps p
	 ./xlpex8
	 ./xmipex1
	 ./xmipex2 $(EXDATA)/mexample.mps
	 ./xmipex3
	 ./xmipex4 $(EXDATA)/p0033.mps l
	 ./xmiqpex1
	 ./xnetex1
	 ./xnetex2 $(EXDATA)/infnet.net
	 ./xqcpex1
	 ./xqpex1
	 ./xqpex2 $(EXDATA)/qpex.lp o
	 ./xsteel -r
	 ./xsteel -c
	 ./xdiet -r $(EXDATA)/diet.dat
	 ./xfixnet
	 ./xfoodmanu
	 ./xpopulate $(EXDATA)/location.lp
	 ./xtuneset $(EXDATA)/p0033.mps
	 ./xbendersatsp 1 $(EXDATA)/atsp.dat
	 ./xsocpex1
	 ./xqcpdual
	 ./xadpreex1
	 ./xadmipex1 $(EXDATA)/mexample.mps
	 ./xadmipex2 $(EXDATA)/p0033.mps
	 ./xadmipex3 $(EXDATA)/sosex3.lp
	 ./xadmipex4
	 ./xadmipex5
	 ./xadmipex6 $(EXDATA)/mexample.mps
	 ./xadmipex7 $(EXDATA)/mexample.mps

execute_cpp: $(CPP_EX)
	 ./blend
	 ./cutstock
	 ./etsp
	 ./facility
	 ./fixcost1
	 ./foodmanufact
	 ./iloadmipex1 $(EXDATA)/mexample.mps
	 ./iloadmipex2 $(EXDATA)/p0033.mps
	 ./iloadmipex3 $(EXDATA)/sosex3.lp
	 ./iloadmipex4
	 ./iloadmipex5
	 ./iloadmipex6 $(EXDATA)/mexample.mps
	 ./ilodiet
	 ./ilogoalex1 $(EXDATA)/mexample.mps
	 ./ilogoalex2
	 ./ilogoalex3 $(EXDATA)/mexample.mps
	 ./iloindefqpex1
	 ./iloglobalqpex1 $(EXDATA)/nonconvexqp.lp g
	 ./iloglobalqpex1 $(EXDATA)/nonconvexmiqp.lp g
	 ./ilolpex1 -r
	 ./ilolpex2 $(EXDATA)/example.mps p
	 ./ilolpex3
	 ./ilolpex4
	 ./ilolpex6
	 ./ilolpex7 $(EXDATA)/afiro.mps p
	 ./ilomipex1
	 ./ilomipex2 $(EXDATA)/mexample.mps
	 ./ilomipex3
	 ./ilomipex4 $(EXDATA)/p0033.mps l
	 ./ilomiqpex1
	 ./ilopopulate $(EXDATA)/location.lp
	 ./iloqcpex1
	 ./iloqpex1
	 ./iloqpex2 $(EXDATA)/qpex.lp o
	 ./iloqpex3
	 ./ilotuneset $(EXDATA)/p0033.mps
	 ./ilobendersatsp 1 $(EXDATA)/atsp.dat
	 ./ilosocpex1
	 ./iloqcpdual
	 ./inout1
	 ./inout3
	 ./mixblend
	 ./rates
	 ./ilosteel
	 ./transport 1
	 ./warehouse

execute_java: $(JAVA_EX)
	 $(JAVA) Goalex1 $(EXDATA)/mexample.mps
	 $(JAVA) Goalex2
	 $(JAVA) Goalex3 $(EXDATA)/mexample.mps
	 $(JAVA) IndefQPex1
	 $(JAVA) GlobalQPex1 $(EXDATA)/nonconvexqp.lp g
	 $(JAVA) GlobalQPex1 $(EXDATA)/nonconvexmiqp.lp g
	 $(JAVA) LPex1 -r
	 $(JAVA) LPex2 $(EXDATA)/example.mps p
	 $(JAVA) LPex3
	 $(JAVA) LPex4
	 $(JAVA) LPex6
	 $(JAVA) LPex7 $(EXDATA)/afiro.mps p
	 $(JAVA) MIPex1
	 $(JAVA) MIPex2 $(EXDATA)/mexample.mps
	 $(JAVA) MIPex3
	 $(JAVA) MIPex4 $(EXDATA)/p0033.mps l
	 $(JAVA) MIQPex1
	 $(JAVA) QCPex1
	 $(JAVA) QPex1
	 $(JAVA) QPex2 $(EXDATA)/qpex.lp o
	 $(JAVA) QPex3 
	 $(JAVA) Blend
	 $(JAVA) CplexServer
	 $(JAVA) CutStock
	 $(JAVA) Diet
	 $(JAVA) Etsp
	 $(JAVA) Facility
	 $(JAVA) FixCost1
	 $(JAVA) FoodManufact
	 $(JAVA) InOut1
	 $(JAVA) InOut3
	 $(JAVA) MixBlend
	 $(JAVA) Populate $(EXDATA)/location.lp
	 $(JAVA) Rates
	 $(JAVA) Steel
	 $(JAVA) Transport 1
	 $(JAVA) TuneSet $(EXDATA)/p0033.mps
	 $(JAVA) BendersATSP 1 $(EXDATA)/atsp.dat
	 $(JAVA) SocpEx1
	 $(JAVA) QCPDual
	 $(JAVA) Warehouse
	 $(JAVA) AdMIPex1 $(EXDATA)/mexample.mps
	 $(JAVA) AdMIPex2 $(EXDATA)/p0033.mps
	 $(JAVA) AdMIPex3 $(EXDATA)/sosex3.lp
	 $(JAVA) AdMIPex4
	 $(JAVA) AdMIPex5
	 $(JAVA) AdMIPex6 $(EXDATA)/mexample.mps

# ------------------------------------------------------------

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CX_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
#
# The examples
#
indefqpex1: indefqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o indefqpex1 indefqpex1.o $(CLNFLAGS)
indefqpex1.o: $(EXSRCC)/indefqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/indefqpex1.c -o indefqpex1.o

globalqpex1: globalqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o globalqpex1 globalqpex1.o $(CLNFLAGS)
globalqpex1.o: $(EXSRCC)/globalqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/globalqpex1.c -o globalqpex1.o

globalmiqpex1: globalmiqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o globalmiqpex1 globalmiqpex1.o $(CLNFLAGS)
globalmiqpex1.o: $(EXSRCC)/globalmiqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/globalmiqpex1.c -o globalmiqpex1.o

lpex1: lpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex1 lpex1.o $(CLNFLAGS)
lpex1.o: $(EXSRCC)/lpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex1.c -o lpex1.o

lpex2: lpex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex2 lpex2.o $(CLNFLAGS)
lpex2.o: $(EXSRCC)/lpex2.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex2.c -o lpex2.o

lpex3: lpex3.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex3 lpex3.o $(CLNFLAGS)
lpex3.o: $(EXSRCC)/lpex3.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex3.c -o lpex3.o

lpex4: lpex4.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex4 lpex4.o $(CLNFLAGS)
lpex4.o: $(EXSRCC)/lpex4.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex4.c -o lpex4.o

lpex5: lpex5.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex5 lpex5.o $(CLNFLAGS)
lpex5.o: $(EXSRCC)/lpex5.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex5.c -o lpex5.o

lpex6: lpex6.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex6 lpex6.o $(CLNFLAGS)
lpex6.o: $(EXSRCC)/lpex6.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex6.c -o lpex6.o

lpex7: lpex7.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex7 lpex7.o $(CLNFLAGS)
lpex7.o: $(EXSRCC)/lpex7.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex7.c -o lpex7.o

lpex8: lpex8.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o lpex8 lpex8.o $(CLNFLAGS)
lpex8.o: $(EXSRCC)/lpex8.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/lpex8.c -o lpex8.o

mipex1: mipex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o mipex1 mipex1.o $(CLNFLAGS)
mipex1.o: $(EXSRCC)/mipex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/mipex1.c -o mipex1.o

mipex2: mipex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o mipex2 mipex2.o $(CLNFLAGS)
mipex2.o: $(EXSRCC)/mipex2.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/mipex2.c -o mipex2.o

mipex3: mipex3.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o mipex3 mipex3.o $(CLNFLAGS)
mipex3.o: $(EXSRCC)/mipex3.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/mipex3.c -o mipex3.o

mipex4: mipex4.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o mipex4 mipex4.o $(CLNFLAGS)
mipex4.o: $(EXSRCC)/mipex4.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/mipex4.c -o mipex4.o

miqpex1: miqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o miqpex1 miqpex1.o $(CLNFLAGS)
miqpex1.o: $(EXSRCC)/miqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/miqpex1.c -o miqpex1.o

netex1: netex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o netex1 netex1.o $(CLNFLAGS)
netex1.o: $(EXSRCC)/netex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/netex1.c -o netex1.o

netex2: netex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o netex2 netex2.o $(CLNFLAGS)
netex2.o: $(EXSRCC)/netex2.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/netex2.c -o netex2.o

qcpex1: qcpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o qcpex1 qcpex1.o $(CLNFLAGS)
qcpex1.o: $(EXSRCC)/qcpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/qcpex1.c -o qcpex1.o

qpex1: qpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o qpex1 qpex1.o $(CLNFLAGS)
qpex1.o: $(EXSRCC)/qpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/qpex1.c -o qpex1.o

qpex2: qpex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o qpex2 qpex2.o $(CLNFLAGS)
qpex2.o: $(EXSRCC)/qpex2.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/qpex2.c -o qpex2.o

steel: steel.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o steel steel.o $(CLNFLAGS)
steel.o: $(EXSRCC)/steel.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/steel.c -o steel.o

diet: diet.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o diet diet.o $(CLNFLAGS)
diet.o: $(EXSRCC)/diet.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/diet.c -o diet.o

fixnet: fixnet.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o fixnet fixnet.o $(CLNFLAGS)
fixnet.o: $(EXSRCC)/fixnet.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/fixnet.c -o fixnet.o

foodmanu: foodmanu.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o foodmanu foodmanu.o $(CLNFLAGS)
foodmanu.o: $(EXSRCC)/foodmanu.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/foodmanu.c -o foodmanu.o

adpreex1: adpreex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o adpreex1 adpreex1.o $(CLNFLAGS)
adpreex1.o: $(EXSRCC)/adpreex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/adpreex1.c -o adpreex1.o

admipex1: admipex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex1 admipex1.o $(CLNFLAGS)
admipex1.o: $(EXSRCC)/admipex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex1.c -o admipex1.o

admipex2: admipex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex2 admipex2.o $(CLNFLAGS)
admipex2.o: $(EXSRCC)/admipex2.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex2.c -o admipex2.o

admipex3: admipex3.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex3 admipex3.o $(CLNFLAGS)
admipex3.o: $(EXSRCC)/admipex3.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex3.c -o admipex3.o

admipex4: admipex4.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex4 admipex4.o $(CLNFLAGS)
admipex4.o: $(EXSRCC)/admipex4.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex4.c -o admipex4.o

admipex5: admipex5.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex5 admipex5.o $(CLNFLAGS)
admipex5.o: $(EXSRCC)/admipex5.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex5.c -o admipex5.o

admipex6: admipex6.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex6 admipex6.o $(CLNFLAGS)
admipex6.o: $(EXSRCC)/admipex6.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex6.c -o admipex6.o

admipex7: admipex7.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o admipex7 admipex7.o $(CLNFLAGS)
admipex7.o: $(EXSRCC)/admipex7.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/admipex7.c -o admipex7.o

populate: populate.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o populate populate.o $(CLNFLAGS)
populate.o: $(EXSRCC)/populate.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/populate.c -o populate.o

tuneset: tuneset.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o tuneset tuneset.o $(CLNFLAGS)
tuneset.o: $(EXSRCC)/tuneset.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/tuneset.c -o tuneset.o

bendersatsp: bendersatsp.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o bendersatsp bendersatsp.o $(CLNFLAGS)
bendersatsp.o: $(EXSRCC)/bendersatsp.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/bendersatsp.c -o bendersatsp.o

socpex1: socpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o socpex1 socpex1.o $(CLNFLAGS)
socpex1.o: $(EXSRCC)/socpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/socpex1.c -o socpex1.o

qcpdual: qcpdual.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o qcpdual qcpdual.o $(CLNFLAGS)
qcpdual.o: $(EXSRCC)/qcpdual.c
	$(CC) -c $(CFLAGS) $(EXSRCC)/qcpdual.c -o qcpdual.o

xindefqpex1: xindefqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xindefqpex1 xindefqpex1.o $(CLNFLAGS)
xindefqpex1.o: $(EXSRCCX)/xindefqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xindefqpex1.c -o xindefqpex1.o

xglobalqpex1: xglobalqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xglobalqpex1 xglobalqpex1.o $(CLNFLAGS)
xglobalqpex1.o: $(EXSRCCX)/xglobalqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xglobalqpex1.c -o xglobalqpex1.o

xglobalmiqpex1: xglobalmiqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xglobalmiqpex1 xglobalmiqpex1.o $(CLNFLAGS)
xglobalmiqpex1.o: $(EXSRCCX)/xglobalmiqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xglobalmiqpex1.c -o xglobalmiqpex1.o

xlpex1: xlpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex1 xlpex1.o $(CLNFLAGS)
xlpex1.o: $(EXSRCCX)/xlpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex1.c -o xlpex1.o

xlpex2: xlpex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex2 xlpex2.o $(CLNFLAGS)
xlpex2.o: $(EXSRCCX)/xlpex2.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex2.c -o xlpex2.o

xlpex3: xlpex3.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex3 xlpex3.o $(CLNFLAGS)
xlpex3.o: $(EXSRCCX)/xlpex3.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex3.c -o xlpex3.o

xlpex4: xlpex4.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex4 xlpex4.o $(CLNFLAGS)
xlpex4.o: $(EXSRCCX)/xlpex4.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex4.c -o xlpex4.o

xlpex5: xlpex5.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex5 xlpex5.o $(CLNFLAGS)
xlpex5.o: $(EXSRCCX)/xlpex5.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex5.c -o xlpex5.o

xlpex6: xlpex6.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex6 xlpex6.o $(CLNFLAGS)
xlpex6.o: $(EXSRCCX)/xlpex6.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex6.c -o xlpex6.o

xlpex7: xlpex7.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex7 xlpex7.o $(CLNFLAGS)
xlpex7.o: $(EXSRCCX)/xlpex7.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex7.c -o xlpex7.o

xlpex8: xlpex8.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xlpex8 xlpex8.o $(CLNFLAGS)
xlpex8.o: $(EXSRCCX)/xlpex8.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xlpex8.c -o xlpex8.o

xmipex1: xmipex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xmipex1 xmipex1.o $(CLNFLAGS)
xmipex1.o: $(EXSRCCX)/xmipex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xmipex1.c -o xmipex1.o

xmipex2: xmipex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xmipex2 xmipex2.o $(CLNFLAGS)
xmipex2.o: $(EXSRCCX)/xmipex2.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xmipex2.c -o xmipex2.o

xmipex3: xmipex3.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xmipex3 xmipex3.o $(CLNFLAGS)
xmipex3.o: $(EXSRCCX)/xmipex3.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xmipex3.c -o xmipex3.o

xmipex4: xmipex4.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xmipex4 xmipex4.o $(CLNFLAGS)
xmipex4.o: $(EXSRCCX)/xmipex4.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xmipex4.c -o xmipex4.o

xmiqpex1: xmiqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xmiqpex1 xmiqpex1.o $(CLNFLAGS)
xmiqpex1.o: $(EXSRCCX)/xmiqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xmiqpex1.c -o xmiqpex1.o

xnetex1: xnetex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xnetex1 xnetex1.o $(CLNFLAGS)
xnetex1.o: $(EXSRCCX)/xnetex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xnetex1.c -o xnetex1.o

xnetex2: xnetex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xnetex2 xnetex2.o $(CLNFLAGS)
xnetex2.o: $(EXSRCCX)/xnetex2.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xnetex2.c -o xnetex2.o

xqcpex1: xqcpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xqcpex1 xqcpex1.o $(CLNFLAGS)
xqcpex1.o: $(EXSRCCX)/xqcpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xqcpex1.c -o xqcpex1.o

xqpex1: xqpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xqpex1 xqpex1.o $(CLNFLAGS)
xqpex1.o: $(EXSRCCX)/xqpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xqpex1.c -o xqpex1.o

xqpex2: xqpex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xqpex2 xqpex2.o $(CLNFLAGS)
xqpex2.o: $(EXSRCCX)/xqpex2.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xqpex2.c -o xqpex2.o

xsteel: xsteel.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xsteel xsteel.o $(CLNFLAGS)
xsteel.o: $(EXSRCCX)/xsteel.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xsteel.c -o xsteel.o

xdiet: xdiet.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xdiet xdiet.o $(CLNFLAGS)
xdiet.o: $(EXSRCCX)/xdiet.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xdiet.c -o xdiet.o

xfixnet: xfixnet.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xfixnet xfixnet.o $(CLNFLAGS)
xfixnet.o: $(EXSRCCX)/xfixnet.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xfixnet.c -o xfixnet.o

xfoodmanu: xfoodmanu.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xfoodmanu xfoodmanu.o $(CLNFLAGS)
xfoodmanu.o: $(EXSRCCX)/xfoodmanu.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xfoodmanu.c -o xfoodmanu.o

xadpreex1: xadpreex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadpreex1 xadpreex1.o $(CLNFLAGS)
xadpreex1.o: $(EXSRCCX)/xadpreex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadpreex1.c -o xadpreex1.o

xadmipex1: xadmipex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex1 xadmipex1.o $(CLNFLAGS)
xadmipex1.o: $(EXSRCCX)/xadmipex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex1.c -o xadmipex1.o

xadmipex2: xadmipex2.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex2 xadmipex2.o $(CLNFLAGS)
xadmipex2.o: $(EXSRCCX)/xadmipex2.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex2.c -o xadmipex2.o

xadmipex3: xadmipex3.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex3 xadmipex3.o $(CLNFLAGS)
xadmipex3.o: $(EXSRCCX)/xadmipex3.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex3.c -o xadmipex3.o

xadmipex4: xadmipex4.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex4 xadmipex4.o $(CLNFLAGS)
xadmipex4.o: $(EXSRCCX)/xadmipex4.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex4.c -o xadmipex4.o

xadmipex5: xadmipex5.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex5 xadmipex5.o $(CLNFLAGS)
xadmipex5.o: $(EXSRCCX)/xadmipex5.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex5.c -o xadmipex5.o

xadmipex6: xadmipex6.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex6 xadmipex6.o $(CLNFLAGS)
xadmipex6.o: $(EXSRCCX)/xadmipex6.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex6.c -o xadmipex6.o

xadmipex7: xadmipex7.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xadmipex7 xadmipex7.o $(CLNFLAGS)
xadmipex7.o: $(EXSRCCX)/xadmipex7.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xadmipex7.c -o xadmipex7.o

xpopulate: xpopulate.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xpopulate xpopulate.o $(CLNFLAGS)
xpopulate.o: $(EXSRCCX)/xpopulate.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xpopulate.c -o xpopulate.o

xtuneset: xtuneset.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xtuneset xtuneset.o $(CLNFLAGS)
xtuneset.o: $(EXSRCCX)/xtuneset.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xtuneset.c -o xtuneset.o

xbendersatsp: xbendersatsp.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xbendersatsp xbendersatsp.o $(CLNFLAGS)
xbendersatsp.o: $(EXSRCCX)/xbendersatsp.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xbendersatsp.c -o xbendersatsp.o

xsocpex1: xsocpex1.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xsocpex1 xsocpex1.o $(CLNFLAGS)
xsocpex1.o: $(EXSRCCX)/xsocpex1.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xsocpex1.c -o xsocpex1.o

xqcpdual: xqcpdual.o
	$(CC) $(CFLAGS) $(CLNDIRS) -o xqcpdual xqcpdual.o $(CLNFLAGS)
xqcpdual.o: $(EXSRCCX)/xqcpdual.c
	$(CC) -c $(CFLAGS) $(EXSRCCX)/xqcpdual.c -o xqcpdual.o

blend: blend.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o blend blend.o $(CCLNFLAGS)
blend.o: $(EXSRCCPP)/blend.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/blend.cpp -o blend.o

cutstock: cutstock.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o cutstock cutstock.o $(CCLNFLAGS)
cutstock.o: $(EXSRCCPP)/cutstock.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/cutstock.cpp -o cutstock.o

etsp: etsp.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o etsp etsp.o $(CCLNFLAGS)
etsp.o: $(EXSRCCPP)/etsp.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/etsp.cpp -o etsp.o

facility: facility.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o facility facility.o $(CCLNFLAGS)
facility.o: $(EXSRCCPP)/facility.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/facility.cpp -o facility.o

fixcost1: fixcost1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o fixcost1 fixcost1.o $(CCLNFLAGS)
fixcost1.o: $(EXSRCCPP)/fixcost1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/fixcost1.cpp -o fixcost1.o

foodmanufact: foodmanufact.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o foodmanufact foodmanufact.o $(CCLNFLAGS)
foodmanufact.o: $(EXSRCCPP)/foodmanufact.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/foodmanufact.cpp -o foodmanufact.o

iloadmipex1: iloadmipex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloadmipex1 iloadmipex1.o $(CCLNFLAGS)
iloadmipex1.o: $(EXSRCCPP)/iloadmipex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloadmipex1.cpp -o iloadmipex1.o

iloadmipex2: iloadmipex2.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloadmipex2 iloadmipex2.o $(CCLNFLAGS)
iloadmipex2.o: $(EXSRCCPP)/iloadmipex2.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloadmipex2.cpp -o iloadmipex2.o

iloadmipex3: iloadmipex3.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloadmipex3 iloadmipex3.o $(CCLNFLAGS)
iloadmipex3.o: $(EXSRCCPP)/iloadmipex3.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloadmipex3.cpp -o iloadmipex3.o

iloadmipex4: iloadmipex4.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloadmipex4 iloadmipex4.o $(CCLNFLAGS)
iloadmipex4.o: $(EXSRCCPP)/iloadmipex4.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloadmipex4.cpp -o iloadmipex4.o

iloadmipex5: iloadmipex5.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloadmipex5 iloadmipex5.o $(CCLNFLAGS)
iloadmipex5.o: $(EXSRCCPP)/iloadmipex5.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloadmipex5.cpp -o iloadmipex5.o

iloadmipex6: iloadmipex6.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloadmipex6 iloadmipex6.o $(CCLNFLAGS)
iloadmipex6.o: $(EXSRCCPP)/iloadmipex6.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloadmipex6.cpp -o iloadmipex6.o

ilodiet: ilodiet.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilodiet ilodiet.o $(CCLNFLAGS)
ilodiet.o: $(EXSRCCPP)/ilodiet.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilodiet.cpp -o ilodiet.o

iloindefqpex1: iloindefqpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloindefqpex1 iloindefqpex1.o $(CCLNFLAGS)
iloindefqpex1.o: $(EXSRCCPP)/iloindefqpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloindefqpex1.cpp -o iloindefqpex1.o

iloglobalqpex1: iloglobalqpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloglobalqpex1 iloglobalqpex1.o $(CCLNFLAGS)
iloglobalqpex1.o: $(EXSRCCPP)/iloglobalqpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloglobalqpex1.cpp -o iloglobalqpex1.o

ilolpex1: ilolpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilolpex1 ilolpex1.o $(CCLNFLAGS)
ilolpex1.o: $(EXSRCCPP)/ilolpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilolpex1.cpp -o ilolpex1.o

ilolpex2: ilolpex2.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilolpex2 ilolpex2.o $(CCLNFLAGS)
ilolpex2.o: $(EXSRCCPP)/ilolpex2.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilolpex2.cpp -o ilolpex2.o

ilolpex3: ilolpex3.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilolpex3 ilolpex3.o $(CCLNFLAGS)
ilolpex3.o: $(EXSRCCPP)/ilolpex3.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilolpex3.cpp -o ilolpex3.o

ilolpex4: ilolpex4.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilolpex4 ilolpex4.o $(CCLNFLAGS)
ilolpex4.o: $(EXSRCCPP)/ilolpex4.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilolpex4.cpp -o ilolpex4.o

ilolpex6: ilolpex6.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilolpex6 ilolpex6.o $(CCLNFLAGS)
ilolpex6.o: $(EXSRCCPP)/ilolpex6.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilolpex6.cpp -o ilolpex6.o

ilolpex7: ilolpex7.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilolpex7 ilolpex7.o $(CCLNFLAGS)
ilolpex7.o: $(EXSRCCPP)/ilolpex7.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilolpex7.cpp -o ilolpex7.o

ilomipex1: ilomipex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilomipex1 ilomipex1.o $(CCLNFLAGS)
ilomipex1.o: $(EXSRCCPP)/ilomipex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilomipex1.cpp -o ilomipex1.o

ilomipex2: ilomipex2.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilomipex2 ilomipex2.o $(CCLNFLAGS)
ilomipex2.o: $(EXSRCCPP)/ilomipex2.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilomipex2.cpp -o ilomipex2.o

ilomipex3: ilomipex3.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilomipex3 ilomipex3.o $(CCLNFLAGS)
ilomipex3.o: $(EXSRCCPP)/ilomipex3.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilomipex3.cpp -o ilomipex3.o

ilomipex4: ilomipex4.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilomipex4 ilomipex4.o $(CCLNFLAGS)
ilomipex4.o: $(EXSRCCPP)/ilomipex4.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilomipex4.cpp -o ilomipex4.o

inout1: inout1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o inout1 inout1.o $(CCLNFLAGS)
inout1.o: $(EXSRCCPP)/inout1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/inout1.cpp -o inout1.o

inout3: inout3.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o inout3 inout3.o $(CCLNFLAGS)
inout3.o: $(EXSRCCPP)/inout3.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/inout3.cpp -o inout3.o

ilomiqpex1: ilomiqpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilomiqpex1 ilomiqpex1.o $(CCLNFLAGS)
ilomiqpex1.o: $(EXSRCCPP)/ilomiqpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilomiqpex1.cpp -o ilomiqpex1.o

ilogoalex1: ilogoalex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilogoalex1 ilogoalex1.o $(CCLNFLAGS)
ilogoalex1.o: $(EXSRCCPP)/ilogoalex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilogoalex1.cpp -o ilogoalex1.o

ilogoalex2: ilogoalex2.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilogoalex2 ilogoalex2.o $(CCLNFLAGS)
ilogoalex2.o: $(EXSRCCPP)/ilogoalex2.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilogoalex2.cpp -o ilogoalex2.o

ilogoalex3: ilogoalex3.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilogoalex3 ilogoalex3.o $(CCLNFLAGS)
ilogoalex3.o: $(EXSRCCPP)/ilogoalex3.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilogoalex3.cpp -o ilogoalex3.o

iloqcpex1: iloqcpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloqcpex1 iloqcpex1.o $(CCLNFLAGS)
iloqcpex1.o: $(EXSRCCPP)/iloqcpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloqcpex1.cpp -o iloqcpex1.o

iloqpex1: iloqpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloqpex1 iloqpex1.o $(CCLNFLAGS)
iloqpex1.o: $(EXSRCCPP)/iloqpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloqpex1.cpp -o iloqpex1.o

iloqpex2: iloqpex2.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloqpex2 iloqpex2.o $(CCLNFLAGS)
iloqpex2.o: $(EXSRCCPP)/iloqpex2.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloqpex2.cpp -o iloqpex2.o

iloqpex3: iloqpex3.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloqpex3 iloqpex3.o $(CCLNFLAGS)
iloqpex3.o: $(EXSRCCPP)/iloqpex3.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloqpex3.cpp -o iloqpex3.o

mixblend: mixblend.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o mixblend mixblend.o $(CCLNFLAGS)
mixblend.o: $(EXSRCCPP)/mixblend.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/mixblend.cpp -o mixblend.o

rates: rates.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o rates rates.o $(CCLNFLAGS)
rates.o: $(EXSRCCPP)/rates.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/rates.cpp -o rates.o

transport: transport.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o transport transport.o $(CCLNFLAGS)
transport.o: $(EXSRCCPP)/transport.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/transport.cpp -o transport.o

warehouse: warehouse.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o warehouse warehouse.o $(CCLNFLAGS)
warehouse.o: $(EXSRCCPP)/warehouse.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/warehouse.cpp -o warehouse.o

ilopopulate: ilopopulate.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilopopulate ilopopulate.o $(CCLNFLAGS)
ilopopulate.o: $(EXSRCCPP)/ilopopulate.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilopopulate.cpp -o ilopopulate.o

ilotuneset: ilotuneset.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilotuneset ilotuneset.o $(CCLNFLAGS)
ilotuneset.o: $(EXSRCCPP)/ilotuneset.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilotuneset.cpp -o ilotuneset.o

ilobendersatsp: ilobendersatsp.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilobendersatsp ilobendersatsp.o $(CCLNFLAGS)
ilobendersatsp.o: $(EXSRCCPP)/ilobendersatsp.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilobendersatsp.cpp -o ilobendersatsp.o

ilosocpex1: ilosocpex1.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilosocpex1 ilosocpex1.o $(CCLNFLAGS)
ilosocpex1.o: $(EXSRCCPP)/ilosocpex1.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilosocpex1.cpp -o ilosocpex1.o

iloqcpdual: iloqcpdual.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o iloqcpdual iloqcpdual.o $(CCLNFLAGS)
iloqcpdual.o: $(EXSRCCPP)/iloqcpdual.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/iloqcpdual.cpp -o iloqcpdual.o

ilosteel: ilosteel.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o ilosteel ilosteel.o $(CCLNFLAGS)
ilosteel.o: $(EXSRCCPP)/ilosteel.cpp
	$(CCC) -c $(CCFLAGS) $(EXSRCCPP)/ilosteel.cpp -o ilosteel.o

IndefQPex1.class: $(EXSRCJAVA)/IndefQPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/IndefQPex1.java 

GlobalQPex1.class: $(EXSRCJAVA)/GlobalQPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/GlobalQPex1.java 

LPex1.class: $(EXSRCJAVA)/LPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/LPex1.java 

LPex2.class: $(EXSRCJAVA)/LPex2.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/LPex2.java 

LPex3.class: $(EXSRCJAVA)/LPex3.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/LPex3.java 

LPex4.class: $(EXSRCJAVA)/LPex4.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/LPex4.java 

LPex6.class: $(EXSRCJAVA)/LPex6.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/LPex6.java 

LPex7.class: $(EXSRCJAVA)/LPex7.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/LPex7.java 

MIPex1.class: $(EXSRCJAVA)/MIPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/MIPex1.java 

MIPex2.class: $(EXSRCJAVA)/MIPex2.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/MIPex2.java 

MIPex3.class: $(EXSRCJAVA)/MIPex3.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/MIPex3.java 

MIPex4.class: $(EXSRCJAVA)/MIPex4.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/MIPex4.java 

MIQPex1.class: $(EXSRCJAVA)/MIQPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/MIQPex1.java 

Goalex1.class: $(EXSRCJAVA)/Goalex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Goalex1.java 

Goalex2.class: $(EXSRCJAVA)/Goalex2.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Goalex2.java 

Goalex3.class: $(EXSRCJAVA)/Goalex3.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Goalex3.java 

AdMIPex1.class: $(EXSRCJAVA)/AdMIPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/AdMIPex1.java

AdMIPex2.class: $(EXSRCJAVA)/AdMIPex2.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/AdMIPex2.java

AdMIPex3.class: $(EXSRCJAVA)/AdMIPex3.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/AdMIPex3.java

AdMIPex4.class: $(EXSRCJAVA)/AdMIPex4.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/AdMIPex4.java

AdMIPex5.class: $(EXSRCJAVA)/AdMIPex5.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/AdMIPex5.java

AdMIPex6.class: $(EXSRCJAVA)/AdMIPex6.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/AdMIPex6.java

QCPex1.class: $(EXSRCJAVA)/QCPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/QCPex1.java 

QPex1.class: $(EXSRCJAVA)/QPex1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/QPex1.java 

QPex2.class: $(EXSRCJAVA)/QPex2.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/QPex2.java 

QPex3.class: $(EXSRCJAVA)/QPex3.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/QPex3.java 

Diet.class: $(EXSRCJAVA)/Diet.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Diet.java \
	                         $(EXSRCJAVA)/InputDataReader.java 

Etsp.class: $(EXSRCJAVA)/Etsp.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Etsp.java \
	                         $(EXSRCJAVA)/InputDataReader.java 

Blend.class: $(EXSRCJAVA)/Blend.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Blend.java

MixBlend.class: $(EXSRCJAVA)/MixBlend.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/MixBlend.java

CplexServer.class: $(EXSRCJAVA)/CplexServer.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/CplexServer.java

CutStock.class: $(EXSRCJAVA)/CutStock.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/InputDataReader.java \
                                 $(EXSRCJAVA)/CutStock.java

Facility.class: $(EXSRCJAVA)/Facility.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/InputDataReader.java \
                                 $(EXSRCJAVA)/Facility.java

FixCost1.class: $(EXSRCJAVA)/FixCost1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/FixCost1.java

FoodManufact.class: $(EXSRCJAVA)/FoodManufact.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/FoodManufact.java

InOut1.class: $(EXSRCJAVA)/InOut1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/InOut1.java

InOut3.class: $(EXSRCJAVA)/InOut3.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/InOut3.java

Populate.class: $(EXSRCJAVA)/Populate.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Populate.java

TuneSet.class: $(EXSRCJAVA)/TuneSet.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/TuneSet.java

BendersATSP.class: $(EXSRCJAVA)/BendersATSP.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/BendersATSP.java \
	                         $(EXSRCJAVA)/InputDataReader.java 

SocpEx1.class: $(EXSRCJAVA)/SocpEx1.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/SocpEx1.java

QCPDual.class: $(EXSRCJAVA)/QCPDual.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/QCPDual.java

Rates.class: $(EXSRCJAVA)/Rates.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/InputDataReader.java \
                                 $(EXSRCJAVA)/Rates.java

Steel.class: $(EXSRCJAVA)/Steel.java $(EXSRCJAVA)/InputDataReader.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/InputDataReader.java \
                                 $(EXSRCJAVA)/Steel.java

Transport.class: $(EXSRCJAVA)/Transport.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Transport.java

Warehouse.class: $(EXSRCJAVA)/Warehouse.java
	$(JAVAC) $(JCFLAGS) -d . $(EXSRCJAVA)/Warehouse.java

# Local Variables:
# mode: makefile
# End:

#------------------------------------------------------------
#
# My test
#
#------------------------------------------------------------
# test_game:
# 	g++ -g -c CN_Constants.cpp -std=c++11
# 	g++ -g -c CN_Node.cpp -std=c++11
# 	g++ -g -c CN_Edge.cpp -std=c++11
# 	# g++ -g -c CN_Graph.cpp -std=c++11
# 	# g++ -g -c CN_WidgetGraph.cpp -std=c++11
# 	# g++ -g -c CN_CplexConverter.cpp -std=c++11
# 	# $(CCC) -c $(CCFLAGS) $(CCLNDIRS) CN_Solver.cpp -std=c++11
# 	# $(CCC) -c $(CCFLAGS) CN_CreditNet.cpp -std=c++11
# 	# $(CCC) $(CCFLAGS) $(CCLNDIRS) -o game test_game.cpp *.o $(CCLNFLAGS) -std=c++11	
# 	# ./game $(folder) $(samps)
	
test_game:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	g++ -c $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CPPFLAGS) -o game test_game.cpp *.o $(CCLNFLAGS) -std=c++11	
	./game $(folder) $(samps)

test_IR:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	g++ -c $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CPPFLAGS) -o test_IR test.cpp *.o $(CCLNFLAGS) -std=c++11	
	./test_IR json 1

simuMech_IRSRC:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o simu2 test.cpp *.o $(CCLNFLAGS) -std=c++11
	# time ./simu 1 4 1 MIN_SRC_COST> out.1AMT.4IR.1CAP_src
	# time ./simu 1 7 1 MIN_SRC_COST> out.1AMT.7IR.1CAP_src
	time ./simu2 1 4 1 MIN_CREDIT_COST> out.1AMT.4IR.1CAP_credit_high
	# time ./simu1 1 4 1 MIN_CREDIT_SRC > out.1AMT.4IR.1CAP_credit_srcmin_high
	# time ./simu1 1 7 1 MIN_CREDIT_SRC> out.1AMT.7IR.1CAP_credit_srcmin
	# time ./simu 1 1 1 MIN_SUMIR_COST> out.1AMT.1IR.1CAP_ircost
	# time ./simu 1 4 1 MIN_SUMIR_COST> out.1AMT.4IR.1CAP_ircost_irmin
	# time ./simu 1 7 1 MIN_SUMIR_COST> out.1AMT.7IR.1CAP_ircost_irmin
	time ./simu2 1 4 1 MIN_DEGREE_COST> out.1AMT.4IR.1CAP_degree_high
	# time ./simu1 1 4 1 MIN_DEGREE_SRC> out.1AMT.4IR.1CAP_degree_srcmin_high
	# time ./simu1 1 7 1 MIN_DEGREE_SRC> out.1AMT.7IR.1CAP_degree_srcmin
	# time ./simu2 10 4 10 > out.10AMT.4IR.10CAP_iravg1
	# time ./simu2 8 4 10 > out.8AMT.4IR.10CAP_iravg1
	# time ./simu2 15 4 10 > out.15AMT.4IR.10CAP_iravg1
	# time ./simu2 20 4 10 > out.20AMT.4IR.10CAP_iravg1


testC:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o testCollateral test_collateral.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json
	# ./testCollateral json 1 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out

simC:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateral collateralSim.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json
	# ./simCollateral json 30 1 > 1IR
	# ./simCollateral json 30 2 > 2IR
	# ./simCollateral json 30 3 > 3IR
	# time ./simCollateral $(folder) $(samps)
	./simCollateral $(folder) $(samps)

	# ./simCollateral json 1 50 20 0.1> out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out

highsimCtest:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateraltesthigh1 highcsim.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./simCollateraltesthigh1 json 1 1 > shortTest0High3
	./simCollateraltesthigh1 json 1 3 > shortTest2High3
	./simCollateraltesthigh1 json 1 2 > shortTest1High3
	# ./simCollateraltesthigh json 1 2 > outTest1High

	# ./simCollateraltesthigh json 1 7 > outTest7High

lowsimCtest:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateraltestlow1 lowcsim.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./simCollateraltestlow1 json 1 1 > shortTest0Low3
	./simCollateraltestlow1 json 1 3 > shortTest2Low3
	./simCollateraltestlow1 json 1 2 > shortTest1Low3
	# ./simCollateraltesthigh json 1 2 > outTest1High

	# ./simCollateraltesthigh json 1 7 > outTest7High

highFFR:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o highFFR highFFR.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	# ./highFFR json 1 1 > outTest0High
	./highFFR json 1 2 > outTestHighFFR_all
	# ./highFFR json 1 2 > outTest1HighFFR

highFFR2:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o highFFR2 highFFR2.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	# ./highFFR json 1 1 > outTest0High
	./highFFR2 json 1 2 > outTestHighFFR_high
	# ./highFFR json 1 2 > outTest1HighFFR

forcesimCtest:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateraltestforce testSimDefaults.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 	1 3 > 2IR
	# ./simCollateraltesthigh json 1 1 > outTest0High
	# ./simCollateraltestforce json 1 2 > outForceTest
	./simCollateraltestforce json 1 2 > outTestrun1_high
	# > outForced2_unreg
	# > outForced1_new
	./simCollateraltestforce json 1 3 > outTestrun2_high

forcesimCtestHigh:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateraltestforceHigh2 testSimDefaults.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	# ./simCollateraltestforceHigh json 1 3 > outTest2High
	./simCollateraltestforceHigh2 json 1 1 7 > outForced0_mid3
	./simCollateraltestforceHigh2 json 1 2 7 > outForced1_mid3
	./simCollateraltestforceHigh2 json 1 3 7 > outForced2_mid3
	./simCollateraltestforceHigh2 json 1 1 6 > outForced0_low3
	./simCollateraltestforceHigh2 json 1 2 6 > outForced1_low3
	./simCollateraltestforceHigh2 json 1 3 6 > outForced2_low3
	./simCollateraltestforceHigh2 json 1 1 8 > outForced0_high3
	./simCollateraltestforceHigh2 json 1 2 8 > outForced1_high3
	./simCollateraltestforceHigh2 json 1 3 8 > outForced2_high3

	# 	./simCollateraltestforceHigh2 json 1 1 7 > outForced0_mid2
	# ./simCollateraltestforceHigh2 json 1 2 7 > outForced1_mid2
	# ./simCollateraltestforceHigh2 json 1 3 7 > outForced2_mid2_asset1
	# ./simCollateraltestforceHigh2 json 1 1 6 > outForced0_low2_asset1
	# ./simCollateraltestforceHigh2 json 1 2 6 > outForced1_low2_asset1
	# ./simCollateraltestforceHigh2 json 1 3 6 > outForced2_low2_asset1
	# ./simCollateraltestforceHigh2 json 1 1 8 > outForced0_high2_asset1
	# ./simCollateraltestforceHigh2 json 1 2 8 > outForced1_high2_asset1
	# ./simCollateraltestforceHigh2 json 1 3 8 > outForced2_high2_asset1
simCtest:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateraltest2 csimtest.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./simCollateraltest2 json 1 1 > adjTest0
	./simCollateraltest2 json 1 3 > adjTest2
	./simCollateraltest2 json 1 2 > adjTest1
	# ./simCollateraltest1 json 1 3 > outTest
	# ./simCollateraltest json 1 2 > outTest1_1

	# ./simCollateraltest json 30 1 > 0R
	# ./simCollateraltest json 30 2 > 1R
	# ./simCollateraltest json 30 3 > 2R
	# ./simCollateraltest json 30 4 > 3R
	# valgrind --tool=memcheck --leak-check=yes ./simCollateraltest json 1 2

	 # > outTest
	# ./simCollateral json 1 50 20 0.1> out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out
	# ./testCollateral json 5 5 10 > out

simCtestP:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) creditNetTest.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCollateraltestP csimtest_payments.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./simCollateraltestP json 1 1 > adjTest0
	./simCollateraltestP json 1 3 > adjTest2
	./simCollateraltestP json 1 2 > adjTest1





simReturns:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simReturns testReturns.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./simReturns json 1 1 > outReturns

simCtestMany:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -g -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) -g $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o simCtestMany csimtest.cpp *.o $(CCLNFLAGS) -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./simCtestMany json 1 7 > outTest7
	
test:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o test test.cpp *.o $(CCLNFLAGS) -std=c++11
	./test json 1 5
	# time ./test 1 4 1 MIN_SRC_COST
	# time ./test 1 7 1 MIN_CREDIT_COST> out1
	
test_route:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -g -c CN_Node.cpp -std=c++11
	g++ -g -c CN_Edge.cpp -std=c++11
	g++ -g -c CN_Graph.cpp -std=c++11
	g++ -g -c CN_WidgetGraph.cpp -std=c++11
	g++ -g -c CN_CplexConverter.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_Solver.cpp -std=c++11
	$(CCC) -c $(CCFLAGS) $(CPPFLAGS) CN_CreditNet.cpp -std=c++11
	$(CCC) $(CCFLAGS) $(CPPFLAGS) $(CCLNDIRS) -o test_route test_route.cpp *.o $(CCLNFLAGS) -std=c++11
	time ./test_route json 1
	# time ./test 1 4 1 MIN_SRC_COST
	# time ./test 1 7 1 MIN_CREDIT_COST> out1

testFunc:
	g++ -g -c CN_Constants.cpp -std=c++11
	g++ -o testfunc testfunction.cpp *.o -std=c++11
	# ./testCollateral json'
	# ./simCollateraltest json 1 2 > 1IR
	# ./simCollateraltest json 1 1 > 0IR
	# ./simCollateraltest json 1 3 > 2IR
	./testfunc