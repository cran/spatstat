
RATFOR = /home/adrian/bin/ratfor77
#RATFOR = /usr/local/bin/ratfor
CPP = /usr/bin/cpp 

##########################################################
# Sources actually written by humans:

RAT_SRC = dppll.r inxypOld.r 

C_DOMINIC = dinfty.c dwpure.c 

C_MISC = raster.h areadiff.c closepair.c connect.c corrections.c \
       discarea.c distances.c distmapbin.c distseg.c \
       exactdist.c  exactPdist.c  \
       massdisthack.c poly2im.c trigraf.c utils.c xyseg.c 

C_MH = methas.h dist2.h areaint.c badgey.c dgs.c \
	diggra.c dist2.c fexitc.c getcif.c geyer.c \
	lookup.c methas.c stfcr.c \
	straush.c straushm.c strauss.c straussm.c

C_SRC = $(C_DOMINIC) $(C_MISC) $(C_MH)

CC_SRC = PerfectStrauss.cc

C_TEMPLATES = Kborder.cpp Kborder.h

HUMAN = $(C_TEMPLATES) $(RAT_SRC) $(C_SRC) $(CC_SRC)  Makefile

##########################################################
# Source to be generated automatically:

RAT_FOR = dppll.f inxypOld.f 

K_C = Kborder.c Kwborder.c

GENERATED = $(K_C) $(RAT_FOR) 


######################################################
###########  TARGETS ################################

target:	$(GENERATED)
	@echo -- Done -------

tar:
	tar cvf src.tar $(HUMAN)

clean:
	rm $(GENERATED)
	-rm src.tar

#######################################################
#########    RULES   ##################################


.r.f:
	$(RATFOR) -o $@ $?


Kborder.c:	Kborder.cpp
	$(CPP) -C -ansi -o $@ Kborder.cpp
	cat Kborder.h $@ > tmp
	mv tmp $@

Kwborder.c:	Kborder.cpp
	$(CPP) -C -ansi -DWEIGHTED -o $@ Kborder.cpp
	cat Kborder.h $@ > tmpw
	mv tmpw $@


