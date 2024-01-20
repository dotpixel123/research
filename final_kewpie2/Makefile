# paths
DSRC = ./Librairies/Sources
DINC = ./Librairies/Headers
DOBJ = ./Librairies/Objects
DEXE = ./Executable
EXEN = KEWPIE2.exe

# flags
FLAGS = -Wall -ansi -O3 -I$(DINC)
LIBS = -lm -lgsl -lgslcblas

# commands for compilation
CCL = g++ -o
CC = g++ $(FLAGS) -c -g

# objects
OBJECTS = $(DOBJ)/Element.o $(DOBJ)/Leveldensity.o $(DOBJ)/Readinput.o \
$(DOBJ)/InitMassArray.o $(DOBJ)/Writeoutput.o $(DOBJ)/Fission.o $(DOBJ)/Potential.o \
$(DOBJ)/Transcoef.o $(DOBJ)/Evaporation.o $(DOBJ)/Decay.o $(DOBJ)/Particle.o \
$(DOBJ)/Cascade.o $(DOBJ)/Synthesis.o $(DOBJ)/Gamma.o $(DOBJ)/Fusion.o $(DOBJ)/main.o

main: $(DEXE)/$(EXEN)

$(DEXE)/$(EXEN): $(OBJECTS)
	$(CCL) $@ $(OBJECTS) $(LIBS)

# precompiled headers
$(DOBJ)/%.o: $(DSRC)/%.cpp
	$(CC) -o $@ $< $(LIBS)

clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/*.exe

run:
	cd $(DEXE); \
	./$(EXEN)

cleansave:
	rm -rf ./*~ $(DSRC)/*.cpp~ $(DINC)/*.h~
