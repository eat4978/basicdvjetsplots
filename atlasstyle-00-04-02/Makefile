
ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)


OBJ = AtlasLabels.o AtlasStyle.o AtlasUtils.o TestLabel.o

ATLASOBJ = AtlasExample.o $(OBJS)
BASICOBJ = BasicExample.o $(OBJS)

all : atlastest basictest

.SUFFIXES: .C .o

.C.o :
	g++ $(ROOTINC) -c $<

atlastest : $(OBJ)
	g++  -o $@ $(OBJ) $(ROOTLIB)

basictest : $(OBJ)
	g++  -o $@ $(OBJ) $(ROOTLIB)

archive : clean
	cd .. ; tar -czf atlasstyle-00-04-00.tgz atlasstyle-00-04-00		  

clean : 
	rm -rf *.o *~
