FILE=rasterizer
CL=clipping
TE=texture
ALL=all

########
#   Directories
S_DIR=src
B_DIR=build

########
#   Output
EXEC1=$(B_DIR)/$(FILE)
EXEC2=$(B_DIR)/$(FILE)_$(CL)
EXEC3=$(B_DIR)/$(FILE)_$(TE)
EXEC4=$(B_DIR)/$(FILE)_$(ALL)

# default build settings
CC_OPTS=-c -pipe -Wall -Wno-switch -ggdb -g3 
LN_OPTS=-Ofast
CC=g++

########
#       SDL options
SDL_CFLAGS := $(shell sdl-config --cflags)
GLM_CFLAGS := -I$(GLMDIR)
SDL_LDFLAGS := $(shell sdl-config --libs)

########
#   This is the default action
all:Build


########
#   Object list
#
OBJ1 = $(B_DIR)/$(FILE).o
OBJ2 = $(B_DIR)/$(FILE)_$(CL).o
OBJ3 = $(B_DIR)/$(FILE)_$(TE).o
OBJ4 = $(B_DIR)/$(FILE)_$(ALL).o


########
#   Objects
$(B_DIR)/$(FILE).o : $(S_DIR)/$(FILE).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

$(B_DIR)/$(FILE)_$(CL).o : $(S_DIR)/$(FILE)_$(CL).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE)_$(CL).o $(S_DIR)/$(FILE)_$(CL).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

$(B_DIR)/$(FILE)_$(TE).o : $(S_DIR)/$(FILE)_$(TE).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel_texture.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE)_$(TE).o $(S_DIR)/$(FILE)_$(TE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

$(B_DIR)/$(FILE)_$(ALL).o : $(S_DIR)/$(FILE)_$(ALL).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel_texture.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE)_$(ALL).o $(S_DIR)/$(FILE)_$(ALL).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)



########
#   Main build rule     
Build : $(OBJ1) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC1) $(OBJ1) $(SDL_LDFLAGS)

clipping : $(OBJ2) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC2) $(OBJ2) $(SDL_LDFLAGS)

texture : $(OBJ3) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC3) $(OBJ3) $(SDL_LDFLAGS)

all : $(OBJ4) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC4) $(OBJ4) $(SDL_LDFLAGS)

run : $(EXEC4)
	./$(EXEC4)



clean:
	rm -f $(B_DIR)/* 
