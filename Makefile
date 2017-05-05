FILE=rasterizer
CL=clipping

########
#   Directories
S_DIR=src
B_DIR=build

########
#   Output
EXEC1=$(B_DIR)/$(FILE)
EXEC2=$(B_DIR)/$(FILE)_$(CL)

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


########
#   Objects
$(B_DIR)/$(FILE).o : $(S_DIR)/$(FILE).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE).o $(S_DIR)/$(FILE).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)

$(B_DIR)/$(FILE)_$(CL).o : $(S_DIR)/$(FILE)_$(CL).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) -o $(B_DIR)/$(FILE)_$(CL).o $(S_DIR)/$(FILE)_$(CL).cpp $(SDL_CFLAGS) $(GLM_CFLAGS)


########
#   Main build rule     
Build : $(OBJ1) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC1) $(OBJ1) $(SDL_LDFLAGS)

clipping : $(OBJ2) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC2) $(OBJ2) $(SDL_LDFLAGS)



clean:
	rm -f $(B_DIR)/* 
