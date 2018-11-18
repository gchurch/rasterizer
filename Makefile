FILE=rasterizer

########
#   Directories
S_DIR=src
B_DIR=build

########
#   Output
EXEC1=$(B_DIR)/$(FILE)

# default build settings
CC_OPTS=-c -pipe -Wall -Wno-switch -ggdb -g3 -Ofast
LN_OPTS=
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


########
#   Objects
$(OBJ1) : $(S_DIR)/$(FILE).cpp $(S_DIR)/SDLauxiliary.h $(S_DIR)/TestModel.h
	$(CC) $(CC_OPTS) $(S_DIR)/$(FILE).cpp -o $(OBJ1) $(SDL_CFLAGS) $(GLM_CFLAGS)


########
#   Main build rule     
Build : $(OBJ1) Makefile
	$(CC) $(LN_OPTS) -o $(EXEC1) $(OBJ1) $(SDL_LDFLAGS)

run : $(EXEC1)
	./$(EXEC1)

clean:
	rm -f $(B_DIR)/* 
