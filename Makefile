LDLIBS=-lm
CFLAGS=-g -Wall
all: draw-ellipse xigel
draw-ellipse: draw-ellipse.c cel.c inverse.c
xigel: xigel.c cel.c inverse.c
