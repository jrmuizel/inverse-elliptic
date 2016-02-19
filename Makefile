LDLIBS=-lm
CFLAGS=-g
all: draw-ellipse xigel
draw-ellipse: draw-ellipse.c cel.c inverse.c
xigel: xigel.c cel.c inverse.c
