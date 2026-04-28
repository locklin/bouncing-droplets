CC = gcc
CFLAGS = -O2 -Wall -Ilib
LIBS = lib/libraylib.a -lGL -lm -lpthread -ldl -lrt -lX11 -lXrandr -lXinerama -lXcursor -lXi

all: lib/libraylib.a walker oza harmonic billiard schrodinger

core.o: core.c core.h
	$(CC) $(CFLAGS) -c -o $@ core.c

walker: main.c core.o core.h
	$(CC) $(CFLAGS) -o $@ main.c core.o $(LIBS)

oza: oza.c core.o core.h
	$(CC) $(CFLAGS) -o $@ oza.c core.o $(LIBS)

harmonic: harmonic.c core.o core.h
	$(CC) $(CFLAGS) -o $@ harmonic.c core.o $(LIBS)

billiard: billiard.c core.o core.h
	$(CC) $(CFLAGS) -o $@ billiard.c core.o $(LIBS)

schrodinger: schrodinger.c core.o core.h
	$(CC) $(CFLAGS) -o $@ schrodinger.c core.o $(LIBS)

lib/libraylib.a:
	cd lib/raylib-src && $(CC) -c -O2 -DPLATFORM_DESKTOP \
		-DGRAPHICS_API_OPENGL_33 -D_GLFW_X11 \
		rcore.c rshapes.c rtextures.c rtext.c rmodels.c utils.c raudio.c rglfw.c \
		-I. -Iexternal/glfw/include && \
	ar rcs ../libraylib.a rcore.o rshapes.o rtextures.o rtext.o rmodels.o utils.o raudio.o rglfw.o && \
	rm -f *.o

clean:
	rm -f walker oza harmonic billiard schrodinger core.o

clean-all: clean
	rm -f lib/libraylib.a

.PHONY: all clean clean-all
