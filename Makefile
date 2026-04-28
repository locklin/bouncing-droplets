CC = gcc
CFLAGS = -O2 -Wall -Ilib
LIBS = lib/libraylib.a -lGL -lm -lpthread -ldl -lrt -lX11 -lXrandr -lXinerama -lXcursor -lXi

all: lib/libraylib.a walker oza harmonic billiard schrodinger

walker: main.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

oza: oza.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

harmonic: harmonic.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

billiard: billiard.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

schrodinger: schrodinger.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

lib/libraylib.a:
	cd lib/raylib-src && $(CC) -c -O2 -DPLATFORM_DESKTOP \
		-DGRAPHICS_API_OPENGL_33 -D_GLFW_X11 \
		rcore.c rshapes.c rtextures.c rtext.c rmodels.c utils.c raudio.c rglfw.c \
		-I. -Iexternal/glfw/include && \
	ar rcs ../libraylib.a rcore.o rshapes.o rtextures.o rtext.o rmodels.o utils.o raudio.o rglfw.o && \
	rm -f *.o

clean:
	rm -f walker oza harmonic billiard schrodinger

clean-all: clean
	rm -f lib/libraylib.a

.PHONY: all clean clean-all
