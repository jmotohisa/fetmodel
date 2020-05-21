SUBDIRS = src

if WITH_LIBCTL
LIBCTL = libctl
SUBDIRS += libctl
else
LIBCTL =
endif

if WITH_PYTHON
SUBDIRS += python
endif

ACLOCAL_AMFLAGS =  -I m4

EXTRA_DIST = AUTHORS COPYING ChangeLog INSTALL NEWS README TODO

libfetmodel:
	cd src && $(MAKE)
