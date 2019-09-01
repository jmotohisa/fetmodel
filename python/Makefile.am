BUILT_SOURCES = pyfet-python.c __init__.py
EXTRA_DIST = $(BUILT_SOURCES)
CLEANFILES = $(BUILT_SOURCES) pyfet.py .coverage

AM_CPPFLAGS = -I$(top_srcdir)/src               \
              -I$(top_builddir) # for config.h

LIBFETMODEL = $(top_builddir)/src/libfetmodel.la

_pyfet_la_SOURCES = pyfet.c pyfet-python.c
_pyfet_la_LIBADD = $(LIBFETMODEL) $(PYTHON_LIBS)
_pyfet_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
# _pyfet_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@ -export-symbols-regex initmyext
_pyfet_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

pkgpyexec_LTLIBRARIES = _pyfet.la

SWIG ?= swig
# SWIG_VERSION = $($(SWIG) -version | grep Version | awk '{print $$3}')
SWIG_VERSION = $(shell $(SWIG) -version | grep Version | awk '{print $$3}')
PYFET_SWIG_SRC = pyfet.i numpy.i

PYFET_SWIG_HDRS =             \
 $(top_srcdir)/src/ccm.h         \
 $(top_srcdir)/src/density1d.h   \
 $(top_srcdir)/src/density2d.h   \
 $(top_srcdir)/src/ballistic.h   \
 $(top_srcdir)/src/capacitor.h


pyfet-python.c: $(PYFET_SWIG_SRC) $(PYFET_SWIG_HDRS)
	$(SWIG) -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -nofastunpack -python -o $@ $(srcdir)/pyfet.i

pyfet.py: pyfet-python.c

__init__.py: pyfet.py
	cp $< $@
	echo "__version__ = '$(shell git describe --tags | sed 's/^v//')'" >> $@
	if [[ "${SWIG_VERSION}" = 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' $@; \
		sed -i.bak 's/    import _pyfet/from . import _pyfet/' $@; \
	fi

INIT_PY = __init__.py

######################################################################
# specification of python source files to be byte-compiled at installation
######################################################################
HL_IFACE =

pkgpython_PYTHON = __init__.py $(HL_IFACE)


######################################################################
# finally, specification of what gets installed in the pyfet python
# module directory of the python site-packages installation
######################################################################
PY_PKG_FILES = $(INIT_PY) $(HL_IFACE) .libs/_pyfet.so

pyfet: _pyfet.la __init__.py $(HL_IFACE)
	mkdir -p pyfet
	cp $(PY_PKG_FILES) pyfet

all-local: pyfet

clean-local:
	rm -rf pyfet __init__.py.bak
