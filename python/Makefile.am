BUILT_SOURCES = fetmodel-python.c __init__.py
EXTRA_DIST = $(BUILT_SOURCES)
CLEANFILES = $(BUILT_SOURCES) fetmodel.py .coverage

AM_CPPFLAGS = -I$(top_srcdir)/src               \
              -I$(top_builddir) # for config.h

LIBFETMODEL = $(top_builddir)/src/libfetmodel.la

_fetmodel_la_SOURCES = fetmodel.c fetmodel-python.c
_fetmodel_la_LIBADD = $(LIBFETMODEL) $(PYTHON_LIBS)
_fetmodel_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@
# _fetmodel_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@ -export-symbols-regex initmyext
_fetmodel_la_CPPFLAGS = $(PYTHON_INCLUDES) $(AM_CPPFLAGS)

pkgpyexec_LTLIBRARIES = _fetmodel.la

SWIG ?= swig
# SWIG_VERSION = $($(SWIG) -version | grep Version | awk '{print $$3}')
SWIG_VERSION = $(shell $(SWIG) -version | grep Version | awk '{print $$3}')
FETMODEL_SWIG_SRC = fetmodel.i numpy.i

FETMODEL_SWIG_HDRS =             \
 $(top_srcdir)/src/ccm.h         \
 $(top_srcdir)/src/density1d.h   \
 $(top_srcdir)/src/density2d.h   \
 $(top_srcdir)/src/ballistic.h   \
 $(top_srcdir)/src/capacitor.h


fetmodel-python.c: $(FETMODEL_SWIG_SRC) $(FETMODEL_SWIG_HDRS)
	$(SWIG) -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -nofastunpack -python -o $@ $(srcdir)/fetmodel.i

fetmodel.py: fetmodel-python.c

__init__.py: fetmodel.py
	cp $< $@
	echo "__version__ = '$(shell git describe --tags | sed 's/^v//')'" >> $@
	if [[ "${SWIG_VERSION}" = 3.0.12 ]]; then \
		sed -i.bak '/^if _swig_python_version_info >= (2, 7, 0):/,/^else:/d' $@; \
		sed -i.bak 's/    import _fetmodel/from . import _fetmodel/' $@; \
	fi

INIT_PY = __init__.py

######################################################################
# specification of python source files to be byte-compiled at installation
######################################################################
HL_IFACE =

pkgpython_PYTHON = __init__.py $(HL_IFACE)


######################################################################
# finally, specification of what gets installed in the fetmodel python
# module directory of the python site-packages installation
######################################################################
PY_PKG_FILES = $(INIT_PY) $(HL_IFACE) .libs/_fetmodel.so

fetmodel: _fetmodel.la __init__.py $(HL_IFACE)
	mkdir -p fetmodel
	cp $(PY_PKG_FILES) fetmodel

all-local: fetmodel

clean-local:
	rm -rf fetmodelx __init__.py.bak
