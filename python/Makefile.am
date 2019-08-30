LIBFETMODEL = $(top_builddir)/src/libfetmodel.la

SWIG ?= swig

PYFET_SWIG_SRC = pyfet.i numpy.i

pyfet-python.c: $(PYFET_SWIG_SRC) $(HPPFILES)
        $(SWIG) -Wextra $(AM_CPPFLAGS) -outdir $(builddir) -nofastunpack -python -o $@ $(srcdir)/pyfet.i

_pyfet_la_SOURCES = pyfet.c
_pyfet_la_LIBADD = $(PYTHON_LIBS) $(top_builddir)/lib/libhollow.la
_pyfet_la_CPPFLAGS = $(PYTHON_CFLAGS) $(PYTHON_INCLUDES) $(AM_CPPFLAGS)
_pyfet_la_LDFLAGS = -module -version-info @SHARED_VERSION_INFO@ -export-symbols-regex initmyext
