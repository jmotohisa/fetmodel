bin_PROGRAMS = fetmodel

SPECIFICATION_FILE = fetmodel.scm
SPECIFICATION_FILE_EXTRA = 

EXTRA_DIST = $(SPECIFICATION_FILE_EXTRA)

pkgdata_DATA = $(SPECIFICATION_FILE) $(SPECIFICATION_FILE_EXTRA)

LIBFETMODEL = $(top_builddir)/src/libfetmodel.la

MY_INCLUDES= 
MY_SOURCES = cfet_libctl.c ballisticFET1d_libctl.c

MY_LIBS    =  # extra libs you need go here
MY_LDFLAGS =         # extra -L flags go here
MY_CPPFLAGS =        # extra -I include flags go here
MY_DEFS =            # extra -D define flags go here

# what is printed out when invoking your program with --version:
VERSION_STRING = "fetmodel @PACKAGE_VERSION@, copyright (c) 2004-2019 by Junichi Motohisa."
# The following variables should be detected and set by autoconf:

# libctl install. dir., e.g. /usr/local/share/libctl
LIBCTL_DIR = @LIBCTL_DIR@

# gen-ctl-io program
GEN_CTL_IO = @GEN_CTL_IO@

# location of libctl headers and library for this project
LIBCTL_CPPFLAGS = -I$(top_srcdir)/utils -I$(top_builddir)/src
LIBCTL = -lctl # change to -lctl for your project

##############################################################################
#   don't (normally) edit below except to replace "example" with your name   #
##############################################################################

CTL_DEFS = -DCTL_SCM='"'$(LIBCTL_DIR)/base/ctl.scm'"' \
       -DINCLUDE_SCM='"'$(LIBCTL_DIR)/base/include.scm'"' \
          -DSPEC_SCM='"'$(pkgdatadir)/$(SPECIFICATION_FILE)'"' \
    -DVERSION_STRING='"'$(VERSION_STRING)'"'

fetmodel_SOURCES = $(MY_SOURCES)
nodist_fetmodel_SOURCES = main.c ctl-io.h ctl-io.c
# do not use geom.c 
BUILT_SOURCES = $(nodist_fetmodel_SOURCES)
fetmodel_LDADD = $(LIBFETMODEL) $(MY_LIBS) $(LIBCTL)
fetmodel_LDFLAGS = $(MY_LDFLAGS)
fetmodel_CPPFLAGS = $(MY_CPPFLAGS) $(MY_DEFS) $(CTL_DEFS) $(LIBCTL_CPPFLAGS)

main.c: $(LIBCTL_DIR)/base/main.c
	cp -f $(LIBCTL_DIR)/base/main.c $@

# geom.c: $(LIBCTL_DIR)/utils/geom.c
# 	cp -f $(LIBCTL_DIR)/utils/geom.c $@

ctl-io.c: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm
	$(GEN_CTL_IO) --code -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

ctl-io.h: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm
	$(GEN_CTL_IO) --header -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

clean-local:
	rm -f ctl-io.* main.* geom.*

