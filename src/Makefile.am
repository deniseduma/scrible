bin_PROGRAMS = hashfilter
hashfilter_SOURCES = Pool4 classIO.cpp hashfilter.cpp parserLexer.o \
	parserYacc.cpp prov


# set the include path found by configure
INCLUDES= $(all_includes)

# the library search path.
hashfilter_LDFLAGS = $(all_libraries) -fopenmp
noinst_HEADERS = classIO.h conf.h const.h hashChaining.h operation.h parser.h
