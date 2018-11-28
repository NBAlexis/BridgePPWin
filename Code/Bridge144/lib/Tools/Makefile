# rebuild evalexpr_parser

.SUFFIXES: .cc .h .hh .yy .o

CXX=		g++
CXX_FLAGS=	-g -Wall
LD=		$(CXX)
LD_FLAGS=
YACC=		bison
#YACC_FLAGS=	-v
YACC_FLAGS=
RM=		rm -f

EXTRA_INCLUDES=
EXTRA_LIBS=	-lm

TARGET=	evalexpr_parser.h evalexpr_parser.cpp
TARGET_EXTRA= stack.hh position.hh location.hh

all: $(TARGET)

evalexpr_parser.cpp: evalexpr_parser.yy
	$(YACC) $(YACC_FLAGS) --defines=evalexpr_parser.h -o $@ $<

evalexpr_parser.h: evalexpr_parser.cpp


realclean: 
	$(RM) $(TARGET) $(TARGET_EXTRA)

