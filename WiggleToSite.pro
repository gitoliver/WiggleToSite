TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/bead_residues.cpp \
    src/io.cpp \
    src/main.cpp \
    src/wiggletosite.cpp

HEADERS += \
    includes/bead_residues.h \
    includes/io.h \
    includes/wiggleTosite.hpp
