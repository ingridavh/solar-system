TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    celestialbody.cpp \
    solarsystem.cpp \
    vec3.cpp \
    solver.cpp

HEADERS += \
    celestialbody.h \
    solarsystem.h \
    vec3.h \
    solver.h \
    catch.hpp

DISTFILES += \
    ../project3_benchmarks.txt \
    ../build-SolarSystem-Desktop-Debug/positions.txt

