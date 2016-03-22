#-------------------------------------------------
#
# Project created by QtCreator 2011-04-18T13:50:29
#
#-------------------------------------------------

QT       += core gui

TARGET = psurve_v1
TEMPLATE = app

QMAKE_CXXFLAGS_RELEASE += -O3


SOURCES += main.cpp\
        mainwindow.cpp \
    konstruktory.cpp \
    newitem.cpp \
    agree.cpp \
    countDispersion_BS.cpp \
    countDispersion_B.cpp \
    countDispersion_bulk_B.cpp \
    drawStructure.cpp \
    countAmplitudes_B.cpp \
    countAmplitudes_bulk_B.cpp \
    countAmplitudes_BS.cpp \
    countBCDet_B.cpp \
    countBCDet_BS.cpp \
    countBCDet_bulk_B.cpp \
    linpack_z.cpp \
    blas1_z.cpp

HEADERS  += mainwindow.h \
    newitem.hpp \
    agree.hpp \
    spectrometer.hpp \
    medium.hpp \
    linpack_z.hpp \
    blas1_z.hpp

FORMS    += mainwindow.ui \
    newitem.ui \
    agree.ui

LIBS += -lkdeui \
    -lgsl \
    -lgslcblas \
