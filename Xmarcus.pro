#-------------------------------------------------
#
# Project created by QtCreator 2018-07-05T13:14:07
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Xmarcus
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++14
CONFIG += staticlib

SOURCES += \
        main.cpp \
        xmarcus.cpp \
        ecdata.cpp \
        qcustomplot.cpp \
        doubledelegate.cpp \
        species.cpp \
        reaction.cpp \
        numericalparams.cpp \
        speciestabledelegate.cpp \
        choosereaction.cpp \
        reactionstabledelegate.cpp \
        workersimcv.cpp \
        simwindow.cpp

HEADERS += \
        xmarcus.h \
        ecdata.h \
        qcustomplot.h \
        doubledelegate.h \
        species.h \
        reaction.h \
        numericalparams.h \
        speciestabledelegate.h \
        choosereaction.h \
        reactionstabledelegate.h \
        workersimcv.h \
        simwindow.h

FORMS += \
        xmarcus.ui \
        numericalparams.ui \
        choosereaction.ui \
        simwindow.ui

RESOURCES += \
        xmarcus.qrc

macx {
INCLUDEPATH += /opt/local/include
LIBS += -L/opt/local/lib -lgsl -lgslcblas -lm
}
win32 {
INCLUDEPATH += c:/msys64/mingw64/include
LIBS += -Lc:/msys64/mingw64/lib -lgsl -lgslcblas -lm
}
unix {
INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm
}
