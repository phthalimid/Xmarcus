#include "xmarcus.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Xmarcus window;
    //window.setAttribute(Qt::WA_QuitOnClose);
    window.show();

    return app.exec();
}
