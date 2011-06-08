#ifndef kbhit_h
#define kbhit_h

#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/time.h>
#include <termios.h>

static int kbhit(int fd);

#endif // kbhit_h
