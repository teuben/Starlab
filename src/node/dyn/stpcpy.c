#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

char *stpcpy(char *restrict dest, const char *restrict src) {
  while (*dest++ = *src++);
  return dest-1;
}

#ifdef __cplusplus
}
#endif
