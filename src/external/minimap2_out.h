#include <stdio.h>

extern FILE *rosella_minimap2_out;
int rosella_minimap2_puts(const char *line);

#undef stdout
#define stdout rosella_minimap2_out
#undef puts
#define puts(line) rosella_minimap2_puts(line)
#undef printf
#define printf(...) fprintf(rosella_minimap2_out, __VA_ARGS__)
