#include <stdio.h>
#include <unistd.h>

FILE *rosella_minimap2_out;

int rosella_minimap2_main(int argc, char *argv[]);

int rosella_minimap2_puts(const char *line)
{
	if (fputs(line, rosella_minimap2_out) == EOF) return EOF;
	return fputc('\n', rosella_minimap2_out);
}

int rosella_minimap2(int argc, char *argv[], int out)
{
	int status;
	rosella_minimap2_out = fdopen(out, "wb");
	if (rosella_minimap2_out == NULL) {
		close(out);
		return -1;
	}
	status = rosella_minimap2_main(argc, argv);
	if (fclose(rosella_minimap2_out) == EOF && status == 0) status = 1;
	rosella_minimap2_out = NULL;
	return status;
}
